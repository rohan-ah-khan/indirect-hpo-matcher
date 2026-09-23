#!/usr/bin/env python3
"""Match patient HPO terms to indirect candidate terms."""

from __future__ import annotations

import math
import re
from collections import defaultdict, deque
from pathlib import Path
from typing import AbstractSet, Iterable, Mapping

import hpotk
import pandas as pd

# Used to restrict candidate discovery to nearby phenotype terms.
PHENOTYPIC_ABNORMALITY = "HP:0000118"  # Main HPO phenotype root.
MAX_DISTANCE = 3  # Maximum patient-to-candidate edges.

# Reduce direct-relation and sibling scores as graph distance increases.
SCORE_STEP = 0.15  # Score reduction per edge.

# Used in all term-label comparisons.
TEXT_STOPWORDS = {"a", "an", "and", "in", "of", "or", "the", "to", "with"}

# Reduce labels to their underlying clinical feature.
_SEED_MODIFIERS = frozenset(
    """
    acute adult antenatal bilateral childhood chronic congenital developmental
    distal early episodic fetal focal generalized generalised infantile juvenile
    late mild moderate neonatal nonprogressive onset pediatric perinatal
    persistent postlingual postnatal prelingual prenatal profound progressive
    proximal recurrent segmental severe subacute transient unilateral young
    """.split()
)

# Used with the modifiers above to remove specimen and measurement wording.
_SEED_MEASUREMENT_WORDS = frozenset(
    """
    anti antibodies antibody activities activity circulating concentration
    content csf cultured decreased diminished elevated fibroblast fibroblasts
    hepatic increased level levels mitochondrial plasma positive positivity
    reduced serum tissue urinary urine
    """.split()
)

# Remove opposing magnitude or direction wording.
_AXIS_PREFIXES = (
    "hyper", "hypo", "macro", "micro", "megalo", "mega", "brady", "tachy",
    "oligo", "poly", "an", "a",
)
_AXIS_SUFFIXES = ("cytosis", "cytopenia", "penia", "philia", "plegia", "paresis")

# Used to require distinctive wording from terms passing the wording-supported sibling gate.
_FINDING_WORDS = frozenset(
    """
    abnormal abnormality abnormalities anomaly anomalies malformation
    morphology defect defects dysfunction dysplasia deformity
    deficiency insufficiency failure disease disorder syndrome
    increased decreased elevated reduced diminished excess excessive
    high low level levels concentration content amount activity
    presence absence absent ratio proportion number count
    aplasia hypoplasia hypoplastic hyperplasia duplication
    """.split()
)

# Initial shared-ancestor and final-evidence sibling checks.
SIBLING_SHARED_SPECIFICITY_MIN = 0.20  # Minimum shared-ancestor specificity.
MIN_SIBLING_EVIDENCE_SCORE = 0.15  # Minimum final ordinary-sibling score.

# Direct siblings supported by similar clinical wording.
WORDING_SUPPORT_CLINICAL_MIN = 0.50  # Minimum label similarity.
WORDING_SUPPORT_SEMANTIC_MIN = 0.45  # Minimum annotation similarity.

# Differently worded direct siblings with exceptionally similar annotations.
PROFILE_RESCUE_SEMANTIC_MIN = 0.90  # Minimum annotation similarity.
PROFILE_RESCUE_SPECIFICITY_MIN = 0.35  # Minimum shared-ancestor specificity.

# Decide whether a candidate's genes should be added.
MAX_GENES_PER_CANDIDATE_TERM = 2000  # Absolute direct-gene limit.
TARGET_EXPANSION_MIN_IC = 0.22  # Informativeness required for larger gene sets.
TARGET_EXPANSION_SMALL_GENES = 100  # Gene sets at or below this size bypass the IC limit.

# Used to suppress broad patient expansions and high-level candidate terms.
BROAD_SOURCE_MAX_ROOT_DISTANCE = 2  # Distance that identifies a broad patient term.
BROAD_SOURCE_MAX_INDIRECT_GENES = 2000  # Limit for a broad source expansion or candidate branch.


def _ancestors(data: dict, term_id: str) -> set[str]:
    """Return the term and all of its ancestors."""
    # Reuse earlier traversals because the same terms are visited during propagation and scoring.
    cache = data["ancestor_cache"]
    if term_id in cache:
        return cache[term_id]

    # Include the starting term so its direct annotations are retained during propagation.
    result = {term_id}
    queue = deque([term_id])

    # Follow each term's direct parents until every broader ancestor has been visited.
    while queue:
        current = queue.popleft()
        for parent in data["parents"].get(current, set()):
            if parent not in result:
                result.add(parent)
                queue.append(parent)

    # Store the completed set so subsequent lookups do not traverse the graph again.
    cache[term_id] = result
    return result


def _limited_distances(start: str, edges: Mapping[str, set[str]], maximum: int) -> dict[str, int]:
    """Return shortest edge distances up to a maximum."""
    # Breadth-first traversal visits terms in increasing distance from the starting term.
    result = {start: 0}
    queue = deque([start])
    while queue:
        current = queue.popleft()

        # Do not follow another edge after reaching the requested search distance.
        if result[current] >= maximum:
            continue
        for neighbour in edges.get(current, set()):
            distance = result[current] + 1

            # Keep the shortest route when HPO contains multiple paths to the same term.
            if neighbour not in result or distance < result[neighbour]:
                result[neighbour] = distance
                queue.append(neighbour)
    return result


def _information_content(term_id: str, propagated: Mapping[str, set[str]], total: int) -> float:
    """Calculate information content on a zero-to-one scale."""
    # Propagated annotations include items assigned to the term or any of its descendants.
    count = max(1, len(propagated.get(term_id, set())))
    if total <= 1:
        return 1.0

    # Rare, specific annotations approach one; common, broad annotations approach zero.
    return min(1.0, max(0.0, math.log(total / count) / math.log(total)))


def _target_informativeness(data: dict, term_id: str) -> float:
    """Average the term's propagated gene and disease information content."""
    # Using both indexes prevents a term from appearing specific in only one annotation source.
    return (
        _information_content(term_id, data["propagated_genes"], data["total_gene_count"])
        + _information_content(term_id, data["propagated_diseases"], data["total_disease_count"])
    ) / 2.0


def _worth_expanding(data: dict, term_id: str) -> bool:
    """Return whether a candidate's direct gene group is useful to add."""
    # Larger groups must be informative; small groups are allowed despite a lower IC value.
    return (
        _target_informativeness(data, term_id) >= TARGET_EXPANSION_MIN_IC
        or len(data["direct_genes"].get(term_id, ())) <= TARGET_EXPANSION_SMALL_GENES
    )


def _shared_specificity(data: dict, term_id: str) -> float:
    """Combine annotation information content and ontology breadth."""
    # Shared ancestors are evaluated repeatedly while comparing nearby terms.
    cache = data["specificity_cache"]
    if term_id in cache:
        return cache[term_id]

    # Annotation specificity is high when relatively few genes and diseases map below the term.
    annotation_specificity = _target_informativeness(data, term_id)

    # Ontology specificity is high when the term has relatively few descendants.
    descendant_count = len(_limited_distances(term_id, data["children"], len(data["terms"]))) - 1
    descendant_specificity = 1.0 - (math.log(descendant_count + 1) / math.log(max(2, len(data["terms"]))))
    # Annotation specificity contains separate gene and disease signals.
    cache[term_id] = (2.0 * annotation_specificity + descendant_specificity) / 3.0
    return cache[term_id]


def _labels(data: dict, term_id: str) -> tuple[str, ...]:
    """Return a term's primary name followed by its synonyms."""
    term = data["terms"][term_id]
    return (term["name"],) + term["synonyms"]


def _similar_token(left: str, right: str) -> bool:
    """Match identical words or words with a near-complete shared stem."""
    # Exact words require no approximate stem comparison.
    if left == right:
        return True

    # Short partial words produce too many accidental matches.
    if len(left) < 5 or len(right) < 5:
        return False

    # Count the common characters from the beginning of both words.
    common = 0
    for a, b in zip(left, right):
        if a != b:
            break
        common += 1
    # Five characters avoids accidental short roots such as hydro-.
    return common >= 5 and common >= min(len(left), len(right)) - 2


def _clinical_tokens(label: str) -> set[str]:
    """Label tokens minus qualifier and assay vocabulary."""
    # Removing age, severity, laterality, specimen, and measurement words reveals the finding.
    tokens = set(re.findall(r"[a-z0-9]+", label.lower()))
    return tokens - TEXT_STOPWORDS - _SEED_MODIFIERS - _SEED_MEASUREMENT_WORDS


def _axis_stem(token: str) -> str:
    """Remove a fused direction or magnitude affix."""
    # Opposing prefixes can describe different directions of the same measured feature.
    for prefix in _AXIS_PREFIXES:
        if token.startswith(prefix) and len(token) - len(prefix) >= 5:
            return token[len(prefix):]

    # Axis suffixes are treated similarly when the remaining stem is long enough.
    for suffix in _AXIS_SUFFIXES:
        if token.endswith(suffix) and len(token) - len(suffix) >= 4:
            # Use a non-word stem so suffix-derived and prefix-derived axes cannot collide.
            return token[: -len(suffix)] + "\x00"
    return token


def _same_clinical_axis(data: dict, patient_hpo: str, candidate_hpo: str) -> bool:
    """Return whether two terms reduce to the same clinical core."""
    # Compare every primary-name and synonym representation of the patient term.
    patient_cores = []
    for label in _labels(data, patient_hpo):
        patient_cores.append(frozenset(_clinical_tokens(label)))

    # Build the equivalent set of clinical cores for the candidate term.
    candidate_cores = []
    for label in _labels(data, candidate_hpo):
        candidate_cores.append(frozenset(_clinical_tokens(label)))

    for patient_core in patient_cores:
        if not patient_core:
            continue

        # Remove fused direction or magnitude affixes before the second comparison.
        patient_axis = frozenset(map(_axis_stem, patient_core))
        for candidate_core in candidate_cores:
            if not candidate_core:
                continue

            # A complete token match after qualifier removal identifies the same feature.
            if _token_retention(patient_core, candidate_core) >= 0.999:
                return True
            candidate_axis = frozenset(map(_axis_stem, candidate_core))

            # Also accept matching stems when an axis affix was removed from either label.
            if patient_axis == candidate_axis:
                return True
    return False


def _token_retention(left_tokens: AbstractSet[str], right_tokens: AbstractSet[str]) -> float:
    """Symmetric fraction of tokens on each side matched on the other side."""
    if not left_tokens or not right_tokens:
        return 0.0

    # Count how many left-side words have an exact or near-stem match on the right.
    left_matches = 0
    for left_token in left_tokens:
        for right_token in right_tokens:
            if _similar_token(left_token, right_token):
                left_matches += 1
                break

    # Repeat in the opposite direction so extra words on either side reduce similarity.
    right_matches = 0
    for right_token in right_tokens:
        for left_token in left_tokens:
            if _similar_token(right_token, left_token):
                right_matches += 1
                break

    left_retention = left_matches / len(left_tokens)
    right_retention = right_matches / len(right_tokens)

    # The lower direction prevents a short label from being contained in a dissimilar long one.
    return min(left_retention, right_retention)


def _clinical_similarity(data: dict, left: str, right: str, drop: AbstractSet[str] = frozenset(),) -> float:
    """Return the best token similarity across names and synonyms."""
    best_similarity = 0.0

    # Test every name/synonym pairing because equivalent HPO wording may be in a synonym.
    for left_label in _labels(data, left):
        left_tokens = _clinical_tokens(left_label) - drop
        for right_label in _labels(data, right):
            right_tokens = _clinical_tokens(right_label) - drop
            similarity = _token_retention(left_tokens, right_tokens)
            best_similarity = max(best_similarity, similarity)
    return best_similarity


def _semantic_similarity(data: dict, patient_hpo: str, candidate_hpo: str, shared_ancestor: str,) -> float:
    """Lin-style shared-information similarity, averaged over genes and diseases."""
    values = []

    # Calculate the same ontology-profile similarity independently for genes and diseases.
    for propagated_annotations, total_annotations in (
        (data["propagated_genes"], data["total_gene_count"]),
        (data["propagated_diseases"], data["total_disease_count"]),
    ):
        # The shared ancestor represents the annotation information common to both terms.
        shared_ic = _information_content(shared_ancestor, propagated_annotations, total_annotations)

        # Compare shared information with the combined information in the two individual terms.
        denominator = _information_content(
            patient_hpo, propagated_annotations, total_annotations
        ) + _information_content(
            candidate_hpo, propagated_annotations, total_annotations
        )
        values.append(2.0 * shared_ic / denominator if denominator else 0.0)

    # Average both annotation sources and keep the result on a zero-to-one scale.
    return max(0.0, min(1.0, sum(values) / len(values)))


def _relationship_path(data: dict, patient_hpo: str, candidate_hpo: str) -> tuple[int, str, int, int] | None:
    """Return the best path through a shared ancestor."""
    # Find each term's reachable ancestors and their upward distances.
    patient_ancestors = _limited_distances(patient_hpo, data["parents"], MAX_DISTANCE)
    candidate_ancestors = _limited_distances(candidate_hpo, data["parents"], MAX_DISTANCE)

    # A relationship path goes up from both terms until they meet at a shared ancestor.
    candidates = []
    shared_ancestors = patient_ancestors.keys() & candidate_ancestors.keys()
    for shared_ancestor in shared_ancestors:
        patient_distance = patient_ancestors[shared_ancestor]
        candidate_distance = candidate_ancestors[shared_ancestor]
        total = patient_distance + candidate_distance

        # Exclude exact identity and paths beyond the permitted matching distance.
        if not 0 < total <= MAX_DISTANCE:
            continue
        candidates.append((total, shared_ancestor, patient_distance, candidate_distance))
    if not candidates:
        return None

    # Retain only paths with the fewest total ontology edges.
    shortest_distance = MAX_DISTANCE + 1
    for path in candidates:
        shortest_distance = min(shortest_distance, path[0])

    shortest_paths = []
    for path in candidates:
        if path[0] == shortest_distance:
            shortest_paths.append(path)
    if len(shortest_paths) == 1:
        return shortest_paths[0]

    # On ties, prefer direct lineage and then the more specific shared ancestor.
    best_path = shortest_paths[0]
    for path in shortest_paths[1:]:
        # False sorts before True to favor direct lineage; negation favors higher specificity.
        path_order = (
            path[2] > 0 and path[3] > 0,
            -_shared_specificity(data, path[1]),
            path[1],
        )
        best_order = (
            best_path[2] > 0 and best_path[3] > 0,
            -_shared_specificity(data, best_path[1]),
            best_path[1],
        )
        # The shared-ancestor HPO ID is the final deterministic tie-breaker.
        if path_order < best_order:
            best_path = path
    return best_path


def score_term_pair(data: dict, patient_hpo: str, candidate_hpo: str) -> float:
    """Score a normalized patient and candidate HPO term; zero rejects."""
    # Classify the closest relationship and reject terms without an allowed path.
    path = _relationship_path(data, patient_hpo, candidate_hpo)
    if path is None:
        return 0.0
    distance, shared_ancestor, patient_distance, candidate_distance = path

    # Direct ancestor and descendant relationships receive distance-based fixed scores.
    ancestor_score = round(1.0 - SCORE_STEP * distance, 6)
    if patient_distance > 0 and candidate_distance == 0:  # Candidate is broader.
        return ancestor_score
    if patient_distance == 0 and candidate_distance > 0:  # Candidate is narrower.
        # Narrower candidates receive an additional penalty because they assume more detail (children).
        return round(ancestor_score - SCORE_STEP, 6)

    # Cross-branch terms require a sufficiently specific shared ancestor.
    shared_specificity = _shared_specificity(data, shared_ancestor)
    if shared_specificity < SIBLING_SHARED_SPECIFICITY_MIN:
        return 0.0

    # Penalize longer cross-branch paths before applying the sibling evidence gates.
    distance_factor = 1.0 - SCORE_STEP * (distance - 1)

    # Same-axis terms may differ only by onset, severity, laterality, or direction.
    if _same_clinical_axis(data, patient_hpo, candidate_hpo):
        return round(distance_factor * shared_specificity, 6)

    # Differently worded cross-branch relationships are limited to direct siblings.
    if distance != 2:  # Other cross-branch terms must be direct siblings.
        return 0.0

    # Ordinary siblings are compared using both their wording and annotation profiles.
    clinical = _clinical_similarity(data, patient_hpo, candidate_hpo)
    semantic = _semantic_similarity(data, patient_hpo, candidate_hpo, shared_ancestor)

    # Very strong annotation agreement can rescue clinically different wording.
    profile_supported = (
        semantic >= PROFILE_RESCUE_SEMANTIC_MIN
        and shared_specificity >= PROFILE_RESCUE_SPECIFICITY_MIN
    )
    if not profile_supported:
        # Otherwise require wording, annotation, and non-generic clinical-word support.
        wording_supported = (
            clinical >= WORDING_SUPPORT_CLINICAL_MIN
            and semantic >= WORDING_SUPPORT_SEMANTIC_MIN
            and _clinical_similarity(
                data, patient_hpo, candidate_hpo, drop=_FINDING_WORDS
            ) > 0.0
        )
        if not wording_supported:
            return 0.0

    # Blend distance, ancestor specificity, clinical wording, and annotation similarity.
    raw_score = distance_factor * shared_specificity * math.sqrt(
        (clinical + semantic) * (1.0 + clinical) / 4.0
    )

    # A zero return means the relationship should not be expanded into the report.
    return round(raw_score, 6) if raw_score >= MIN_SIBLING_EVIDENCE_SCORE else 0.0


def _nearby_terms(data: dict, patient_hpo: str) -> set[str]:
    """Return terms within the allowed up-then-down graph distance."""
    # Include the source term initially; the matcher later removes this exact match.
    nearby = {patient_hpo}

    # Move up to each reachable ancestor, then use the remaining distance to move down.
    for shared, up in _limited_distances(patient_hpo, data["parents"], MAX_DISTANCE).items():
        remaining = MAX_DISTANCE - up
        nearby.update(_limited_distances(shared, data["children"], remaining))
    return nearby


def load_references(ontology_path: str | Path, gene_annotations: pd.DataFrame) -> dict:
    """Load the ontology and annotation indexes."""
    # HPO Toolkit parses hp.json into term objects and a graph of ontology relationships.
    ontology = hpotk.load_ontology(Path(ontology_path))

    # Keep only the ontology fields needed by matching and annotation propagation.
    terms = {}
    alternative_ids = {}
    parents: dict = {}
    children: dict = defaultdict(set)

    # ontology.terms iterates over every current HPO term loaded from the JSON file.
    for term in ontology.terms:
        # HPO Toolkit identifiers are objects; .value provides the HP:####### string.
        term_id = term.identifier.value
        synonyms = []
        if term.synonyms is not None:
            for synonym in term.synonyms:
                synonyms.append(synonym.name)
        terms[term_id] = {"name": term.name, "synonyms": tuple(synonyms)}

        # Resolve legacy alternative IDs to the current primary HPO ID.
        for alternative_id in term.alt_term_ids:
            alternative_ids[alternative_id.value] = term_id

        # get_parents returns direct parents only; reverse these edges to derive children.
        term_parents = set()
        for parent in ontology.graph.get_parents(term.identifier):
            term_parents.add(parent.value)
            children[parent.value].add(term_id)
        parents[term_id] = term_parents

    # Store the graph, term metadata, and caches together for the scoring functions.
    data: dict = {
        "terms": terms,
        "alternative_ids": alternative_ids,
        "parents": parents,
        "children": children,
        "ancestor_cache": {},
        "specificity_cache": {},
    }

    # Normalize the three annotation fields used to build gene and disease indexes.
    annotations = gene_annotations[["gene_symbol", "disease_id", "hpo_id"]]
    annotations = annotations.fillna("").astype(str)
    direct_genes: dict = defaultdict(set)
    gene_terms: dict = defaultdict(set)
    disease_terms: dict = defaultdict(set)

    # Record direct term-to-gene links and reverse gene/disease-to-term links.
    for row in annotations.itertuples(index=False):
        gene = row.gene_symbol.strip()
        disease = row.disease_id.strip() or "-"
        term_id = row.hpo_id.strip()

        # Convert an older HPO ID to its current ID before checking that the term exists.
        term_id = alternative_ids.get(term_id, term_id)
        if not gene or gene == "-" or term_id not in terms:
            continue
        direct_genes[term_id].add(gene)
        gene_terms[gene].add(term_id)
        disease_terms[disease].add(term_id)

    # Propagated indexes associate an annotation with its term and every broader ancestor.
    propagated_genes: dict = defaultdict(set)
    propagated_diseases: dict = defaultdict(set)
    for term_annotations, propagated in ((gene_terms, propagated_genes), (disease_terms, propagated_diseases),):
        for item, term_ids in term_annotations.items():
            for term_id in term_ids:
                for ancestor in _ancestors(data, term_id):
                    propagated[ancestor].add(item)

    # These indexes and totals support broadness checks and information-content scoring.
    data.update(
        {
            "direct_genes": direct_genes,
            "propagated_genes": propagated_genes,
            "propagated_diseases": propagated_diseases,
            "total_gene_count": max(1, len(gene_terms)),
            "total_disease_count": max(1, len(disease_terms)),
        }
    )

    return data


def get_indirect_hpo_gene_matches(patient_hpo_ids: Iterable[str], data: dict) -> pd.DataFrame:
    """Return one row per indirect gene and candidate-HPO pair."""
    # Keep the best score if multiple patient terms produce the same gene/HPO pair.
    pair_scores: dict[tuple[str, str], float] = {}

    for hpo_id in patient_hpo_ids:
        # Resolve older IDs and ignore unknown terms or terms outside Phenotypic abnormality.
        patient_hpo = data["alternative_ids"].get(hpo_id, hpo_id)
        if (patient_hpo not in data["terms"] or PHENOTYPIC_ABNORMALITY not in _ancestors(data, patient_hpo)):
            continue

        # Score eligible nearby terms before deciding whether the whole expansion is too broad.
        candidate_matches = []
        combined_genes = set()
        for candidate_hpo in _nearby_terms(data, patient_hpo):
            candidate_genes = data["direct_genes"].get(candidate_hpo)

            # Exclude exact matches, unannotated terms, and uninformative gene expansions.
            if (
                candidate_hpo == patient_hpo
                or not candidate_genes
                or len(candidate_genes) > MAX_GENES_PER_CANDIDATE_TERM
                or not _worth_expanding(data, candidate_hpo)
            ):
                continue

            # Identify high-level candidate terms whose branches cover over 2,000 genes.
            candidate_root_distance = _limited_distances(
                candidate_hpo, data["parents"], BROAD_SOURCE_MAX_ROOT_DISTANCE
            ).get(PHENOTYPIC_ABNORMALITY)
            broad_candidate = (
                candidate_root_distance is not None
                and len(data["propagated_genes"].get(candidate_hpo, set()))
                > BROAD_SOURCE_MAX_INDIRECT_GENES
            )

            # A positive score means the relationship passed its applicable matching gates.
            score = score_term_pair(data, patient_hpo, candidate_hpo)
            if score <= 0.0:
                continue

            # Include every accepted term when measuring source breadth, but do not report
            # a high-level candidate itself. This preserves the original source safeguard.
            combined_genes.update(candidate_genes)
            if not broad_candidate:
                candidate_matches.append((candidate_hpo, candidate_genes, score))

        # Patient terms close to the root are suppressed when their combined expansion is huge.
        root_distance = _limited_distances(
            patient_hpo, data["parents"], BROAD_SOURCE_MAX_ROOT_DISTANCE
        ).get(PHENOTYPIC_ABNORMALITY)
        if (
            root_distance is not None
            and len(combined_genes) > BROAD_SOURCE_MAX_INDIRECT_GENES
        ):
            continue

        # Expand each accepted candidate term into its directly annotated gene/HPO pairs.
        for candidate_hpo, candidate_genes, score in candidate_matches:
            for gene in candidate_genes:
                pair = (gene, candidate_hpo)
                pair_scores[pair] = max(pair_scores.get(pair, 0.0), score)

    # Format indirect matches with the same core columns used by the exact-match table.
    rows = []
    for pair, score in sorted(pair_scores.items()):
        gene, term_id = pair
        rows.append({
            "gene_symbol": gene,
            "hpo_id": term_id,
            "hpo_name": data["terms"][term_id]["name"],
            "HPO Match Score": round(score, 3),
        })
    return pd.DataFrame(rows, columns=["gene_symbol", "hpo_id", "hpo_name", "HPO Match Score"])
