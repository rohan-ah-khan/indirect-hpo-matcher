#!/usr/bin/env python3
"""Create a gene-level table of exact and indirect HPO matches."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from hpo_indirect_matcher import get_indirect_hpo_gene_matches, load_references


SCRIPT_DIR = Path(__file__).resolve().parent
HPO_ID_PATTERN = re.compile(r"HP:\d{7}", re.IGNORECASE)


def read_patient_hpo_ids(hpo_file: Path) -> list[str]:
    """Extract unique HPO IDs from a patient_hpo_ids line or plain text."""
    try:
        text = hpo_file.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError:
        text = hpo_file.read_text(encoding="ISO-8859-1")

    patient_lines = [
        line for line in text.splitlines() if "patient_hpo_ids" in line.lower()
    ]
    hpo_text = "\n".join(patient_lines) if patient_lines else text
    hpo_ids = [match.upper() for match in HPO_ID_PATTERN.findall(hpo_text)]
    hpo_ids = list(dict.fromkeys(hpo_ids))
    if not hpo_ids:
        raise ValueError(f"No HPO IDs were found in {hpo_file}")
    return hpo_ids


def format_score(score: float) -> str:
    """Format a score without unnecessary trailing zeroes."""
    value = f"{float(score):.3f}".rstrip("0").rstrip(".")
    return value if "." in value else f"{value}.0"


def join_unique(values) -> str:
    """Join values in their existing order without duplicates."""
    return "; ".join(dict.fromkeys(str(value) for value in values))


def make_matches(
    patient_hpo_ids: list[str],
    hpo_data: dict,
    gene_annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Combine direct matches scored at 1.0 with indirect matches."""
    patient_hpo_ids = [
        hpo_data["alternative_ids"].get(hpo_id, hpo_id)
        for hpo_id in patient_hpo_ids
    ]
    patient_hpo_ids = list(dict.fromkeys(patient_hpo_ids))

    unknown = [hpo_id for hpo_id in patient_hpo_ids if hpo_id not in hpo_data["terms"]]
    if unknown:
        print(
            f"Warning: ignoring unknown HPO IDs: {', '.join(unknown)}",
            file=sys.stderr,
        )
    patient_hpo_ids = [
        hpo_id for hpo_id in patient_hpo_ids if hpo_id in hpo_data["terms"]
    ]
    if not patient_hpo_ids:
        raise ValueError("None of the input HPO IDs occur in the supplied ontology")

    direct = gene_annotations[gene_annotations["hpo_id"].isin(patient_hpo_ids)][
        ["gene_symbol", "hpo_id", "hpo_name"]
    ].copy()
    direct["gene_symbol"] = direct["gene_symbol"].fillna("").str.strip()
    direct = direct[~direct["gene_symbol"].isin(["", "-"])]
    direct["HPO Match Score"] = 1.0

    indirect = get_indirect_hpo_gene_matches(patient_hpo_ids, hpo_data)
    matches = pd.concat([direct, indirect], ignore_index=True)
    matches["HPO Match Score"] = pd.to_numeric(
        matches["HPO Match Score"], errors="coerce"
    )
    matches = matches.sort_values("HPO Match Score", ascending=False)
    return matches.drop_duplicates(["gene_symbol", "hpo_id"])


def group_matches(matches: pd.DataFrame, gene_map_file: Path) -> pd.DataFrame:
    """Create the one-row-per-gene format used by the report pipelines."""
    gene_map = pd.read_csv(gene_map_file, dtype=str)
    required = {"hgnc_symbol", "ensembl_gene_id"}
    missing = required - set(gene_map.columns)
    if missing:
        raise ValueError(
            f"{gene_map_file} is missing required columns: {', '.join(sorted(missing))}"
        )

    gene_ids = dict(zip(gene_map["hgnc_symbol"], gene_map["ensembl_gene_id"]))
    matches = matches.copy()
    matches["Gene ID"] = matches["gene_symbol"].map(gene_ids)
    matches["Features"] = matches.apply(
        lambda row: f'{row["hpo_name"]} ({format_score(row["HPO Match Score"])})',
        axis=1,
    )
    matches["HPO Match Score"] = matches["HPO Match Score"].map(format_score)

    grouped = matches.groupby(
        ["gene_symbol", "Gene ID"], as_index=False, dropna=False
    ).agg(
        **{
            "Number of occurrences": ("hpo_id", "nunique"),
            "Features": ("Features", join_unique),
            "HPO IDs": ("hpo_id", join_unique),
            "HPO Match Score": ("HPO Match Score", lambda scores: "; ".join(scores)),
        }
    )
    return grouped.rename(columns={"gene_symbol": "Gene Symbol"})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write exact and indirect HPO gene matches for a list of patient HPO IDs."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Text file containing patient HPO IDs, preferably a patient_hpo_ids line",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output matched TSV",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=SCRIPT_DIR / "references",
        help="Directory containing hp.json and genes_to_phenotype.txt",
    )
    parser.add_argument(
        "--gene-map",
        type=Path,
        default=SCRIPT_DIR / "ensembl_to_NCBI_ID.csv",
        help="CSV containing hgnc_symbol and ensembl_gene_id",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotations_file = args.reference_dir / "genes_to_phenotype.txt"
    ontology_file = args.reference_dir / "hp.json"

    patient_hpo_ids = read_patient_hpo_ids(args.input)
    gene_annotations = pd.read_csv(annotations_file, sep="\t", dtype=str)
    hpo_data = load_references(ontology_file, gene_annotations)
    matches = make_matches(patient_hpo_ids, hpo_data, gene_annotations)
    grouped = group_matches(matches, args.gene_map)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    grouped.to_csv(args.output, sep="\t", index=False)
    print(
        f"Wrote {len(grouped)} genes from {len(patient_hpo_ids)} patient HPO IDs "
        f"to {args.output}"
    )


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, KeyError, OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error
