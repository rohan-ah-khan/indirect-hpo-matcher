# Indirect HPO matcher

This repository provides a standalone version of the HPO matching used by the
short- and long-read annotation pipelines. It takes a list of proband HPO IDs and
creates a gene-level TSV containing both exact and indirectly related phenotype
annotations. Snakemake and a pipeline-generated HPO file are not required.

Exact gene–phenotype annotations receive a score of `1.0`. The matching
logic in `hpo_indirect_matcher.py` adds nearby ontology terms that pass its
relationship, wording, annotation-profile, and broadness checks. Indirect scores
are rule-based indicators of matching support rather than probabilities.

## Repository contents

- `run_hpo_matcher.py`: command-line wrapper.
- `hpo_indirect_matcher.py`: ontology matching and scoring logic.
- `download_hpo_references.sh`: downloads the two required HPO reference files.
- `ensembl_to_NCBI_ID.csv`: maps HGNC symbols to Ensembl gene IDs.
- `environment.yml`: Python environment definition.
- `examples/`: a five-term input and its generated output.

## Create the environment in this repository

From the repository root, you can create the Conda environment under `./env`:

```bash
conda env create --prefix "$PWD/env" --file environment.yml
conda activate "$PWD/env"
```

`hpo-toolkit` is available through PyPI, so `environment.yml` installs it using `pip` inside this Conda environment.

## Download the HPO references

```bash
bash download_hpo_references.sh
```

This creates `references/` in the repository and downloads:

- `references/hp.json`, containing the HPO terms and ontology relationships.
- `references/genes_to_phenotype.txt`, containing gene–phenotype annotations.

The matching thresholds were calibrated with HPO release `2026-06-23` and checked
against `2026-09-01`. A substantially different release should be validated before
routine use.

## Input format

Create a text file containing the proband HPO IDs. The preferred format is as follows:

```text
# proband_hpo_ids=HP:0007328,HP:0000970,HP:0001250
```

The leading `#` is optional. If a `proband_hpo_ids` line is present, IDs elsewhere
in the file are ignored. Otherwise, every unique `HP:#######` ID in the file is
used. Alternative HPO IDs are converted to their current primary IDs, while IDs
absent from the supplied ontology are reported and ignored.

## Run the matcher

With the default reference and gene-map locations:

```bash
python run_hpo_matcher.py \
  --input proband_hpo_ids.txt \
  --output proband.matched.tsv
```

Alternate locations can be specified when needed:

```bash
python run_hpo_matcher.py \
  --input proband_hpo_ids.txt \
  --output proband.matched.tsv \
  --reference-dir /path/to/references \
  --gene-map /path/to/ensembl_to_NCBI_ID.csv
```

The reference directory must contain `hp.json` and `genes_to_phenotype.txt`. The
gene map must contain `hgnc_symbol` and `ensembl_gene_id` columns.

## Output format

The output contains one row per gene:

| Column | Description |
| --- | --- |
| `Gene Symbol` | HGNC gene symbol from the HPO annotation file. |
| `Gene ID` | Ensembl gene ID from `ensembl_to_NCBI_ID.csv`. |
| `Number of occurrences` | Number of unique matched HPO terms for the gene. |
| `Features` | Matched term names with their scores. |
| `HPO IDs` | Matched HPO IDs in the same order as `Features`. |
| `HPO Match Score` | Scores in the same order as the matched terms. |

If the same gene–HPO pair is reached from more than one proband term, only its
highest score is retained. Direct and indirect matches are then combined before
the table is grouped by gene.

## Included example

`examples/example_proband_hpo_ids.txt` contains five example proband terms:

- Facial flushing after alcohol intake (`HP:0001033`)
- Trimethylaminuria (`HP:0003614`)
- Fractured radius (`HP:0003978`)
- Mitral annular calcification (`HP:0005136`)
- Susceptibility to chickenpox (`HP:0005360`)

The checked-in output was generated with HPO release `2026-09-01`. After
downloading the references, run the example with:

```bash
python run_hpo_matcher.py \
  --input examples/example_proband_hpo_ids.txt \
  --output examples/example.matched.tsv
```

Results may change slightly when the ontology or gene annotations are updated.
