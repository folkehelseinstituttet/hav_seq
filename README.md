# HAV sequence analysis

This repository contains scripts for Sanger HAV sequence analysis. The current
supported pipeline runs BLAST and Nextclade classification, builds phylogenetic
trees, and generates an HTML batch report.

**WGS analysis is not currently supported.** The project wrapper rejects
`--mode wgs`; WGS is a planned future addition, and its assembly step has not
yet been implemented.

## How to run

For routine institute use, start from `wrapper.sh` in the separate
[`ngs_scripts/hav_sequencing`](https://github.com/folkehelseinstituttet/ngs_scripts/tree/main/hav_sequencing)
directory. That wrapper stages inputs and reference data from shared storage,
updates or clones `hav_seq`, calls `scripts/hav_wrapper.sh`, then copies
analysis results and the updated sequence database back to shared storage.

Its current argument format is:

```bash
bash ~/ngs_scripts/hav_sequencing/wrapper.sh --mode sanger <batch_name> <year>
```

The batch name and year are positional arguments. The outer wrapper provides
runtime variables such as `TMP_DIR` and `HAV_SEQ_REPO` and stages the input
files expected by `hav_wrapper.sh`.

The project wrapper can also be run in an already prepared Linux environment:

```bash
bash scripts/hav_wrapper.sh --mode sanger <batch_name> <year>
```

It expects staged inputs in `$TMP_DIR`, including the Sanger FASTA files,
local HAV sequence FASTA (`local_dataset/2PA.fa`), and metadata files used by
the pipeline. A samplesheet can be supplied with `--samplesheet`; otherwise,
one is generated from FASTA filenames. The samplesheet generator expects names
like `LabWareID_SekvensID-HAV.fa` or `.fasta`.

The Sanger pipeline prepares a batch FASTA, builds the local BLAST database,
prepares metadata, runs BLAST and Nextclade, builds per-sequence, batch and
outbreak phylogenies, and renders `batch_report.html`. Analysis outputs are
written to `$HOME/<batch_name>_results`; the outer wrapper copies them to
shared storage.

## Conda environments

The Conda environment files currently present are:

| Environment | Definition | Use |
|---|---|---|
| `BLAST` | `envs/BLAST.yml` | Local BLAST database creation and batch BLAST. |
| `NEXTCLADE` | `envs/NEXTCLADE.yml` | Downloading Nextclade datasets and running classification. |
| `iqtree-mafft` | `envs/iqtree-mafft.yml` | MAFFT alignments and IQ-TREE phylogenetic analysis. |
| `SEQKIT` | `envs/SEQKIT.yml` | Selecting reference sequences during tree construction. |
| `R_shared` | `envs/R_shared.yml` | R scripts, SNP/tree analysis, and HTML report rendering. |

There is no `HAVDEV` environment in the current setup. Create the named
environments from the project root, for example:

```bash
conda env create -f envs/BLAST.yml
conda env create -f envs/NEXTCLADE.yml
conda env create -f envs/iqtree-mafft.yml
conda env create -f envs/SEQKIT.yml
conda env create -f envs/R_shared.yml
```

The analysis scripts invoke these environments with `conda run -n <name>`.
The environment files are Linux-oriented; use the institute's supported Linux
or WSL environment.

## Current limitations and future work

- **WGS:** `--mode wgs` is rejected by `hav_wrapper.sh`. The old assembly
  script, `scripts/run_pipeline.sh`, is not present and has not been
  implemented. WGS support will require that assembly stage, consensus
  collection, and an end-to-end validation of downstream classification.
- **Primer trimming:** optional primer trimming has been removed from the
  active pipeline. The former `scripts/trim_primers_by_reference.py` helper is
  not present. `run_all_analyses.sh` documents this as not yet implemented and
  runs analyses on the untrimmed input FASTA. A replacement can be added later
  if primer trimming is needed.

These steps are intentionally disabled rather than allowed to fail later due
to missing scripts.

## Script reference

| File | Purpose |
|---|---|
| `scripts/hav_wrapper.sh` | Sanger entry point; prepares the batch FASTA and orchestrates analysis. |
| `scripts/run_all_analyses.sh` | Runs BLAST, Nextclade, phylogenetic analysis, and report generation. |
| `scripts/generate_samplesheet.R` | Creates a samplesheet from Sanger FASTA filenames. |
| `scripts/prepare_input_fasta.R` | Collects sample sequences into a batch multi-FASTA. |
| `scripts/prepare_metadata.R` | Combines request and HAV metadata into `metadata.tsv` (id, genotype, variant, sample_date). `sample_date` comes from `SAMPLED_DATE` in `HAV_lw_uttrekk.tsv`; for `Requests.xlsx` only `Year` is available and is used instead. |
| `scripts/make_blast_db.sh` | Prepares and builds the local nucleotide BLAST database. |
| `scripts/blast_batch.sh` | Searches batch sequences against the local BLAST database. |
| `scripts/build_per_seq_trees.sh` | Selects nearby sequences, aligns them, and builds per-sequence IQ-TREE trees. |
| `scripts/build_batch_analysis.sh` | Builds batch phylogenies and SNP-distance results. |
| `scripts/propose_variant_names.R` | Proposes variant names for batch sequences; writes `proposed_variants.tsv` (genotype, BLAST hit and its sample date, SNPs, proposed variant, notes). |
| `scripts/build_outbreak_trees.R` | Generates outbreak/variant trees using batch and database sequences. |
| `scripts/batch_report.Rmd` | Renders the HTML report (including sample dates for BLAST hits and variants) and `summary.tsv`. |
| `scripts/R-database.R` | Utility for importing, querying, and exporting external sequences in SQLite. |

## Runtime data

Batch inputs, sequence database exports, metadata, and generated results are
runtime data and are not included in this repository. The institute wrapper
stages the operational inputs from shared storage. Keep confidential sequence
data, metadata, credentials, and batch outputs out of Git.
