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
cluster phylogenies, and renders `batch_report.html`. Analysis outputs are
written to `$HOME/<batch_name>_results`; the outer wrapper copies them to
shared storage.

## Analysis choices and interpretation

The values below describe the current implementation. Some are set directly in
the scripts rather than exposed as wrapper options, so changing them may
require updating the relevant script.

### Per-sequence trees

- The pipeline requests up to **30 reference neighbors per query**, in addition
  to the query sequence itself. A per-sequence tree therefore has up to 31
  sequences; fewer may be available.
- Neighbor selection starts with BLAST hits, ranked by lowest mismatch count
  and then highest bitscore, with at most **2 hits per lineage**. Up to two
  thirds of the neighbor slots are initially reserved for BLAST (20 of 30 at
  the current setting).
- Remaining slots are filled from the Nextclade community tree: leaves below
  the query's nearest Nextclade node are ranked by accumulated branch
  mutations. If there are too few descendants, sibling subtrees are considered.
  If no community-tree neighbors are available, remaining slots are filled
  with additional eligible BLAST hits.
- A tree is attempted only when at least **3 sequences including the query**
  are available. The MAFFT alignment is trimmed to positions covered by the
  query. IQ-TREE uses GTR+G; trees with at least 4 sequences use 1,000 UFBoot
  replicates, while 3-sequence trees are built without bootstrap support.

### Batch tree and SNP distances

- The batch tree and `batch_snp_distances.tsv` use the batch query sequences
  extracted from per-sequence BLAST-only alignments. Reference sequences are
  not included: the batch tree compares the batch sequences with one another.
  Sequences without a usable per-sequence alignment are omitted, and at least
  2 sequences are required.
- The batch sequences are realigned together with MAFFT before tree and SNP
  analysis. Pairwise SNP values are raw mismatch counts, not percentages.
  Columns containing gaps, `n`, or `?` in either sequence are omitted from the
  batch SNP count.
- Per-sequence SNP matrices and heatmaps use the query plus its selected BLAST
  neighbors (normally up to 20 BLAST references at the 30-neighbor setting).
  They do not use the Nextclade community-tree neighbors. Pairwise values are
  mismatch counts over the aligned positions that are not gaps in either
  sequence. For older outputs without a BLAST-only alignment, the report falls
  back to the full per-sequence alignment.

### Cluster identification trees

- The report groups a query's best-hit variant with alternative database
  variants that are within **3 BLAST mismatches** of that query. Variant names
  that share a base name (for example, `V6`, `V6a`, and `V6b`) are treated as
  related; shared variants can connect groups transitively.
- The report's variant overview includes groups with at least **2 batch and
  database sequences combined**. The tree builder first requires at least
  **3 batch and database records** for a group, and then needs at least 3
  successfully loaded sequences in total to build a tree. The tree includes the
  group's batch sequences and matching local-database sequences, plus up to
  **2 outgroup sequences** from other detected groups when available. There is
  no fixed maximum tree size.

### Proposed variant names

The nearest BLAST hit is selected by the lowest mismatch count, breaking ties
by highest bitscore. The naming rules use the **calendar year when the script
runs**, not the dataset date:

- No BLAST hit, more than **2 mismatches**, or a hit assigned to another year:
  propose the next `NOR-<current-year>-V<n>` number.
- **0 mismatches** to a variant from the current year: keep that database
  variant name.
- **1-2 mismatches** to a variant from the current year: use its base name with
  the suffix after the highest existing letter (`b` if there is no suffix;
  suffix `a` is skipped).
- Identical sequences in the batch share a proposed name. New numbers are
  assigned in the order sequences appear in the batch FASTA.

These are analysis and naming rules, not epidemiological decision thresholds.
Cluster membership and proposed names should be reviewed with the sequence
results and relevant metadata before being used as conclusions about
transmission or as final database assignments. BLAST mismatch counts also
depend on the alignment and reference database used; they are not normalized
pairwise distances.

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
| `scripts/build_outbreak_trees.R` | Generates cluster phylogenetic trees using batch and database sequences; saves them under `cluster_trees/`. |
| `scripts/batch_report.Rmd` | Renders the HTML report (including sample dates for BLAST hits and variants) and `summary.tsv`. |
| `scripts/R-database.R` | Utility for importing, querying, and exporting external sequences in SQLite. |

## Runtime data

Batch inputs, sequence database exports, metadata, and generated results are
runtime data and are not included in this repository. The institute wrapper
stages the operational inputs from shared storage. Keep confidential sequence
data, metadata, credentials, and batch outputs out of Git.
