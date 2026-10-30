#!/usr/bin/env bash
# scripts/run_all_analyses.sh
#
# Unified wrapper: runs BLAST + NextClade analyses on a batch and builds per-sequence phylogenetic trees.
#
# Called by hav_wrapper.sh with:
#   bash scripts/run_all_analyses.sh <batch_dir> <dataset_date> <batch_fa> <dataset_dir> <out_base> <year> [options]
#
# Arguments:
#   <batch_dir>     Batch directory path
#   <dataset_date>  Dataset version date (e.g., 2026-04-10)
#   <batch_fa>      Path to batch FASTA file
#   <dataset_dir>   Path to local_datasets/<dataset_date> directory
#   <out_base>      Output base directory for results
#   <year>          Year folder (for context)
#   [options]       Optional Nextclade settings (see below)
#
# NextClade alignment options:
#   --min-match-length <int>       Minimum seed match length for NextClade (default: 20)
#
# Analyses run (in order):
#   1. BLAST     — batch sequences vs the local HAV BLAST database
#                  (closest-match / genotype confirmation)
#   2. NextClade — against the public-derived hav-vp1-2b-lineages dataset
#                  (assigns genotype `clade` and QC metrics)
#   3. Per-sequence trees — builds IQ-TREE phylogenetic trees for each query
#                  (using nearest neighbors from BLAST and NextClade)
#   4. Batch analysis — batch phylogenetic trees and SNP matrices
#   5. Variant names + outbreak trees — propose_variant_names.R proposes variant names,
#                  then trees combine batch + database sequences
#   6. Batch report — generates interactive HTML report with visualizations
#
# Writes to <out_base>/:
#   blast_results.tsv
#   lineages/nextclade.tsv  (+ aligned fasta, json, ndjson)
#   lineages/nextclade.auspice.json
#   trees/<seqName>/*  (per-sequence tree artifacts)
#   trees/batch/*  (batch tree artifacts)
#   snp_matrices/batch_snp_distances.tsv
#   proposed_variants.tsv  (variant naming table: genotype, BLAST hit + sample date, SNPs, proposed variant, notes)
#   outbreak_trees/  (PNG images of outbreak lineage trees)
#   summary.tsv  (per-sequence summary, written by the report)
#   batch_report.html  (interactive HTML report)
#
# Prerequisites:
#   - BLAST database built:   bash scripts/make_blast_db.sh [dataset_date]
#   - Nextclade dataset staged by hav_wrapper.sh
#   - Conda environments: BLAST, NEXTCLADE, iqtree-mafft, SEQKIT, R_shared

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_all_analyses.sh <batch_dir> <dataset_date> <batch_fa> <dataset_dir> <out_base> <year> [options]

Required arguments:
  <batch_dir>      Batch directory path
  <dataset_date>   Dataset date (e.g., 2026-04-10)
  <batch_fa>       Path to batch FASTA file
  <dataset_dir>    Path to local_dataset/<date>
  <out_base>       Output base directory
  <year>           Year folder

Optional flags:
  --min-match-length <int>       Minimum seed match length for Nextclade (default: 20)
USAGE
}


if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

BATCH_DIR="$1"
DATASET_DATE="$2"
BATCH_FA="$3"
DATASET_DIR="$4"
OUT_BASE="$5"
YEAR="$6"

# Skip the 6 positional arguments before processing flags
shift 6

MIN_MATCH_LENGTH=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --min-match-length)
      MIN_MATCH_LENGTH="${2:?--min-match-length requires a value}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

THREADS=4


# ── Validation ────────────────────────────────────────────────────────────────
if [[ ! -f "$BATCH_FA" ]]; then
  echo "ERROR: Batch FASTA not found: $BATCH_FA" >&2
  exit 1
fi
if [[ ! -d "$DATASET_DIR" ]]; then
  echo "ERROR: Local dataset directory not found: $DATASET_DIR" >&2
  exit 1
fi
if [[ ! -d "$LINEAGES_DATASET" ]]; then
  echo "ERROR: Lineages dataset not found at $LINEAGES_DATASET" >&2
  exit 1
fi

mkdir -p "$OUT_BASE"

echo "════════════════════════════════════════════════════════════════"
echo "HAV Batch Analysis: $BATCH_NAME  (dataset: $DATASET_DATE)"
echo "  Batch     : $BATCH_NAME"
echo "  Year      : $YEAR"
echo "  Batch Dir       : $BATCH_DIR"
echo "  Dataset   : $DATASET_DATE"
echo "  Threads   : $THREADS"
echo "  Output    : $OUT_BASE"
echo "  Lineages dataset : $LINEAGES_DATASET"
echo "════════════════════════════════════════════════════════════════"
echo ""


# Primer trimming is intentionally disabled: the replacement
# trim_primers_by_reference.py script has not been implemented yet.
ANALYSIS_FA="$BATCH_FA"

# ── Analysis 1: BLAST ─────────────────────────────────────────────────────────
echo "▶ Step 1/5: BLAST"
echo "─────────────────────────────────────────────────────────────────"

# blast_batch.sh writes directly to OUT_BASE
# Pass arguments directly to preserve spaces in paths

conda run -n BLAST bash "$HAV_SEQ_REPO/scripts/blast_batch.sh" "$ANALYSIS_FA" "$DATASET_DATE" "$DATASET_DIR" "$OUT_BASE"
echo ""

# ── Analysis 2: NextClade (lineages dataset) ──────────────────────────────────
echo "▶ Step 2/6: NextClade — clade assignment and QC"
echo "─────────────────────────────────────────────────────────────────"
echo "  Lineages dataset → $LINEAGES_DATASET"
LINEAGES_OUT="$OUT_BASE/lineages"
mkdir -p "$LINEAGES_OUT"

if [[ ! -d "$LINEAGES_DATASET" ]]; then
  echo "ERROR: Lineages dataset not found at $LINEAGES_DATASET" >&2
  exit 1
fi
# Pass arguments directly to preserve spaces in paths
conda run -n NEXTCLADE nextclade run \
  --input-dataset    "$LINEAGES_DATASET" \
  --alignment-preset high-diversity \
  --output-all       "$LINEAGES_OUT" \
  --jobs             "$THREADS" \
  --min-match-length "$MIN_MATCH_LENGTH" \
  "$ANALYSIS_FA"

echo "  Results → $LINEAGES_OUT/nextclade.tsv"
N_SEQS=$(tail -n +2 "$LINEAGES_OUT/nextclade.tsv" | wc -l)
echo "  Sequences processed: $N_SEQS"
echo ""


# ── Analysis 3: Per-sequence trees ────────────────────────────────────────────
echo "▶ Step 3/6: Per-sequence trees (IQ-TREE with nearest neighbors)"
echo "─────────────────────────────────────────────────────────────────"

# Pass arguments directly to preserve spaces in paths
conda run -n iqtree-mafft bash "$HAV_SEQ_REPO/scripts/build_per_seq_trees.sh" "$BATCH_DIR" "$DATASET_DATE" 30 "$OUT_BASE" "$DATASET_DIR" "$BATCH_FA"
echo "  Results → $OUT_BASE/trees/"
echo ""


# ── Analysis 4: Batch-specific phylogenetics ─────────────────────────────────
echo "▶ Step 4/6: Batch-specific trees and SNP matrices"
echo "─────────────────────────────────────────────────────────────────"
# Pass arguments directly to preserve spaces in paths
conda run -n iqtree-mafft bash "$HAV_SEQ_REPO/scripts/build_batch_analysis.sh" "$BATCH_DIR" "$DATASET_DATE" "$OUT_BASE"
echo "  Results → $OUT_BASE/trees/batch/"
echo "  SNP matrix → $OUT_BASE/snp_matrices/batch_snp_distances.tsv"
echo ""


# ── Analysis 5: Outbreak-specific phylogenetics ──────────────────────────────
echo "▶ Step 5/6: Outbreak lineage trees (batch + database sequences)"
echo "─────────────────────────────────────────────────────────────────"

MAFFT_BIN="$(
  conda run -n iqtree-mafft which mafft 2>/dev/null |
  tr -d '\r' |
  grep -v '^[[:space:]]*$' |
  tail -n 1
)"

IQTREE_BIN="$(
  conda run -n iqtree-mafft which iqtree 2>/dev/null |
  tr -d '\r' |
  grep -v '^[[:space:]]*$' |
  tail -n 1
)"

printf 'MAFFT=[%s]\n' "$MAFFT_BIN"
printf 'IQTREE=[%s]\n' "$IQTREE_BIN"

[[ -x "$MAFFT_BIN" ]] && echo "MAFFT OK" || echo "MAFFT NOT EXECUTABLE"
[[ -x "$IQTREE_BIN" ]] && echo "IQ-TREE OK" || echo "IQ-TREE NOT EXECUTABLE"


# Proposed variant names (shared by outbreak trees and the report)
conda run -n R_shared Rscript \
  "$HAV_SEQ_REPO/scripts/propose_variant_names.R" \
  "$OUT_BASE" \
  "$DATASET_DATE" \
  "$DATASET_DIR" \
  "$BATCH_FA" || {
    echo "WARNING: Variant name proposal failed"
    echo "         Trees and report will lack proposed variant names"
  }

conda run -n R_shared Rscript \
  "$HAV_SEQ_REPO/scripts/build_outbreak_trees.R" \
  "$OUT_BASE" \
  "$DATASET_DATE" \
  "$DATASET_DIR" \
  "$MAFFT_BIN" \
  "$IQTREE_BIN" || {
    echo "WARNING: Outbreak tree generation failed"
    echo "         Report will not include outbreak lineage trees"
  }

echo "  Results → $OUT_BASE/outbreak_trees/"
echo ""

# ── Analysis 6: Batch report ──────────────────────────────────────────────────
echo "▶ Step 6/6: Generate batch report (HTML)"
echo "─────────────────────────────────────────────────────────────────"

REPORT_OUTPUT="$OUT_BASE/batch_report.html"
# Pass arguments directly to preserve spaces in paths
conda run -n R_shared Rscript -e "
rmarkdown::render(
  '$HAV_SEQ_REPO/scripts/batch_report.Rmd',
  output_file = '$REPORT_OUTPUT',
  params = list(
    batch_dir    = '$BATCH_DIR',
    batch_name   = '$BATCH_NAME',
    batch_fa     = '$BATCH_FA',    
    dataset_date = '$DATASET_DATE',
    n_neighbors  = 30,
    out_base     = '$OUT_BASE',
    dataset_dir  = '$DATASET_DIR',
    final_output_dir = '$SMB_DIR'
  ),
  quiet = TRUE
)
" || { echo "ERROR: Batch report generation failed" >&2; exit 1; }

echo "  Report written: $REPORT_OUTPUT"
echo ""


# ── Summary ───────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════"
echo "All analyses complete for $BATCH_NAME"
echo ""
echo "  BLAST results            : $OUT_BASE/blast_results.tsv"
echo "  NC lineages (FHI)        : $LINEAGES_OUT/nextclade.tsv"
echo "  Per-sequence trees       : $OUT_BASE/trees/"
echo "  Batch phylogenetic tree  : $OUT_BASE/trees/batch/tree.treefile"
echo "  Batch SNP matrix         : $OUT_BASE/snp_matrices/batch_snp_distances.tsv"
echo "  Outbreak lineage trees   : $OUT_BASE/outbreak_trees/"
echo "  Batch report             : $REPORT_OUTPUT"
echo "════════════════════════════════════════════════════════════════"
