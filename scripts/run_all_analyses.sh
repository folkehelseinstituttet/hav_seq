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
#   [options]       Optional flags (see below)
#
# Optional Sanger-specific options:
#   --sanger                       Enable primer trimming pre-step
#   --primer-names <csv>           Comma-separated primer names to use
#   --primers-file <path>          Optional custom primer file (FASTA/CSV/TSV)
#   --trimmed-fasta <path>         Output path for trimmed FASTA
#   --resolved-primers-tsv <path>  Output path for resolved primer coordinates
#
# NextClade alignment options:
#   --min-match-length <int>       Minimum seed match length for NextClade (default: 20)
#
# Analyses run (in order):
#   1. BLAST     — batch sequences vs the local HAV BLAST database
#                  (closest-match / genotype confirmation)
#   2. NextClade — against the public-derived hav-vp1-2b-lineages dataset
#                  (assigns genotype `clade` + lineage information)
#   3. Per-sequence trees — builds IQ-TREE phylogenetic trees for each query
#                  (using nearest neighbors from BLAST and NextClade)
#   4. Batch report — generates interactive HTML report with visualizations
#
# Writes to <out_base>/:
#   blast_results.tsv
#   lineages/nextclade.tsv  (+ aligned fasta, json, ndjson)
#   lineages/nextclade.auspice.json
#   trees/<seqName>/*  (per-sequence tree artifacts)
#   batch_report.html  (interactive HTML report)
#
# Prerequisites:
#   - BLAST database built:   bash scripts/make_blast_db.sh [dataset_date]
#   - Lineages dataset built: bash scripts/build_vp1b_lineage_dataset.sh
#   - HAVDEV conda environment active (nextclade, blastn, seqkit)

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_all_analyses.sh <batch_dir> <dataset_date> <batch_fa> <dataset_dir> <out_base> <year> [options]

Required arguments:
  <batch_dir>      Batch directory path
  <dataset_date>   Dataset date (e.g., 2026-04-10)
  <batch_fa>       Path to batch FASTA file
  <dataset_dir>    Path to local_datasets/<date>
  <out_base>       Output base directory
  <year>           Year folder

Optional flags:
  --sanger                       Enable Sanger primer trimming pre-step
  --primer-names <csv>           Comma-separated primer names
  --primers-file <path>          Optional custom primer file (FASTA/CSV/TSV)
  --trimmed-fasta <path>         Optional trimmed FASTA output path
  --resolved-primers-tsv <path>  Optional resolved primer coordinate TSV path
  --min-match-length <int>       Minimum seed match length for NextClade (default: 20)
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

IS_SANGER=0
PRIMER_NAMES=""
PRIMERS_FILE=""
TRIMMED_FASTA_OVERRIDE=""
RESOLVED_PRIMERS_TSV=""
MIN_MATCH_LENGTH=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sanger)
      IS_SANGER=1
      shift
      ;;
    --primer-names)
      PRIMER_NAMES="${2:-}"
      shift 2
      ;;
    --primers-file)
      PRIMERS_FILE="${2:-}"
      shift 2
      ;;
    --trimmed-fasta)
      TRIMMED_FASTA_OVERRIDE="${2:-}"
      shift 2
      ;;
    --resolved-primers-tsv)
      RESOLVED_PRIMERS_TSV="${2:-}"
      shift 2
      ;;
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

if [[ "$IS_SANGER" -eq 1 && -z "$PRIMER_NAMES" ]]; then
  echo "ERROR: --primer-names is required when --sanger is enabled." >&2
  exit 1
fi

PROJECT_DIR="$(pwd)"
LINEAGES_DATASET="$PROJECT_DIR/data/nextclade_datasets/hav-vp1-2b-lineages"

BATCH_NAME=$(basename "$BATCH_DIR")
THREADS=4

cd "$PROJECT_DIR"

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
echo "  Primer names : $PRIMER_NAMES"
echo "  Primers file : $PRIMERS_FILE"
echo "════════════════════════════════════════════════════════════════"
echo ""


# ── Optional pre-step: coordinate-based primer trimming ──────────────────────
ANALYSIS_FA="$BATCH_FA"
if [[ "$IS_SANGER" -eq 1 ]]; then
  echo "▶ Pre-step: Sanger primer trimming (reference coordinates)"
  echo "─────────────────────────────────────────────────────────────────"

  TRIMMED_FASTA="$PROJECT_DIR/$BATCH_DIR/${BATCH_NAME}.trimmed.fa"
  if [[ -n "$TRIMMED_FASTA_OVERRIDE" ]]; then
    TRIMMED_FASTA="$TRIMMED_FASTA_OVERRIDE"
  fi

  PRIMER_TRIM_REPORT="$OUT_BASE/primer_trimming.tsv"
  if [[ -n "$RESOLVED_PRIMERS_TSV" ]]; then
    RESOLVED_OUT="$RESOLVED_PRIMERS_TSV"
  else
    RESOLVED_OUT="$OUT_BASE/resolved_primers.tsv"
  fi

  TRIM_CMD=(
    python3 scripts/trim_primers_by_reference.py
    --input-fasta "$BATCH_FA"
    --output-fasta "$TRIMMED_FASTA"
    --reference-fasta "$LINEAGES_DATASET/reference.fasta"
    --primer-names "$PRIMER_NAMES"
    --builtin-primers "$PROJECT_DIR/data/Primers/PCR_primers/primers.fa"
    --report-tsv "$PRIMER_TRIM_REPORT"
    --resolved-primers-tsv "$RESOLVED_OUT"
    --threads "$THREADS"
  )

  if [[ -n "$PRIMERS_FILE" ]]; then
    TRIM_CMD+=(--custom-primers "$PRIMERS_FILE")
  fi

  "${TRIM_CMD[@]}"

  ANALYSIS_FA="$TRIMMED_FASTA"
  echo "  Trimmed FASTA → $ANALYSIS_FA"
  echo "  Trimming report → $PRIMER_TRIM_REPORT"
  echo ""
fi

# ── Analysis 1: BLAST ─────────────────────────────────────────────────────────
echo "▶ Step 1/5: BLAST"
echo "─────────────────────────────────────────────────────────────────"

# blast_batch.sh writes directly to OUT_BASE
# Pass arguments directly to preserve spaces in paths
bash scripts/blast_batch.sh "$ANALYSIS_FA" "$DATASET_DATE" "$DATASET_DIR" "$OUT_BASE"
echo ""

# ── Analysis 2: NextClade (lineages dataset) ──────────────────────────────────
echo "▶ Step 2/5: NextClade — lineages dataset (clade + lineage_phylo)"
echo "─────────────────────────────────────────────────────────────────"
LINEAGES_OUT="$OUT_BASE/lineages"
mkdir -p "$LINEAGES_OUT"

if [[ ! -d "$LINEAGES_DATASET" ]]; then
  echo "ERROR: Lineages dataset not found at $LINEAGES_DATASET" >&2
  exit 1
fi
# Pass arguments directly to preserve spaces in paths
nextclade run \
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
echo "▶ Step 3/5: Per-sequence trees (IQ-TREE with nearest neighbors)"
echo "─────────────────────────────────────────────────────────────────"
# Pass arguments directly to preserve spaces in paths
bash scripts/build_per_seq_trees.sh "$BATCH_DIR" "$DATASET_DATE" 30 "$OUT_BASE" "$DATASET_DIR" "$BATCH_FA"
echo "  Results → $OUT_BASE/trees/"
echo ""

# ── Analysis 4: Batch-specific phylogenetics ─────────────────────────────────
echo "▶ Step 4/5: Batch-specific trees and SNP matrices"
echo "─────────────────────────────────────────────────────────────────"
# Pass arguments directly to preserve spaces in paths
bash scripts/build_batch_analysis.sh "$BATCH_DIR" "$DATASET_DATE" "$OUT_BASE"
echo "  Results → $OUT_BASE/trees/batch/"
echo "  SNP matrix → $OUT_BASE/snp_matrices/batch_snp_distances.tsv"
echo ""

# ── Analysis 5: Batch report ──────────────────────────────────────────────────
echo "▶ Step 5/5: Generate batch report (HTML)"
echo "─────────────────────────────────────────────────────────────────"

REPORT_OUTPUT="$OUT_BASE/batch_report.html"
# Pass arguments directly to preserve spaces in paths
"$HOME/.conda/R_shared/bin/Rscript" -e "
rmarkdown::render(
  'scripts/batch_report.Rmd',
  output_file = '$REPORT_OUTPUT',
  params = list(
    batch_dir    = '$BATCH_DIR',    
    batch_fa     = '$BATCH_FA',    
    dataset_date = '$DATASET_DATE',
    n_neighbors  = 30,
    out_base     = '$OUT_BASE',
    dataset_dir  = '$DATASET_DIR'
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
echo "  Batch report             : $REPORT_OUTPUT"
echo "════════════════════════════════════════════════════════════════"


