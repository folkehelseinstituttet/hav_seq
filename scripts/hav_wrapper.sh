#!/usr/bin/env bash
# scripts/hav_wrapper.sh
#
# Master wrapper for the HAV batch analysis pipeline.
#
# Usage (from project root, HAVDEV conda active):
#   bash scripts/hav_wrapper.sh --mode <sanger|wgs> <batch_name> <year> [OPTIONS]
#
# Modes:
#   sanger   Short-fragment Sanger sequences (VP1/2A or VP1/2B junction)
#   wgs      Whole-genome Nanopore tiling-PCR sequences (requires ViroConstrictor)
#
# Required arguments:
#   <batch_name>  Name of this batch (e.g., HAGEN1-0706-1 or HAV-2026_01)
#   <year>        Year folder where batch is stored
#
# Required for both modes:
#   --samplesheet <file>    TSV with columns:
#                           Sanger: Sample, LabWareID (auto-generated if not provided;
#                                   file names must be: LabWareID_SekvensID-HAV.fa(sta))
#                           WGS:    Sample, Virus, MatchRef, Segmented, Primers,
#                                   Reference, Features, MinCov, Mismatch, InputDir
#
# Common options:
#   --dataset-date <date>   Local dataset version to use (default: 2026-04-10)
#   --threads <n>           CPU threads (default: 4)
#   --n-neighbors <n>       Phylo-tree neighbors per query sequence (default: 30)
#   --skip-trees            Skip build_per_seq_trees.sh
#   --skip-report           Skip batch_report.Rmd rendering
#   -h, --help
#
# Sanger-specific options:
#   --primer-names <csv>    Comma-separated primer names — enables primer trimming
#   --primers-file <path>   Custom primer FASTA/CSV/TSV (optional supplement)
#
# WGS-specific options:
#   --skip-assembly         Skip ViroConstrictor; batch FASTA must already exist
#                           (--samplesheet still required unless FASTA exists)
#
# Workflow steps by mode:
#   sanger:
#     0. generate_samplesheet (auto if --samplesheet not provided)
#     1. prepare_input_fasta.R   — fasta_dir + samplesheet → <batch>.fasta
#     2. run_all_analyses.sh     — BLAST + NextClade (with optional primer trim)
#     3. build_per_seq_trees.sh  — per-sequence ML trees
#     4. batch_report.Rmd        — HTML report
#
#   wgs:
#     1. run_pipeline.sh         — NanoPlot QC + Chopper + ViroConstrictor assembly
#     2. collect consensus        — ViroConstrictor all_consensus.fasta → <batch>.fasta
#     3. run_all_analyses.sh     — BLAST + NextClade
#     4. build_per_seq_trees.sh  — per-sequence ML trees
#     5. batch_report.Rmd        — HTML report
#
# Prerequisites:
#   - HAVDEV conda environment active (nextclade, blastn, mafft, iqtree3)
#   - $HOME/.conda/R_shared  (for R analyses; built by setup_havdev_env.sh)
#   - Local dataset built: bash scripts/make_blast_db.sh <dataset_date>
#   - Community dataset present: data/nextclade_datasets/hav-vp1-2b-lineages/

set -euo pipefail


#Legge til aktivering av conda miljø i starten av script??
#Hva med mount til N: ? Løser seg på NGS4. 

# ── Constants ─────────────────────────────────────────────────────────────────

DEFAULT_THREADS=4
DEFAULT_N_NEIGHBORS=30
BATCH_DIR="$TMP_DIR"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'USAGE'
Usage:
  bash scripts/hav_wrapper.sh --mode <sanger|wgs> <batch_name> <year> [OPTIONS]

Modes:
  sanger   Short-fragment Sanger sequences (VP1/2A or VP1/2B junction)
  wgs      Whole-genome Nanopore tiling-PCR sequences

Required argument:
  <batch_name>             Name of this batch
  <year>                   Year of this batch 

Common options:
  --dataset-date <date>    Dataset version (default: 2026-04-10)
  --threads <n>            CPU threads (default: 4)
  --n-neighbors <n>        Tree neighbors per query (default: 30)
  --skip-trees             Skip phylogenetic tree building
  --skip-report            Skip HTML report generation
  -h, --help

Required (both modes):
  --samplesheet <path>     TSV with Sample + other columns
                           Sanger: auto-generated from FASTA filenames if not provided
                                   Expected format: LabWareID_SekvensID-HAV.fa(sta)
                           WGS:    InputDir = dir with raw FASTQ files

Sanger options:
  --primer-names <csv>     Enable primer trimming; comma-separated primer names
  --primers-file <path>    Custom primer file (FASTA/CSV/TSV)

WGS options:
  --skip-assembly          Skip ViroConstrictor; batch FASTA must already exist

Examples:
  # Sanger batch — auto-generate samplesheet from FASTA filenames
  bash scripts/hav_analyse.sh --mode sanger data/Batch-1

  # Sanger batch with custom samplesheet and primer trimming
  bash scripts/hav_analyse.sh --mode sanger data/Batch-1 \
    --samplesheet data/Batch-1/samplesheet.tsv \
    --primer-names "HAV_6.1_codehop,HAV_10_codehop,HAV_8.2_codehop,HAV_11_codehop"

  # WGS batch — full pipeline from raw Nanopore reads
  bash scripts/hav_analyse.sh --mode wgs data/HAV-2026_01 \
    --samplesheet data/HAV-2026_01/samplesheet.tsv

  # WGS batch — assembly already done, run analyses only
  bash scripts/hav_analyse.sh --mode wgs data/HAV-2026_01 \
    --samplesheet data/HAV-2026_01/samplesheet.tsv --skip-assembly
USAGE
}

# ── Argument parsing ──────────────────────────────────────────────────────────
MODE=""
BATCH_NAME=""
YEAR=""
THREADS="$DEFAULT_THREADS"
N_NEIGHBORS="$DEFAULT_N_NEIGHBORS"
SKIP_TREES=0
SKIP_REPORT=0
PRIMER_NAMES=""
PRIMERS_FILE=""
SAMPLESHEET=""
SKIP_ASSEMBLY=0
BATCH_FA=""

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:?--mode requires a value (sanger or wgs)}"
      shift 2 ;;
    --dataset-date)
      DATASET_DATE="${2:?--dataset-date requires a value}"
      shift 2 ;;
    --threads)
      THREADS="${2:?--threads requires a value}"
      shift 2 ;;
    --n-neighbors)
      N_NEIGHBORS="${2:?--n-neighbors requires a value}"
      shift 2 ;;
    --skip-trees)
      SKIP_TREES=1; shift ;;
    --skip-report)
      SKIP_REPORT=1; shift ;;
    --primer-names)
      PRIMER_NAMES="${2:?--primer-names requires a value}"
      shift 2 ;;
    --primers-file)
      PRIMERS_FILE="${2:?--primers-file requires a value}"
      shift 2 ;;
    --samplesheet)
      SAMPLESHEET="${2:?--samplesheet requires a value}"
      shift 2 ;;
    --skip-assembly)
      SKIP_ASSEMBLY=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
    *)
      POSITIONAL+=("$1"); shift ;;
  esac
done

[[ ${#POSITIONAL[@]} -ge 1 ]] && BATCH_NAME="${POSITIONAL[0]}"
[[ ${#POSITIONAL[@]} -ge 2 ]] && YEAR="${POSITIONAL[1]}"


# ── Resolve paths ─────────────────────────────────────────────────────────────
SCRIPT_DIR="$HAV_SEQ_REPO"

# For Sanger mode, samplesheet will be auto-generated if not provided
# For WGS mode, it will be validated later
shopt -s nullglob
SAMPLESHEET_ARRAY=("$BATCH_DIR"/*samplesheet.tsv)
shopt -u nullglob
if [[ ${#SAMPLESHEET_ARRAY[@]} -gt 0 ]]; then
  SAMPLESHEET="${SAMPLESHEET_ARRAY[0]}"
fi

FASTA_DIR="$TMP_DIR/Fasta"
BATCH_FA="$FASTA_DIR/$BATCH_NAME.fasta"
OUT_BASE="$HOME/$BATCH_NAME"
DATASET_DIR="$TMP_DIR/local_dataset"


# ── Validation ────────────────────────────────────────────────────────────────
if [[ -z "$MODE" ]]; then
  echo "ERROR: --mode <sanger|wgs> is required" >&2; usage; exit 1
fi
if [[ "$MODE" != "sanger" && "$MODE" != "wgs" ]]; then
  echo "ERROR: --mode must be 'sanger' or 'wgs', got: '$MODE'" >&2; exit 1
fi
if [[ -z "$BATCH_DIR" ]]; then
  echo "ERROR: <batch_dir> is required" >&2; usage; exit 1
fi
# Samplesheet is required for both modes, with exceptions:
# - WGS + --skip-assembly when the batch FASTA already exists
# - Sanger: can be auto-generated from FASTA filenames in $FASTA_DIR
SAMPLESHEET_OPTIONAL=0
if [[ "$MODE" == "wgs" && "$SKIP_ASSEMBLY" -eq 1 && -f "$BATCH_FA" ]]; then
  SAMPLESHEET_OPTIONAL=1
fi
if [[ "$MODE" == "sanger" ]]; then
  SAMPLESHEET_OPTIONAL=1
fi
if [[ -z "$SAMPLESHEET" && "$SAMPLESHEET_OPTIONAL" -eq 0 ]]; then
  echo "ERROR: --samplesheet <path> is required" >&2
  if [[ "$MODE" == "wgs" ]]; then
    echo "       Use --skip-assembly to reuse an existing batch FASTA without a sheet." >&2
  fi
  exit 1
fi
if [[ -n "$SAMPLESHEET" && ! -f "$SAMPLESHEET" ]]; then
  echo "ERROR: Samplesheet not found: $SAMPLESHEET" >&2; exit 1
fi



# Sjekk at HAV_lw_uttrekk.tsv finnes
if [[ ! -f "$TMP_DIR/HAV_lw_uttrekk.tsv" ]]; then
    echo "ERROR: Mangler nødvendig fil: $TMP_DIR/HAV_lw_uttrekk.tsv" >&2
    echo "Kopier fil  fra V:\Prod\FromSecure\LW_Datauttrekk 
                      til $TMP_DIR" >&2
    exit 1
fi

# Ensure output directory exists
mkdir -p "$OUT_BASE"

# ── Logging ───────────────────────────────────────────────────────────────────
LOGFILE="/home/ngs/hav_sequencing_wrapper.log"
ERRORLOG="/home/ngs/hav_sequencing_wrapper.error.log"

# Opprett loggfil hvis den ikke finnes
if [[ ! -f "$LOGFILE" ]]; then
    touch "$LOGFILE" || {
        echo "ERROR: Could not create logfile: $LOGFILE" >&2
        exit 1
    }
fi
# Opprett error loggfil hvis den ikke finnes
if [[ ! -f "$ERRORLOG" ]]; then
    touch "$ERRORLOG" || {
        echo "ERROR: Could not create error logfile: $ERRORLOG" >&2
        exit 1
    }
fi

exec > >(tee -a "$LOGFILE") \
     2> >(tee -a "$LOGFILE" >> "$ERRORLOG")


step() {
  echo ""
  echo "▶ $*"
  echo "─────────────────────────────────────────────────────────────────"
}

echo "════════════════════════════════════════════════════════════════"
echo "HAV Analysis Pipeline"
echo "  Batch     : $BATCH_NAME"
echo "  Year      : $YEAR"
echo "  Batch Dir       : $BATCH_DIR"
echo "  Batch fasta dir : $FASTA_DIR"
#echo "  Samplesheet     : $SAMPLESHEET"
echo "  Mode      : $MODE"
#echo "  Dataset   : $DATASET_DATE"
echo "  Threads   : $THREADS"
echo "  Neighbors : $N_NEIGHBORS"
echo "  Started   : $(date)"
echo "  Log       : $LOGFILE"
echo "  Error Log : $ERRORLOG"
echo "  Output    : $OUT_BASE"
echo "════════════════════════════════════════════════════════════════"



# ════════════════════════════════════════════════════════════════════════════
# BUILD BLAST DATABASE
# ════════════════════════════════════════════════════════════════════════════

step "Build BLAST database"
conda run -n BLAST bash "$HAV_SEQ_REPO/scripts/make_blast_db.sh" "$TMP_DIR/local_dataset/2PA.fa" "$TMP_DIR/local_dataset" || exit 1


# ════════════════════════════════════════════════════════════════════════════
# PREPARE METADATA
# ════════════════════════════════════════════════════════════════════════════

step "Preparing metadata"
conda run -n R_shared Rscript "$HAV_SEQ_REPO/scripts/prepare_metadata.R" "$TMP_DIR/Requests.xlsx" "$TMP_DIR/HAV_lw_uttrekk.tsv" "$TMP_DIR/metadata.tsv" || exit 1

METADATA="$TMP_DIR/metadata.tsv"







# ════════════════════════════════════════════════════════════════════════════
# WGS BRANCH
# ════════════════════════════════════════════════════════════════════════════
#conda activate viroconstrictor
if [[ "$MODE" == "wgs" ]]; then

conda run -n NEXTCLADE nextclade dataset get \
  --name="community/masphl-bioinformatics/hav/whole-genome" \
  --output-dir="/mnt/tempdata/hav_input/nextclade_hav_whole_genome"

  if [[ "$SKIP_ASSEMBLY" -eq 0 ]]; then
    step "WGS 1/2: Nanopore QC + filtering + ViroConstrictor assembly"
    bash "$HAV_SEQ_REPO/scripts/run_pipeline.sh" "$SAMPLESHEET" "$BATCH_ABS" "$THREADS"
    echo "  Samplesheet : $SAMPLESHEET"
  else
    echo ""
    echo "  Skipping assembly (--skip-assembly)"
  fi

  step "WGS 2/2: Collecting ViroConstrictor consensus sequences → batch FASTA"
  VC_CONSENSUS="$BATCH_ABS/viroconstrictor/results/combined/all_samples/all_consensus.fasta"
  if [[ ! -f "$VC_CONSENSUS" ]]; then
    echo "ERROR: ViroConstrictor consensus not found:" >&2
    echo "  $VC_CONSENSUS" >&2
    echo "  Check that assembly completed, or run without --skip-assembly." >&2
    exit 1
  fi
  # ViroConstrictor writes headers like ">barcode84 HAV LC128713.1_IB mincov=20".
  # Keep only the first token (the sample name) so the sample is referred to
  # identically by every downstream tool: BLAST truncates qseqid at the first
  # space anyway, while Nextclade and the MSA/heatmap keep the full header —
  # stripping here makes all of them agree on "barcode84".
  awk '/^>/{print $1; next} {print}' "$VC_CONSENSUS" > "$BATCH_FA"
  N_SEQS=$(grep -c "^>" "$BATCH_FA")
  echo "  Collected $N_SEQS consensus sequences → $BATCH_FA (headers trimmed to sample name)"

fi
#conda deactivate

# ════════════════════════════════════════════════════════════════════════════
# SANGER BRANCH
# ════════════════════════════════════════════════════════════════════════════
if [[ "$MODE" == "sanger" ]]; then

conda run -n NEXTCLADE nextclade dataset get \
  --name="community/masphl-bioinformatics/hav/vp1-2b-junction" \
  --output-dir="/mnt/tempdata/hav_input/nextclade_hav_vp1_2b"

  # Auto-generate samplesheet from FASTA filenames if not provided
  if [[ -z "$SAMPLESHEET" ]] || [[ ! -f "$SAMPLESHEET" ]]; then
    step "Sanger 1/2: Generate samplesheet from FASTA directory"
    AUTO_SAMPLESHEET="$BATCH_DIR/auto_samplesheet.tsv"
    conda run -n R_shared Rscript "$HAV_SEQ_REPO/scripts/generate_samplesheet.R" "$FASTA_DIR" "$AUTO_SAMPLESHEET" || exit 1
    SAMPLESHEET="$AUTO_SAMPLESHEET"
  fi

  step "Sanger 2/2: Prepare batch FASTA from samplesheet"
  conda run -n R_shared Rscript "$HAV_SEQ_REPO/scripts/prepare_input_fasta.R" "$FASTA_DIR" "$SAMPLESHEET" "$BATCH_FA"
  if [[ ! -f "$BATCH_FA" ]]; then
    echo "ERROR: prepare_input_fasta.R completed but $BATCH_FA was not created." >&2
    exit 1
  fi

fi


# ════════════════════════════════════════════════════════════════════════════
# COMMON STEPS — BLAST + NextClade → trees → report
# ════════════════════════════════════════════════════════════════════════════
step "Running all analyses"

if [[ ! -f "$BATCH_FA" ]]; then
  echo "ERROR: Batch FASTA not found: $BATCH_FA" >&2
  exit 1
fi

# Build arguments array for optional flags
OPTIONAL_ARGS=()
if [[ "$MODE" == "sanger" && -n "$PRIMER_NAMES" ]]; then
  OPTIONAL_ARGS+=(--sanger --primer-names "$PRIMER_NAMES")
  [[ -n "$PRIMERS_FILE" ]] && OPTIONAL_ARGS+=(--primers-file "$PRIMERS_FILE")
fi

export DATASET_DATE=$(find "$TMP_DIR/local_dataset" -maxdepth 1 -type d -name "????-??-??" | xargs basename)

echo "Running all analyses with the following parameters:"
echo "  BATCH_DIR     = $BATCH_DIR"
echo "  DATASET_DATE  = $DATASET_DATE"
echo "  BATCH_FA      = $BATCH_FA"
echo "  DATASET_DIR   = $DATASET_DIR"
echo "  OUT_BASE      = $OUT_BASE"
echo "  YEAR          = $YEAR"
echo "  OPTIONAL_ARGS = ${OPTIONAL_ARGS[@]}"

# Pass arguments directly to preserve spaces in paths
bash "$HAV_SEQ_REPO/scripts/run_all_analyses.sh" "$BATCH_DIR" "$DATASET_DATE" "$BATCH_FA" "$DATASET_DIR" "$OUT_BASE" "$YEAR" "${OPTIONAL_ARGS[@]}"

# ════════════════════════════════════════════════════════════════════════════
# ADD FASTA TO DATABASE
# ════════════════════════════════════════════════════════════════════════════

step "Add FASTA to database"
#bash "$HAV_SEQ_REPO/add_fasta_to_db.sh" "$BATCH_FA" "$OUT_BASE" || exit 1
if [[ -f "$BATCH_FA" ]]; then
  bash "$HAV_SEQ_REPO/scripts/add_fasta_to_db.sh" "$BATCH_FA" "$OUT_BASE" || exit 1
fi