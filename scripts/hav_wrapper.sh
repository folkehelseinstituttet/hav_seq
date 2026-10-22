#!/usr/bin/env bash
# scripts/hav_wrapper.sh
#
# Master wrapper for the HAV batch analysis pipeline.
#
# Usage:
#   bash scripts/hav_wrapper.sh --mode sanger <batch_name> <year> [OPTIONS]
#
# Modes:
#   sanger   Short-fragment Sanger sequences (VP1/2A or VP1/2B junction)
#   wgs      Not implemented yet
#
# Required arguments:
#   <batch_name>  Name of this batch (e.g., HAGEN1-0706-1 or HAV-2026_01)
#   <year>        Year folder where batch is stored
#
# Required for Sanger:
#   --samplesheet <file>    TSV with columns:
#                           Sanger: Sample, LabWareID (auto-generated if not provided;
#                                   file names must be: LabWareID_SekvensID-HAV.fa(sta))
#
# Common options:
#   --dataset-date <date>   Local dataset version to use (default: 2026-04-10)
#   --threads <n>           CPU threads (default: 4)
#   -h, --help
#
# Workflow steps by mode:
#   sanger:
#     0. generate_samplesheet (auto if --samplesheet not provided)
#     1. prepare_input_fasta.R   — fasta_dir + samplesheet → <batch>.fasta
#     2. run_all_analyses.sh     — BLAST + NextClade
#     3. build_per_seq_trees.sh  — per-sequence ML trees
#     4. batch_report.Rmd        — HTML report
#
# Prerequisites:
#   - Conda environments in envs/ (BLAST, NEXTCLADE, iqtree-mafft, SEQKIT, R_shared)
#   - Local dataset built: bash scripts/make_blast_db.sh <dataset_date>
#   - Institute inputs staged by ngs_scripts/hav_sequencing/wrapper.sh

set -euo pipefail


#Legge til aktivering av conda miljø i starten av script??
#Hva med mount til N: ? Løser seg på NGS4. 

# ── Constants ─────────────────────────────────────────────────────────────────

DEFAULT_THREADS=4
BATCH_DIR="$TMP_DIR"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'USAGE'
Usage:
  bash scripts/hav_wrapper.sh --mode sanger <batch_name> <year> [OPTIONS]

Modes:
  sanger   Short-fragment Sanger sequences (VP1/2A or VP1/2B junction)
  wgs      Not implemented yet

Required argument:
  <batch_name>             Name of this batch
  <year>                   Year of this batch 

Common options:
  --dataset-date <date>    Dataset version (default: 2026-04-10)
  --threads <n>            CPU threads (default: 4)
  -h, --help

Required:
  --samplesheet <path>     TSV with Sample + other columns
                           Sanger: auto-generated from FASTA filenames if not provided
                                   Expected format: LabWareID_SekvensID-HAV.fa(sta)

Examples:
  # Sanger batch — auto-generate samplesheet from FASTA filenames
  bash scripts/hav_wrapper.sh --mode sanger HAGEN1-0706-1 2026
USAGE
}

# ── Argument parsing ──────────────────────────────────────────────────────────
MODE=""
BATCH_NAME=""
YEAR=""
THREADS="$DEFAULT_THREADS"
SAMPLESHEET=""
BATCH_FA=""

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:?--mode requires a value (sanger)}"
      shift 2 ;;
    --dataset-date)
      DATASET_DATE="${2:?--dataset-date requires a value}"
      shift 2 ;;
    --threads)
      THREADS="${2:?--threads requires a value}"
      shift 2 ;;
    --samplesheet)
      SAMPLESHEET="${2:?--samplesheet requires a value}"
      shift 2 ;;
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

# For Sanger mode, samplesheet will be auto-generated if not provided.
shopt -s nullglob
SAMPLESHEET_ARRAY=("$BATCH_DIR"/*samplesheet.tsv)
shopt -u nullglob
if [[ ${#SAMPLESHEET_ARRAY[@]} -gt 0 ]]; then
  SAMPLESHEET="${SAMPLESHEET_ARRAY[0]}"
fi

FASTA_DIR="$TMP_DIR/Fasta"
BATCH_FA="$FASTA_DIR/$BATCH_NAME.fasta"
OUT_BASE="$HOME/$BATCH_NAME"_results
DATASET_DIR="$TMP_DIR/local_dataset"


# ── Validation ────────────────────────────────────────────────────────────────
if [[ -z "$MODE" ]]; then
  echo "ERROR: --mode sanger is required" >&2; usage; exit 1
fi
if [[ "$MODE" != "sanger" && "$MODE" != "wgs" ]]; then
  echo "ERROR: --mode must be 'sanger' or 'wgs', got: '$MODE'" >&2; exit 1
fi
if [[ "$MODE" == "wgs" ]]; then
  echo "ERROR: WGS analysis is not implemented yet. Only --mode sanger is supported." >&2
  exit 1
fi
if [[ -z "$BATCH_DIR" ]]; then
  echo "ERROR: <batch_dir> is required" >&2; usage; exit 1
fi
# Sanger samplesheets can be auto-generated from FASTA filenames in $FASTA_DIR.
SAMPLESHEET_OPTIONAL=0
if [[ "$MODE" == "sanger" ]]; then
  SAMPLESHEET_OPTIONAL=1
fi
if [[ -z "$SAMPLESHEET" && "$SAMPLESHEET_OPTIONAL" -eq 0 ]]; then
  echo "ERROR: --samplesheet <path> is required" >&2
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







# TODO(WGS): Implement run_pipeline.sh and consensus collection before
# restoring the WGS path. WGS mode is rejected during validation above.

# ════════════════════════════════════════════════════════════════════════════
# SANGER BRANCH
# ════════════════════════════════════════════════════════════════════════════
if [[ "$MODE" == "sanger" ]]; then

conda run -n NEXTCLADE nextclade dataset get \
  --name="community/masphl-bioinformatics/hav/vp1-2b-junction" \
  --output-dir="/mnt/tempdata/hav_input/nextclade_hav_vp1_2b"

export LINEAGES_DATASET="$TMP_DIR/nextclade_hav_vp1_2b"

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

export DATASET_DATE=$(find "$TMP_DIR/local_dataset" -maxdepth 1 -type d -name "????-??-??" | xargs basename)

echo "Running all analyses with the following parameters:"
echo "  BATCH_DIR     = $BATCH_DIR"
echo "  DATASET_DATE  = $DATASET_DATE"
echo "  BATCH_FA      = $BATCH_FA"
echo "  DATASET_DIR   = $DATASET_DIR"
echo "  OUT_BASE      = $OUT_BASE"
echo "  YEAR          = $YEAR"

# Pass arguments directly to preserve spaces in paths
bash "$HAV_SEQ_REPO/scripts/run_all_analyses.sh" "$BATCH_DIR" "$DATASET_DATE" "$BATCH_FA" "$DATASET_DIR" "$OUT_BASE" "$YEAR"

# ════════════════════════════════════════════════════════════════════════════
# ADD FASTA TO DATABASE
# ════════════════════════════════════════════════════════════════════════════

step "Add FASTA to database"

# Opprett mappe for ny database
HAV_DB_DIR="$HOME/hav_database"
mkdir -p "$HAV_DB_DIR"

# Dagens dato
TODAY=$(date +%F)

# Flytt gammel database til arkivnavn
cp "$DATASET_DIR/2PA.fa" "$HAV_DB_DIR/2PA_pre_${TODAY}.fa"

# Kopier gammel database tilbake som utgangspunkt for ny
cp "$HAV_DB_DIR/2PA_pre_${TODAY}.fa" "$HAV_DB_DIR/2PA.fa"

# Legg til nye sekvenser
# BATCH_FA inneholder multifasta med alle prøvene i batchen
cat "$BATCH_FA" >> "$HAV_DB_DIR/2PA.fa"
echo "  New sequences added to $HAV_DB_DIR/2PA.fa"
