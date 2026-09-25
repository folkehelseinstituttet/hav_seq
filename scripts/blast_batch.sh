#!/usr/bin/env bash
# scripts/blast_batch.sh
#
# BLASTs a batch FASTA of query sequences against the local HAV BLAST database
# built by make_blast_db.sh.
#
# Run from the project root with HAVDEV conda active:
#   conda activate /path/to/hav_dev/.conda/HAVDEV
#   bash scripts/blast_batch.sh <batch_fasta> <dataset_date> <dataset_dir> <out_base>
#
# Arguments:
#   batch_fasta    Path to the multi-FASTA file (e.g. data/Batch-1/Batch-1.fa)
#   dataset_date   Dataset date folder (default: 2026-04-10)
#   dataset_dir    Path to local dataset directory (e.g. /mnt/n/.../local_datasets/2026-04-10)
#   out_base       Output directory where blast_results.tsv will be written
#
# Output columns (tab-separated):
#   qseqid    query sequence ID
#   sseqid    subject (database hit) sequence ID
#   pident    % identity
#   length    alignment length
#   mismatch  number of mismatches (SNPs)
#   gapopen   number of gap-opening events
#   qstart    query alignment start
#   qend      query alignment end
#   sstart    subject alignment start
#   send      subject alignment end
#   evalue    E-value
#   bitscore  bit score
#   stitle    full FASTA header of hit (includes ID)
#   qlen      query sequence length
#   slen      subject sequence length
#   nident    number of identical positions
#   gaps      total number of gap characters
#   qcovs     % of query sequence covered by alignment
#
# All hits with e-value <= 1e-5 are reported (one best HSP per subject).
# Join with metadata_corrected.tsv on sseqid == id for genotype/lineage etc.
#
# To retrieve the full sequence for a hit, use the deduplicated FASTA:
#   grep -A2 '^><sseqid>' data/local_datasets/<date>/blast_db/input_dedup.fa
#   seqkit grep -p <sseqid> data/local_datasets/<date>/blast_db/input_dedup.fa

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
BATCH_FA="${1:?Usage: bash scripts/blast_batch.sh <batch_fasta> <dataset_date> <dataset_dir> <out_base>}"
DATASET_DATE="${2:?dataset_date is required}"
DATASET_DIR="${3:?dataset_dir is required}"
OUT_BASE="${4:?out_base is required}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DB_PATH="$DATASET_DIR/$DATASET_DATE/blast_db/hav"
THREADS=4
EVALUE="1e-5"

# Normalise batch FASTA to a project-root-relative path
if [[ "$BATCH_FA" == /* ]]; then
    # strip absolute PROJECT_DIR prefix if present
    BATCH_FA="${BATCH_FA#"$PROJECT_DIR/"}"
fi

# Output goes directly to OUT_BASE
mkdir -p "$OUT_BASE"
OUT_TSV="$OUT_BASE/blast_results.tsv"

# ── Checks ────────────────────────────────────────────────────────────────────
if [[ ! -f "$BATCH_FA" ]]; then
    echo "ERROR: batch FASTA not found: $PROJECT_DIR/$BATCH_FA" >&2
    exit 1
fi

if [[ ! -f "${DB_PATH}.nhr" && ! -f "${DB_PATH}.ndb" ]]; then
    echo "ERROR: BLAST database not found at $DB_PATH" >&2
    echo "  Expected to find: ${DB_PATH}.nhr or ${DB_PATH}.ndb" >&2
    echo "  DATASET_DIR: $DATASET_DIR" >&2
    echo "  Build with: bash scripts/make_blast_db.sh $DATASET_DATE" >&2
    exit 1
fi

if ! command -v blastn &> /dev/null; then
    echo "ERROR: blastn not found. Activate the HAVDEV conda environment:" >&2
    echo "  conda activate /path/to/hav_dev/.conda/HAVDEV" >&2
    exit 1
fi


# ── Run BLAST ─────────────────────────────────────────────────────────────────
echo "── Running blastn ────────────────────────────────────────────────────────"
echo "  Query  : $BATCH_FA"
echo "  DB     : $DB_PATH"
echo "  Output : $OUT_TSV"
echo ""

# Write header
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    qseqid sseqid pident length mismatch gapopen \
    qstart qend sstart send evalue bitscore \
    stitle qlen slen nident gaps qcovs \
    > "$OUT_TSV"

# Opprett en midlertidig mappe i din home-direktory (uten mellomrom)
mkdir -p "$HOME/blast_tmp_db"

# Lag en symbolsk lenke direkte til hav-filene i den ønskede datasettversjonen
# (unngår mellomrom i "/mnt/n/.../HAV genteknologi/..." som blastn kan streve med)
ln -sf "$DATASET_DIR/$DATASET_DATE/blast_db"/hav.* "$HOME/blast_tmp_db/"

# Sett den nye stien som DB_PATH (uten mellomrom)
DB_PATH="$HOME/blast_tmp_db/hav"

# Build BLAST command as array to properly handle paths with spaces
declare -a BLAST_CMD=(
    blastn
    -query       "$BATCH_FA"
    -db          "$DB_PATH"
    -out         "$OUT_TSV.tmp"
    -evalue      "$EVALUE"
    -max_hsps    1
    -num_threads "$THREADS"
    -outfmt      "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen nident gaps qcovs"
)

# Execute BLAST command
"${BLAST_CMD[@]}"

cat "$OUT_TSV.tmp" >> "$OUT_TSV"
rm  "$OUT_TSV.tmp"

# Dataset-specific exclusion: ID 534954 is known reverse-complemented in
# dataset 2026-04-10 and should be ignored in downstream analyses.
if [[ "$DATASET_DATE" == "2026-04-10" ]]; then
    awk -F '\t' 'NR==1 || $2 != "534954"' "$OUT_TSV" > "$OUT_TSV.filtered"
    mv "$OUT_TSV.filtered" "$OUT_TSV"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
N_HITS=$(tail -n +2 "$OUT_TSV" | wc -l)
N_QUERIES=$(grep -c "^>" "$BATCH_FA")

echo "── Done ──────────────────────────────────────────────────────────────────"
echo "  Queries  : $N_QUERIES"
echo "  Hits     : $N_HITS (e-value <= $EVALUE)"
echo "  Output   : $PROJECT_DIR/$OUT_TSV"
echo ""
echo "  To retrieve a hit sequence (from deduplicated FASTA):"
echo "    grep -A2 '^><sseqid>' $DATASET_DIR/$DATASET_DATE/blast_db/input_dedup.fa"
  echo ""
  echo "  To join with metadata:"
  echo "    data/local_datasets/$DATASET_DATE/metadata_corrected.tsv  (join on sseqid == id)"
