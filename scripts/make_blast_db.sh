#!/usr/bin/env bash
# scripts/make_blast_db.sh
#
# Builds a BLAST nucleotide database from the prepared input FASTA for a
# given dataset date.  The database uses -parse_seqids so that individual
# sequences can be retrieved later with blastdbcmd.
#
# Run from the project root with HAVDEV conda active:
#   conda activate /path/to/hav_dev/.conda/HAVDEV
#   bash scripts/make_blast_db.sh [dataset_date]
#
# Output:
#   data/local_datasets/<date>/blast_db/hav.{nhr,nin,nsq,...}
#
# To retrieve a sequence by ID after running BLAST, use the deduplicated FASTA:
#   grep -A2 '^><id>' data/local_datasets/<date>/blast_db/input_dedup.fa
#   seqkit grep -p <id> data/local_datasets/<date>/blast_db/input_dedup.fa

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
DATASET_DATE=$(date +%Y-%m-%d)
SCRIPT_DIR="$HAV_SEQ_REPO"
BASE_REL="$2"
INPUT_FA="$1"
DB_DIR="$BASE_REL"
DB_PATH="$DB_DIR/$DATASET_DATE"

cd "$BASE_REL"

# ── Checks ────────────────────────────────────────────────────────────────────
if [[ ! -f "$INPUT_FA" ]]; then
    echo "ERROR: input FASTA not found: $BASE_REL/$INPUT_FA" >&2
    exit 1
fi

if ! command -v makeblastdb &> /dev/null; then
    echo "ERROR: makeblastdb not found. Activate the HAVDEV conda environment:" >&2
    echo "  conda activate /path/to/hav_dev/.conda/HAVDEV" >&2
    exit 1
fi

# ── Build database ────────────────────────────────────────────────────────────
mkdir -p "$DB_PATH"

# Deduplicate FASTA: keep first occurrence of each sequence ID.
# input.fa may contain entries that cause makeblastdb to fail or warn:
#   - duplicate IDs (e.g. >2PA)
#   - non-ASCII or control characters in headers
#   - empty sequences
# This step deduplicates and sanitizes headers before building the database.
DEDUP_FA="$DB_PATH/input_dedup.fa"
python3 - "$INPUT_FA" "$DEDUP_FA" << 'PYEOF'
import sys, re
in_fa, out_fa = sys.argv[1], sys.argv[2]
seen = set(); cur_id = cur_lines = cur_seq_len = None
kept = skipped_dup = skipped_bad = 0

def is_valid(id_str, seq_lines):
    """Reject sequences with no data, non-ASCII IDs, or control chars in header."""
    seq = "".join(l.strip() for l in seq_lines)
    if not seq:
        return False
    # header must be ASCII printable only
    try:
        id_str.encode("ascii")
    except UnicodeEncodeError:
        return False
    if any(ord(c) < 32 for c in id_str):
        return False
    return True

def sanitize_header(header_line):
    """Replace spaces with underscores, strip non-ASCII from the header."""
    h = header_line.rstrip("\n")
    # encode to ASCII, replacing non-ASCII with empty
    h = h.encode("ascii", errors="ignore").decode("ascii")
    # replace control chars and spaces with underscore
    h = re.sub(r"[\x00-\x1f\s]+", "_", h[1:])
    return ">" + h + "\n"

with open(in_fa, encoding="utf-8", errors="replace") as f, open(out_fa, "w") as out:
    for line in f:
        if line.startswith(">"):
            if cur_id is not None:
                raw_id = cur_id
                if raw_id not in seen and is_valid(cur_id, cur_lines[1:]):
                    seen.add(raw_id)
                    out.write(sanitize_header(cur_lines[0]))
                    out.writelines(cur_lines[1:])
                    kept += 1
                elif raw_id in seen:
                    skipped_dup += 1
                else:
                    skipped_bad += 1
            cur_id = line[1:].split()[0].strip()
            cur_lines = [line]
        elif cur_id is not None:
            cur_lines.append(line)
    if cur_id is not None:
        if cur_id not in seen and is_valid(cur_id, cur_lines[1:]):
            seen.add(cur_id)
            out.write(sanitize_header(cur_lines[0]))
            out.writelines(cur_lines[1:])
            kept += 1
        elif cur_id in seen:
            skipped_dup += 1
        else:
            skipped_bad += 1

print(f"  Preprocessing: kept {kept}, skipped {skipped_dup} duplicates, "
      f"{skipped_bad} bad headers/empty sequences", file=sys.stderr)
PYEOF

echo "── Building BLAST database ───────────────────────────────────────────────"
echo "  Input : $INPUT_FA"
echo "  Output: $DB_PATH"
echo ""

BLAST_DIR=/tmp/blastdb

ln -sfn "$DB_PATH" "$BLAST_DIR" #link til sti pga mellomrom i filbane
ln -sfn "$DEDUP_FA" "$BLAST_DIR/input_dedup.fa" #link til sti pga mellomrom i filbane

mkdir -p "$DB_PATH/blast_db"

makeblastdb \
    -in      "$BLAST_DIR/input_dedup.fa" \
    -dbtype  nucl \
    -out     "$BLAST_DIR/blast_db/hav" \
    -title   "$DATASET_DATE"

mv "$DEDUP_FA" "$DB_PATH/blast_db"

echo ""
echo "── Done ──────────────────────────────────────────────────────────────────"
echo "  Sequences in DB : $(grep -c '^>' "$DEDUP_FA")"
echo "  Database files  :"
ls -lh "$DB_PATH/blast_db"
echo ""
echo "  To retrieve a sequence by ID (run from project root):"
echo "    grep -A999 '^><id>' $DEDUP_FA | tail -n +2 | grep -B999 '^>' | head -n -1"
echo "  Or with seqkit:"
echo "    seqkit grep -p <id> $DEDUP_FA"
