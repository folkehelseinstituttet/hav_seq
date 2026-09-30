#!/usr/bin/env Rscript
#
# Reads a samplesheet TSV (Sample column only) and collects one consensus
# sequence per sample from a specified FASTA directory into a single
# multi-FASTA output file.
#
# Samplesheet format (TSV, header required):
#   Sample — sequence ID used in the output FASTA header
#            (InputDir column no longer required)
#
# FASTA directory: contains .fa, .fasta, or .TXT files for each sample.
#                  Files are matched by sample name (wildcard search).
#
# For .TXT files (raw Sanger contig text): all lines are concatenated
# as sequence.
# For .fa / .fasta files: the first sequence in the file is taken; the Sample
#   column overrides whatever header is in the file.
#
# Usage:
#   Rscript scripts/prepare_input_fasta.R <fasta_dir> <samplesheet.tsv> <output.fasta>

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(seqinr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript scripts/prepare_input_fasta.R <fasta_dir> <samplesheet.tsv> <output.fasta>",
    call. = FALSE
  )
}

fasta_dir        <- args[1]
samplesheet_path <- args[2]
output_fa        <- args[3]

# Ensure output file has .fasta extension
if (!grepl("\\.fasta$", output_fa)) {
  output_fa <- sub("\\.fa$|\\.txt$", ".fasta", output_fa, ignore.case = TRUE)
  if (!grepl("\\.fasta$", output_fa)) {
    output_fa <- paste0(output_fa, ".fasta")
  }
}

if (file.exists(output_fa) && file.info(output_fa)$size > 0) {
  message("Output FASTA already exists and is non-empty: ", output_fa)
  quit(save = "no", status = 0)
}

# ── Validate FASTA directory ──────────────────────────────────────────────────
if (!dir.exists(fasta_dir)) {
  stop("FASTA directory not found: ", fasta_dir, call. = FALSE)
}

# ── Read & validate samplesheet ───────────────────────────────────────────────
if (!file.exists(samplesheet_path)) {
  stop("Samplesheet not found: ", samplesheet_path, call. = FALSE)
}

df <- read_tsv(samplesheet_path, show_col_types = FALSE)

if (!"Sample" %in% colnames(df)) {
  stop("Samplesheet missing required column: Sample", call. = FALSE)
}

df$Sample <- trimws(df$Sample)

if (any(is.na(df$Sample) | df$Sample == "")) {
  stop("Empty Sample values detected in samplesheet", call. = FALSE)
}
if (any(duplicated(df$Sample))) {
  dups <- unique(df$Sample[duplicated(df$Sample)])
  stop("Duplicate Sample IDs: ", paste(dups, collapse = ", "), call. = FALSE)
}

message("Samples: ", nrow(df))

# ── Read one sequence per sample ──────────────────────────────────────────────
# Returns list(name, seq) or stops with an informative error.
read_sample_seq <- function(sample_name, input_dir) {
  # First token only in FASTA headers (tools truncate at whitespace)
  clean_name <- str_replace(sample_name, "\\s+.*$", "")

  # Look for files matching sample name:  .fa
  fasta_files <- list.files(input_dir, pattern = paste0("^", sample_name, ".*\\.fa$"),
                            full.names = TRUE, ignore.case = TRUE)
 
  if (length(fasta_files) > 0) {
    seqs <- read.fasta(fasta_files[[1]], seqtype = "DNA",
                       as.string = TRUE, forceDNAtolower = FALSE)
    if (length(seqs) == 0) {
      stop("No sequences in FASTA for sample '", sample_name, "': ",
           fasta_files[[1]], call. = FALSE)
    }
    # Use Sample name from samplesheet, not the embedded FASTA header
    return(list(name = clean_name, seq = as.character(seqs[[1]])))
  }

  stop("No .fa, file found for sample '",
       sample_name, "' in: ", input_dir, call. = FALSE)
}

all_seqs <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  all_seqs[[i]] <- read_sample_seq(df$Sample[i], fasta_dir)
}

# ── Write multi-FASTA ─────────────────────────────────────────────────────────
out_dir <- dirname(output_fa)
if (out_dir == "") out_dir <- "."
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fasta_lines <- unlist(lapply(
  all_seqs, function(s) c(paste0(">", s$name), s$seq)
))
writeLines(fasta_lines, output_fa)

message("Written ", length(all_seqs), " sequence(s) to: ", output_fa)
