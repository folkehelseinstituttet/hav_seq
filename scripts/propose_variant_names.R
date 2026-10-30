#!/usr/bin/env Rscript
#
# propose_variant_names.R
# ────────────────────────────────────────────────────────────────────────────
# Propose a variant name for every batch sequence. Runs after BLAST and before
# the outbreak trees, so trees and the report (batch_report.Rmd) share the same
# names.
#
# Usage: Rscript propose_variant_names.R <out_base> <dataset_date> <dataset_dir> <batch_fa>
#
# Output: <out_base>/proposed_variants.tsv (same information as the report's
#         "Variant naming" table, plus dup_of)
#
# Rules (current year = year of today's date):
#   - Hit from another year, no hit, or >2 SNPs: next free NOR-<year>-V<n>
#   - 1-2 SNPs from a same-year variant: same stem + first free letter (b, c, ...)
#   - 0 SNPs, same year: the database variant
#   - Identical batch sequences share one name
# ────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript propose_variant_names.R <out_base> <dataset_date> <dataset_dir> <batch_fa>\n")
  quit(save = "no", status = 1)
}
out_base     <- args[1]
dataset_date <- args[2]
dataset_dir  <- args[3]
batch_fa     <- args[4]

meta_path  <- file.path(dirname(dataset_dir), "metadata.tsv")
blast_path <- file.path(out_base, "blast_results.tsv")
out_path   <- file.path(out_base, "proposed_variants.tsv")

for (p in c(meta_path, blast_path, batch_fa)) {
  if (!file.exists(p)) {
    cat("ERROR: file not found:", p, "\n")
    quit(save = "no", status = 1)
  }
}

variant_year <- as.integer(format(Sys.Date(), "%Y"))

meta <- read_tsv(meta_path, show_col_types = FALSE) %>%
  mutate(id = as.character(id), variant = as.character(variant))

blast <- read_tsv(
  blast_path,
  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
  col_types = "ccddddiiiidd",
  show_col_types = FALSE
)
if (dataset_date == "2026-04-10") blast <- blast %>% filter(sseqid != "534954")

# Same definition of nearest hit as batch_report.Rmd: fewest SNPs, then bitscore
best_hit <- blast %>%
  left_join(meta %>% transmute(sseqid = id, blast_variant = variant), by = "sseqid") %>%
  arrange(qseqid, mismatch, desc(bitscore)) %>%
  group_by(qseqid) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(qseqid, sseqid, blast_variant, blast_mismatch = mismatch)

# Batch FASTA (order defines numbering order)
lines <- trimws(readLines(batch_fa, warn = FALSE))
lines <- lines[nzchar(lines)]
is_hdr <- startsWith(lines, ">")
grp <- cumsum(is_hdr)
fasta <- tibble(
  seqName = sub("^>", "", lines[is_hdr]),
  sequence = vapply(split(lines[!is_hdr], grp[!is_hdr]), paste, character(1), collapse = "")
) %>%
  mutate(
    qseqid = word(seqName, 1),
    seq_key = str_remove_all(str_to_upper(sequence), "[-\\s]")
  ) %>%
  distinct(qseqid, .keep_all = TRUE) %>%
  left_join(best_hit, by = "qseqid") %>%
  mutate(
    hit_year = as.integer(str_extract(blast_variant, "NOR-([0-9]+)", group = 1)),
    is_different_year = !is.na(hit_year) & hit_year != variant_year
  )

all_db_variants <- unique(meta$variant[!is.na(meta$variant) & meta$variant != ""])

# Next free number for the current year
cur <- all_db_variants[grepl(paste0("^NOR-", variant_year, "-V"), all_db_variants)]
nums <- suppressWarnings(as.numeric(gsub("^NOR-[0-9]+-V([0-9]+).*", "\\1", cur)))
next_num <- if (length(nums) > 0 && any(!is.na(nums))) max(nums, na.rm = TRUE) + 1 else 1

# Next free suffix letter for a stem, skipping 'a'
suggest_suffix <- function(base_variant, existing_variants) {
  variants <- existing_variants[grepl(paste0("^", base_variant, "[a-z]?$"), existing_variants)]
  suffixes <- str_extract(variants, "[a-z]$")
  suffixes <- suffixes[!is.na(suffixes)]
  if (length(suffixes) == 0) return("b")
  highest <- sort(suffixes, decreasing = TRUE)[1]
  rawToChar(as.raw(as.integer(charToRaw(highest)) + 1L))
}

used_variants <- all_db_variants
n <- nrow(fasta)
proposed <- character(n)
display  <- character(n)
dup_of   <- rep(NA_character_, n)

for (i in seq_len(n)) {
  key <- fasta$seq_key[i]
  if (!is.na(key) && key != "") {
    j <- which(fasta$seq_key[seq_len(i - 1)] == key)[1]
    if (!is.na(j)) {
      proposed[i] <- proposed[j]
      display[i]  <- display[j]
      dup_of[i]   <- fasta$seqName[j]
      next
    }
  }
  hv <- fasta$blast_variant[i]
  mm <- fasta$blast_mismatch[i]
  if (is.na(hv) || is.na(mm) || mm > 2 || fasta$is_different_year[i]) {
    proposed[i] <- paste0("NOR-", variant_year, "-V", next_num)
    used_variants <- c(used_variants, proposed[i])
    next_num <- next_num + 1
    display[i] <- proposed[i]
  } else if (mm == 0) {
    proposed[i] <- hv
    display[i] <- hv
  } else {
    base <- str_remove(hv, "[a-z]$")
    proposed[i] <- paste0(base, suggest_suffix(base, used_variants))
    used_variants <- c(used_variants, proposed[i])
    display[i] <- paste0(proposed[i], " (", mm, " SNP diff)")
  }
}

# Genotype (NextClade clade), as in the report's variant naming table
nc_path <- file.path(out_base, "lineages", "nextclade.tsv")
nc_clade <- if (file.exists(nc_path)) {
  read_tsv(nc_path, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    transmute(qseqid = word(seqName, 1), genotype = clade) %>%
    distinct(qseqid, .keep_all = TRUE)
} else {
  tibble(qseqid = character(), genotype = character())
}

hit_dates <- meta %>%
  transmute(sseqid = id,
            blast_hit_sample_date = if ("sample_date" %in% names(meta)) as.character(sample_date) else NA_character_)

result <- fasta %>%
  mutate(proposed_variant = proposed, proposed_display = display, dup_of = dup_of) %>%
  left_join(nc_clade, by = "qseqid") %>%
  left_join(hit_dates, by = "sseqid") %>%
  mutate(
    notes = case_when(
      !is.na(dup_of) ~ paste0("Identical to ", dup_of, " in this batch"),
      !is.na(blast_mismatch) & blast_mismatch == 0 & is_different_year ~
        paste0("Identical to ", blast_variant, " (but from different year)"),
      !is.na(blast_mismatch) & blast_mismatch == 0 ~ paste0("Identical to ", blast_variant),
      !is.na(blast_mismatch) & blast_mismatch <= 2 ~ paste0("Close to ", blast_variant),
      is.na(blast_variant) ~ "No BLAST hit",
      TRUE ~ paste0(">2 SNPs from nearest variant ", blast_variant)
    )
  ) %>%
  select(qseqid, seqName, genotype, blast_hit = sseqid, blast_hit_sample_date,
         blast_variant, blast_mismatch, is_different_year,
         proposed_variant, proposed_display, dup_of, notes)

write_tsv(result, out_path, na = "")
cat(sprintf("Wrote %d proposed variant name(s) to %s\n", nrow(result), out_path))
