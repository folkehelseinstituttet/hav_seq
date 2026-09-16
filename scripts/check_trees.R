#!/usr/bin/env Rscript
# Check which sequence names are expected vs which tree directories exist

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript check_trees.R <out_base> <nc_lineages.tsv> <batch_dir>\n")
  cat("  out_base: e.g., /mnt/n/.../output\n")
  cat("  nc_lineages.tsv: e.g., /mnt/n/.../output/lineages/nextclade.tsv\n")
  cat("  batch_dir: e.g., /mnt/n/.../TestBatchHAV\n")
  quit(status = 1)
}

out_base <- args[1]
nc_file <- args[2]
batch_dir <- args[3]

trees_dir <- file.path(out_base, "trees")

# Read NextClade results
nc <- if (file.exists(nc_file)) {
  read_tsv(nc_file, show_col_types = FALSE) %>% select(seqName)
} else {
  tibble(seqName = character())
}

# List actual tree directories
if (dir.exists(trees_dir)) {
  tree_dirs <- list.dirs(trees_dir, recursive = FALSE, full.names = FALSE)
} else {
  tree_dirs <- character()
}

# Compare
cat("NextClade sequences: ", length(nc$seqName), "\n")
cat("Tree directories: ", length(tree_dirs), "\n\n")

# Check for each sequence
results <- tibble(
  seqName = nc$seqName,
  dir_exists = seqName %in% tree_dirs,
  treefile = NA_character_,
  aligned_fasta = NA_character_
)

for (i in seq_len(nrow(results))) {
  dir_path <- file.path(trees_dir, results$seqName[i])
  treefile_path <- file.path(dir_path, "tree.treefile")
  aln_path <- file.path(dir_path, paste0(results$seqName[i], ".aligned_trimmed.fasta"))
  
  results$treefile[i] <- if_else(file.exists(treefile_path), "✓", "✗")
  results$aligned_fasta[i] <- if_else(file.exists(aln_path), "✓", "✗")
}

print(results, n = Inf)

cat("\n=== Summary ===\n")
cat("Directories missing: ", sum(!results$dir_exists), " / ", nrow(results), "\n")
cat("Tree files missing: ", sum(results$treefile == "✗"), " / ", nrow(results), "\n")
cat("Aligned FASTA missing: ", sum(results$aligned_fasta == "✗"), " / ", nrow(results), "\n")

# Show mismatches
if (any(!results$dir_exists)) {
  cat("\nSequence names in NextClade but NOT in tree directories:\n")
  cat(paste(" - ", results$seqName[!results$dir_exists], "\n", sep = ""))
}

# Show tree dirs not in NextClade
extra_dirs <- setdiff(tree_dirs, nc$seqName)
if (length(extra_dirs) > 0) {
  cat("\nTree directories NOT in NextClade:\n")
  cat(paste(" - ", extra_dirs, "\n", sep = ""))
}
