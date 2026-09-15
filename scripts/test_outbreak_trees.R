#!/usr/bin/env Rscript
#
# test_outbreak_trees.R - Diagnostic script to check why trees aren't generating
#

cat("=== Outbreak Trees Diagnostic ===\n\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript test_outbreak_trees.R <output_dir>\n")
  cat("  Example: Rscript test_outbreak_trees.R ~/output\n")
  quit(save = "no", status = 1)
}

output_dir <- args[1]
cat(sprintf("Output dir: %s\n", output_dir))
cat(sprintf("Exists: %s\n\n", dir.exists(output_dir)))

# Check key files
cat("=== Checking required files ===\n\n")

files_to_check <- list(
  "NextClade results" = file.path(output_dir, "lineages", "nextclade.tsv"),
  "Dataset metadata" = file.path(output_dir, "..", "..", "data", "local_datasets", "2026-08-11", "metadata_corrected.tsv"),
  "Per-sequence trees dir" = file.path(output_dir, "trees")
)

for (name in names(files_to_check)) {
  path <- files_to_check[[name]]
  exists <- file.exists(path) || dir.exists(path)
  cat(sprintf("[%s] %s\n  Path: %s\n  Exists: %s\n\n", 
              if (exists) "✓" else "✗", name, path, exists))
}

# Try to load packages
cat("=== Checking required packages ===\n\n")

packages <- c("tidyverse", "ape", "phangorn", "ggtree")
for (pkg in packages) {
  result <- require(pkg, quietly = TRUE, character.only = TRUE)
  cat(sprintf("[%s] %s\n\n", if (result) "✓" else "✗", pkg))
}

# Try to read NextClade results
cat("=== Loading NextClade results ===\n\n")

nc_file <- file.path(output_dir, "lineages", "nextclade.tsv")
if (file.exists(nc_file)) {
  tryCatch({
    nc_data <- readr::read_tsv(nc_file, show_col_types = FALSE)
    cat(sprintf("✓ Loaded %d sequences\n\n", nrow(nc_data)))
    cat("Columns:\n")
    cat(paste("  -", names(nc_data), "\n"), sep = "")
    cat("\n")
    
    # Check for lineage_phylo column
    if ("lineage_phylo" %in% names(nc_data)) {
      n_with_lineage <- sum(!is.na(nc_data$lineage_phylo))
      cat(sprintf("✓ lineage_phylo column found: %d/%d sequences have lineage assignments\n\n", 
                  n_with_lineage, nrow(nc_data)))
    } else {
      cat("✗ lineage_phylo column NOT found!\n")
      cat("Available columns: ", paste(names(nc_data), collapse = ", "), "\n\n")
    }
    
    # Show first few sequences
    cat("First 5 sequences:\n")
    print(nc_data %>% select(seqName, lineage_phylo, clade) %>% slice(1:5))
    cat("\n")
    
  }, error = function(e) {
    cat(sprintf("✗ Error loading NextClade results: %s\n\n", e$message))
  })
} else {
  cat(sprintf("✗ NextClade results file not found: %s\n\n", nc_file))
}

# Check trees directory
cat("=== Checking per-sequence trees ===\n\n")

trees_dir <- file.path(output_dir, "trees")
if (dir.exists(trees_dir)) {
  tree_seqs <- list.dirs(trees_dir, full.names = FALSE, recursive = FALSE)
  tree_seqs <- setdiff(tree_seqs, "batch")
  cat(sprintf("Found %d per-sequence tree directories:\n", length(tree_seqs)))
  cat(paste("  -", head(tree_seqs, 10), "\n"), sep = "")
  if (length(tree_seqs) > 10) cat(sprintf("  ... and %d more\n", length(tree_seqs) - 10))
  cat("\n")
  
  # Check for alignment files
  n_aln_files <- 0
  for (seq in head(tree_seqs, 3)) {
    aln_file <- file.path(trees_dir, seq, "aligned_trimmed_blast_only.fa")
    if (file.exists(aln_file)) {
      n_aln_files <- n_aln_files + 1
    }
  }
  cat(sprintf("Checked first 3: %d have aligned_trimmed_blast_only.fa\n\n", n_aln_files))
} else {
  cat(sprintf("✗ Per-sequence trees directory not found: %s\n\n", trees_dir))
}

cat("=== Diagnostic complete ===\n")
