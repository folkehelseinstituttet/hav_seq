#!/usr/bin/env Rscript
#
# build_outbreak_trees.R
# ────────────────────────────────────────────────────────────────────────────
# Generate and visualize phylogenetic trees for each lineage found in batch.
# Creates PNG images combining batch sequences + database sequences per lineage.
#
# Usage: Rscript build_outbreak_trees.R <output_dir> <dataset_date> <dataset_dir>
#
# Output:
#   - <output_dir>/outbreak_trees/  (directory with PNG images per lineage)
#   - <output_dir>/.outbreak_trees.log  (detailed log of execution)
#
# ────────────────────────────────────────────────────────────────────────────

# === SETUP ===================================================================

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript build_outbreak_trees.R <output_dir> <dataset_date> <dataset_dir>\n")
  quit(save = "no", status = 1)
}

output_dir <- args[1]
dataset_date <- args[2]
dataset_dir <- args[3]

# Setup logging
log_file <- file.path(output_dir, ".outbreak_trees.log")
log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

log_msg("=== Starting outbreak tree generation ===")
log_msg("Output dir: %s", output_dir)
log_msg("Dataset date: %s", dataset_date)
log_msg("Dataset dir: %s", dataset_dir)

# === LOAD PACKAGES ===========================================================

log_msg("\n--- Loading R packages ---")
required_packages <- c("tidyverse", "ape", "phangorn", "ggtree")
for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    log_msg("ERROR: Cannot load required package '%s'", pkg)
    quit(save = "no", status = 1)
  }
  log_msg("✓ Loaded %s", pkg)
}

# === VALIDATE PATHS ===========================================================

log_msg("\n--- Validating paths ---")
if (!dir.exists(output_dir)) {
  log_msg("ERROR: Output directory not found: %s", output_dir)
  quit(save = "no", status = 1)
}
log_msg("✓ Output directory exists")

# Setup paths
trees_dir <- file.path(output_dir, "trees")
lineages_file <- file.path(output_dir, "lineages", "nextclade.tsv")
dataset_meta <- file.path(dataset_dir, "metadata_corrected.tsv")
outbreak_trees_dir <- file.path(output_dir, "outbreak_trees")

if (!dir.exists(trees_dir)) {
  log_msg("ERROR: Per-sequence trees directory not found: %s", trees_dir)
  quit(save = "no", status = 1)
}
log_msg("✓ Per-sequence trees directory exists")

if (!file.exists(lineages_file)) {
  log_msg("ERROR: NextClade results not found: %s", lineages_file)
  quit(save = "no", status = 1)
}
log_msg("✓ NextClade results file exists")

if (!file.exists(dataset_meta)) {
  log_msg("ERROR: Dataset metadata not found: %s", dataset_meta)
  quit(save = "no", status = 1)
}
log_msg("✓ Dataset metadata file exists")

# Create output directory
dir.create(outbreak_trees_dir, showWarnings = FALSE, recursive = TRUE)
log_msg("✓ Created/verified outbreak_trees directory")

# === LOAD DATA ===============================================================

log_msg("\n--- Loading batch lineage assignments ---")

# Read BLAST results (raw BLAST hits for batch sequences against database)
blast_file <- file.path(output_dir, "blast_results.tsv")

if (!file.exists(blast_file)) {
  log_msg("ERROR: BLAST results not found: %s", blast_file)
  quit(save = "no", status = 1)
}

tryCatch({
  # BLAST output columns (standard -outfmt 6 format)
  blast_results <- read_tsv(
    blast_file, 
    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
    col_types = "ccddddiiiidd",
    show_col_types = FALSE
  )
  log_msg("✓ Loaded BLAST results")
  log_msg("  Rows: %d (BLAST hits)", nrow(blast_results))
  
  # Verify columns are present
  if (!"qseqid" %in% names(blast_results) || !"sseqid" %in% names(blast_results)) {
    log_msg("ERROR: BLAST results missing qseqid or sseqid columns")
    quit(save = "no", status = 1)
  }
}, error = function(e) {
  log_msg("ERROR loading BLAST results: %s", conditionMessage(e))
  quit(save = "no", status = 1)
})

log_msg("\n--- Loading database metadata ---")
tryCatch({
  meta <- read_tsv(dataset_meta, show_col_types = FALSE)
  log_msg("✓ Loaded %d database sequences", nrow(meta))
  log_msg("  Columns: %s", paste(names(meta), collapse = ", "))
  
  # Check if lineage column exists (this is the variant in BLAST database metadata)
  if (!"lineage" %in% names(meta)) {
    log_msg("ERROR: 'lineage' column not found in database metadata")
    log_msg("  Available columns: %s", paste(names(meta), collapse = ", "))
    quit(save = "no", status = 1)
  }
  
  # Keep id and lineage columns (lineage = variants like NOR-2026-V6b)
  meta <- meta %>%
    select(id, lineage, genotype, date) %>%
    filter(!is.na(lineage))
  
  log_msg("✓ Extracted %d database sequences with lineage (variant) assignments", nrow(meta))
}, error = function(e) {
  log_msg("ERROR loading database metadata: %s", conditionMessage(e))
  quit(save = "no", status = 1)
})

# === HELPER FUNCTIONS ========================================================

norm_id <- function(x) {
  # Normalize IDs for matching, with encoding protection
  # Convert to UTF-8 first to avoid encoding issues
  x <- iconv(x, to = "UTF-8", sub = "")
  x <- tolower(x)
  # Remove hyphens and underscores
  x <- gsub("[-_]", "", x, perl = TRUE)
  x
}

# === ASSIGN BATCH SEQUENCES TO VARIANTS ======================================

log_msg("\n--- Assigning batch sequences to variants ---")

# Extract base variant (remove trailing a/b/c/d suffixes)
# NOR-2024-V6a → NOR-2024-V6
# NOR-2025-V8 → NOR-2025-V8 (no suffix)
extract_base_variant <- function(variant) {
  ifelse(is.na(variant), NA_character_, sub("[a-z]$", "", variant))
}

# Assign each batch sequence to the base variant of its best BLAST match
batch_lineages <- blast_results %>%
  mutate(sseqid_norm = norm_id(sseqid)) %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(
    meta %>%
      mutate(id_norm = norm_id(id)) %>%
      select(id_norm, lineage),
    by = c("sseqid_norm" = "id_norm")
  ) %>%
  select(seqName = qseqid, subject_id = sseqid, variant = lineage) %>%
  filter(!is.na(variant)) %>%
  # Group by base variant (remove sub-variant suffixes)
  mutate(base_variant = extract_base_variant(variant))

n_assigned <- nrow(batch_lineages)
log_msg("✓ Assigned %d batch sequences to variants via BLAST", n_assigned)

if (n_assigned == 0) {
  log_msg("WARNING: No batch sequences could be assigned to variants. Exiting.")
  quit(save = "no", status = 0)
}

batch_lineages_clean <- batch_lineages %>%
  filter(!is.na(base_variant))

unique_lineages <- batch_lineages_clean %>%
  distinct(base_variant) %>%
  pull(base_variant) %>%
  sort()

n_unique_variants <- length(unique_lineages)
log_msg("✓ Found %d unique base variants: %s", n_unique_variants, 
        paste(unique_lineages, collapse = ", "))

# Only generate trees if >= 3 unique base variants
if (n_unique_variants < 3) {
  log_msg("WARNING: Only %d variant(s) found. Skipping tree generation (need ≥3).", n_unique_variants)
  quit(save = "no", status = 0)
}

# === GENERATE TREES ===========================================================

log_msg("\n--- Generating phylogenetic trees ---")

n_trees_generated <- 0

for (outbreak_variant in unique_lineages) {
  log_msg("\nProcessing variant: %s", outbreak_variant)
  
  # Get batch sequences with this base variant
  batch_seqs_list <- batch_lineages_clean %>%
    filter(base_variant == outbreak_variant) %>%
    pull(seqName)
  
  n_batch <- length(batch_seqs_list)
  log_msg("  Batch sequences: %d", n_batch)
  
  # Get database sequences with this base variant
  # Match both the variant and any sub-variants (V6, V6a, V6b all match NOR-2024-V6)
  db_seqs_for_variant <- meta %>%
    mutate(base_var = extract_base_variant(lineage)) %>%
    filter(base_var == outbreak_variant) %>%
    select(id, lineage, genotype, date)
  
  n_db <- nrow(db_seqs_for_variant)
  log_msg("  Database sequences: %d", n_db)
  
  # Collect all sequences (both batch and database)
  all_outbreak_seqs <- list()
  
  # Add batch sequences
  for (seq_name in batch_seqs_list) {
    seq_aln_file <- file.path(trees_dir, seq_name, "aligned_trimmed_blast_only.fa")
    if (!file.exists(seq_aln_file)) {
      log_msg("    ✗ Alignment file not found: %s", basename(seq_aln_file))
      next
    }
    
    tryCatch({
      seqs <- ape::read.dna(seq_aln_file, format = "fasta", as.character = TRUE)
      if (nrow(seqs) > 0) {
        for (j in seq_len(nrow(seqs))) {
          seq_id <- rownames(seqs)[j]
          if (j == 1) seq_id <- paste0(seq_name, "_BATCH")
          all_outbreak_seqs[[seq_id]] <- seqs[j, ]
        }
        log_msg("    ✓ Loaded %d sequence(s) from %s", nrow(seqs), seq_name)
      }
    }, error = function(e) {
      log_msg("    ✗ Error reading %s: %s", seq_name, e$message)
    })
  }
  
  # Add database sequences (read FASTA files from database directory)
  if (nrow(db_seqs_for_variant) > 0) {
    # Try to find database FASTA file
    # Look for input.fa, input_dedup.fa, or references.fa
    db_fasta_candidates <- c(
      file.path(dataset_dir, "input_dedup.fa"),
      file.path(dataset_dir, "input.fa"),
      file.path(dataset_dir, "reference.fasta"),
      file.path(dataset_dir, "..", "References", "IA_IB_IIA_IIB_IIIA_IIIB_references.fa")
    )
    
    db_fasta_file <- NULL
    for (candidate in db_fasta_candidates) {
      if (file.exists(candidate)) {
        db_fasta_file <- candidate
        log_msg("    Found database FASTA: %s", basename(candidate))
        break
      }
    }
    
    if (!is.null(db_fasta_file)) {
      tryCatch({
        # Try to read with error handling for encoding issues
        db_seqs_all <- NULL
        tryCatch({
          db_seqs_all <- ape::read.dna(db_fasta_file, format = "fasta", as.character = TRUE)
        }, error = function(e) {
          # If encoding error, try reading with Biostrings as fallback
          log_msg("    ⚠ Encoding issue with ape::read.dna, attempting alternative method")
          return(NULL)
        })
        
        if (!is.null(db_seqs_all)) {
          log_msg("    Loaded %d total sequences from database FASTA", nrow(db_seqs_all))
          
          # Match database sequences by their ID
          n_matched <- 0
          for (i in seq_len(nrow(db_seqs_for_variant))) {
            db_id <- db_seqs_for_variant$id[i]
            db_id_norm <- norm_id(db_id)
            
            # Find matching sequence in FASTA file
            # Try exact match first, then normalized match
            seq_idx <- which(tolower(rownames(db_seqs_all)) == tolower(db_id))
            if (length(seq_idx) == 0) {
              seq_idx <- which(norm_id(rownames(db_seqs_all)) == db_id_norm)
            }
            
            if (length(seq_idx) > 0) {
              all_outbreak_seqs[[db_id]] <- db_seqs_all[seq_idx[1], ]
              n_matched <- n_matched + 1
            }
          }
          
          log_msg("    ✓ Loaded %d database sequence(s)", n_matched)
        } else {
          log_msg("    ⚠ Could not load database sequences due to encoding issues")
        }
      }, error = function(e) {
        log_msg("    ✗ Error reading database sequences: %s", e$message)
      })
    } else {
      log_msg("    ⚠ Database FASTA file not found in expected locations")
    }
  }
  
  n_total <- length(all_outbreak_seqs)
  log_msg("  Total sequences for tree: %d", n_total)
  
  # Only generate trees with sufficient sequences for meaningful phylogeny
  if (n_total < 3) {
    log_msg("  ⚠ SKIP: Only %d sequence(s), need ≥3 for phylogenetic analysis", n_total)
    next
  }
  
  # Generate tree
  tryCatch({
    # Convert to matrix
    outbreak_seqs_mat <- do.call(rbind, all_outbreak_seqs)
    
    # Write temporary alignment to output directory (for better control)
    aln_file_tmp <- file.path(outbreak_trees_dir, sprintf(".tmp_tree_%s.fa", gsub("[^A-Za-z0-9.-]", "_", outbreak_variant)))
    cat_lines <- character()
    for (j in seq_len(nrow(outbreak_seqs_mat))) {
      cat_lines <- c(cat_lines, paste0(">", rownames(outbreak_seqs_mat)[j]))
      cat_lines <- c(cat_lines, paste(outbreak_seqs_mat[j, ], collapse = ""))
    }
    writeLines(cat_lines, aln_file_tmp)
    log_msg("    Alignment written to: %s", basename(aln_file_tmp))
    
    # Build tree with IQ-TREE
    tree_out_prefix <- sub("\\.fa$", "", aln_file_tmp)
    
    # Capture both stdout and stderr to temporary files for debugging
    stdout_log <- file.path(outbreak_trees_dir, sprintf(".iqtree_stdout_%s.log", gsub("[^A-Za-z0-9.-]", "_", outbreak_variant)))
    stderr_log <- file.path(outbreak_trees_dir, sprintf(".iqtree_stderr_%s.log", gsub("[^A-Za-z0-9.-]", "_", outbreak_variant)))
    
    iqtree_cmd <- sprintf("iqtree -s %s -m JC -nt AUTO -fast -redo 2>&1 | tee %s", aln_file_tmp, stdout_log)
    iqtree_exit <- system(iqtree_cmd)
    
    log_msg("    IQ-TREE exit code: %d", iqtree_exit)
    
    # IQ-TREE appends .treefile to the INPUT filename, not the stem
    # So if input is file.fa, it creates file.fa.treefile
    tree_file <- paste0(aln_file_tmp, ".treefile")
    
    log_msg("    Looking for: %s", basename(tree_file))
    log_msg("    Exists: %s", file.exists(tree_file))
    
    if (!file.exists(tree_file)) {
      # Show IQ-TREE output for debugging
      if (file.exists(stdout_log)) {
        log_msg("    IQ-TREE stdout:")
        iqtree_out <- readLines(stdout_log)
        for (line in iqtree_out) {
          log_msg("      %s", line)
        }
      }
      log_msg("  ✗ IQ-TREE tree file not found")
      unlink(paste0(tree_out_prefix, ".*"), force = TRUE)
      next
    }
    
    # Read and root tree
    outbreak_tree <- read.tree(tree_file)
    outbreak_tree <- phangorn::midpoint(outbreak_tree)
    log_msg("  ✓ Tree built with %d tips", length(outbreak_tree$tip.label))
    
    # Validate tree structure before visualization
    if (length(outbreak_tree$tip.label) == 0 || any(is.na(outbreak_tree$tip.label))) {
      log_msg("  ✗ Invalid tree structure (empty or missing tip labels)")
      unlink(paste0(tree_out_prefix, ".*"), force = TRUE)
      next
    }
    
    # Prepare metadata for tree tips - simple source indicator
    outbreak_tip_info <- tibble(label = outbreak_tree$tip.label) %>%
      mutate(
        source = if_else(grepl("_BATCH$", label, ignore.case = TRUE), "Batch", "Database")
      )
    
    # Plot tree with improved error handling
    p <- NULL
    tryCatch({
      p <- ggtree(outbreak_tree, layout = "rectangular", branch.length = "branch.length") %<+% outbreak_tip_info +
        geom_tippoint(aes(color = source, shape = source), size = 3) +
        geom_tiplab(aes(label = label), size = 2.5, hjust = -0.05) +
        scale_color_manual(
          values = c("Batch" = "#e74c3c", "Database" = "#3498db"),
          name = "Source"
        ) +
        scale_shape_manual(
          values = c("Batch" = 19, "Database" = 17),
          guide = "none"
        ) +
        theme_tree2() +
        theme(legend.position = "bottom") +
        hexpand(.6) +
        labs(
          title = sprintf("Variant %s — Phylogenetic tree", outbreak_variant),
          subtitle = sprintf("Batch vs database sequences (n=%d total)", n_total)
        )
    }, error = function(e) {
      log_msg("  ✗ Error in ggtree visualization: %s", e$message)
    })
    
    if (is.null(p)) {
      log_msg("  ✗ Failed to create tree visualization")
      unlink(paste0(tree_out_prefix, ".*"), force = TRUE)
      next
    }
    
    # Save PNG with error handling
    tryCatch({
      png_file <- file.path(outbreak_trees_dir, sprintf("%s.png", outbreak_variant))
      ggsave(png_file, plot = p, width = 12, height = 8, dpi = 150)
      log_msg("  ✓ Saved tree PNG: %s", basename(png_file))
      n_trees_generated <- n_trees_generated + 1
    }, error = function(e) {
      log_msg("  ✗ Error saving PNG: %s", e$message)
    })
    
    # Cleanup temp files
    unlink(paste0(tree_out_prefix, ".*"), force = TRUE)
    unlink(stdout_log, force = TRUE)
    unlink(stderr_log, force = TRUE)
    
  }, error = function(e) {
    log_msg("  ✗ Error: %s", e$message)
  })
}

# === COMPLETION ==============================================================

log_msg("\n=== Completion ===")
log_msg("Generated %d outbreak trees", n_trees_generated)
log_msg("Saved to: %s", outbreak_trees_dir)
log_msg("Log file: %s", log_file)
