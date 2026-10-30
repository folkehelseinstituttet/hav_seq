#!/usr/bin/env Rscript
#
# build_outbreak_trees.R
# ────────────────────────────────────────────────────────────────────────────
# Generate and visualize phylogenetic trees for each variant found in batch.
# Creates PNG images combining batch sequences + database sequences per variant.
#
# Usage: Rscript build_outbreak_trees.R <output_dir> <dataset_date> <dataset_dir>
#
# Output:
#   - <output_dir>/outbreak_trees/  (directory with PNG images per variant)
#   - <output_dir>/.outbreak_trees.log  (detailed log of execution)
#
# ────────────────────────────────────────────────────────────────────────────

# === SETUP ===================================================================

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat(
    paste0(
      "Usage: Rscript build_outbreak_trees.R ",
      "<output_dir> <dataset_date> <dataset_dir> ",
      "<mafft_bin> <iqtree_bin>\n"
    )
  )
  quit(save = "no", status = 1)
}

output_dir <- args[1]
dataset_date <- args[2]
dataset_dir <- args[3]
mafft_bin <- args[4]
iqtree_bin <- args[5]

dataset_version_dir <- file.path(dataset_dir, dataset_date)
blast_db_dir <- file.path(dataset_version_dir, "blast_db")
dataset_root <- dirname(dataset_dir)
dataset_meta <- file.path(dataset_root, "metadata.tsv")

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
log_msg("Dataset version dir: %s", dataset_version_dir)
log_msg("BLAST DB dir: %s", blast_db_dir)
log_msg("Dataset root: %s", dataset_root)
log_msg("Dataset metadata: %s", dataset_meta)

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

log_msg("\n--- Loading batch variant assignments ---")

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
  
  # Check if variant column exists (this is the variant in BLAST database metadata)
  if (!"variant" %in% names(meta)) {
    log_msg("ERROR: 'variant' column not found in database metadata")
    log_msg("  Available columns: %s", paste(names(meta), collapse = ", "))
    quit(save = "no", status = 1)
  }
  
  # Keep id and variant columns (variant = variants like NOR-2026-V6b)
  meta <- meta %>%
    select(id, genotype, variant) %>%
    filter(!is.na(variant))
  
  log_msg("✓ Extracted %d database sequences with variant assignments", nrow(meta))
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

# Robust FASTA reader (readLines + iconv) — avoids ape::read.dna crashing on
# non-UTF8 bytes in headers (common with Norwegian names in local database FASTA).
# Returns a named character vector: name = header (without ">"), value = ungapped sequence.
read_fasta_raw <- function(path) {
  raw <- readLines(path, warn = FALSE)
  raw <- iconv(raw, sub = "")
  raw <- raw[!is.na(raw)]
  
  ids <- character()
  seqs <- character()
  cur_id <- NA_character_
  cur_seq <- character()
  
  flush <- function() {
    if (!is.na(cur_id)) {
      ids[[length(ids) + 1]] <<- cur_id
      seqs[[length(seqs) + 1]] <<- gsub("-", "", paste(cur_seq, collapse = ""))
    }
  }
  
  for (ln in raw) {
    if (startsWith(ln, ">")) {
      flush()
      cur_id <- sub("^>", "", ln)
      cur_seq <- character()
    } else if (!is.na(cur_id)) {
      cur_seq <- c(cur_seq, trimws(ln))
    }
  }
  flush()
  
  setNames(seqs, ids)
}

# === ASSIGN BATCH SEQUENCES TO VARIANTS ======================================

log_msg("\n--- Assigning batch sequences to variants ---")

# Extract base variant (remove trailing a/b/c/d suffixes)
# NOR-2024-V6a → NOR-2024-V6
# NOR-2025-V8 → NOR-2025-V8 (no suffix)
extract_base_variant <- function(variant) {
  variant <- iconv(
    as.character(variant),
    from = "",
    to = "UTF-8",
    sub = ""
  )

  variant <- trimws(variant)
  variant[is.na(variant) | variant == ""] <- NA_character_

  sub("[a-z]$", "", variant, perl = TRUE)
}

# Assign each batch sequence to a group. Sequences are linked when they share a base
# variant (V6 ~ V6a ~ V6b) through their best BLAST hit or any hit within <=3 SNPs,
# directly or transitively. This mirrors the grouping in batch_report.Rmd.
batch_hits <- blast_results %>%
  arrange(qseqid, mismatch, desc(bitscore)) %>%
  mutate(sseqid_norm = norm_id(sseqid)) %>%
  left_join(
    meta %>%
      mutate(id_norm = norm_id(id)) %>%
      select(id_norm, variant),
    by = c("sseqid_norm" = "id_norm")
  ) %>%
  filter(!is.na(variant), nzchar(trimws(variant))) %>%
  group_by(qseqid) %>%
  filter(row_number() == 1 | mismatch <= 3) %>%
  mutate(is_best = row_number() == 1) %>%
  ungroup() %>%
  mutate(base_variant = extract_base_variant(variant)) %>%
  filter(!is.na(base_variant))

batch_variants <- batch_hits %>%
  filter(is_best) %>%
  select(seqName = qseqid, subject_id = sseqid, variant, base_variant)

n_assigned <- nrow(batch_variants)
log_msg("✓ Assigned %d batch sequences to variants via BLAST", n_assigned)

if (n_assigned == 0) {
  log_msg("WARNING: No batch sequences could be assigned to variants. Exiting.")
  quit(save = "no", status = 0)
}

# Connected components of base variants (union via shared batch sequences)
seq_base <- batch_hits %>% distinct(seqName = qseqid, base_variant)
bases_all <- sort(unique(seq_base$base_variant))
comp <- setNames(seq_along(bases_all), bases_all)
repeat {
  seq_min <- tapply(comp[seq_base$base_variant], seq_base$seqName, min)
  cand <- as.integer(seq_min[seq_base$seqName])
  upd <- tapply(cand, seq_base$base_variant, min)
  new_comp <- comp
  new_comp[names(upd)] <- pmin(comp[names(upd)], as.integer(upd))
  if (all(new_comp == comp)) break
  comp <- new_comp
}

# group_key = sorted base variants joined by "_" (also the PNG/FASTA file name)
group_bases <- split(names(comp), unname(comp))
group_keys <- vapply(group_bases, function(b) paste(sort(b), collapse = "_"), character(1))
group_bases <- setNames(lapply(group_bases, sort), group_keys)
base_to_key <- setNames(rep(group_keys, lengths(group_bases)), unlist(group_bases))

batch_variants_clean <- batch_variants %>%
  mutate(group_key = unname(base_to_key[base_variant]))

unique_variant <- sort(unique(batch_variants_clean$group_key))

n_unique_variants <- length(unique_variant)
log_msg("✓ Found %d variant group(s): %s", n_unique_variants,
        paste(unique_variant, collapse = ", "))

# Proposed variant names (from propose_variant_names.R), used in tree labels
proposed_file <- file.path(output_dir, "proposed_variants.tsv")
proposed_names <- if (file.exists(proposed_file)) {
  read_tsv(proposed_file, col_types = cols(.default = col_character()), show_col_types = FALSE)
} else {
  log_msg("⚠ proposed_variants.tsv not found; labels show nearest-hit variant")
  tibble(qseqid = character(), proposed_variant = character())
}

if (n_unique_variants == 0) {
  log_msg("WARNING: No variants found. Skipping tree generation.")
  quit(save = "no", status = 0)
}

# === GENERATE TREES ===========================================================

log_msg("\n--- Generating phylogenetic trees ---")

# Locate and read the database FASTA once (shared across all variants)
db_fasta_candidates <- c(
  file.path(blast_db_dir, "input_dedup.fa"),
  file.path(dataset_version_dir, "input_dedup.fa"),
  file.path(dataset_dir, "input_dedup.fa"),
  file.path(dataset_dir, "input.fa"),
  file.path(dataset_dir, "reference.fasta"),
  file.path(
    dataset_dir,
    "..",
    "References",
    "IA_IB_IIA_IIB_IIIA_IIIB_references.fa"
  )
)

db_seqs_raw <- NULL
db_fasta_file <- NULL
for (candidate in db_fasta_candidates) {
  if (file.exists(candidate)) {
    db_fasta_file <- candidate
    log_msg("✓ Found database FASTA: %s", basename(candidate))
    break
  }
}

if (!is.null(db_fasta_file)) {
  db_seqs_raw <- tryCatch(read_fasta_raw(db_fasta_file), error = function(e) {
    log_msg("✗ Error reading database FASTA: %s", conditionMessage(e))
    NULL
  })
  if (!is.null(db_seqs_raw)) {
    log_msg("✓ Loaded %d sequence(s) from database FASTA", length(db_seqs_raw))
  }
} else {
  log_msg("⚠ Database FASTA file not found in expected locations")
}

# Lookup a database sequence by ID (exact then normalized match)
lookup_db_seq <- function(id) {
  if (is.null(db_seqs_raw)) return(NULL)
  idx <- which(tolower(names(db_seqs_raw)) == tolower(id))
  if (length(idx) == 0) {
    idx <- which(norm_id(names(db_seqs_raw)) == norm_id(id))
  }
  if (length(idx) == 0) return(NULL)
  db_seqs_raw[[idx[1]]]
}


log_msg("MAFFT executable: %s", mafft_bin)
log_msg("IQ-TREE executable: %s", iqtree_bin)

if (
  is.na(mafft_bin) ||
  !nzchar(mafft_bin) ||
  !file.exists(mafft_bin)
) {
  log_msg("ERROR: MAFFT executable not found: %s", mafft_bin)
  quit(save = "no", status = 1)
}

if (file.access(mafft_bin, mode = 1) != 0) {
  log_msg("ERROR: MAFFT is not executable: %s", mafft_bin)
  quit(save = "no", status = 1)
}

if (
  is.na(iqtree_bin) ||
  !nzchar(iqtree_bin) ||
  !file.exists(iqtree_bin)
) {
  log_msg("ERROR: IQ-TREE executable not found: %s", iqtree_bin)
  quit(save = "no", status = 1)
}

if (file.access(iqtree_bin, mode = 1) != 0) {
  log_msg("ERROR: IQ-TREE is not executable: %s", iqtree_bin)
  quit(save = "no", status = 1)
}

log_msg("✓ MAFFT executable validated")
log_msg("✓ IQ-TREE executable validated")

n_trees_generated <- 0
unlink(file.path(outbreak_trees_dir, "tree_index.tsv"))
for (outbreak_variant in unique_variant) {
  log_msg("\nProcessing variant: %s", outbreak_variant)
  
  # Get batch sequences with this base variant
  batch_seqs_list <- batch_variants_clean %>%
    filter(group_key == outbreak_variant) %>%
    pull(seqName)
  
  n_batch <- length(batch_seqs_list)
  log_msg("  Batch sequences: %d", n_batch)
  batch_hit_variant <- setNames(batch_variants_clean$variant, batch_variants_clean$seqName)
  batch_label_variant <- batch_hit_variant
  if (nrow(proposed_names) > 0) {
    prop <- setNames(proposed_names$proposed_variant, proposed_names$qseqid)
    batch_label_variant <- ifelse(
      !is.na(prop[names(batch_hit_variant)]) & nzchar(prop[names(batch_hit_variant)]),
      prop[names(batch_hit_variant)],
      batch_hit_variant
    )
    names(batch_label_variant) <- names(batch_hit_variant)
  }

  # Tree/file name: group stems plus any NEW proposed variant names of the batch
  # sequences (names that do not exist in the database yet)
  new_names <- unique(unname(batch_label_variant[batch_seqs_list]))
  new_names <- new_names[!is.na(new_names) & nzchar(new_names) &
                           !(new_names %in% meta$variant)]
  tree_label <- paste(sort(unique(c(group_bases[[outbreak_variant]], new_names))), collapse = "_")
  log_msg("  Tree name: %s", tree_label)
  
  # Get database sequences with any base variant in this group
  # (V6, V6a, V6b all match NOR-2024-V6; linked bases are included too)
  db_seqs_for_variant <- meta %>%
    mutate(base_var = extract_base_variant(variant)) %>%
    filter(base_var %in% group_bases[[outbreak_variant]]) %>%
    select(id, variant, genotype)
  
  n_db <- nrow(db_seqs_for_variant)
  log_msg("  Database sequences: %d", n_db)
  
  # Skip early based on real batch+database counts, before loading the (much larger)
  # per-sequence-tree alignment files, which would otherwise inflate n_total artificially
  if (n_batch + n_db < 3) {
    log_msg("  ⚠ SKIP: Only %d real sequence(s) (batch=%d, db=%d), need ≥3 for phylogenetic analysis", n_batch + n_db, n_batch, n_db)
    next
  }
  
  # Collect all sequences as raw (ungapped) strings, to be aligned together with MAFFT
  all_outbreak_seqs <- list()
  
  # Add batch sequences (query only — row 1 of that sample's per-sequence alignment file)
  for (seq_name in batch_seqs_list) {
    seq_aln_file <- file.path(trees_dir, seq_name, "aligned_trimmed_blast_only.fa")
    if (!file.exists(seq_aln_file)) {
      log_msg("    ✗ Alignment file not found: %s", basename(seq_aln_file))
      next
    }
    
    tryCatch({
      seqs <- ape::read.dna(seq_aln_file, format = "fasta", as.character = TRUE)
      if (nrow(seqs) > 0) {
        query_seq <- paste(toupper(seqs[1, ]), collapse = "")
        query_seq <- gsub("-", "", query_seq)
        all_outbreak_seqs[[paste0(seq_name, "_BATCH")]] <- query_seq
        log_msg("    ✓ Loaded query sequence from %s", seq_name)
      }
    }, error = function(e) {
      log_msg("    ✗ Error reading %s: %s", seq_name, e$message)
    })
  }
  
  # Add database sequences matching this variant (from the shared, pre-loaded database FASTA)
  n_matched <- 0
  if (nrow(db_seqs_for_variant) > 0 && !is.null(db_seqs_raw)) {
    for (i in seq_len(nrow(db_seqs_for_variant))) {
      db_id <- db_seqs_for_variant$id[i]
      db_seq <- lookup_db_seq(db_id)
      if (!is.null(db_seq)) {
        all_outbreak_seqs[[db_id]] <- db_seq
        n_matched <- n_matched + 1
      }
    }
    log_msg("  ✓ Matched %d/%d database sequence(s) in FASTA", n_matched, nrow(db_seqs_for_variant))
  } else if (nrow(db_seqs_for_variant) > 0) {
    log_msg("  ⚠ Could not load database sequences (FASTA not available)")
  }
  
  n_total <- length(all_outbreak_seqs)
  log_msg("  Total sequences for tree: %d (batch + matched database)", n_total)
  
  # Add 1-2 outgroup sequences from OTHER variants for phylogenetic context
  other_variants <- unique_variant[unique_variant != outbreak_variant]
  if (length(other_variants) > 0 && !is.null(db_seqs_raw)) {
    n_outgroup <- min(2, length(other_variants))
    outgroup_variants <- sample(other_variants, n_outgroup)
    
    for (other_var in outgroup_variants) {
      other_id <- meta %>%
        mutate(base_var = extract_base_variant(variant)) %>%
        filter(base_var %in% group_bases[[other_var]]) %>%
        slice(1) %>%
        pull(id)
      
      if (length(other_id) > 0) {
        other_seq <- lookup_db_seq(other_id[1])
        if (!is.null(other_seq)) {
          outgroup_label <- sprintf("%s_OUTGROUP_%s", other_id[1], other_var)
          all_outbreak_seqs[[outgroup_label]] <- other_seq
          log_msg("    ✓ Added outgroup from %s: %s", other_var, other_id[1])
        }
      }
    }
  }
  
  n_total <- length(all_outbreak_seqs)
  log_msg("  Final tree composition: %d sequences (main variant + outgroups)", n_total)
  
  # Only generate trees with sufficient sequences for meaningful phylogeny
  if (n_total < 3) {
    log_msg("  ⚠ SKIP: Only %d sequence(s), need ≥3 for phylogenetic analysis", n_total)
    next
  }
  
  # Generate tree
  tryCatch({
    variant_slug <- gsub("[^A-Za-z0-9.-]", "_", tree_label)
    
    # Write raw (unaligned) sequences, then align them together with MAFFT
    raw_file_tmp <- file.path(outbreak_trees_dir, sprintf("raw_%s.fasta", variant_slug))
    cat_lines <- character()
    for (seq_id in names(all_outbreak_seqs)) {
      cat_lines <- c(cat_lines, paste0(">", seq_id))
      cat_lines <- c(cat_lines, all_outbreak_seqs[[seq_id]])
    }
    writeLines(cat_lines, raw_file_tmp)
    
    aln_file_tmp <- file.path(outbreak_trees_dir, sprintf("%s.fasta", variant_slug))
    mafft_stderr <- file.path(
  outbreak_trees_dir,
  sprintf(".mafft_stderr_%s.log", variant_slug)
)

mafft_exit <- system2(
  command = mafft_bin,
  args = c(
    "--auto",
    "--quiet",
    shQuote(raw_file_tmp)
  ),
  stdout = aln_file_tmp,
  stderr = mafft_stderr
)

if (
  mafft_exit != 0 ||
  !file.exists(aln_file_tmp) ||
  file.size(aln_file_tmp) == 0
) {
  log_msg("  ✗ MAFFT alignment failed (exit code %d)", mafft_exit)

  if (file.exists(mafft_stderr)) {
    mafft_output <- readLines(mafft_stderr, warn = FALSE)
    for (line in mafft_output) {
      log_msg("      %s", line)
    }
  }

  next
}

unlink(raw_file_tmp, force = TRUE)
unlink(mafft_stderr, force = TRUE)
    
    
    # Build tree with IQ-TREE
    tree_out_prefix <- sub("\\.fasta$", "", aln_file_tmp)
    
    # Capture both stdout and stderr to temporary files for debugging
    stdout_log <- file.path(outbreak_trees_dir, sprintf(".iqtree_stdout_%s.log", variant_slug))
    stderr_log <- file.path(outbreak_trees_dir, sprintf(".iqtree_stderr_%s.log", variant_slug))
    
   # IQ-TREE appends .treefile to the full input filename
tree_file <- paste0(aln_file_tmp, ".treefile")

stdout_log <- file.path(
  outbreak_trees_dir,
  sprintf(".iqtree_stdout_%s.log", variant_slug)
)

stderr_log <- file.path(
  outbreak_trees_dir,
  sprintf(".iqtree_stderr_%s.log", variant_slug)
)

iqtree_exit <- system2(
  command = iqtree_bin,
  args = c(
    "-s", aln_file_tmp,
    "-m", "JC",
    "-nt", "AUTO",
    "-fast",
    "-redo"
  ),
  stdout = stdout_log,
  stderr = stderr_log
)

log_msg("    IQ-TREE exit code: %d", iqtree_exit)
log_msg("    Looking for: %s", basename(tree_file))
log_msg("    Exists: %s", file.exists(tree_file))

if (iqtree_exit != 0 || !file.exists(tree_file)) {
  if (file.exists(stdout_log)) {
    log_msg("    IQ-TREE stdout:")
    iqtree_out <- readLines(stdout_log, warn = FALSE)

    for (line in iqtree_out) {
      log_msg("      %s", line)
    }
  }

  if (file.exists(stderr_log)) {
    log_msg("    IQ-TREE stderr:")
    iqtree_err <- readLines(stderr_log, warn = FALSE)

    for (line in iqtree_err) {
      log_msg("      %s", line)
    }
  }

  log_msg("  ✗ IQ-TREE failed or tree file was not created")
  next
}
  # Read the IQ-TREE result into R
outbreak_tree <- tryCatch(
  {
    ape::read.tree(tree_file)
  },
  error = function(e) {
    log_msg(
      "  ✗ Could not read IQ-TREE tree file %s: %s",
      basename(tree_file),
      conditionMessage(e)
    )
    NULL
  }
)

if (is.null(outbreak_tree)) {
  next
}

log_msg(
  "  ✓ Loaded phylogenetic tree with %d tips",
  length(outbreak_tree$tip.label)
)

# Validate tree structure before visualization
if (
  length(outbreak_tree$tip.label) == 0 ||
  any(is.na(outbreak_tree$tip.label)) ||
  any(!nzchar(outbreak_tree$tip.label))
) {
  log_msg("  ✗ Invalid tree structure: empty or missing tip labels")
  next
}  
    # Validate tree structure before visualization
    if (length(outbreak_tree$tip.label) == 0 || any(is.na(outbreak_tree$tip.label))) {
      log_msg("  ✗ Invalid tree structure (empty or missing tip labels)")
      unlink(paste0(tree_out_prefix, ".*"), force = TRUE)
      next
    }
    
    # Prepare rich metadata for tree tips: ID, variant, genotype, source
    # Handle batch sequences, database sequences for this variant, and outgroups
    outbreak_tip_info <- tibble(label = outbreak_tree$tip.label) %>%
      mutate(
        is_batch = grepl("_BATCH$", label, ignore.case = TRUE),
        is_outgroup = grepl("_OUTGROUP_", label, ignore.case = TRUE),
        # Extract base ID (remove _BATCH or _OUTGROUP_* suffixes for matching)
        base_id = case_when(
          is_batch ~ sub("_BATCH$", "", label),
          is_outgroup ~ sub("_OUTGROUP_.*$", "", label),
          TRUE ~ label
        ),
        source = case_when(
          is_batch ~ "Batch (red)",
          is_outgroup ~ "Outgroup (gray)",
          TRUE ~ "Database"
        )
      ) %>%
      left_join(
        db_seqs_for_variant %>%
          mutate(genotype = replace_na(genotype, "Unknown")),
        by = c("base_id" = "id")
      ) %>%
      mutate(
        # Display label: ID [Variant|Genotype]
        display_label = case_when(
          is_batch ~ sprintf("%s [%s]", base_id, coalesce(unname(batch_label_variant[base_id]), "")),
          is_outgroup ~ sprintf("%s [OUTGROUP]", base_id),
          !is.na(variant) ~ sprintf("%s [%s %s]", base_id, variant, replace_na(genotype, "Unknown")),
          TRUE ~ sprintf("%s [%s]", base_id, replace_na(genotype, "Unknown"))
        )
      )
    
    # Plot tree with enhanced visualization
    p <- NULL
    tryCatch({
      p <- ggtree(outbreak_tree, layout = "rectangular", branch.length = "branch.length") %<+% outbreak_tip_info +
        geom_tippoint(aes(color = source), size = 3) +
        geom_tiplab(aes(label = display_label), size = 2.2, hjust = -0.05) +
        scale_color_manual(
          values = c(
            "Batch (red)" = "#e74c3c",
            "Database" = "#3498db",
            "Outgroup (gray)" = "#95a5a6"
          ),
          name = "Source"
        ) +
        theme_tree2() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(size = 8)
        ) +
        hexpand(.8) +
        labs(
          title = sprintf("Variant group %s — Phylogenetic outbreak tree",
                          gsub("_", " / ", tree_label, fixed = TRUE)),
          subtitle = sprintf("Batch (n=%d) vs database (n=%d) + outgroups, total=%d",
                            n_batch, nrow(db_seqs_for_variant), n_total)
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
      png_file <- file.path(outbreak_trees_dir, sprintf("%s.png", tree_label))
      ggsave(png_file, plot = p, width = 12, height = 8, dpi = 150)
      # Index so the report can map tree file -> variant group
      write_tsv(tibble(tree_label = tree_label, group_key = outbreak_variant),
                file.path(outbreak_trees_dir, "tree_index.tsv"), append = TRUE)
      log_msg("  ✓ Saved tree PNG: %s", basename(png_file))
      n_trees_generated <- n_trees_generated + 1
    }, error = function(e) {
      log_msg("  ✗ Error saving PNG: %s", e$message)
    })
    
    # Preserve the MAFFT alignment, but remove IQ-TREE intermediate files
iqtree_prefix <- aln_file_tmp

unlink(paste0(iqtree_prefix, ".bionj"), force = TRUE)
unlink(paste0(iqtree_prefix, ".ckp.gz"), force = TRUE)
unlink(paste0(iqtree_prefix, ".iqtree"), force = TRUE)
unlink(paste0(iqtree_prefix, ".log"), force = TRUE)
unlink(paste0(iqtree_prefix, ".mldist"), force = TRUE)
#unlink(paste0(iqtree_prefix, ".treefile"), force = TRUE)
unlink(paste0(iqtree_prefix, ".uniqueseq.phy"), force = TRUE)

unlink(stdout_log, force = TRUE)
unlink(stderr_log, force = TRUE)

log_msg("  ✓ Alignment preserved: %s", basename(aln_file_tmp))
    
  }, error = function(e) {
    log_msg("  ✗ Error: %s", e$message)
  })
}

# === COMPLETION ==============================================================

log_msg("\n=== Completion ===")
log_msg("Generated %d outbreak trees", n_trees_generated)
log_msg("Saved to: %s", outbreak_trees_dir)
log_msg("Log file: %s", log_file)
