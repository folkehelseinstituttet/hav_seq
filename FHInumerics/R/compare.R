# Compare query sequences to a reference database:
# search -> closest neighbours -> MAFFT alignment -> trim to the query region ->
# SNP matrix (ape) -> neighbour-joining tree (ape, midpoint rooted).
#
# Used by both the Shiny app (app.R) and the command line (cli/fhinumerics),
# so both give identical results. Everything runs on this computer, inside R or
# the local Python (see R/python.R); the method of each step is set in config.yaml.

DEFAULT_SETTINGS <- list(
  n_neighbours = 20,      # closest database sequences kept per query
  min_coverage = 90,      # % of the query a hit must cover
  search = "decipher",    # decipher (R) | blast (rBLAST; needs BLAST+ installed)
  aligner = "decipher",   # decipher (R) | mafft (pymafft, Python)
  bootstrap = TRUE,       # 100 bootstraps for the tree
  threads = 4             # search and alignment
)

ALIGNERS <- c("MAFFT (--auto, pymafft)" = "mafft", "DECIPHER AlignSeqs (R only)" = "decipher")

run_compare <- function(queries, db, cfg, settings = list(), out_dir = NULL,
                        progress = function(step, detail = NULL) NULL, log = message) {
  settings <- utils::modifyList(DEFAULT_SETTINGS, settings)
  if (length(queries) > 50) stop("At most 50 query sequences per run.", call. = FALSE)
  if (!settings$aligner %in% ALIGNERS) stop("Unknown aligner '", settings$aligner, "'.", call. = FALSE)
  started <- Sys.time()

  if (is.null(out_dir)) {
    label <- gsub("[^A-Za-z0-9._-]+", "_", names(queries)[1])
    out_dir <- file.path(cfg$output_dir, paste0(format(started, "%Y%m%d-%H%M%S"), "_", label))
  }
  work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  log("Output folder: ", out_dir)
  write_fasta(queries, file.path(out_dir, "query.fasta"))

  progress("Preparing database")
  ref <- prepare_db(db, cfg, settings$search, log = log)
  meta <- read_metadata(db)

  # Internal IDs (Q1, Q2, ...) keep tool input free of odd characters.
  q_ids <- data.frame(internal = paste0("Q", seq_along(queries)), id = names(queries))
  q_internal <- queries
  names(q_internal) <- q_ids$internal
  write_fasta(q_internal, file.path(work, "query_internal.fasta"))

  # ---- 1. Search ----------------------------------------------------------------
  progress(paste0("Searching database (", settings$search, ")"))
  hits <- search_db(q_internal, ref, settings$search, settings$threads)
  utils::write.table(hits, file.path(work, "search_raw.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  hits$query <- q_ids$id[match(hits$qseqid, q_ids$internal)]
  hits$subject <- ref$ids$id[match(hits$sseqid, ref$ids$internal)]

  # ---- 2. Closest neighbours ----------------------------------------------------
  # Ties are broken by database order, so the same input always gives the same neighbours.
  db_order <- as.integer(sub("^D", "", hits$sseqid))
  hits <- hits[order(hits$qseqid, -hits$pident, -hits$score, db_order), ]
  hits$passes_coverage <- hits$qcovhsp >= settings$min_coverage
  kept <- do.call(rbind, lapply(split(hits[hits$passes_coverage, ], hits$qseqid[hits$passes_coverage]),
                                function(h) utils::head(h, settings$n_neighbours)))
  if (is.null(kept) || nrow(kept) == 0) {
    stop("No database sequence covers at least ", settings$min_coverage,
         "% of any query. Check the sequence, or lower the coverage setting.", call. = FALSE)
  }
  hits$selected <- paste(hits$qseqid, hits$sseqid) %in% paste(kept$qseqid, kept$sseqid)
  no_hits <- setdiff(q_ids$internal, kept$qseqid)
  if (length(no_hits)) {
    log("WARNING: no neighbours found for ", paste(q_ids$id[q_ids$internal %in% no_hits], collapse = ", "),
        " (no hit covers >= ", settings$min_coverage, "% of it)")
  }
  neighbour_ids <- unique(kept$sseqid)
  log("Neighbours selected: ", length(neighbour_ids), " database sequences")

  # How many database sequences are as close as the best hit? With many identical
  # sequences in a database, more can tie than fit in the tree.
  passing <- hits[hits$passes_coverage, ]
  ties <- do.call(rbind, lapply(split(passing, passing$qseqid), function(h) {
    best <- h[h$pident == h$pident[1] & h$length == h$length[1], ]
    data.frame(query = q_ids$id[match(h$qseqid[1], q_ids$internal)], best_pident = h$pident[1],
               best_mismatches = h$mismatch[1], n_equally_close = nrow(best),
               n_shown = sum(best$sseqid %in% kept$sseqid))
  }))
  for (i in seq_len(NROW(ties))) {
    if (ties$n_equally_close[i] > ties$n_shown[i]) {
      log("NOTE: ", ties$n_equally_close[i], " database sequences are equally close to ", ties$query[i],
          " (", ties$best_pident[i], "% identity); ", ties$n_shown[i], " are shown. See search_hits.tsv for all.")
    }
  }

  # Queries that match the database on the minus strand are reverse-complemented.
  best_strand <- tapply(kept$sstrand, kept$qseqid, `[`, 1)
  flipped <- names(best_strand)[best_strand == "minus"]
  for (q in flipped) q_internal[[q]] <- reverse_complement(q_internal[[q]])
  if (length(flipped)) log("Reverse-complemented: ", paste(q_ids$id[match(flipped, q_ids$internal)], collapse = ", "))

  # ---- 3. Alignment -------------------------------------------------------------
  progress(paste0("Aligning sequences (", names(ALIGNERS)[ALIGNERS == settings$aligner], ")"),
           paste(length(q_internal) + length(neighbour_ids), "sequences"))
  combined <- c(q_internal, ref$seqs[match(neighbour_ids, ref$ids$internal)])
  names(combined) <- c(names(q_internal), neighbour_ids)
  write_fasta(combined, file.path(work, "combined.fasta"))
  if (settings$aligner == "mafft") {
    run_python(cfg, c("align", file.path(work, "combined.fasta"), file.path(work, "aligned.fasta"),
                      "--strategy", "auto", "--threads", settings$threads), log = log)
  } else {
    aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(combined), processors = settings$threads, verbose = FALSE)
    write_fasta(stats::setNames(as.character(aln), names(aln)), file.path(work, "aligned.fasta"))
  }
  aln <- read_alignment(file.path(work, "aligned.fasta"))

  # Trim to the region every query covers, so all comparisons use the same positions.
  trimmed <- trim_to_queries(aln, names(q_internal))
  if (!is.null(trimmed$dropped)) {
    log("WARNING: dropped (covers < 50% of the compared region): ",
        paste(ref$ids$id[match(trimmed$dropped, ref$ids$internal)], collapse = ", "))
  }
  aln_t <- trimmed$aln
  log("Compared region: ", nchar(aln_t[1]), " alignment columns")
  aln_file <- file.path(work, "aligned_trimmed.fasta")
  write_fasta(aln_t, aln_file)

  # ---- 4. SNP distances ---------------------------------------------------------
  # ape::dist.dna(model = "N", pairwise.deletion = TRUE) counts differences at positions
  # where both sequences have A/C/G/T: the same numbers as snp-dists.
  progress("Counting SNP differences (ape)")
  dna <- ape::read.dna(aln_file, format = "fasta")
  snp <- as.matrix(ape::dist.dna(dna, model = "N", pairwise.deletion = TRUE))
  storage.mode(snp) <- "integer"
  sites <- comparable_sites(aln_t)

  # ---- 5. Tree ------------------------------------------------------------------
  tree <- NULL
  if (length(aln_t) >= 3) {
    boot <- settings$bootstrap && length(aln_t) >= 4
    progress("Building tree (neighbour-joining)", if (boot) "with 100 bootstraps" else NULL)
    tree <- nj_tree(dna, bootstrap = boot)
    ape::write.tree(tree, file.path(work, "nj.nwk"))
    tree <- phangorn::midpoint(tree, node.labels = "support")
  } else {
    log("WARNING: fewer than 3 sequences; no tree built")
  }

  # ---- 6. Results with real names -------------------------------------------------
  progress("Writing results")
  all_ids <- rbind(q_ids, ref$ids[ref$ids$internal %in% names(aln_t), ])
  real <- function(x) all_ids$id[match(x, all_ids$internal)]

  dimnames(snp) <- list(real(rownames(snp)), real(colnames(snp)))
  dimnames(sites) <- list(real(rownames(sites)), real(colnames(sites)))
  if (!is.null(tree)) tree$tip.label <- real(tree$tip.label)

  tips <- data.frame(id = real(names(aln_t)), is_query = names(aln_t) %in% q_ids$internal)
  tip_meta <- merge(tips, meta, by = "id", all.x = TRUE, sort = FALSE)

  hits_out <- merge(hits[, c("query", "subject", "pident", "length", "mismatch", "gapopen",
                             "qcovhsp", "qlen", "slen", "score", "sstrand", "selected")],
                    meta, by.x = "subject", by.y = "id", all.x = TRUE, sort = FALSE)
  hits_out <- hits_out[order(hits_out$query, -hits_out$pident, -hits_out$score),
                       c("query", setdiff(names(hits_out), "query"))]

  utils::write.table(hits_out, file.path(out_dir, "search_hits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(cbind(id = rownames(snp), as.data.frame(snp, check.names = FALSE)),
                     file.path(out_dir, "snp_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(cbind(id = rownames(sites), as.data.frame(sites, check.names = FALSE)),
                     file.path(out_dir, "compared_sites.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(tip_meta, file.path(out_dir, "tree_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  aln_real <- aln_t
  names(aln_real) <- real(names(aln_t))
  write_fasta(aln_real, file.path(out_dir, "alignment.fasta"))
  if (!is.null(tree)) ape::write.tree(tree, file.path(out_dir, "tree.nwk"))

  manifest <- list(
    created = format(started, "%Y-%m-%d %H:%M:%S"),
    user = Sys.info()[["user"]],
    computer = Sys.info()[["nodename"]],
    app_version = app_version(),
    database = list(name = db$name, fasta = db$fasta, metadata = db$metadata,
                    fasta_md5 = ref$checksum, n_sequences = length(ref$seqs)),
    queries = as.list(names(queries)),
    reverse_complemented = as.list(q_ids$id[match(flipped, q_ids$internal)]),
    settings = settings,
    neighbours = length(neighbour_ids),
    compared_columns = nchar(aln_t[1]),
    tools = tool_versions(cfg, settings),
    method = method_text(settings)
  )
  yaml::write_yaml(manifest, file.path(out_dir, "run_manifest.yaml"))
  log("Finished in ", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), " s")

  utils::write.table(ties, file.path(out_dir, "closest_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(out_dir = out_dir, hits = hits_out, ties = ties, snp = snp, sites = sites, tree = tree,
       tip_meta = tip_meta, queries = names(queries), manifest = manifest,
       display_columns = display_columns(db, meta))
}

# Neighbour-joining tree (ape): BIONJ on TN93 distances (pairwise deletion), with
# bootstrap support (100 replicates, % of replicates) as node labels.
nj_tree <- function(dna, bootstrap = TRUE) {
  build <- function(x) ape::bionj(ape::dist.dna(x, model = "TN93", pairwise.deletion = TRUE))
  tree <- build(dna)
  if (bootstrap) {
    set.seed(1)
    tree$node.label <- ape::boot.phylo(tree, dna, build, B = 100, quiet = TRUE)
  }
  tree$edge.length[tree$edge.length < 0] <- 0
  tree
}

# Keep the columns between the latest query start and the earliest query end;
# then drop all-gap columns and sequences covering < 50% of that region.
trim_to_queries <- function(aln, query_ids) {
  chars <- do.call(rbind, strsplit(aln, ""))
  rownames(chars) <- names(aln)
  spans <- t(vapply(query_ids, function(q) range(which(chars[q, ] != "-")), numeric(2)))
  from <- max(spans[, 1])
  to <- min(spans[, 2])
  if (from >= to) {
    stop("The query sequences cover different parts of the genome and do not overlap. ",
         "Run them separately.", call. = FALSE)
  }
  chars <- chars[, from:to, drop = FALSE]
  chars <- chars[, colSums(chars != "-") > 0, drop = FALSE]
  acgt <- matrix(chars %in% c("A", "C", "G", "T"), nrow = nrow(chars))
  covered <- rowMeans(acgt)
  drop <- setdiff(rownames(chars)[covered < 0.5], query_ids)
  chars <- chars[!rownames(chars) %in% drop, , drop = FALSE]
  out <- apply(chars, 1, paste, collapse = "")
  list(aln = out, dropped = if (length(drop)) drop else NULL)
}

# Number of positions where both sequences have A/C/G/T (the denominator of the SNP counts).
comparable_sites <- function(aln) {
  ok <- do.call(rbind, strsplit(aln, "")) %in% c("A", "C", "G", "T")
  ok <- matrix(ok, nrow = length(aln))
  sites <- tcrossprod(ok * 1)
  dimnames(sites) <- list(names(aln), names(aln))
  sites
}

method_text <- function(settings) {
  paste(
    switch(settings$search,
      decipher = "Search: DECIPHER IndexSeqs (k = 11) + SearchIndex, both strands, up to 500 hits per query, each aligned with AlignPairs;",
      blast = "Search: blastn via rBLAST (max_target_seqs 500, max_hsps 1);"),
    "per query, hits covering >= min_coverage % of the query, ranked by percent identity then score",
    "(ties by database order), top n_neighbours kept.",
    switch(settings$aligner,
      mafft = "Alignment: MAFFT --auto (pymafft).",
      decipher = "Alignment: DECIPHER AlignSeqs."),
    "Alignment trimmed to the columns covered by every query.",
    "SNPs: ape::dist.dna(model = 'N', pairwise.deletion = TRUE), i.e. differences at positions where",
    "both sequences have A/C/G/T (as snp-dists).",
    "Tree: ape::bionj on TN93 distances (pairwise deletion), negative branch lengths set to 0,",
    "+ 100 bootstraps (ape::boot.phylo, seed 1) when bootstrap is on and >= 4 sequences;",
    "midpoint rooted (phangorn::midpoint)."
  )
}

tool_versions <- function(cfg, settings) {
  pkg <- function(p) as.character(utils::packageVersion(p))
  out <- list(R = R.version.string, ape = pkg("ape"), phangorn = pkg("phangorn"),
              DECIPHER = pkg("DECIPHER"), Biostrings = pkg("Biostrings"))
  if (settings$search == "blast") {
    out$rBLAST <- pkg("rBLAST")
    out$blastn <- trimws(system2(Sys.which("blastn"), "-version", stdout = TRUE, stderr = TRUE)[1])
  }
  if (settings$aligner == "mafft") out$python <- python_versions(cfg)
  out
}

app_version <- function() {
  root <- getOption("fhinumerics.root", ".")
  v <- tryCatch(readLines(file.path(root, "VERSION"), warn = FALSE)[1], error = function(e) "unknown")
  git <- suppressWarnings(tryCatch(
    system2("git", c("-C", root, "describe", "--always", "--dirty", "--tags"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)))
  if (length(git)) paste0(v, " (", git[1], ")") else v
}
