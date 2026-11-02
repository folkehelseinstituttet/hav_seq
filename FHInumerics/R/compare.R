# Compare query sequences to a reference database:
# search -> closest neighbours -> MAFFT alignment -> trim to the query region ->
# SNP matrix (ape) -> neighbour-joining tree (ape, midpoint rooted).
# With several queries, the steps after the search run for all queries together and
# for each query alone (its own neighbours, compared over its own length).
# Sets (decision 0013) can be the database searched, or be compared with each other
# without a search (compare_within()); the steps from trimming on are the same
# (compare_alignment()).
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
  threads = 4,            # search and alignment
  exclude_self = FALSE    # drop hits to the database entry with the query's own name; or a
                          # named list: query name -> the database ID to leave out (sets)
)

# Name of the result for all queries together (a FASTA name cannot contain spaces).
ALL_QUERIES <- "All queries together"

ALIGNERS <- c("MAFFT (--auto, pymafft)" = "mafft", "DECIPHER AlignSeqs (R only)" = "decipher")

# Three modes (decision 0013), recorded as `mode` in run_manifest.yaml:
#   database      queries searched against a database from config.yaml
#   set_database  queries searched against sets: `db` from set_db() (R/sets.R)
#   within_set    no queries and no search: the sequences of sets compared with each other,
#                 `within` from set_within() (R/sets.R); `queries` and `db` are not used
# `input`: where the queries came from (a set, R/sets.R set_queries()), for the manifest.
run_compare <- function(queries, db, cfg, settings = list(), out_dir = NULL,
                        progress = function(step, detail = NULL) NULL, log = message, input = NULL,
                        within = NULL) {
  settings <- utils::modifyList(DEFAULT_SETTINGS, settings)
  if (!is.null(within)) return(compare_within(within, cfg, settings, out_dir, progress, log))
  mode <- if (is.null(db$from_sets)) "database" else "set_database"
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
  if (mode == "set_database") {
    log("Searching only ", db$name, " (", length(ref$seqs), " sequences)")
    turned <- db$from_sets$turned
    if (length(turned)) log("Turned (reverse complement) to the orientation of the other sequences: ", paste(turned, collapse = ", "))
  }

  # Internal IDs (Q1, Q2, ...) keep tool input free of odd characters. A query with the
  # name of a database sequence gets "_query" added, so the two do not share a name in
  # the tree and SNP matrix.
  q_ids <- data.frame(internal = paste0("Q", seq_along(queries)), id = names(queries), given = names(queries))
  clash <- q_ids$id %in% ref$ids$id
  q_ids$id[clash] <- paste0(q_ids$id[clash], "_query")
  if (any(clash)) {
    log("Renamed (a database sequence has the same name): ",
        paste(q_ids$given[clash], "->", q_ids$id[clash], collapse = ", "))
  }
  q_internal <- queries
  names(q_internal) <- q_ids$internal
  write_fasta(q_internal, file.path(work, "query_internal.fasta"))

  # ---- 1. Search ----------------------------------------------------------------
  progress(paste0("Searching database (", settings$search, ")"))
  hits <- search_db(q_internal, ref, settings$search, settings$threads)
  utils::write.table(hits, file.path(work, "search_raw.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  hits$query <- q_ids$id[match(hits$qseqid, q_ids$internal)]
  hits$subject <- ref$ids$id[match(hits$sseqid, ref$ids$internal)]
  # Queries taken from the database itself (Database viewer, sets) would find their own
  # entry: with exclude_self TRUE each query's own name, or for the queries named in the
  # list, the database ID given there.
  own <- if (isTRUE(settings$exclude_self)) q_ids$given else if (is.list(settings$exclude_self)) {
    unname(unlist(settings$exclude_self)[q_ids$given])
  } else rep(NA_character_, nrow(q_ids))
  self <- own[match(hits$qseqid, q_ids$internal)] == hits$subject
  self[is.na(self)] <- FALSE
  if (any(self)) log("Left out each query's own database entry: ", paste(unique(hits$subject[self]), collapse = ", "))
  hits <- hits[!self, ]

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

  # ---- 3-5. Alignment, SNPs and tree, for all queries together and for each alone ----
  # With one query there is one result, in the run folder. With several, the result for
  # all queries together is in the run folder and each query's own in per_query/.
  sets <- as.list(q_ids$internal)
  names(sets) <- q_ids$id
  folders <- out_dir
  if (length(queries) > 1) {
    sets <- c(stats::setNames(list(q_ids$internal), ALL_QUERIES), sets)
    folders <- c(out_dir, file.path(out_dir, "per_query", sprintf("%02d_%s", seq_along(queries),
                                                                  gsub("[^A-Za-z0-9._-]+", "_", q_ids$id))))
  }
  run_set <- function(i) {
    say <- if (length(sets) > 1) function(...) log(names(sets)[i], ": ", ...) else log
    compare_set(sets[[i]], names(sets)[i], q_internal, q_ids, kept, ref, meta, cfg, settings, folders[i],
                progress = function(step) progress(step, if (length(sets) > 1) names(sets)[i]), log = say)
  }
  views <- if (length(sets) == 1) list(run_set(1)) else lapply(seq_along(sets), function(i) {
    tryCatch(run_set(i), error = function(e) {
      log("WARNING: no result for ", names(sets)[i], ": ", conditionMessage(e))
      list(label = names(sets)[i], error = conditionMessage(e))
    })
  })
  names(views) <- names(sets)
  failed <- vapply(views, function(v) !is.null(v$error), TRUE)
  if (all(failed)) stop(views[[1]]$error, call. = FALSE)

  # ---- 6. Search results and manifest ----------------------------------------------
  progress("Writing results")
  hits_out <- merge(hits[, c("query", "subject", "pident", "length", "mismatch", "gapopen",
                             "qcovhsp", "qlen", "slen", "score", "sstrand", "selected")],
                    meta, by.x = "subject", by.y = "id", all.x = TRUE, sort = FALSE)
  hits_out <- hits_out[order(hits_out$query, -hits_out$pident, -hits_out$score),
                       c("query", setdiff(names(hits_out), "query"))]
  utils::write.table(hits_out, file.path(out_dir, "search_hits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  manifest <- list(
    created = format(started, "%Y-%m-%d %H:%M:%S"),
    user = Sys.info()[["user"]],
    computer = Sys.info()[["nodename"]],
    app_version = app_version(),
    mode = mode,
    database = c(list(name = db$name, fasta = db$fasta, metadata = db$metadata,
                      fasta_md5 = ref$checksum, n_sequences = length(ref$seqs)), db$from_sets$manifest),
    queries = as.list(names(queries)),
    input = input,
    renamed = as.list(stats::setNames(q_ids$id[clash], q_ids$given[clash])),
    reverse_complemented = as.list(q_ids$id[match(flipped, q_ids$internal)]),
    settings = settings,
    results = lapply(seq_along(views), function(i) {
      v <- views[[i]]
      c(list(queries = names(sets)[i], folder = if (i == 1) "." else file.path("per_query", basename(folders[i]))),
        if (is.null(v$error)) list(neighbours = v$neighbours, compared_columns = v$compared_columns,
                                   compared_region = v$compared_region, dropped = as.list(v$dropped))
        else list(error = v$error))
    }),
    tools = tool_versions(cfg, settings),
    method = method_text(settings, mode)
  )
  manifest <- Filter(Negate(is.null), manifest)   # no `input` unless the queries came from a set
  yaml::write_yaml(manifest, file.path(out_dir, "run_manifest.yaml"))
  # Sets searched: a copy of them (as given, before turning) in the run folder.
  if (mode == "set_database") set_write(db$from_sets$set, file.path(out_dir, "set"))
  log("Finished in ", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), " s")

  utils::write.table(ties, file.path(out_dir, "closest_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(out_dir = out_dir, hits = hits_out, ties = ties, queries = q_ids$id, views = views, manifest = manifest,
       display_columns = display_columns(db, meta), mode = mode, set = db$from_sets$set)
}

# Alignment, SNP matrix and tree for one set of queries (internal IDs `q`) and their
# kept neighbours; writes them to `out_dir`. Returns the results with real names.
compare_set <- function(q, label, q_internal, q_ids, kept, ref, meta, cfg, settings, out_dir,
                        progress = function(step) NULL, log = message) {
  work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  neighbour_ids <- unique(kept$sseqid[kept$qseqid %in% q])
  if (!length(neighbour_ids)) {
    stop("No database sequence covers at least ", settings$min_coverage, "% of ", label, ".", call. = FALSE)
  }
  log("Neighbours selected: ", length(neighbour_ids), " database sequences")

  # ---- 3. Alignment -------------------------------------------------------------
  progress(paste0("Aligning sequences (", names(ALIGNERS)[ALIGNERS == settings$aligner], ", ",
                  length(q) + length(neighbour_ids), " sequences)"))
  combined <- c(q_internal[q], ref$seqs[match(neighbour_ids, ref$ids$internal)])
  names(combined) <- c(q, neighbour_ids)
  write_fasta(combined, file.path(work, "combined.fasta"))
  if (settings$aligner == "mafft") {
    run_python(cfg, c("align", file.path(work, "combined.fasta"), file.path(work, "aligned.fasta"),
                      "--strategy", "auto", "--threads", settings$threads), log = log)
  } else {
    aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(combined), processors = settings$threads, verbose = FALSE)
    write_fasta(stats::setNames(as.character(aln), names(aln)), file.path(work, "aligned.fasta"))
  }
  aln <- read_alignment(file.path(work, "aligned.fasta"))
  all_ids <- rbind(q_ids[, c("internal", "id")], ref$ids[ref$ids$internal %in% names(aln), c("internal", "id")])
  c(compare_alignment(aln, q, all_ids, meta, settings, out_dir, trim_to_queries(aln, q), progress, log),
    list(label = label, neighbours = length(neighbour_ids)))
}

# Trimming, SNP matrix and tree of one alignment (internal names; `ids` maps them to real
# names), and its result files in `out_dir`. `q`: the queries (red; each gets its own base
# numbers for the hover box); `number`: whose base numbers are kept (default the queries).
# `trimmed`: from trim_to_queries() or trim_to_region().
compare_alignment <- function(aln, q, ids, meta, settings, out_dir, trimmed, progress = function(step) NULL,
                              log = message, number = q) {
  work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  real <- function(x) ids$id[match(x, ids$internal)]
  if (!is.null(trimmed$dropped)) {
    log("WARNING: dropped (covers < 50% of the compared region): ", paste(real(trimmed$dropped), collapse = ", "))
  }
  aln_t <- trimmed$aln
  log("Compared region: ", nchar(aln_t[1]), " alignment columns")
  aln_file <- file.path(work, "aligned_trimmed.fasta")
  write_fasta(aln_t, aln_file)

  # Each query's own base number for every column left after trimming (for the hover
  # text that lists differing positions).
  positions <- lapply(number, function(id) cumsum(strsplit(aln[[id]], "")[[1]] != "-")[trimmed$cols])

  # ---- 4. SNP distances ---------------------------------------------------------
  # Ambiguity codes count only when they cannot match; N and gaps never (R/snps.R).
  progress("Counting SNP differences (phangorn)")
  snp <- snp_matrix(aln_t)
  sites <- comparable_sites(aln_t)
  dna <- ape::read.dna(aln_file, format = "fasta")

  # ---- 5. Tree ------------------------------------------------------------------
  tree <- NULL
  if (length(aln_t) >= 3) {
    boot <- settings$bootstrap && length(aln_t) >= 4
    progress(paste0("Building tree (neighbour-joining", if (boot) ", 100 bootstraps", ")"))
    tree <- nj_tree(dna, bootstrap = boot)
    ape::write.tree(tree, file.path(work, "nj.nwk"))
    tree <- phangorn::midpoint(tree, node.labels = "support")
  } else {
    log("WARNING: fewer than 3 sequences; no tree built")
  }

  # ---- Results with real names --------------------------------------------------
  dimnames(snp) <- list(real(rownames(snp)), real(colnames(snp)))
  dimnames(sites) <- list(real(rownames(sites)), real(colnames(sites)))
  if (!is.null(tree)) tree$tip.label <- real(tree$tip.label)

  tips <- data.frame(id = real(names(aln_t)), is_query = names(aln_t) %in% q)
  tip_meta <- merge(tips, meta, by = "id", all.x = TRUE, sort = FALSE)

  utils::write.table(cbind(id = rownames(snp), as.data.frame(snp, check.names = FALSE)),
                     file.path(out_dir, "snp_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(cbind(id = rownames(sites), as.data.frame(sites, check.names = FALSE)),
                     file.path(out_dir, "compared_sites.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(tip_meta, file.path(out_dir, "tree_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  aln_real <- aln_t
  names(aln_real) <- real(names(aln_t))
  write_fasta(aln_real, file.path(out_dir, "alignment.fasta"))
  if (!is.null(tree)) ape::write.tree(tree, file.path(out_dir, "tree.nwk"))
  # The whole alignment (before trimming) for the Alignment page: compared_cols are its
  # columns that were compared, `dropped` the sequences left out by trimming.
  aln_full <- aln
  names(aln_full) <- real(names(aln))
  write_fasta(aln_full, file.path(out_dir, "alignment_untrimmed.fasta"))

  names(positions) <- real(number)
  list(queries = real(q), snp = snp, sites = sites, tree = tree, tip_meta = tip_meta,
       aln = aln_real, positions = positions, aln_full = aln_full, compared_cols = trimmed$cols,
       aln_meta = merge(data.frame(id = names(aln_full), is_query = names(aln) %in% q), meta,
                        by = "id", all.x = TRUE, sort = FALSE),
       compared_region = sprintf("%d-%d", min(trimmed$cols), max(trimmed$cols)),
       dropped = if (is.null(trimmed$dropped)) character(0) else real(trimmed$dropped),
       limited_by = real(trimmed$limited_by), compared_columns = nchar(aln_t[1]), out_dir = out_dir)
}

# ---- Compare the sequences of a set with each other (decision 0013 (b)) ----------------
SET_COMPARE_MAX <- 300   # sequences per comparison within sets

# No search: the set's alignment (aligned first, as "Align the set", when it is missing or
# out of date), trimmed to the columns every sequence covers or to settings$region, then
# SNPs and tree as in runs. `w` is from set_within() (R/sets.R). The run folder gets a copy
# of the set compared (set/, a set folder with the alignment used).
compare_within <- function(w, cfg, settings, out_dir = NULL, progress = function(step, detail = NULL) NULL,
                           log = message) {
  started <- Sys.time()
  s <- w$set
  keys <- s$members$key
  n <- length(keys)
  if (n < 2) stop("At least 2 sequences are needed to compare a set.", call. = FALSE)
  if (n > SET_COMPARE_MAX) {
    stop(sprintf("%s has %d sequences; at most %d can be compared with each other. On the Sets page, click up to %d rows and press Compare the chosen with each other.",
                 w$what %||% "The set", n, SET_COMPARE_MAX, SET_COMPARE_MAX), call. = FALSE)
  }
  if (is.null(out_dir)) {
    out_dir <- file.path(cfg$output_dir, paste0(format(started, "%Y%m%d-%H%M%S"), "_", gsub("[^A-Za-z0-9._-]+", "_", s$name)))
  }
  dir.create(file.path(out_dir, "work"), recursive = TRUE, showWarnings = FALSE)
  log("Output folder: ", out_dir)
  log("Comparing the ", n, " sequences of ", s$name, " with each other (no database search)")
  names_in <- set_names(s$members)
  ids <- data.frame(internal = paste0("S", seq_len(n)), id = names_in)
  red <- ids$internal[s$members$source %in% c("query", "file")]

  progress("Preparing the alignment")
  al <- set_alignment(s)
  aligned_now <- is.null(al) || length(al$missing) > 0
  if (aligned_now) {
    log("The set has no alignment of all these sequences: aligning them now, as Align the set does")
    progress(sprintf("Aligning the sequences (DECIPHER, %d sequences)", n))
    s$alignment <- set_align(stats::setNames(unname(s$seqs[keys]), keys), settings$threads)[c("aln", "turned", "made", "from", "method", "versions")]
    s$history <- c(s$history, sprintf("%s aligned (%d sequences) for a comparison within the set", s$alignment$made, n))
    al <- set_alignment(s)
  } else {
    log("Using the set's alignment (", al$from, ", ", al$made, ")")
  }
  if (length(al$turned)) log("Turned (reverse complement) for the alignment: ", paste(al$turned, collapse = ", "))
  aln <- stats::setNames(unname(al$aln[names_in]), ids$internal)
  write_fasta(aln, file.path(out_dir, "work", "aligned.fasta"))

  region <- w$region %||% settings$region
  trimmed <- if (is.null(region)) {
    trim_to_queries(aln, ids$internal, no_overlap = "These sequences do not overlap; pick a region in the Alignment page.")
  } else {
    r <- parse_region(region, nchar(aln[[1]]))
    log("Region: alignment columns ", region)
    trim_to_region(aln, r[1], r[2])
  }
  settings$region <- region
  meta <- set_table(s)
  meta <- meta[, setdiff(names(meta), "database_md5"), drop = FALSE]
  view <- compare_alignment(aln, red, ids, meta, settings, out_dir, trimmed,
                            progress = function(step) progress(step), log = log, number = ids$internal)
  view$label <- s$name
  view$within <- TRUE
  view$sequences <- n

  progress("Writing results")
  set_write(s, file.path(out_dir, "set"))
  turned <- al$turned
  manifest <- list(
    created = format(started, "%Y-%m-%d %H:%M:%S"),
    user = Sys.info()[["user"]],
    computer = Sys.info()[["nodename"]],
    app_version = app_version(),
    mode = "within_set",
    input = w$input,
    alignment = list(file = "set/aligned.fasta", md5 = unname(tools::md5sum(file.path(out_dir, "set", "aligned.fasta"))),
                     made_in_this_run = aligned_now, from = al$from, made = al$made, method = al$method,
                     versions = al$versions, reverse_complemented = as.list(turned)),
    queries = as.list(view$queries),
    reverse_complemented = as.list(turned),
    settings = settings[intersect(c("bootstrap", "threads", "region"), names(settings))],
    results = list(list(queries = s$name, folder = ".", sequences = n, compared_columns = view$compared_columns,
                        compared_region = view$compared_region, dropped = as.list(view$dropped),
                        limited_by = if (is.null(region)) as.list(view$limited_by))),
    tools = tool_versions(cfg, settings),
    method = method_text(settings, "within_set")
  )
  manifest$results[[1]] <- Filter(Negate(is.null), manifest$results[[1]])
  yaml::write_yaml(manifest, file.path(out_dir, "run_manifest.yaml"))
  log("Finished in ", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), " s")
  empty <- data.frame(query = character(0))
  list(out_dir = out_dir, hits = empty, ties = empty, queries = view$queries, views = stats::setNames(list(view), s$name),
       manifest = manifest, display_columns = setdiff(names(meta), "id"), mode = "within_set", set = s,
       aligned_now = aligned_now)
}

# "FROM-TO" (alignment columns) -> c(from, to), checked against the alignment length.
parse_region <- function(region, ncol) {
  r <- suppressWarnings(as.integer(strsplit(trimws(region), "\\s*-\\s*")[[1]]))
  if (length(r) != 2 || anyNA(r) || r[1] < 1 || r[2] > ncol || r[1] >= r[2]) {
    stop(sprintf("The region '%s' is not two alignment columns FROM-TO between 1 and %d (e.g. 2001-3000).", region, ncol),
         call. = FALSE)
  }
  r
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
# Returns the trimmed alignment, the dropped sequences, the kept columns' numbers and the
# queries whose start or end limits the region (`limited_by`).
trim_to_queries <- function(aln, query_ids, no_overlap = paste(
  "The query sequences cover different parts of the genome and do not overlap,",
  "so they cannot be compared together. Each query's own result is still made.")) {
  chars <- do.call(rbind, strsplit(aln, ""))
  rownames(chars) <- names(aln)
  spans <- t(vapply(query_ids, function(q) range(which(chars[q, ] != "-")), numeric(2)))
  from <- max(spans[, 1])
  to <- min(spans[, 2])
  if (from >= to) stop(no_overlap, call. = FALSE)
  limited_by <- query_ids[spans[, 1] == from | spans[, 2] == to]
  trimmed <- trim_columns(chars, from:to, keep = query_ids)
  trimmed$limited_by <- limited_by
  trimmed
}

# Keep alignment columns from..to (a region picked by the user); then drop all-gap columns
# and every sequence covering < 50% of the region (with A/C/G/T).
trim_to_region <- function(aln, from, to) {
  chars <- do.call(rbind, strsplit(aln, ""))
  rownames(chars) <- names(aln)
  out <- trim_columns(chars, from:to)
  if (length(out$aln) < 2) {
    stop(sprintf(paste("Fewer than 2 sequences have bases in at least half of columns %d-%d.",
                       "Pick a larger region in the Alignment page, or another part of the genome."), from, to),
         call. = FALSE)
  }
  out
}

# Columns `cols` of a character matrix; drop all-gap columns, then the sequences (except
# `keep`) with A/C/G/T in < 50% of what is left.
trim_columns <- function(chars, cols, keep = character(0)) {
  chars <- chars[, cols, drop = FALSE]
  cols <- cols[colSums(chars != "-") > 0]
  chars <- chars[, colSums(chars != "-") > 0, drop = FALSE]
  acgt <- matrix(chars %in% c("A", "C", "G", "T"), nrow = nrow(chars))
  covered <- rowMeans(acgt)
  drop <- setdiff(rownames(chars)[covered < 0.5], keep)
  chars <- chars[!rownames(chars) %in% drop, , drop = FALSE]
  out <- apply(chars, 1, paste, collapse = "")
  list(aln = out, dropped = if (length(drop)) drop else NULL, cols = cols)
}

method_text <- function(settings, mode = "database") {
  tree <- paste(
    "Tree: ape::bionj on TN93 distances (pairwise deletion), negative branch lengths set to 0,",
    "+ 100 bootstraps (ape::boot.phylo, seed 1) when bootstrap is on and >= 4 sequences;",
    "midpoint rooted (phangorn::midpoint).")
  snps <- paste(
    "SNPs: phangorn::dist.hamming(ratio = FALSE, exclude = 'none'): a position counts when the two",
    "bases cannot be the same (IUPAC ambiguity codes count only when they cannot match; N and gaps",
    "never count).")
  if (mode == "within_set") {
    return(paste(
      "Compare within sets: no database search. The sequences of the set(s) are compared with each other",
      "in the set's alignment (when it has all of them), else aligned in the run: DECIPHER OrientNucleotides",
      "(reverse complement only, against 6 HAV genomes, one per subgenotype), then DECIPHER AlignSeqs.",
      if (is.null(settings$region)) "Alignment trimmed to the columns covered by every sequence."
      else paste0("Alignment trimmed to the region (alignment columns ", settings$region, "); sequences with A/C/G/T",
                  " in < 50% of it left out."),
      snps, tree))
  }
  text <- paste(
    switch(settings$search,
      decipher = "Search: DECIPHER IndexSeqs (k = 11) + SearchIndex, both strands, up to 500 hits per query, each aligned with AlignPairs;",
      blast = "Search: blastn via rBLAST (max_target_seqs 500, max_hsps 1);"),
    "per query, hits covering >= min_coverage % of the query, ranked by percent identity then score",
    "(ties by database order), top n_neighbours kept.",
    if (isTRUE(settings$exclude_self)) "Hits to the database entry with the query's own name were left out.",
    if (is.list(settings$exclude_self)) "For the queries listed under exclude_self, hits to the database entry given there (their own) were left out.",
    switch(settings$aligner,
      mafft = "Alignment: MAFFT --auto (pymafft).",
      decipher = "Alignment: DECIPHER AlignSeqs."),
    "Alignment trimmed to the columns covered by every query.",
    "With several queries, alignment, SNPs and tree are made for all queries together (all their",
    "neighbours) and for each query alone (its own neighbours, trimmed to its own length).",
    snps, tree
  )
  if (mode != "set_database") return(text)
  paste("Database: the sequences of the set(s) only, turned to one orientation (DECIPHER OrientNucleotides,",
        "reverse complement only, against 6 HAV genomes, one per subgenotype), in name order.", text)
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
