# Alignment viewer (Alignment page, decision 0011): what each cell of an alignment shows
# when compared with an anchor (a sequence, or the majority consensus).
#
# Every cell gets one class, from the same bit rule as the SNP count (R/snps.R):
#   0 same as the anchor (".")           4 indel: a gap where the anchor has a base, or a
#   1 mismatch (A/C/G/T vs A/C/G/T)        base where the anchor has a gap (never a SNP)
#   2 ambiguity code that fits           5 N (or another unknown letter): no data
#   3 ambiguity code that cannot match   6 missing start or end: no data
#                                        7 the anchor has no base here (N, missing end)
# A gap where the anchor also has a gap is 0. Classes 1 + 3 are exactly the SNP columns
# (snp_columns()); tests/run_tests.R checks it.
# The anchor row itself: 8 base, 2 ambiguity code, 4 gap, 5 N, 6 missing end.
# The browser (app/www/msa.js) only draws these classes; it does not compare bases.
#
# Alignments are kept as raw matrices (one byte per cell), so a whole database
# (~1,800 x ~8,300) fits in ~15 MB.

MSA_GAP <- charToRaw("-")
MSA_N <- charToRaw("N")
MSA_CONSENSUS <- "Majority consensus"
# Most important class first: what a pixel shows when it covers several columns.
MSA_PRIORITY <- c(1L, 3L, 4L, 2L, 7L, 5L, 6L, 0L)

# Bits of each byte (IUPAC_BITS, R/snps.R); anything else is unknown, like N.
msa_bits_lut <- function() {
  lut <- rep(15L, 256)
  lut[as.integer(charToRaw(paste(names(IUPAC_BITS), collapse = ""))) + 1L] <- unname(IUPAC_BITS)
  lut
}
msa_bits <- function(x) msa_bits_lut()[as.integer(x) + 1L]

# Aligned sequences (named character vector, equal lengths) -> raw matrix, upper case,
# "." read as a gap.
msa_matrix <- function(aln) {
  aln <- chartr(".", "-", toupper(aln))
  if (length(unique(nchar(aln))) != 1) {
    stop("This file is not aligned: the sequences have different lengths. ",
         "Choose an aligned FASTA file (all sequences the same length, with gaps as -).", call. = FALSE)
  }
  m <- matrix(unlist(lapply(aln, charToRaw), use.names = FALSE), nrow = length(aln), byrow = TRUE)
  rownames(m) <- names(aln)
  m
}

# Leading and trailing gaps of one sequence (raw vector): its missing ends.
msa_ends <- function(x) {
  base <- which(x != MSA_GAP)
  cols <- seq_along(x)
  if (!length(base)) return(rep(TRUE, length(x)))
  cols < base[1] | cols > base[length(base)]
}

# Each sequence's own base number at each column (cumulative count of non-gaps). At a gap,
# the number of the base before it (0 before the first base).
msa_positions <- function(x) cumsum(x != MSA_GAP)

# The column of base number `pos` of a sequence (raw vector), or NA.
msa_column_of <- function(x, pos) {
  col <- which(x != MSA_GAP & msa_positions(x) == pos)
  if (length(col)) col[1] else NA_integer_
}

# Classes of one sequence x against the anchor a (raw vectors; see the top of this file).
# `a_has_base`: the anchor has a base other than N there (FALSE also at its gaps and
# missing ends); `a_gap`: the anchor has a gap inside it.
msa_row_classes <- function(x, a_bits, a_acgt, a_has_base, a_gap, a) {
  x_bits <- msa_bits(x)
  x_acgt <- x_bits %in% c(1L, 2L, 4L, 8L)
  fits <- bitwAnd(x_bits, a_bits) != 0L
  cls <- ifelse(fits, ifelse(x_acgt & x == a, 0L, 2L), ifelse(x_acgt & a_acgt, 1L, 3L))
  is_gap <- x == MSA_GAP
  cls[is_gap] <- 4L
  cls[a_gap] <- ifelse(is_gap[a_gap], 0L, 4L)
  cls[!a_has_base & !a_gap] <- 7L
  cls[x_bits == 15L & !is_gap] <- 5L
  cls[msa_ends(x)] <- 6L
  cls
}

# Classes of every row (raw matrix m) against the anchor (raw vector). Returns a raw matrix.
msa_classes <- function(m, anchor) {
  a_bits <- msa_bits(anchor)
  a_bits[msa_ends(anchor)] <- 15L
  a_acgt <- a_bits %in% c(1L, 2L, 4L, 8L)
  a_has_base <- a_bits != 15L
  a_gap <- anchor == MSA_GAP & !msa_ends(anchor)
  out <- matrix(as.raw(0), nrow(m), ncol(m), dimnames = dimnames(m))
  for (i in seq_len(nrow(m))) out[i, ] <- as.raw(msa_row_classes(m[i, ], a_bits, a_acgt, a_has_base, a_gap, anchor))
  out
}

# Display classes of the anchor row itself.
msa_anchor_classes <- function(anchor) {
  bits <- msa_bits(anchor)
  cls <- ifelse(bits %in% c(1L, 2L, 4L, 8L), 8L, 2L)
  cls[anchor == MSA_GAP] <- 4L
  cls[bits == 15L & anchor != MSA_GAP] <- 5L
  cls[msa_ends(anchor)] <- 6L
  cls
}

# Majority consensus (raw vector). Per column, A/C/G/T votes: ambiguity codes split their
# vote (R = 1/2 A + 1/2 G); N, gaps and missing ends do not vote. A column where more
# sequences have a gap inside them than a base is "-"; a column without votes is N.
# The most common base wins. When it is not clear - a tie, or with a threshold (50-100%)
# the most common base has less than that share - the result is the ambiguity code of
# the fewest top bases that together reach it (ties = "iupac") or N (ties = "n").
msa_consensus <- function(m, ties = c("iupac", "n"), threshold = NULL) {
  ties <- match.arg(ties)
  nc <- ncol(m)
  # Gaps inside a sequence: between its first and last base.
  nongap <- m != MSA_GAP
  has_base <- rowSums(nongap) > 0
  first <- ifelse(has_base, max.col(nongap, "first"), nc + 1L)
  last <- ifelse(has_base, max.col(nongap, "last"), 0L)
  column <- matrix(rep(seq_len(nc), each = nrow(m)), nrow(m))
  gaps_inside <- colSums(!nongap & column > first & column < last)
  # Votes: count each letter per column, then share its vote among the bases it can be.
  votes <- matrix(0, 4, nc)
  n_base <- numeric(nc)
  lut <- msa_bits_lut()
  for (code in setdiff(unique(as.vector(m)), MSA_GAP)) {
    bits <- lut[as.integer(code) + 1L]
    if (bits == 15L) next
    count <- colSums(m == code)
    n_base <- n_base + count
    for (b in 1:4) if (bitwAnd(bits, bitwShiftL(1L, b - 1L))) votes[b, ] <- votes[b, ] + count / sum(intToBits(bits)[1:4] == 1)
  }
  votes <- round(votes, 9)   # split votes (1/3) must tie exactly
  code <- stats::setNames(names(IUPAC_BITS)[1:15], IUPAC_BITS[1:15])
  need <- if (is.null(threshold)) NA else threshold / 100
  out <- vapply(seq_len(nc), function(j) {
    v <- votes[, j]
    total <- sum(v)
    if (gaps_inside[j] > n_base[j]) return("-")
    if (total == 0) return("N")
    o <- order(-v)
    k <- if (is.na(need)) sum(v == v[o[1]]) else which(cumsum(v[o]) / total >= need - 1e-9)[1]
    top <- o[v[o] >= v[o[k]]]   # bases tied with the k-th are in too
    if (length(top) == 1) return(c("A", "C", "G", "T")[top])
    if (ties == "n") return("N")
    code[[as.character(sum(bitwShiftL(1L, top - 1L)))]]
  }, "")
  charToRaw(paste(out, collapse = ""))
}

# ---- Loading -------------------------------------------------------------------------

# An alignment for the viewer: the matrix, plus what the window shows about each row.
#   meta      data.frame with id (+ any metadata columns) for the hover box
#   queries   names shown in red; dropped: left out of the SNP table (grey)
#   compared  columns compared in the run (others greyed), or NULL (all)
#   revcomp   names of queries numbered on their reverse complement
#   order     row order (tree order) or NULL (sorted by differences to the anchor)
#   snp       SNP matrix of the run (for "query and its closest"), or NULL
msa_load <- function(aln, label, meta = NULL, queries = character(0), dropped = character(0),
                     compared = NULL, revcomp = character(0), order = NULL, snp = NULL, note = NULL) {
  m <- msa_matrix(aln)
  if (anyDuplicated(rownames(m))) stop("Duplicate sequence names in the alignment.", call. = FALSE)
  if (is.null(meta)) meta <- data.frame(id = rownames(m))
  keep <- if (is.null(order)) rownames(m) else c(intersect(order, rownames(m)), setdiff(rownames(m), order))
  m <- m[keep, , drop = FALSE]
  list(label = label, m = m, cols = seq_len(ncol(m)), meta = meta, queries = intersect(queries, keep),
       dropped = intersect(dropped, keep), revcomp = intersect(revcomp, keep),
       compared = if (is.null(compared)) NULL else seq_len(ncol(m)) %in% compared,
       tree_order = !is.null(order), snp = snp, note = note)
}

# Only some rows (no re-alignment). Columns that are a gap in every chosen row go; `cols`
# keeps each column's number in the full alignment.
msa_subset <- function(a, rows) {
  rows <- intersect(rownames(a$m), rows)
  m <- a$m[rows, , drop = FALSE]
  keep <- colSums(m != MSA_GAP) > 0
  a$m <- m[, keep, drop = FALSE]
  a$cols <- a$cols[keep]
  if (!is.null(a$compared)) a$compared <- a$compared[keep]
  a$queries <- intersect(a$queries, rows)
  a$dropped <- intersect(a$dropped, rows)
  a$revcomp <- intersect(a$revcomp, rows)
  a
}

# A query and its n closest sequences by SNPs (snp matrix of the run). Ties at the n-th
# place are all included. Returns list(ids, n_tied_extra).
msa_closest <- function(snp, query, n) {
  s <- snp[query, setdiff(colnames(snp), query)]
  s <- s[order(s)]
  if (!length(s)) return(list(ids = query, extra = 0L))
  limit <- s[min(n, length(s))]
  ids <- names(s)[s <= limit]
  list(ids = c(query, ids), extra = max(0L, length(ids) - as.integer(n)))
}

# Cells of class k (one or several) in a raw class matrix or vector. (== is ~20 times
# faster than %in% on raw values; this matters for a whole database.)
msa_is <- function(cls, k) {
  hit <- cls == as.raw(k[1])
  for (x in k[-1]) hit <- hit | cls == as.raw(x)
  hit
}

# The alignment compared with an anchor: a row name, or NULL for the majority consensus.
# Rows are in tree order (run results) or by differences to the anchor; the anchor row is
# pinned at the top. Per row, over the compared columns (run results) or the whole
# alignment: SNPs (classes 1 + 3; = the SNP table; NA for dropped rows), indel events (a
# run of indel columns counts once) and columns, ambiguity codes and no-data columns.
# Per column: how many rows differ from the anchor (SNP or indel), for the overview.
msa_view <- function(a, anchor = NULL, ties = "iupac", threshold = NULL) {
  m <- a$m
  if (!is.null(anchor) && !anchor %in% rownames(m)) anchor <- NULL
  anchor_seq <- if (is.null(anchor)) msa_consensus(m, ties, threshold) else m[anchor, ]
  cls <- msa_classes(m, anchor_seq)
  counted <- if (is.null(a$compared)) rep(TRUE, ncol(m)) else a$compared
  diffs <- rowSums(msa_is(cls, c(1, 3)) & rep(counted, each = nrow(m)))
  names(diffs) <- rownames(m)
  diffs[a$dropped] <- NA
  indel <- msa_is(cls, 4) & rep(counted, each = nrow(m))
  starts <- indel & cbind(TRUE, indel[, -ncol(m), drop = FALSE] == FALSE)
  stats <- list(indel_events = unname(rowSums(starts)), indel_cols = unname(rowSums(indel)),
                ambiguity = unname(rowSums(msa_is(cls, c(2, 3)) & rep(counted, each = nrow(m)))),
                nodata = unname(rowSums(msa_is(cls, c(5, 6)) & rep(counted, each = nrow(m)))),
                n_counted = sum(counted))
  others <- setdiff(rownames(m), anchor)
  if (!a$tree_order) others <- others[order(diffs[others], match(others, rownames(m)))]
  rows <- match(others, rownames(m))

  # Columns that differ from the anchor in any (other) row, and the anchor's own gaps.
  anchor_gap <- anchor_seq == MSA_GAP & !msa_ends(anchor_seq)
  c_other <- cls[rows, , drop = FALSE]
  variable <- colSums(msa_is(c_other, 1:4)) > 0
  ndiff <- unname(colSums(msa_is(c_other, c(1, 3, 4))))

  list(anchor = anchor, anchor_name = anchor %||% MSA_CONSENSUS, anchor_seq = anchor_seq,
       anchor_classes = msa_anchor_classes(anchor_seq), rows = rows, classes = cls, diffs = diffs,
       stats = stats, variable = variable, anchor_gap = anchor_gap, ndiff = ndiff)
}

# First column (after `from`) where a row differs from the anchor (classes 1 and 3), within
# the compared columns when there are any. NA when there is none.
msa_first_difference <- function(a, v, row, from = 0L) {
  hit <- msa_is(v$classes[row, ], c(1, 3)) & (if (is.null(a$compared)) TRUE else a$compared)
  hit <- which(hit & seq_along(hit) > from)
  if (length(hit)) hit[1] else NA_integer_
}

# What the window draws, as plain strings (one per row) for the browser. Dates (for
# sorting and the labels) come from the metadata column that holds dates ("date" first).
msa_payload <- function(a) {
  meta <- a$meta[match(rownames(a$m), a$meta$id), setdiff(names(a$meta), c("id", "is_query")), drop = FALSE]
  date_col <- c(intersect("date", date_columns(meta)), date_columns(meta))[1]
  dates <- if (is.na(date_col)) NULL else {
    d <- parse_partial_date(meta[[date_col]])
    list(text = ifelse(is.na(d), "", as.character(meta[[date_col]])), sort = ifelse(is.na(d), "", format(d)))
  }
  list(label = a$label, names = rownames(a$m), ncol = ncol(a$m), cols = a$cols, dates = dates,
       tree_order = a$tree_order,
       chars = unname(apply(a$m, 1, rawToChar)),
       query = rownames(a$m) %in% a$queries, dropped = rownames(a$m) %in% a$dropped,
       revcomp = rownames(a$m) %in% a$revcomp,
       compared = if (is.null(a$compared)) NULL else paste(as.integer(a$compared), collapse = ""),
       meta = lapply(seq_len(nrow(meta)), function(i) {
         vals <- vapply(meta[i, , drop = FALSE], function(x) if (is.na(x)) "?" else as.character(x), "")
         as.list(vals)
       }),
       note = a$note)
}

msa_view_payload <- function(a, v) {
  digits <- function(x) paste(as.integer(x), collapse = "")
  list(anchor = if (is.null(v$anchor)) NULL else match(v$anchor, rownames(a$m)) - 1L,
       anchor_name = v$anchor_name, anchor_chars = rawToChar(v$anchor_seq),
       anchor_classes = digits(v$anchor_classes), rows = v$rows - 1L,
       classes = unname(apply(v$classes, 1, function(x) rawToChar(as.raw(48L + as.integer(x))))),
       diffs = unname(v$diffs), stats = v$stats, variable = digits(v$variable),
       anchor_gap = digits(v$anchor_gap), ndiff = v$ndiff)
}

# ---- Export --------------------------------------------------------------------------

# Colours of the classes, as in app/www/msa.js (no data = light grey, indel = purple).
MSA_COLOURS <- c(`0` = "#ffffff", `1` = "#e6550d", `2` = "#c6dbef", `3` = "#2166ac", `4` = "#7b3294",
                 `5` = "#e6e9ed", `6` = "#e6e9ed", `7` = "#cfd4da", `8` = "#fff3cd")

# PDF of the shown window: rows (anchor first) and columns (column numbers of a$m). With
# more than 250 columns, each drawn cell covers several columns and shows the most
# important class among them (as on screen); letters only up to 250 columns.
msa_pdf <- function(a, v, rows, cols, file) {
  label_rows <- c(v$anchor_name, rownames(a$m)[rows])
  cls <- rbind(v$anchor_classes[cols], matrix(as.integer(v$classes[rows, cols, drop = FALSE]), length(rows)))
  chars <- rbind(strsplit(rawToChar(v$anchor_seq[cols]), "")[[1]],
                 matrix(strsplit(paste(apply(a$m[rows, cols, drop = FALSE], 1, rawToChar), collapse = ""), "")[[1]],
                        length(rows), byrow = TRUE))
  letters <- length(cols) <= 250
  if (!letters) {
    bin <- ceiling(seq_along(cols) / (length(cols) / 250))
    rank <- match(cls, c(MSA_PRIORITY, 8L))
    cls <- t(apply(matrix(rank, nrow(cls)), 1, function(r) c(MSA_PRIORITY, 8L)[tapply(r, bin, min)]))
  }
  nr <- nrow(cls)
  nc <- ncol(cls)
  grDevices::pdf(file, width = min(50, 3 + nc * 0.11), height = 1.5 + nr * 0.17)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(2, 12, 2, 1), xpd = NA)
  graphics::plot.new()
  graphics::plot.window(c(0, nc), c(nr, 0), xaxs = "i", yaxs = "i")
  x <- rep(seq_len(nc), each = nr)
  y <- rep(seq_len(nr), nc)
  graphics::rect(x - 1, y - 1, x, y, col = MSA_COLOURS[as.character(cls)], border = NA)
  if (!is.null(a$compared)) {
    out <- which(!a$compared[cols])
    if (letters && length(out)) graphics::rect(out - 1, 0, out, nr, col = grDevices::adjustcolor("grey50", 0.35), border = NA)
  }
  if (letters) {
    shown <- ifelse(cls == 0L, ".", ifelse(cls == 6L, "", chars))
    graphics::text(x - 0.5, y - 0.5, shown, cex = 0.55, family = "mono",
                   col = ifelse(cls %in% c(3L, 4L), "white", ifelse(cls %in% c(0L, 5L), "grey45", "black")))
  }
  graphics::text(-0.3, seq_len(nr) - 0.5, label_rows, adj = 1, cex = 0.6,
                 font = ifelse(label_rows %in% a$queries, 2, 1),
                 col = ifelse(label_rows %in% a$queries, "#d7301f", "black"))
  first <- a$cols[cols[1]]
  last <- a$cols[cols[length(cols)]]
  graphics::title(main = sprintf("%s - columns %d-%d, compared with %s", a$label, first, last, v$anchor_name),
                  cex.main = 0.8, font.main = 1)
  invisible(file)
}

# A standalone page (one HTML file, no internet needed) with the viewer as shown: the
# same drawing code (app/www/msa.js) and the classes of the current anchor. Zoom, move,
# sort, hover and the two column toggles work; another anchor needs FHInumerics.
msa_html <- function(a, v, file, sort = NULL, created = format(Sys.time(), "%Y-%m-%d %H:%M")) {
  www <- file.path(getOption("fhinumerics.root", "."), "app", "www")
  read <- function(f) paste(readLines(file.path(www, f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  json <- function(x) gsub("</", "<\\/", jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null"), fixed = TRUE)
  esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE)
  sub_text <- sprintf("%d sequences, %d columns. Compared with %s. Made with FHInumerics %s on %s.",
                      nrow(a$m), ncol(a$m), if (is.null(v$anchor)) "the majority consensus" else v$anchor,
                      app_version(), created)
  page <- c("<!doctype html>", "<html lang=\"en\"><head><meta charset=\"utf-8\">",
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
            paste0("<title>", esc(a$label), "</title>"),
            "<style>body { margin: 0; padding: 14px 20px; font: 14px/1.4 system-ui, 'Segoe UI', sans-serif; color: #1f2328; background: #fff; }",
            ".fhi-tip { background: white; border: 1px solid #bbb; border-radius: 6px; box-shadow: 0 2px 8px rgba(0,0,0,.15); padding: 6px 9px; font-size: 12px; max-width: 420px; z-index: 1000; }",
            "h1 { font-size: 18px; margin: 0 0 2px; } .sub { color: #656d76; font-size: 13px; margin-bottom: 8px; }",
            read("msa.css"), "</style></head><body>",
            paste0("<h1>", esc(a$label), "</h1><div class=\"sub\">", esc(sub_text), "</div>"),
            "<div id=\"msa\"></div>",
            "<script>", read("msa.js"), "</script>",
            "<script>",
            paste0("fhiMsa.load(", json(msa_payload(a)), ", ", json(msa_view_payload(a, v)),
                   ", {standalone: true, sort: ", json(sort %||% NA), "});"),
            "</script></body></html>")
  writeLines(page, file, useBytes = TRUE)
  invisible(file)
}

# Aligned FASTA of rows (row numbers of a$m) and the columns from..to (inclusive).
msa_fasta <- function(a, rows, from = 1L, to = ncol(a$m)) {
  stats::setNames(apply(a$m[rows, from:to, drop = FALSE], 1, rawToChar), rownames(a$m)[rows])
}

# An aligned FASTA file chosen by the user (names up to the first space), checked.
msa_read_fasta <- function(path) {
  aln <- read_alignment(path)
  if (!length(aln)) stop("This file has no sequences in FASTA format (lines starting with >).", call. = FALSE)
  names(aln) <- sub("\\s.*", "", names(aln))
  aln <- gsub("[[:space:]]", "", aln)
  bad <- names(aln)[grepl(paste0("[^", IUPAC_NT, "U.-]"), aln)]
  if (length(bad)) {
    stop("Not a nucleotide alignment: ", paste(utils::head(bad, 5), collapse = ", "),
         if (length(bad) > 5) " and more", " (unexpected letters).", call. = FALSE)
  }
  dups <- unique(names(aln)[duplicated(names(aln))])
  if (length(dups)) stop("Duplicate sequence names: ", paste(dups, collapse = ", "), call. = FALSE)
  chartr("U", "T", aln)
}

# ---- Run folders ------------------------------------------------------------------------

# The results in a run folder (from its run_manifest.yaml): names = labels, values = folders.
msa_run_results <- function(run_dir) {
  file <- file.path(run_dir, "run_manifest.yaml")
  if (!file.exists(file)) {
    stop("This folder is not a FHInumerics run folder (it has no run_manifest.yaml). ",
         "Choose the folder of one run, inside the results folder.", call. = FALSE)
  }
  man <- yaml::read_yaml(file)
  res <- Filter(function(r) is.null(r$error), man$results %||% list(list(queries = man$queries[[1]], folder = ".")))
  folders <- vapply(res, function(r) file.path(run_dir, r$folder), "")
  names(folders) <- vapply(res, function(r) paste(unlist(r$queries), collapse = ", "), "")
  attr(folders, "manifest") <- man
  attr(folders, "regions") <- lapply(res, function(r) r$compared_region)
  attr(folders, "dropped") <- lapply(res, function(r) unlist(r$dropped))
  folders
}

# One result of an earlier run as a viewer alignment. Older runs have only the trimmed
# alignment (alignment.fasta); it is shown with a note.
msa_load_run_folder <- function(folder, label, manifest, region = NULL, dropped = NULL) {
  untrimmed <- file.path(folder, "alignment_untrimmed.fasta")
  note <- NULL
  if (!file.exists(untrimmed)) {
    untrimmed <- file.path(folder, "alignment.fasta")
    region <- NULL
    note <- "This run is older than the Alignment page: only the compared region was saved, so only that is shown."
  }
  aln <- read_alignment(untrimmed)
  meta_file <- file.path(folder, "tree_metadata.tsv")
  meta <- if (file.exists(meta_file)) utils::read.delim(meta_file, check.names = FALSE, colClasses = "character") else NULL
  queries <- if (!is.null(meta)) meta$id[meta$is_query == "TRUE"] else character(0)
  order <- NULL
  if (file.exists(file.path(folder, "tree.nwk"))) order <- tree_order(ape::read.tree(file.path(folder, "tree.nwk")))
  snp_file <- file.path(folder, "snp_matrix.tsv")
  snp <- if (file.exists(snp_file)) as.matrix(utils::read.delim(snp_file, row.names = 1, check.names = FALSE)) else NULL
  compared <- if (length(region)) { r <- as.integer(strsplit(region, "-")[[1]]); r[1]:r[2] } else NULL
  msa_load(aln, label, meta = meta, queries = queries, dropped = dropped %||% character(0), compared = compared,
           revcomp = unlist(manifest$reverse_complemented), order = order, snp = snp, note = note)
}

# Tip labels from top to bottom, as the tree is drawn (and the SNP table is ordered).
tree_order <- function(tree) {
  d <- ggtree::fortify(tree)
  d$label[d$isTip][order(d$y[d$isTip], decreasing = TRUE)]
}

# A result of the current run as a viewer alignment (the untrimmed alignment).
msa_from_result <- function(v, revcomp = character(0)) {
  msa_load(v$aln_full, v$label, meta = v$aln_meta, queries = v$queries, dropped = v$dropped,
           compared = v$compared_cols, revcomp = revcomp,
           order = if (is.null(v$tree)) NULL else tree_order(v$tree), snp = v$snp)
}
