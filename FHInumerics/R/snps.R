# SNP counting that takes IUPAC ambiguity codes into account.
#
# Two sequences differ at a position only when their bases cannot be the same:
#   A vs T            -> 1 SNP
#   K (G/T) vs A      -> 1 SNP  (K cannot be A)
#   M (A/C) vs A      -> 0      (M can be A)
#   N, or a gap       -> 0      (never counted; gaps are ignored for now, see docs)
# The SNP matrix is counted by phangorn::dist.hamming(ratio = FALSE, exclude = "none").
# The same rule written out per position (IUPAC_BITS) gives the number of compared
# positions and the list of differing positions shown in the window;
# tests/run_tests.R checks that both give the same numbers.

# The bases each IUPAC code can be, as bits: A = 1, C = 2, G = 4, T = 8. Two bases can be
# the same when their bits overlap. N and gaps can be anything (all four bits).
IUPAC_BITS <- c(A = 1L, C = 2L, G = 4L, T = 8L, R = 5L, Y = 10L, S = 6L, W = 9L, K = 12L, M = 3L,
                B = 14L, D = 13L, H = 11L, V = 7L, N = 15L, `-` = 15L)

# Alignment (named character vector) -> integer matrix of bits, one row per sequence.
base_bits <- function(aln) {
  chars <- do.call(rbind, strsplit(toupper(aln), ""))
  bits <- matrix(IUPAC_BITS[chars], nrow = length(aln), dimnames = list(names(aln), NULL))
  bits[is.na(bits)] <- 15L   # any other character: unknown, like N
  bits
}

# Number of SNPs between every pair of sequences in the alignment (phangorn).
snp_matrix <- function(aln) {
  chars <- do.call(rbind, strsplit(tolower(aln), ""))
  rownames(chars) <- names(aln)
  snp <- as.matrix(phangorn::dist.hamming(phangorn::phyDat(chars, type = "DNA"), ratio = FALSE, exclude = "none"))
  snp <- round(snp[names(aln), names(aln), drop = FALSE])
  storage.mode(snp) <- "integer"
  snp
}

# Number of positions where both sequences have a base other than N or a gap
# (the denominator of the SNP count).
comparable_sites <- function(aln) {
  known <- (base_bits(aln) != 15L) * 1L
  sites <- tcrossprod(known)
  storage.mode(sites) <- "integer"
  dimnames(sites) <- list(names(aln), names(aln))
  sites
}

# Alignment columns where sequences a and b differ (same rule as snp_matrix()).
snp_columns <- function(bits, a, b) which(bitwAnd(bits[a, ], bits[b, ]) == 0L)

# Hover text for the SNP table: for each pair with at least one query, the number of
# SNPs and compared positions, then what differs, numbered along the query:
#   SNPs                               counted (snp_columns)
#   Gaps (not counted)                 each gap inside a sequence as one line per run of
#                                      gaps; missing start/end of a sequence as one range
#   Ambiguity codes that fit           e.g. A/M: not counted, listed so they can be seen
#   N                                  number of Ns in each sequence (not compared)
# Returns list(text, flagged): named by "<row>\t<column>" (both orders); flagged = the
# pair also differs by gaps, Ns or fitting ambiguity codes (marked in the table).
# `view` is one result of run_compare() (needs aln, positions, queries, snp, sites).
snp_tips <- function(view, reverse_complemented = character(0), max_listed = 20) {
  bits <- base_bits(view$aln)
  chars <- do.call(rbind, strsplit(toupper(view$aln), ""))
  rownames(chars) <- names(view$aln)
  ids <- names(view$aln)
  pairs <- expand.grid(row = ids, col = ids, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$row != pairs$col & (pairs$row %in% view$queries | pairs$col %in% view$queries), ]
  tips <- lapply(seq_len(nrow(pairs)), function(i) {
    pair_tip(view, chars, bits, pairs$row[i], pairs$col[i], reverse_complemented, max_listed)
  })
  keys <- paste(pairs$row, pairs$col, sep = "\t")
  list(text = stats::setNames(vapply(tips, `[[`, "", "text"), keys),
       flagged = stats::setNames(vapply(tips, `[[`, TRUE, "flagged"), keys))
}

# Hover text of one pair: row r, column k (see snp_tips()).
pair_tip <- function(view, chars, bits, r, k, reverse_complemented, max_listed) {
  q <- if (r %in% view$queries) r else k
  pos <- view$positions[[q]]
  a <- chars[r, ]
  b <- chars[k, ]
  cap <- function(lines) c(utils::head(lines, max_listed),
                           if (length(lines) > max_listed) sprintf("... and %d more", length(lines) - max_listed))
  bases <- function(cols) paste0(paste(a[cols], collapse = ""), "/", paste(b[cols], collapse = ""))
  span <- function(cols) if (min(pos[cols]) == max(pos[cols])) pos[cols[1]] else paste0(min(pos[cols]), "-", max(pos[cols]))

  snp <- snp_columns(bits, r, k)

  # Gaps: exactly one of the two has a gap. Gaps before a sequence's first base or after
  # its last are "missing" ends; the rest are gaps inside, one line per run.
  gap <- (a == "-") != (b == "-")
  end_gap <- function(x) cumsum(x != "-") == 0 | rev(cumsum(rev(x != "-")) == 0)
  missing <- character(0)
  for (who in c(r, k)) {
    x <- if (who == r) a else b
    lead <- which(gap & x == "-" & cumsum(x != "-") == 0)
    trail <- which(gap & x == "-" & rev(cumsum(rev(x != "-")) == 0))
    if (length(lead)) missing <- c(missing, paste(span(lead), "- missing at the start of", who))
    if (length(trail)) missing <- c(missing, paste(span(trail), "- missing at the end of", who))
  }
  inside <- gap & !((a == "-" & end_gap(a)) | (b == "-" & end_gap(b)))
  side <- ifelse(inside, ifelse(a == "-", "r", "k"), "")
  runs <- rle(side)
  last <- cumsum(runs$lengths)
  gaps_inside <- unlist(lapply(which(runs$values != ""), function(j) {
    cols <- (last[j] - runs$lengths[j] + 1):last[j]
    in_query <- (runs$values[j] == "r") == (q == r)
    # A gap in the query has no base number of its own: give the base before it.
    paste(if (in_query) paste("after", pos[cols[1]]) else span(cols), "-", bases(cols))
  }))

  # Ambiguity codes (not N) that fit the other base, so they are not SNPs.
  ambiguous <- function(x) !x %in% c("A", "C", "G", "T", "N", "-")
  fits <- which((ambiguous(a) | ambiguous(b)) & !a %in% c("N", "-") & !b %in% c("N", "-") &
                  bitwAnd(bits[r, ], bits[k, ]) != 0L)
  n_count <- c(sum(a == "N"), sum(b == "N"))
  names(n_count) <- c(r, k)
  n_count <- n_count[n_count > 0]

  n <- view$snp[r, k]
  lines <- c(paste(r, "/", k), sprintf("%d SNP%s over %d compared positions", n, if (n == 1) "" else "s", view$sites[r, k]))
  details <- c(
    if (length(snp)) c("SNPs:", cap(paste(pos[snp], "-", vapply(snp, bases, "")))),
    if (length(missing) || length(gaps_inside)) c("Gaps (not counted):", cap(c(missing, gaps_inside))),
    if (length(fits)) c("Ambiguity codes that fit (not counted):", cap(paste(pos[fits], "-", vapply(fits, bases, "")))),
    if (length(n_count)) paste0("N (not compared): ", paste(n_count, "in", names(n_count), collapse = ", ")))
  if (length(details)) {
    lines <- c(lines, paste0("Positions in ", q, ":"), details,
               if (q %in% reverse_complemented) paste("Positions are on the reverse complement of", q))
  }
  list(text = paste(lines, collapse = "\n"),
       flagged = length(missing) + length(gaps_inside) + length(fits) + length(n_count) > 0)
}
