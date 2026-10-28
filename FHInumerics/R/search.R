# Database search: find database sequences similar to each query.
# Both methods return one row per query/database sequence pair, with BLAST-style
# columns (as blastn -outfmt 6 with max_hsps 1):
#   qseqid sseqid pident length mismatch gapopen qstart qend sstart send
#   score qlen slen qcovhsp sstrand
# pident = identical positions / alignment length (gaps included), as BLAST.

HIT_FIELDS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "score", "qlen", "slen", "qcovhsp", "sstrand")
MAX_HITS <- 500   # per query, as blastn -max_target_seqs 500

search_db <- function(queries, ref, method = "decipher", threads = 1) {
  hits <- switch(method,
    decipher = search_decipher(queries, ref, threads),
    blast = search_blast(queries, ref, threads),
    stop("Unknown search method '", method, "' (use decipher or blast).", call. = FALSE))
  if (!NROW(hits)) hits <- data.frame(matrix(ncol = length(HIT_FIELDS), nrow = 0, dimnames = list(NULL, HIT_FIELDS)))
  hits[, HIT_FIELDS]
}

# DECIPHER: k-mer index search (SearchIndex), then pairwise alignment of each hit
# (AlignPairs). Both strands are searched by also searching the reverse complements.
search_decipher <- function(queries, ref, threads) {
  n <- length(queries)
  pattern <- Biostrings::DNAStringSet(queries)
  pattern <- c(pattern, Biostrings::reverseComplement(pattern))
  subject <- db_stringset(ref)
  found <- DECIPHER::SearchIndex(pattern, ref$index, subject, perPatternLimit = MAX_HITS,
                                 processors = threads, verbose = FALSE)
  if (!nrow(found)) return(NULL)
  al <- DECIPHER::AlignPairs(pattern, subject, found, processors = threads, verbose = FALSE)

  minus <- al$Pattern > n
  q <- ifelse(minus, al$Pattern - n, al$Pattern)
  qlen <- nchar(queries)[q]
  hits <- data.frame(
    qseqid = names(queries)[q], sseqid = names(subject)[al$Subject],
    pident = round(100 * al$Matches / al$AlignmentLength, 3), length = al$AlignmentLength,
    mismatch = al$Mismatches,
    gapopen = lengths(al$PatternGapPosition) + lengths(al$SubjectGapPosition),
    qstart = ifelse(minus, qlen - al$PatternEnd + 1, al$PatternStart),
    qend = ifelse(minus, qlen - al$PatternStart + 1, al$PatternEnd),
    sstart = ifelse(minus, al$SubjectEnd, al$SubjectStart),
    send = ifelse(minus, al$SubjectStart, al$SubjectEnd),
    score = al$Score, qlen = qlen, slen = Biostrings::width(subject)[al$Subject],
    qcovhsp = round(100 * (al$PatternEnd - al$PatternStart + 1) / qlen),
    sstrand = ifelse(minus, "minus", "plus")
  )
  # Keep the better strand for each query/database pair, and at most MAX_HITS per query.
  hits <- hits[order(hits$qseqid, -hits$score), ]
  hits <- hits[!duplicated(hits[c("qseqid", "sseqid")]), ]
  do.call(rbind, lapply(split(hits, hits$qseqid), utils::head, MAX_HITS))
}

# rBLAST (BLAST+ installed on this computer): blastn, one HSP per pair.
search_blast <- function(queries, ref, threads) {
  fields <- sub("^score$", "bitscore", HIT_FIELDS)
  bl <- rBLAST::blast(db = ref$prefix, type = "blastn")
  hits <- stats::predict(bl, Biostrings::DNAStringSet(queries), custom_format = paste(fields, collapse = " "),
                         BLAST_args = paste("-max_target_seqs", MAX_HITS, "-max_hsps 1 -num_threads", threads))
  if (!nrow(hits)) return(NULL)
  names(hits)[names(hits) == "bitscore"] <- "score"
  hits
}
