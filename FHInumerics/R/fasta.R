# FASTA reading, writing and validation of user input.

IUPAC_NT <- "ACGTRYKMSWBDHVN"

# Parse FASTA text as pasted by a user (NCBI BLAST style).
# - Text without any ">" header is treated as one unnamed sequence.
# - Whitespace, digits and gap characters are removed (GenBank flat-file copy-paste).
# - Returns a named character vector of upper-case sequences, or stops with a
#   message meant for the user.
parse_fasta_text <- function(text, prefix = "query") {
  lines <- unlist(strsplit(gsub("\r", "", text), "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) stop("No sequence was given.", call. = FALSE)
  if (!startsWith(lines[1], ">")) lines <- c(">", lines)

  is_header <- startsWith(lines, ">")
  record <- cumsum(is_header)
  headers <- sub("^>\\s*", "", lines[is_header])
  names <- vapply(strsplit(headers, "\\s+"), function(x) if (length(x)) x[1] else "", "")
  seqs <- vapply(split(lines[!is_header], record[!is_header])[as.character(seq_along(names))],
                 function(x) if (is.null(x)) "" else paste(x, collapse = ""), "")

  unnamed <- !nzchar(names)
  names[unnamed] <- paste0(prefix, "_", which(unnamed))

  seqs <- toupper(gsub("[[:space:][:digit:]./-]", "", seqs))
  seqs <- chartr("U", "T", seqs)
  names(seqs) <- names

  empty <- names[!nzchar(seqs)]
  if (length(empty)) stop("These entries have no sequence: ", paste(empty, collapse = ", "), call. = FALSE)
  bad <- vapply(seqs, function(s) gsub(paste0("[", IUPAC_NT, "]"), "", s), "")
  if (any(nzchar(bad))) {
    which_bad <- names(bad)[nzchar(bad)]
    chars <- unique(strsplit(paste(bad[nzchar(bad)], collapse = ""), "")[[1]])
    stop("Not a nucleotide sequence: ", paste(which_bad, collapse = ", "),
         " (unexpected characters: ", paste(head(chars, 10), collapse = " "), ")", call. = FALSE)
  }
  dups <- unique(names[duplicated(names)])
  if (length(dups)) stop("Duplicate sequence names: ", paste(dups, collapse = ", "), call. = FALSE)
  seqs
}

read_fasta <- function(path) {
  parse_fasta_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

write_fasta <- function(seqs, path, width = 70) {
  out <- character(0)
  for (i in seq_along(seqs)) {
    s <- seqs[[i]]
    starts <- seq(1, max(nchar(s), 1), by = width)
    out <- c(out, paste0(">", names(seqs)[i]), substring(s, starts, pmin(starts + width - 1, nchar(s))))
  }
  writeLines(out, path)
  invisible(path)
}

# Aligned FASTA (keeps gap characters) -> named character vector.
read_alignment <- function(path) {
  lines <- readLines(path, warn = FALSE)
  is_header <- startsWith(lines, ">")
  record <- cumsum(is_header)
  seqs <- vapply(split(lines[!is_header], record[!is_header]), paste, "", collapse = "")
  names(seqs) <- sub("^>", "", lines[is_header])
  toupper(seqs)
}

reverse_complement <- function(s) {
  paste(rev(strsplit(chartr("ACGTRYKMBDHV", "TGCAYRMKVHDB", s), "")[[1]]), collapse = "")
}
