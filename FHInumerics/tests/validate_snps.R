# Validate the SNP count (R/snps.R) on a whole database. Slower than run_tests.R; run it
# when the SNP method changes or for a new database release.
#   Rscript tests/validate_snps.R <database.fasta> [number of other sequences, default 150]
# Aligns every sequence with an ambiguity code plus a random sample of the others
# (DECIPHER), then checks for every pair:
#   1. phangorn's count (snp_matrix) = an independent count from Biostrings' IUPAC table
#   2. phangorn's count = the old count (ape, A/C/G/T only) + differences that involve
#      an ambiguity code, so pairs without ambiguity codes are unchanged
# and prints how many pairs change and by how much.

file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(file_arg), ".."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript tests/validate_snps.R <database.fasta> [n_other]")
n_other <- if (length(args) > 1) as.integer(args[2]) else 150

seqs <- read_fasta(args[1])
amb <- grepl("[RYSWKMBDHV]", seqs)
codes <- table(unlist(strsplit(gsub("[ACGTN]", "", seqs), "")))
cat(sprintf("%d sequences; %d with ambiguity codes (%s), %d with N\n", length(seqs), sum(amb),
            paste(names(codes), codes, sep = "=", collapse = " "), sum(grepl("N", seqs))))
set.seed(1)
pick <- c(names(seqs)[amb], sample(names(seqs)[!amb], min(n_other, sum(!amb))))
aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs[pick]), verbose = FALSE)
aln <- stats::setNames(as.character(aln), names(aln))
cat(sprintf("Aligned %d sequences, %d columns\n", length(aln), nchar(aln[1])))

new <- snp_matrix(aln)
dna <- ape::as.DNAbin(lapply(strsplit(tolower(aln), ""), identity))
old <- round(as.matrix(ape::dist.dna(dna, "N", pairwise.deletion = TRUE))[names(aln), names(aln)])

# Independent count: two codes differ when their sets of bases (Biostrings) do not overlap.
map <- c(Biostrings::IUPAC_CODE_MAP, `-` = "ACGT")
sets <- lapply(map, function(x) strsplit(x, "")[[1]])
fits <- outer(names(map), names(map), Vectorize(function(a, b) length(intersect(sets[[a]], sets[[b]])) > 0))
chars <- do.call(rbind, strsplit(aln, ""))
code <- matrix(match(chars, names(map)), nrow = nrow(chars))
acgt <- matrix(chars %in% c("A", "C", "G", "T"), nrow = nrow(chars))
n <- length(aln)
indep <- ambiguous <- matrix(0L, n, n)
for (i in seq_len(n - 1)) for (j in (i + 1):n) {
  differ <- !fits[cbind(code[i, ], code[j, ])]
  indep[i, j] <- indep[j, i] <- sum(differ)
  ambiguous[i, j] <- ambiguous[j, i] <- sum(differ & !(acgt[i, ] & acgt[j, ]))
}

ok <- c(all(new == indep), all(new == old + ambiguous))
cat(if (ok[1]) "PASS" else "FAIL", " phangorn count = independent count (all pairs)\n", sep = "")
cat(if (ok[2]) "PASS" else "FAIL", " phangorn count = A/C/G/T-only count + ambiguity-code differences\n", sep = "")
change <- (new - old)[upper.tri(new)]
cat(sprintf("Pairs: %d; changed: %d (%.1f%%); by: %s\n", length(change), sum(change > 0), 100 * mean(change > 0),
            paste0("+", names(table(change)), " (", table(change), ")", collapse = ", ")))
quit(save = "no", status = if (all(ok)) 0 else 1)
