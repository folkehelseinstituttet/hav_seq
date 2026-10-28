# Run the test queries (see tests/README.md) and check the known answers.
#   Rscript tests/run_tests.R [--aligner decipher|mafft]
# Uses only the demo database from config.example.yaml (never your config.yaml),
# writes results to a temporary folder, and prints PASS/FAIL per check.
# Exit status 0 = all checks passed. Run it after every change.

args <- commandArgs(trailingOnly = TRUE)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(file_arg), ".."))
options(fhinumerics.root = root)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)

cfg <- load_config(file.path(root, "config.example.yaml"))
out <- tempfile("fhinumerics-tests-")
cfg$output_dir <- out
cfg$cache_dir <- file.path(out, "cache")
db <- cfg$databases[[1]]
settings <- cfg$settings
if (length(args) == 2 && args[1] == "--aligner") settings$aligner <- args[2]
message("Aligner: ", settings$aligner, "; results in ", out)

failures <- 0
check <- function(what, ok) {
  ok <- tryCatch(ok, error = function(e) FALSE)   # a check that cannot be evaluated fails
  cat(if (isTRUE(ok)) "PASS" else "FAIL", " ", what, "\n", sep = "")
  if (!isTRUE(ok)) failures <<- failures + 1
}
run <- function(fasta) {
  run_compare(read_fasta(file.path(root, "tests", fasta)), db, cfg, settings,
              out_dir = file.path(out, sub("\\.fasta$", "", fasta)), log = function(...) NULL)
}

res <- run("fragments.fasta")
# KY782331.1 is one of the 101 equally close sequences, not necessarily one of the
# 20 in the tree, so it is checked in the search hits.
hit <- res$hits[res$hits$query == "test_KY782331_3snps" & res$hits$subject == "KY782331.1", ]
check("fragment with 3 substitutions: 3 mismatches to KY782331.1", hit$mismatch[1] == 3)
check("fragment with 3 substitutions: closest sequences have 3 SNPs",
      res$ties$best_mismatches[res$ties$query == "test_KY782331_3snps"] == 3)
check("fragment with 3 substitutions: 101 database sequences equally close",
      res$ties$n_equally_close[res$ties$query == "test_KY782331_3snps"] == 101)
check("exact copy: 0 SNPs to AY649695.1",
      res$snp["test_AY649695_1_copy", "AY649695.1"] == 0)

res <- run("genome_revcomp.fasta")
check("whole genome: detected as reverse complement",
      "test_LC128713_1_rc_2snps" %in% unlist(res$manifest$reverse_complemented))
check("whole genome: 2 SNPs to LC128713.1",
      res$snp["test_LC128713_1_rc_2snps", "LC128713.1"] == 2)

cat(if (failures) paste(failures, "check(s) FAILED") else "All checks passed", "\n")
quit(save = "no", status = if (failures) 1 else 0)
