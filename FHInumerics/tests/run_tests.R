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

check("whole genome: hover box of the query names its closest sequence",
      grepl("LC128713.1", tip_details(res, "test_LC128713_1_rc_2snps")[["Closest in tree"]]))
check("tree plot is a plain ggplot (otherwise Shiny reports wrong hover positions)",
      !inherits(plot_tree(res), "ggtree"))
check("hover box of a database sample shows its metadata",
      tip_details(res, "LC128713.1")[["Sample"]] == "LC128713.1" && "country" %in% names(tip_details(res, "LC128713.1")))

# A database sample sent from the Database viewer: its own entry is left out.
ref <- load_db_seqs(db, cfg)
self_query <- ref$seqs["AY649695.1"]
res <- run_compare(self_query, db, cfg, settings, out_dir = file.path(out, "self"), log = function(...) NULL)
check("database sample as query: finds itself without exclude_self",
      res$hits$subject[1] == "AY649695.1" && res$hits$pident[1] == 100)
res <- run_compare(self_query, db, cfg, utils::modifyList(settings, list(exclude_self = TRUE)),
                   out_dir = file.path(out, "self_excluded"), log = function(...) NULL)
check("database sample as query: own entry left out with exclude_self",
      !"AY649695.1" %in% res$hits$subject && isTRUE(res$manifest$settings$exclude_self))

# Database viewer helpers.
check("partial dates: year and month count as their first day",
      identical(parse_partial_date(c("1985", "2019-04", "2014-04-30", "", "2023-13-01")),
                as.Date(c("1985-01-01", "2019-04-01", "2014-04-30", NA, NA))))
ids <- match_sample_ids(">AY649695.1\nnot_an_id, LC128713.1;", ref$ids$id)
check("pasted sample IDs: found and missing",
      identical(ids$found, c("AY649695.1", "LC128713.1")) && identical(ids$missing, "not_an_id"))
meta <- read_metadata(db)
check("date columns found in the demo metadata: date and year",
      identical(date_columns(meta), c("date", "year")))
map <- plot_country_map(meta, "country", split_by = "clade")
check("map: every demo country is placed on the map", length(map$unknown) == 0)
japan <- circle_details(map$circles[map$circles$region == "Japan", ], "clade")
check("map hover box: Japan has 365 samples, 318 of clade IA",
      japan[["Samples"]] == "365" && japan[["clade IA"]] == "318")
hover <- list(domain = list(left = 0, right = 10, bottom = 0, top = 10),
              range = list(left = 0, right = 100, bottom = 100, top = 0), coords_css = list(x = 52, y = 50))
check("hover: a point 2 pixels from the pointer is found, one 20 pixels away is not",
      point_at(5, 5, hover)$i == 1 && is.null(point_at(7, 5, hover)))

cat(if (failures) paste(failures, "check(s) FAILED") else "All checks passed", "\n")
quit(save = "no", status = if (failures) 1 else 0)
