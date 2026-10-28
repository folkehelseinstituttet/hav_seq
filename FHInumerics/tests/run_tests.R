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

# SNP counting with ambiguity codes (R/snps.R): one sequence per case against "ref".
toy <- c(ref     = "ACGTACGTAC",
         snp     = "TCGTACGTAC",   # A/T: differs
         M_fits  = "MCGTACGTAC",   # M = A/C can be A: same
         K_clash = "KCGTACGTAC",   # K = G/T cannot be A: differs
         N       = "NCGTACGTAC",   # N: never counted
         gap     = "-CGTACGTAC",   # gap: never counted
         R_Y     = "ACGTRCGTAC",   # R = A/G vs A: same (counts against Y below)
         Y       = "ACGTYCGTAC")   # Y = C/T vs A: differs; vs R: differs
toy_snp <- snp_matrix(toy)
check("SNPs with ambiguity codes: A/T, K/A and Y/A differ; M/A, N, gap and R/A do not",
      identical(unname(toy_snp["ref", -1]), c(1L, 0L, 1L, 0L, 0L, 0L, 1L)) && toy_snp["R_Y", "Y"] == 1)
check("compared positions: N and gaps are not compared, ambiguity codes are",
      identical(unname(comparable_sites(toy)["ref", c("N", "gap", "M_fits")]), c(9L, 9L, 10L)))
check("SNPs with ambiguity codes: the positions rule gives the same counts as phangorn", {
  bits <- base_bits(toy)
  all(outer(names(toy), names(toy), Vectorize(function(a, b) length(snp_columns(bits, a, b)))) == toy_snp)
})

# Hover box: SNPs, gaps (missing start, gap inside each sequence), an ambiguity code that
# fits and an N, numbered along the query q.
hover_aln <- c(q = "ACGTACG-TACG",
               o = "--GT-MGATNCT",
               p = "ACGTACG-TACG")
hover_view <- list(aln = hover_aln, queries = "q", snp = snp_matrix(hover_aln), sites = comparable_sites(hover_aln),
                   positions = list(q = cumsum(strsplit(hover_aln[["q"]], "")[[1]] != "-")))
hover <- snp_tips(hover_view)
check("hover box: SNP, missing start, gap in either sequence, ambiguity code that fits, N",
      identical(hover$text[["q\to"]], paste(
        "q / o", "1 SNP over 7 compared positions", "Positions in q:", "SNPs:", "11 - G/T",
        "Gaps (not counted):", "1-2 - missing at the start of o", "5 - A/-", "after 7 - -/A",
        "Ambiguity codes that fit (not counted):", "6 - C/M", "N (not compared): 1 in o", sep = "\n")))
check("hover box: bases follow row/column order; marked pairs only where gaps/N/ambiguity occur",
      grepl("\n11 - T/G\n", hover$text[["o\tq"]]) && hover$flagged[["q\to"]] && !hover$flagged[["q\tp"]] &&
        identical(hover$text[["q\tp"]], "q / p\n0 SNPs over 11 compared positions"))

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
      res$views[[ALL_QUERIES]]$snp["test_AY649695_1_copy", "AY649695.1"] == 0)
check("two queries: a result for both together and one for each alone",
      identical(names(res$views), c(ALL_QUERIES, "test_KY782331_3snps", "test_AY649695_1_copy")) &&
        file.exists(file.path(res$out_dir, "per_query", "02_test_AY649695_1_copy", "snp_matrix.tsv")))
own <- res$views[["test_AY649695_1_copy"]]
check("one query alone: only that query, compared over at least as many positions",
      identical(own$queries, "test_AY649695_1_copy") && own$snp["test_AY649695_1_copy", "AY649695.1"] == 0 &&
        own$compared_columns >= res$views[[ALL_QUERIES]]$compared_columns)

res <- run("genome_revcomp.fasta")
check("whole genome: detected as reverse complement",
      "test_LC128713_1_rc_2snps" %in% unlist(res$manifest$reverse_complemented))
check("whole genome: 2 SNPs to LC128713.1",
      res$views[[1]]$snp["test_LC128713_1_rc_2snps", "LC128713.1"] == 2)
check("one query: a single result, in the run folder",
      length(res$views) == 1 && !dir.exists(file.path(res$out_dir, "per_query")))

tip <- snp_tips(res$views[[1]], unlist(res$manifest$reverse_complemented))$text[["test_LC128713_1_rc_2snps\tLC128713.1"]]
snp_lines <- sub(".*\nSNPs:\n(([0-9]+ - [^\n]*\n?)*).*", "\\1", tip)
check("whole genome: SNP hover box lists the 2 substitutions (positions 3001 and 3101)",
      grepl("2 SNPs over", tip) &&
        identical(regmatches(snp_lines, gregexpr("[0-9]+(?= - )", snp_lines, perl = TRUE))[[1]], c("3001", "3101")))
check("whole genome: every SNP hover box lists as many positions as the SNP matrix counts", {
  v <- res$views[[1]]
  bits <- base_bits(v$aln)
  all(vapply(rownames(v$snp), function(a) all(vapply(colnames(v$snp), function(b)
    length(snp_columns(bits, a, b)) == v$snp[a, b], TRUE)), TRUE))
})
check("whole genome: hover box of the query names its closest sequence",
      grepl("LC128713.1", tip_details(res$views[[1]], "test_LC128713_1_rc_2snps")[["Closest in tree"]]))
check("tree plot is a plain ggplot (otherwise Shiny reports wrong hover positions)",
      !inherits(plot_tree(res$views[[1]]), "ggtree"))
check("hover box of a database sample shows its metadata",
      tip_details(res$views[[1]], "LC128713.1")[["Sample"]] == "LC128713.1" &&
        "country" %in% names(tip_details(res$views[[1]], "LC128713.1")))

# A database sample sent from the Database viewer: its own entry is left out.
ref <- load_db_seqs(db, cfg)
self_query <- ref$seqs["AY649695.1"]
res <- run_compare(self_query, db, cfg, settings, out_dir = file.path(out, "self"), log = function(...) NULL)
check("database sample as query: finds itself without exclude_self",
      res$hits$subject[1] == "AY649695.1" && res$hits$pident[1] == 100)
check("query named like a database sequence: renamed, both in the SNP matrix, 0 SNPs apart",
      identical(res$queries, "AY649695.1_query") && res$views[[1]]$snp["AY649695.1_query", "AY649695.1"] == 0 &&
        !anyDuplicated(rownames(res$views[[1]]$snp)))
res <- run_compare(self_query, db, cfg, utils::modifyList(settings, list(exclude_self = TRUE)),
                   out_dir = file.path(out, "self_excluded"), log = function(...) NULL)
check("database sample as query: own entry left out with exclude_self",
      !"AY649695.1" %in% res$hits$subject && isTRUE(res$manifest$settings$exclude_self))

# Two fragments of different genome regions: no result together, one for each alone.
apart <- c(region_a = substr(ref$seqs[["LC128713.1"]], 1001, 1500),
           region_b = substr(ref$seqs[["LC128713.1"]], 5001, 5500))
res <- run_compare(apart, db, cfg, settings, out_dir = file.path(out, "apart"), log = function(...) NULL)
check("queries that do not overlap: no result together, a result for each alone",
      grepl("do not overlap", res$views[[ALL_QUERIES]]$error) &&
        all(vapply(res$views[c("region_a", "region_b")], function(v) is.matrix(v$snp), TRUE)))

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
