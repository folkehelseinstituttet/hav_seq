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

# Alignment viewer (R/msa.R): every cell class against an anchor sequence.
msa_toy <- c(anchor = "ACGTRACGT-ACGA",
             x      = "--CTAANKTA-CGA",   # end, end, mismatch, ., amb.fit vs R, ., N, K fits G, ., base vs anchor gap, gap, ., ., .
             y      = "ACGTCACG-TACG-",   # C vs R cannot match; gap, base vs anchor gap (one indel event); missing end
             z      = "ACGTAACGT-ACGA")   # gap where the anchor has a gap: same
toy_m <- msa_matrix(msa_toy)
toy_cls <- matrix(as.integer(msa_classes(toy_m, toy_m["anchor", ])), 4)
check("alignment viewer: cell classes (same, SNP, ambiguity fits/clashes, indel both ways, N, missing end)",
      identical(toy_cls[2, ], c(6L, 6L, 1L, 0L, 2L, 0L, 5L, 2L, 0L, 4L, 4L, 0L, 0L, 0L)) &&
        identical(toy_cls[3, ], c(0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 4L, 4L, 0L, 0L, 0L, 6L)) &&
        identical(toy_cls[4, ], c(0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)) &&
        identical(msa_anchor_classes(toy_m["anchor", ]), c(8L, 8L, 8L, 8L, 2L, 8L, 8L, 8L, 8L, 4L, 8L, 8L, 8L, 8L)))
check("alignment viewer: mismatch + strong ambiguity cells = SNP columns, for every anchor and row", {
  bits <- base_bits(msa_toy)
  all(vapply(names(msa_toy), function(a) {
    cls <- msa_classes(toy_m, toy_m[a, ])
    all(vapply(names(msa_toy), function(b) identical(which(cls[b, ] %in% as.raw(c(1, 3))), snp_columns(bits, a, b)), TRUE))
  }, TRUE))
})
frag_anchor <- msa_matrix(c(a = "--ACGT--", b = "TTACGTAA"))
check("alignment viewer: a fragment anchor has nothing to compare outside its own length (grey)",
      identical(as.integer(msa_classes(frag_anchor, frag_anchor["a", ])["b", ]), c(7L, 7L, 0L, 0L, 0L, 0L, 7L, 7L)))
toy_view <- msa_view(msa_load(msa_toy, "toy"), "anchor")
check("alignment viewer: differences to the anchor, variable columns, anchor gaps, rows sorted by differences",
      identical(unname(toy_view$diffs), c(0, 1, 1, 0)) &&
        identical(which(toy_view$variable), c(3L, 5L, 8L, 9L, 10L, 11L)) && identical(which(toy_view$anchor_gap), 10L) &&
        identical(rownames(toy_m)[toy_view$rows], c("z", "x", "y")))
check("alignment viewer: per row indel events (a run counts once) and columns, ambiguity codes, no data; per column differing rows",
      identical(toy_view$stats$indel_events, c(0, 1, 1, 0)) && identical(toy_view$stats$indel_cols, c(0, 2, 2, 0)) &&
        identical(toy_view$stats$ambiguity, c(1, 2, 1, 1)) && identical(toy_view$stats$nodata, c(0, 3, 1, 0)) &&
        identical(toy_view$ndiff, c(0, 0, 1, 0, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0)))

# Consensus: majority with split votes; tie -> ambiguity code (or N); gap inside most
# sequences -> "-"; no votes -> N; threshold.
cons_aln <- c(s1 = "AAGNA-AA", s2 = "AGGNA-CA", s3 = "CGRNC-G-", s4 = "CARN-AT-")
# col: 1 A/C tie, 2 A/G tie, 3 G 3 votes (R splits) vs A 1, 4 no votes, 5 A 2 vs C 1 (+ a gap),
#      6 gap inside 3 of 4, 7 four-way tie, 8 A (missing ends do not vote)
cons <- function(...) rawToChar(msa_consensus(msa_matrix(cons_aln), ...))
check("consensus: majority, split votes for ambiguity codes, ties as ambiguity code or N, gaps, no coverage",
      identical(cons(), "MRGNA-NA") && identical(cons(ties = "n"), "NNGNA-NA"))
check("consensus: with a threshold the most common base must reach it, else the code of the top bases",
      identical(cons(threshold = 75), "MRGNM-NA") && identical(cons(threshold = 50), "MRGNA-NA") &&
        identical(cons(threshold = 80, ties = "n"), "NNNNN-NA"))
check("consensus: the default anchor (no anchor sequence chosen)",
      identical(msa_view(msa_load(cons_aln, "c"))$anchor_name, MSA_CONSENSUS))

# Coordinates: own base numbers, "gap after", the column of a base number.
check("alignment viewer: base numbers along a sequence and the column of a base number",
      identical(msa_positions(charToRaw("--AC-GT")), c(0L, 0L, 1L, 2L, 2L, 3L, 4L)) &&
        msa_column_of(charToRaw("--AC-GT"), 3) == 6L && is.na(msa_column_of(charToRaw("--AC"), 5)))

# Subsets keep the original column numbers; closest N includes ties.
sub_a <- msa_subset(msa_load(c(p = "A--CA", q = "A-GCA", r = "AT-CA"), "s", compared = 2:4), c("p", "q"))
check("alignment viewer: a subset drops columns that are a gap in every chosen row, keeps column numbers",
      identical(sub_a$cols, c(1L, 3L, 4L, 5L)) && identical(sub_a$compared, c(FALSE, TRUE, TRUE, FALSE)))
snp_toy <- matrix(c(0, 1, 2, 2, 5,  1, 0, 0, 0, 0,  2, 0, 0, 0, 0,  2, 0, 0, 0, 0,  5, 0, 0, 0, 0), 5,
                  dimnames = list(c("q", "a", "b", "c", "d"), c("q", "a", "b", "c", "d")))
check("closest N: ties at the N-th place are all included",
      identical(msa_closest(snp_toy, "q", 2), list(ids = c("q", "a", "b", "c"), extra = 1L)) &&
        identical(msa_closest(snp_toy, "q", 1)$ids, c("q", "a")))

res <- run("fragments.fasta")
frag <- res
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

# Alignment viewer on the real runs: the untrimmed alignment is kept and saved; with any
# sequence as anchor, its differences over the compared region = the SNP table, and the
# mismatch + strong ambiguity cells = snp_columns() over the whole alignment.
msa_invariant <- function(v) {
  a <- msa_from_result(v)
  bits <- base_bits(v$aln_full)
  all(vapply(setdiff(rownames(a$m), a$dropped), function(anchor) {
    mv <- msa_view(a, anchor)
    rows <- setdiff(rownames(a$m), a$dropped)
    all(mv$diffs[rows] == v$snp[anchor, rows]) &&
      all(vapply(rownames(a$m), function(b) identical(which(mv$classes[b, ] %in% as.raw(c(1, 3))),
                                                      snp_columns(bits, anchor, b)), TRUE))
  }, TRUE))
}
check("alignment viewer: for every anchor and row of every result, differences = SNP table and cells = snp_columns()",
      all(vapply(Filter(function(v) is.null(v$error), c(frag$views, res$views)), msa_invariant, TRUE)))
gv <- res$views[[1]]
check("run folder: untrimmed alignment saved with real names, compared region in the manifest",
      file.exists(file.path(res$out_dir, "alignment_untrimmed.fasta")) &&
        setequal(names(read_alignment(file.path(res$out_dir, "alignment_untrimmed.fasta"))), names(gv$aln_full)) &&
        identical(res$manifest$results[[1]]$compared_region, sprintf("%d-%d", min(gv$compared_cols), max(gv$compared_cols))))
ga <- msa_from_result(gv, unlist(res$manifest$reverse_complemented))
q <- "test_LC128713_1_rc_2snps"
gq <- msa_view(ga, q)
diff_cols <- which(gq$classes["LC128713.1", ] %in% as.raw(c(1, 3)))
check("alignment viewer: reverse-complemented query numbered on its reverse complement, as in the SNP hover box (3001, 3101)",
      q %in% ga$revcomp && identical(msa_positions(ga$m[q, ])[diff_cols], c(3001L, 3101L)) &&
        msa_first_difference(ga, gq, "LC128713.1") == diff_cols[1] &&
        msa_column_of(ga$m[q, ], 3001) == diff_cols[1])
check("alignment viewer: rows in tree order, as the SNP table",
      identical(rownames(ga$m)[seq_along(gv$tree$tip.label)], tree_order(gv$tree)))
folders <- msa_run_results(res$out_dir)
fa <- msa_load_run_folder(folders[[1]], names(folders)[1], attr(folders, "manifest"), attr(folders, "regions")[[1]])
check("earlier run folder: opens the same alignment, compared region, queries and order",
      identical(fa$m, ga$m) && identical(fa$compared, ga$compared) && identical(fa$queries, ga$queries) &&
        identical(fa$revcomp, ga$revcomp))
pdf_file <- file.path(out, "window.pdf")
msa_pdf(ga, gq, gq$rows[1:5], diff_cols[1] + (-20:20), pdf_file)
msa_pdf(ga, gq, gq$rows, seq_len(ncol(ga$m)), file.path(out, "whole.pdf"))
check("alignment viewer: PDF of a window (letters) and of the whole alignment", all(file.size(c(pdf_file, file.path(out, "whole.pdf"))) > 1000))
html_file <- file.path(out, "alignment.html")
msa_html(ga, gq, html_file, sort = "name")
html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
check("alignment viewer: standalone HTML page with the viewer, the data and the current anchor",
      grepl("fhiMsa.load(", html, fixed = TRUE) && grepl("function drawOverview", html, fixed = TRUE) &&
        grepl(q, html, fixed = TRUE) && lengths(regmatches(html, gregexpr("</script>", html, fixed = TRUE))) == 2)
check("alignment viewer: FASTA of the shown region",
      identical(nchar(msa_fasta(ga, 1:2, 101, 200)), stats::setNames(c(100L, 100L), rownames(ga$m)[1:2])))

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
