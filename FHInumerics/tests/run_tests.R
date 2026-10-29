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
genome <- res
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

# Sets (R/sets.R, decision 0012). A second database (3 demo samples) for name clashes.
db2_dir <- file.path(out, "db2")
dir.create(db2_dir)
ids2 <- c("AY649695.1", "LC128713.1", "KY782331.1")
write_fasta(ref$seqs[ids2], file.path(db2_dir, "seqs.fasta"))
utils::write.table(read_metadata(db)[match(ids2, read_metadata(db)$id), c("id", "country")],
                   file.path(db2_dir, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cfg$databases[["Second db"]] <- list(name = "Second db", fasta = file.path(db2_dir, "seqs.fasta"),
                                     metadata = file.path(db2_dir, "meta.tsv"), id_column = "id")
db2 <- cfg$databases[["Second db"]]
meta1 <- read_metadata(db)
st <- set_new("Test set")
st <- set_add(st, set_from_database(db, cfg, c("AY649695.1", "AB279735.1"), meta1))$set
r <- set_add(st, set_from_database(db2, cfg, c("AY649695.1", "LC128713.1"), read_metadata(db2)))
st <- r$set
pasted <- set_candidates(c("AY649695.1", "q1", "q1", "q1"), c(ref$seqs[["AY649695.1"]], "ACGTACGT", "ACGTACGT", "ACGTTCGT"),
                         "query", "pasted")
r2 <- set_add(st, pasted)
st <- r2$set
again <- set_add(st, set_from_database(db, cfg, c("AY649695.1", "AB279735.1", "LC128713.1"), meta1))
check("sets: a name in two databases becomes <name>_<database> for both (short_name, else the name made safe)",
      identical(set_names(r$set$members), c("AY649695.1_demo", "AB279735.1", "AY649695.1_Second_db", "LC128713.1")) &&
        identical(unname(r$renamed), c("AY649695.1_demo", "AY649695.1_Second_db")))
check("sets: a query named like a database sample becomes <name>_query; another sequence with a used name <name>_2",
      identical(set_names(st$members)[5:7], c("AY649695.1_query", "q1", "q1_2")) &&
        identical(r2$already, "q1") && identical(r2$added, c("AY649695.1_query", "q1", "q1_2")))
check("sets: the same database sample is kept once (the same name from another database is a new member)",
      identical(again$already, c("AY649695.1", "AB279735.1")) && identical(again$added, "LC128713.1_demo") &&
        nrow(again$set$members) == 8 && grepl("2 were already in the set", set_add_message(again)))
check("sets: checksum = md5 of sorted 'name<TAB>SEQUENCE' lines (known answer)",
      set_checksum(c(b = "gg-t", a = "ACGT")) == "f0cc91e49490ef264ac6bd412834387e" &&
        set_checksum(c(a = "ACGT", b = "GGT")) == "f0cc91e49490ef264ac6bd412834387e")

# A set made from run results: queries as given (turned back), with a note; database
# samples with the run's database checksum; the set starts with the run's alignment.
ga_run <- msa_from_result(genome$views[[1]], unlist(genome$manifest$reverse_complemented))
ga_run$origin <- msa_run_origin(genome$manifest, basename(genome$out_dir), cfg$databases)
from_run <- set_add(set_new("From run"), set_from_msa(ga_run, added_from = "Run test"))$set
rq <- from_run$members$source == "query"
rc_query <- read_fasta(file.path(root, "tests", "genome_revcomp.fasta"))
check("sets from a run: the query is kept as given (not turned), with a note; samples carry the database checksum",
      identical(unname(from_run$seqs[from_run$members$key[rq]]), unname(rc_query)) &&
        grepl("reverse-complemented in the run", from_run$members$note[rq]) &&
        all(from_run$members$database_md5[!rq] == genome$manifest$database$fasta_md5))
fr_aln <- set_alignment(from_run)
check("sets from a run: start with the run's (untrimmed) alignment, the query turned in it",
      identical(unname(fr_aln$aln[rownames(ga_run$m)]), unname(genome$views[[1]]$aln_full[rownames(ga_run$m)])) &&
        identical(fr_aln$turned, q) && !length(fr_aln$missing))
fr_a <- msa_from_set(from_run)
fr_v <- msa_view(fr_a, q)
fr_cols <- which(fr_v$classes["LC128713.1", ] %in% as.raw(c(1, 3)))
check("sets from a run: in the Alignment page the query is numbered on its reverse complement (3001, 3101)",
      q %in% fr_a$revcomp && identical(msa_positions(fr_a$m[q, ])[fr_cols], c(3001L, 3101L)))
smaller <- set_remove(from_run, from_run$members$key[!rq][1:5])
check("sets: removing members keeps the alignment (their rows taken out); adding makes it out of date",
      nrow(msa_from_set(smaller)$m) == nrow(ga_run$m) - 5 && !length(set_alignment(smaller)$missing) &&
        identical(set_alignment(set_add(smaller, set_from_database(db, cfg, "AB279735.1", meta1))$set)$missing,
                  "AB279735.1"))

# Orientation and alignment: the reverse-complemented genome among database genomes is
# turned (and only it), and numbered on its reverse complement, as in runs.
genomes <- c("LC128713.1", "EU526089.1", "AB279735.1", "MH933766.1", "AY644670.1")
orient_set <- set_add(set_new("Orient"), set_from_database(db, cfg, genomes, meta1))$set
orient_set <- set_add(orient_set, set_candidates(names(rc_query), rc_query, "query", "pasted"))$set
orient_seqs <- stats::setNames(unname(orient_set$seqs), orient_set$members$key)
in_process <- set_align(orient_seqs, threads = 2)
check("set alignment: OrientNucleotides turns the reverse-complemented genome and nothing else",
      identical(in_process$turned, orient_set$members$key[orient_set$members$source == "query"]))
orient_set$alignment <- in_process[c("aln", "turned", "made", "from", "method", "versions")]
os_a <- msa_from_set(orient_set)
os_v <- msa_view(os_a, q)
os_cols <- which(os_v$classes["LC128713.1", ] %in% as.raw(c(1, 3)))
check("set alignment: the turned sequence is numbered on its reverse complement (SNPs at 3001 and 3101)",
      q %in% os_a$revcomp && identical(msa_positions(os_a$m[q, ])[os_cols], c(3001L, 3101L)))
job <- align_start(orient_seqs, threads = 2, dir = file.path(out, "align_job"))
repeat { polled <- align_poll(job); if (polled$done) break; Sys.sleep(0.5) }
check("background alignment: the same alignment and turned sequences as in-process",
      is.null(polled$error) && identical(polled$result$aln, in_process$aln) &&
        identical(polled$result$turned, in_process$turned) && polled$value %||% 1 == 1)
big <- ref$seqs[seq_len(400)]
job <- align_start(big, threads = 1, dir = file.path(out, "cancel_job"))
Sys.sleep(3)
check("background alignment: Cancel stops the process and leaves no result",
      job$proc$is_alive() && align_cancel(job) && !job$proc$is_alive() && !dir.exists(file.path(out, "cancel_job")))
check("background alignment: progress read from DECIPHER's output",
      isTRUE(all.equal(align_progress(paste0("STEP Orienting\nDetermining distance matrix based on shared 10-mers:\n",
                                             strrep("=", 40)))$value, 0.085)) &&
        isTRUE(all.equal(align_progress(paste0("Iteration 2 of 2:\n\nDetermining distance matrix based on alignment:\n",
                                               strrep("=", 80), "\n\nTime difference of 0 secs\n\n",
                                               "Reclustering into groups by similarity:\n", strrep("=", 40)))$value,
                         0.796875)) &&
        align_progress("")$text == "Starting")

# Save and open again: sequences, table, sources with database checksums, checksum, the
# alignment and its turned rows; Replace and Delete touch only the set's files.
orient_set <- set_add(orient_set, set_from_database(db2, cfg, "AY649695.1", read_metadata(db2)))$set
saved_dir <- file.path(out, "sets", set_folder_name("Orient: test?"))
saved <- set_save(orient_set, saved_dir)
loaded <- set_load(saved_dir)
yml <- yaml::read_yaml(file.path(saved_dir, "set.yaml"))
check("set save/open: same sequences, table, alignment and turned rows",
      identical(set_fasta(loaded), set_fasta(orient_set)) && identical(set_table(loaded), set_table(orient_set)) &&
        identical(set_alignment(loaded)$aln, set_alignment(orient_set)$aln) && identical(set_alignment(loaded)$turned, q) &&
        identical(set_alignment(loaded)$missing, "AY649695.1") && !loaded$dirty && basename(saved_dir) == "Orient_ test_")
check("set save/open: set.yaml has the sources with database checksums, the checksum and the alignment",
      identical(yml$checksum, set_checksum(set_fasta(orient_set))) && yml$n_sequences == 7 &&
        identical(vapply(yml$sources, function(s) s$database %||% s$type, ""), c(db$name, "Second db", "query")) &&
        identical(yml$sources[[1]]$database_md5, unname(tools::md5sum(db$fasta))) && yml$sources[[1]]$n == 5 &&
        identical(yml$sources[[2]]$database_md5, unname(tools::md5sum(db2$fasta))) &&
        identical(unlist(yml$alignment$reverse_complemented), q) && nrow(set_list_saved(file.path(out, "sets"))) == 1)
writeLines("mine", file.path(saved_dir, "notes.txt"))
check("sets: Delete removes only the set's files; a folder with other files is not a free place to save",
      identical(set_delete(saved_dir), "notes.txt") && identical(list.files(saved_dir), "notes.txt") &&
        set_folder_state(saved_dir) == "other" && set_folder_state(file.path(out, "sets", "new")) == "free")
aligned_file <- file.path(out, "aligned_in.fasta")
three <- msa_subset(ga_run, rownames(ga_run$m)[1:3])
write_fasta(msa_fasta(three, 1:3), aligned_file)
from_file <- set_add(set_new("F"), set_from_fasta(aligned_file, "aligned_in.fasta"))$set
unaligned <- set_from_fasta(file.path(root, "tests", "fragments.fasta"))
check("sets from a FASTA file: an aligned file brings its alignment, an unaligned one does not",
      identical(unname(set_alignment(from_file)$aln), unname(msa_fasta(three, 1:3))) && is.null(unaligned$aln) &&
        all(unaligned$members$source == "file"))

# A set as the queries of a run (0012): its sample from the searched database leaves out
# its own entry; a pasted query named like a database sample does not; a sample from
# another database does not; the manifest records the set. Same with the command line.
run_set <- set_add(set_new("Run set"), set_from_database(db, cfg, "AY649695.1", meta1))$set
run_set <- set_add(run_set, set_candidates("KY782331.1", ref$seqs[["KY782331.1"]], "query", "pasted"))$set
rq <- set_queries(run_set, db$name)
set_run <- run_compare(rq$seqs, db, cfg, utils::modifyList(settings, list(exclude_self = rq$exclude_self)),
                       out_dir = file.path(out, "set_run"), log = function(...) NULL, input = rq$input)
set_man <- yaml::read_yaml(file.path(set_run$out_dir, "run_manifest.yaml"))
hit_of <- function(q, subj) set_run$hits[set_run$hits$query == q & set_run$hits$subject == subj, ]
check("set as queries: the database sample leaves out its own entry, the pasted query still finds its namesake (100%)",
      identical(names(rq$seqs), c("AY649695.1", "KY782331.1")) && identical(rq$exclude_self, list(AY649695.1 = "AY649695.1")) &&
        !nrow(hit_of("AY649695.1_query", "AY649695.1")) && isTRUE(hit_of("KY782331.1_query", "KY782331.1")$pident == 100))
check("set as queries: the manifest records the set (name, checksum) and whose entry was left out; other runs have no input",
      identical(set_man$input$set, "Run set") && identical(set_man$input$set_checksum, set_checksum(set_fasta(run_set))) &&
        identical(set_man$settings$exclude_self, list(AY649695.1 = "AY649695.1")) &&
        is.null(yaml::read_yaml(file.path(genome$out_dir, "run_manifest.yaml"))$input))
other_db <- set_add(set_new("Other"), set_from_database(db2, cfg, "AY649695.1", read_metadata(db2)))$set
big_set <- set_add(set_new("Big"), set_from_database(db, cfg, ref$ids$id[1:60], meta1))$set
check("set as queries: a sample from another database leaves out nothing; more than 50 -> plain message; chosen rows",
      isFALSE(set_queries(other_db, db$name)$exclude_self) &&
        grepl("has 60 sequences; Outbreak detection takes at most 50", tryCatch(set_queries(big_set, db$name), error = conditionMessage)) &&
        length(set_queries(big_set, db$name, big_set$members$key[1:10])$seqs) == 10)
cli_cfg <- file.path(out, "cli_config.yaml")
yaml::write_yaml(list(output_dir = file.path(out, "cli_results"), cache_dir = cfg$cache_dir,
                      databases = list(list(name = db$name, fasta = db$fasta, metadata = db$metadata, id_column = "id"))), cli_cfg)
invisible(set_save(run_set, file.path(out, "sets", "Run set")))
cli <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                c(file.path(root, "cli", "fhinumerics"), "--query-set", shQuote(file.path(out, "sets", "Run set")),
                                  "--config", shQuote(cli_cfg), "--out", shQuote(file.path(out, "cli_set_run"))),
                                stdout = TRUE, stderr = TRUE))
cli_man <- tryCatch(yaml::read_yaml(file.path(out, "cli_set_run", "run_manifest.yaml")), error = function(e) NULL)
check("command line --query-set: same exclusion and SNP matrix as in the window, set recorded",
      identical(cli_man$settings$exclude_self, list(AY649695.1 = "AY649695.1")) && identical(cli_man$input$set, "Run set") &&
        identical(utils::read.delim(file.path(out, "cli_set_run", "snp_matrix.tsv")),
                  utils::read.delim(file.path(set_run$out_dir, "snp_matrix.tsv"))))

# CSV for Norwegian Excel: UTF-8 with BOM, ";" between columns, decimal comma.
csv_file <- file.path(out, "table.csv")
write_csv2(data.frame(id = c("0012", "b"), share = c("2.5", "10"), n = c(1.5, 2), place = c("Tromsø", NA)), csv_file)
csv_bytes <- readBin(csv_file, "raw", file.size(csv_file))
csv_lines <- readLines(csv_file, encoding = "UTF-8", warn = FALSE)
check("CSV: UTF-8 byte order mark, semicolons, decimal commas, leading zeros and ø kept",
      identical(csv_bytes[1:3], as.raw(c(0xEF, 0xBB, 0xBF))) &&
        identical(sub("^﻿", "", csv_lines[1]), "\"id\";\"share\";\"n\";\"place\"") &&
        identical(csv_lines[2], "\"0012\";\"2,5\";1,5;\"Tromsø\"") && identical(csv_lines[3], "\"b\";\"10\";2;"))

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
