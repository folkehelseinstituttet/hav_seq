# Sets (Sets page, decision 0012): sequences collected from databases, run results, the
# Alignment page and FASTA files, kept together to view and export. The window
# (app/app.R) only shows what these functions return.
#
# A set is a list:
#   name, created, history (one line per change), folder (where it was saved or opened
#   from, or NULL), save_in (the folder Save writes into), dirty (unsaved changes)
#   members   data.frame, one row per sequence, in the order added: SET_FIELDS, then the
#             members' metadata columns (character; NA where a member has none)
#   seqs      sequences by key: no gaps, as given (original orientation)
#   alignment NULL, or list(aln = aligned sequences by key, turned = keys that are
#             reverse-complemented in it, made, from, method, versions)
# Keys are internal ("db\t<database>\t<ID>" for database samples, "s<n>" otherwise); the
# names shown and saved are worked out from the members each time (set_names()). The
# alignment is "out of date" when a member is not in it: its other rows are still shown.

SET_FIELDS <- c("key", "source", "database", "database_md5", "db_short", "original_id", "added_from", "note")
# Columns of the set table (the window, metadata.tsv, the CSV download), before metadata.
SET_TABLE_FIELDS <- c("id", "source", "database", "database_md5", "original_id", "added_from", "note")
# The files Save writes; Delete removes only these.
SET_FILES <- c("sequences.fasta", "metadata.tsv", "set.yaml", "aligned.fasta")
SET_ALIGN_WARN <- 500   # ask before aligning more sequences than this

set_new <- function(name, save_in = NULL) {
  members <- as.data.frame(stats::setNames(replicate(length(SET_FIELDS), character(0), simplify = FALSE), SET_FIELDS))
  list(name = name, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), history = character(0), folder = NULL,
       save_in = save_in, dirty = FALSE, members = members, seqs = character(0), alignment = NULL, next_key = 1L)
}

# The name used for a database in member names (<name>_<database>): its short_name in
# config.yaml, else its name with spaces and other characters made safe.
db_short_name <- function(name, databases = list()) {
  short <- if (!is.na(name) && !is.null(databases[[name]]$short_name)) databases[[name]]$short_name else name
  gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._-]+", "_", short))
}

# Names in the set (decisions 0008, 0012): a name found in two or more databases becomes
# <name>_<database> for all of them; a query or file sequence named like a database sample
# becomes <name>_query; different sequences given the same name are numbered in the order
# added (<name>_2, _3, ...).
set_names <- function(m) {
  if (!nrow(m)) return(character(0))
  id <- m$original_id
  db <- m$source == "database"
  out <- id
  n_db <- tapply(m$database[db], id[db], function(x) length(unique(x)))
  clash <- db & id %in% names(n_db)[n_db > 1]
  out[clash] <- paste0(id[clash], "_", m$db_short[clash])
  query <- !db & id %in% id[db]
  out[query] <- paste0(out[query], "_query")
  k <- stats::ave(seq_along(id), id, db, FUN = seq_along)
  again <- !db & k > 1
  out[again] <- paste0(out[again], "_", k[again])
  make.unique(out, sep = "_")
}

# Sequences to add (from set_from_database(), set_from_msa(), set_from_fasta()): members
# (SET_FIELDS without key, plus metadata), their sequences, and optionally the alignment
# they came with (aligned rows, which rows are reverse-complemented in it, where from).
set_candidates <- function(ids, seqs, source, added_from, database = NA, database_md5 = NA, db_short = NA,
                           note = NA, meta = NULL, aln = NULL, turned = NULL, aln_info = NULL) {
  n <- length(ids)
  members <- data.frame(source = rep(source, length.out = n), database = rep(database, length.out = n),
                        database_md5 = rep(database_md5, length.out = n), db_short = rep(db_short, length.out = n),
                        original_id = ids, added_from = rep(added_from, length.out = n), note = rep(note, length.out = n))
  members[] <- lapply(members, as.character)
  if (!is.null(meta) && ncol(meta)) {
    meta <- meta[, setdiff(names(meta), c("id", "is_query")), drop = FALSE]
    meta[] <- lapply(meta, as.character)
    clash <- names(meta) %in% c(SET_TABLE_FIELDS, SET_FIELDS)
    names(meta)[clash] <- paste0(names(meta)[clash], "_metadata")
    rownames(meta) <- NULL
    members <- cbind(members, meta)
  }
  list(members = members, seqs = unname(as.character(seqs)), aln = if (!is.null(aln)) unname(aln), turned = turned,
       aln_info = aln_info)
}

# Rows of two data frames with different columns (missing columns are NA).
bind_rows_fill <- function(a, b) {
  for (col in setdiff(names(b), names(a))) a[[col]] <- rep(NA_character_, nrow(a))
  for (col in setdiff(names(a), names(b))) b[[col]] <- rep(NA_character_, nrow(b))
  rbind(a, b[, names(a), drop = FALSE])
}

# Add candidates to a set. The same database sample (same database and ID) is kept once;
# a query or file sequence with the same name and the same sequence as a member is kept
# once. A set that was empty starts with the alignment the candidates came with.
# Returns list(set, added (names), already (names as given), renamed (old -> new names of
# members whose name changed, and of new members not named as given)).
set_add <- function(set, cand, when = format(Sys.time(), "%Y-%m-%d %H:%M")) {
  m <- set$members
  was_empty <- !nrow(m)
  before <- stats::setNames(set_names(m), m$key)
  new <- cand$members
  keys <- character(nrow(new))
  keep <- logical(nrow(new))
  for (i in seq_len(nrow(new))) {
    if (new$source[i] == "database") {
      keys[i] <- paste("db", new$database[i], new$original_id[i], sep = "\t")
      keep[i] <- !keys[i] %in% c(m$key, keys[keep])
    } else {
      own <- m$source != "database" & m$original_id == new$original_id[i]
      same <- set$seqs[m$key[own]]
      earlier <- keep & new$source != "database" & new$original_id == new$original_id[i]
      keep[i] <- !cand$seqs[i] %in% c(same, cand$seqs[earlier])
      if (keep[i]) { keys[i] <- paste0("s", set$next_key); set$next_key <- set$next_key + 1L }
    }
  }
  already <- new$original_id[!keep]
  if (any(keep)) {
    add <- cbind(key = keys[keep], new[keep, , drop = FALSE])
    set$members <- if (was_empty) add[, union(SET_FIELDS, names(add)), drop = FALSE] else bind_rows_fill(m, add)
    rownames(set$members) <- NULL
    set$seqs <- c(set$seqs, stats::setNames(cand$seqs[keep], keys[keep]))
    if (was_empty && !is.null(cand$aln)) {
      info <- cand$aln_info %||% list()
      set$alignment <- list(aln = stats::setNames(cand$aln[keep], keys[keep]),
                            turned = keys[keep][(cand$turned %||% rep(FALSE, length(keep)))[keep]],
                            made = info$made %||% when, from = info$from %||% "", method = info$method %||% "",
                            versions = info$versions)
    }
    set$dirty <- TRUE
    set$history <- c(set$history, sprintf("%s added %d from %s", when, sum(keep),
                                          paste(unique(new$added_from[keep]), collapse = ", ")))
  }
  after <- stats::setNames(set_names(set$members), set$members$key)
  old <- intersect(names(before), names(after))
  changed <- old[before[old] != after[old]]
  added <- keys[keep]
  given <- stats::setNames(new$original_id[keep], added)
  renamed <- c(stats::setNames(after[changed], before[changed]),
               stats::setNames(after[added][after[added] != given], given[after[added] != given]))
  list(set = set, added = unname(after[added]), already = already, renamed = renamed)
}

# What a set_add() did, in one message for the window.
set_add_message <- function(r, missing = character(0)) {
  s <- r$set
  few <- function(x) paste0(paste(utils::head(x, 10), collapse = ", "), if (length(x) > 10) sprintf(" and %d more", length(x) - 10))
  paste0(
    sprintf("%d added to the set '%s' (now %d).", length(r$added), s$name, nrow(s$members)),
    if (length(r$already)) sprintf(" %d %s already in the set: %s.", length(r$already),
                                   if (length(r$already) == 1) "was" else "were", few(r$already)),
    if (length(r$renamed)) paste0(" Renamed because of a name clash: ",
                                  few(paste(names(r$renamed), "->", r$renamed)), "."),
    if (length(missing)) paste0(" No sequence in the database for: ", few(missing), "."),
    if (length(r$added) && length(set_alignment_status(s)$missing))
      " The set's alignment is now out of date: press Align the set to include the new sequences.")
}

set_remove <- function(set, keys, when = format(Sys.time(), "%Y-%m-%d %H:%M")) {
  keys <- intersect(keys, set$members$key)
  if (!length(keys)) return(set)
  set$members <- set$members[!set$members$key %in% keys, , drop = FALSE]
  rownames(set$members) <- NULL
  set$seqs <- set$seqs[set$members$key]
  set$dirty <- TRUE
  set$history <- c(set$history, sprintf("%s removed %d", when, length(keys)))
  set
}

# The table of the set (window, metadata.tsv and the CSV download).
set_table <- function(set) {
  m <- set$members
  out <- data.frame(id = set_names(m), m[, setdiff(SET_TABLE_FIELDS, "id"), drop = FALSE],
                    m[, setdiff(names(m), SET_FIELDS), drop = FALSE], check.names = FALSE)
  rownames(out) <- NULL
  out
}

# The sequences by their names in the set (original orientation).
set_fasta <- function(set) stats::setNames(unname(set$seqs[set$members$key]), set_names(set$members))

# Checksum of a set (0012): lines "name<TAB>SEQUENCE\n" sorted by name in byte order
# (sequence upper case, without gaps or spaces), md5 of the bytes.
set_checksum <- function(seqs) {
  nm <- names(seqs)
  lines <- paste0(nm, "\t", toupper(gsub("[-.[:space:]]", "", seqs)))[order(nm, method = "radix")]
  tmp <- tempfile()
  on.exit(unlink(tmp))
  con <- file(tmp, "wb")
  writeLines(lines, con, sep = "\n", useBytes = TRUE)
  close(con)
  unname(tools::md5sum(tmp))
}

# The alignment of the current members, by their names in the set: rows of members that
# are in it (columns that are a gap in every one of them dropped), the names that are
# reverse-complemented in it, and the members that are not in it (added since).
# NULL when the set has no alignment, or none of its members is in it.
set_alignment <- function(set) {
  st <- set_alignment_status(set)
  if (is.null(st)) return(NULL)
  al <- set$alignment
  m <- msa_matrix(al$aln[st$have])
  keep <- colSums(m != MSA_GAP) > 0
  aln <- if (all(keep)) unname(al$aln[st$have]) else unname(apply(m[, keep, drop = FALSE], 1, rawToChar))
  list(aln = stats::setNames(aln, st$names[st$have]), turned = unname(st$names[intersect(al$turned, st$have)]),
       missing = st$missing, made = al$made, from = al$from, method = al$method, versions = al$versions)
}

# Which members are in the set's alignment (keys), the names of those that are not, and
# the members' names by key. NULL when no member is in it.
set_alignment_status <- function(set) {
  al <- set$alignment
  if (is.null(al)) return(NULL)
  keys <- set$members$key
  have <- keys[keys %in% names(al$aln)]
  if (!length(have)) return(NULL)
  nm <- stats::setNames(set_names(set$members), keys)
  list(have = have, names = nm, missing = unname(nm[setdiff(keys, have)]),
       turned = unname(nm[intersect(al$turned, have)]))
}

# One line about the set's alignment, for the Sets and Alignment pages.
set_alignment_text <- function(set) {
  al <- set_alignment_status(set)
  if (is.null(al)) return(if (nrow(set$members) < 2) "Not aligned (add at least 2 sequences)." else "Not aligned yet.")
  paste0(sprintf("Aligned %s (%s).", set$alignment$made, set$alignment$from),
         if (length(al$turned)) sprintf(" %d turned (reverse complement).", length(al$turned)),
         if (length(al$missing)) sprintf(" Out of date: %d added since %s not in it.", length(al$missing),
                                         if (length(al$missing) == 1) "is" else "are"))
}

# The set's alignment as a viewer alignment (R/msa.R), or NULL. Rows remember their
# members, so rows chosen in the Alignment page can be added to another set as they are.
msa_from_set <- function(set) {
  al <- set_alignment(set)
  if (is.null(al)) return(NULL)
  tab <- set_table(set)
  meta <- tab[match(names(al$aln), tab$id), setdiff(names(tab), c("database_md5", "added_from")), drop = FALSE]
  few <- function(x) paste0(paste(utils::head(x, 10), collapse = ", "), if (length(x) > 10) sprintf(" and %d more", length(x) - 10))
  note <- paste0(sprintf("Alignment: %s, %s.", al$from, al$made),
                 if (length(al$turned)) paste0(" Turned (reverse complement) for the alignment, so numbered on the reverse complement: ",
                                               few(al$turned), "."),
                 if (length(al$missing)) paste0(" Not shown: ", length(al$missing), " added after the alignment was made (",
                                                few(al$missing), "). Press Align the set to include them."))
  a <- msa_load(al$aln, set$name, meta = meta, queries = tab$id[tab$source == "query"], revcomp = al$turned, note = note)
  a$origin <- list(kind = "set", members = cbind(id = set_names(set$members), set$members), seqs = set$seqs)
  a
}

# The set (or the members with these keys) as queries for Outbreak detection (0012):
# sequences by their names in the set, as given; exclude_self for the members from the
# searched database (name in the set -> their ID there; others leave out nothing); and
# what run_manifest.yaml records about the set (`input`).
set_queries <- function(set, database, keys = NULL) {
  m <- set$members
  if (length(keys)) m <- m[m$key %in% keys, , drop = FALSE]
  if (!nrow(m)) stop("The set '", set$name, "' has no sequences. Add sequences on the Sets page first.", call. = FALSE)
  if (nrow(m) > 50) {
    stop(sprintf(paste("%s has %d sequences; Outbreak detection takes at most 50. On the Sets page, click up to 50",
                       "rows of the set and press Use in Outbreak detection."),
                 if (length(keys)) sprintf("The choice from the set '%s'", set$name) else sprintf("The set '%s'", set$name),
                 nrow(m)), call. = FALSE)
  }
  names_in_set <- set_names(set$members)[match(m$key, set$members$key)]
  own <- m$source == "database" & m$database %in% database
  list(seqs = stats::setNames(unname(set$seqs[m$key]), names_in_set),
       exclude_self = if (any(own)) as.list(stats::setNames(m$original_id[own], names_in_set[own])) else FALSE,
       n_own = sum(own),
       input = Filter(Negate(is.null), list(
         set = set$name, set_checksum = set_checksum(set_fasta(set)), set_sequences = nrow(set$members),
         set_folder = set$folder, set_saved = !is.null(set$folder) && !isTRUE(set$dirty),
         queries_from_set = if (nrow(m) < nrow(set$members)) "chosen rows" else "all")))
}

# ---- Where new members come from -------------------------------------------------------

# Database samples (in the order given). IDs without a sequence are returned as `missing`.
set_from_database <- function(db, cfg, ids, meta, added_from = "Database viewer") {
  ref <- load_db_seqs(db, cfg)
  found <- intersect(ids, names(ref$seqs))
  cand <- set_candidates(found, ref$seqs[found], "database", added_from, database = db$name,
                         database_md5 = ref$checksum, db_short = db_short_name(db$name, cfg$databases),
                         meta = meta[match(found, meta$id), , drop = FALSE])
  cand$missing <- setdiff(ids, found)
  cand
}

# Where the rows of a run's alignment come from (msa_load_run_folder(), msa_from_result()).
msa_run_origin <- function(manifest, run, databases = list()) {
  aligner <- manifest$settings$aligner %||% "decipher"
  method <- if (aligner == "mafft") "MAFFT --auto (pymafft)" else paste("DECIPHER AlignSeqs", manifest$tools$DECIPHER %||% "")
  list(kind = "run", run = run, database = manifest$database$name, md5 = manifest$database$fasta_md5,
       db_short = db_short_name(manifest$database$name, databases), renamed = unlist(manifest$renamed),
       made = manifest$created, method = trimws(method), versions = manifest$tools)
}

# Rows of a viewer alignment (run result, run folder, FASTA file or set). Queries and
# sequences turned for the alignment are kept as given (turned back); the alignment comes
# along for a set that was empty.
set_from_msa <- function(a, rows = rownames(a$m), added_from = a$label) {
  o <- a$origin %||% list(kind = "file", file = a$label)
  sub <- msa_subset(a, rows)
  rows <- rownames(sub$m)
  aln <- msa_fasta(sub, seq_len(nrow(sub$m)))
  turned <- rows %in% a$revcomp
  as_given <- vapply(seq_along(rows), function(i) {
    s <- gsub("-", "", aln[[i]], fixed = TRUE)
    if (turned[i]) reverse_complement(s) else s
  }, "")
  meta <- if (!is.null(a$meta)) a$meta[match(rows, a$meta$id), , drop = FALSE]
  info <- list(made = o$made, method = o$method, versions = o$versions,
               from = switch(o$kind, run = paste("from the run", o$run), set = "from another set",
                             paste("from the file", o$file %||% a$label)))
  if (o$kind == "set") {
    m <- o$members[match(rows, o$members$id), , drop = FALSE]
    cand <- set_candidates(m$original_id, o$seqs[m$key], m$source, m$added_from, database = m$database,
                           database_md5 = m$database_md5, db_short = m$db_short, note = m$note,
                           meta = m[, setdiff(names(m), c(SET_FIELDS, "id")), drop = FALSE],
                           aln = aln, turned = turned, aln_info = info)
  } else if (o$kind == "run") {
    query <- rows %in% a$queries
    given <- ifelse(rows %in% o$renamed, names(o$renamed)[match(rows, o$renamed)], rows)
    cand <- set_candidates(given, as_given, ifelse(query, "query", "database"), added_from,
                           database = ifelse(query, NA, o$database), database_md5 = ifelse(query, NA, o$md5 %||% NA),
                           db_short = ifelse(query, NA, o$db_short),
                           note = ifelse(query & turned, paste("reverse-complemented in the run", o$run, "(kept as given here)"), NA),
                           meta = meta, aln = aln, turned = turned, aln_info = info)
  } else {
    cand <- set_candidates(rows, as_given, "file", added_from,
                           note = ifelse(turned, "reverse-complemented when aligned (kept as given here)", NA),
                           meta = meta, aln = aln, turned = turned, aln_info = info)
  }
  cand
}

# A FASTA file, aligned or not (names up to the first space). An aligned file (two or more
# sequences, all the same length, as in the Alignment page) comes with its alignment.
set_from_fasta <- function(path, file_name = basename(path)) {
  seqs <- read_fasta(path)
  raw <- gsub("[[:space:]]", "", read_alignment(path))
  raw <- chartr(".", "-", chartr("U", "T", raw))
  aligned <- length(raw) == length(seqs) && length(seqs) > 1 && length(unique(nchar(raw))) == 1 &&
    all(grepl(paste0("^[", IUPAC_NT, "-]+$"), raw))
  aln <- if (aligned) raw
  set_candidates(names(seqs), seqs, "file", paste("File", file_name), aln = aln,
                 turned = if (aligned) rep(FALSE, length(seqs)),
                 aln_info = list(from = paste("from the file", file_name), method = "aligned before it was added"))
}

# ---- Saving and opening ----------------------------------------------------------------------

# A folder name from a set name (characters Windows does not allow become "_").
set_folder_name <- function(name) {
  x <- sub("[. ]+$", "", gsub("[<>:\"/\\\\|?*[:cntrl:]]", "_", trimws(name)))
  if (nzchar(x)) x else "set"
}

# "free" (no folder or an empty one), "set" (a saved set) or "other" (a folder with other files).
set_folder_state <- function(folder) {
  if (file.exists(file.path(folder, "set.yaml"))) return("set")
  if (!dir.exists(folder) || !length(list.files(folder, all.files = TRUE, no.. = TRUE))) "free" else "other"
}

set_yaml <- function(set, al = set_alignment(set)) {
  m <- set$members
  db <- m$source == "database"
  groups <- unique(m[db, c("database", "database_md5", "db_short"), drop = FALSE])
  sources <- c(
    lapply(seq_len(nrow(groups)), function(i) list(
      type = "database", database = groups$database[i], short_name = groups$db_short[i],
      database_md5 = groups$database_md5[i],
      n = sum(db & m$database == groups$database[i] & m$database_md5 %in% groups$database_md5[i]))),
    if (any(m$source == "query")) list(list(type = "query", n = sum(m$source == "query"))),
    lapply(unique(m$added_from[m$source == "file"]), function(f) list(type = "file", file = sub("^File ", "", f),
                                                                       n = sum(m$source == "file" & m$added_from == f))))
  fasta <- set_fasta(set)
  list(name = set$name, made_with = paste("FHInumerics", app_version()), created = set$created,
       saved = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_sequences = nrow(m), checksum = set_checksum(fasta),
       checksum_rule = "md5 of lines 'name<TAB>SEQUENCE' (upper case, no gaps), sorted by name in byte order, each ending in a newline",
       sources = sources, history = as.list(set$history),
       alignment = if (!is.null(al)) list(file = "aligned.fasta", made = al$made, from = al$from, method = al$method,
                                          versions = al$versions, n_sequences = length(al$aln), columns = nchar(al$aln[[1]]),
                                          reverse_complemented = as.list(al$turned), not_in_it = as.list(al$missing)))
}

# Write a set to its folder: sequences.fasta, metadata.tsv, set.yaml and (if aligned)
# aligned.fasta. Returns the set, now saved there.
set_save <- function(set, folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(folder)) stop("Could not make the folder ", folder, ". Choose another folder.", call. = FALSE)
  tab <- set_table(set)
  tab[] <- lapply(tab, function(x) gsub("[\t\r\n]+", " ", x))
  write_fasta(set_fasta(set), file.path(folder, "sequences.fasta"))
  utils::write.table(tab, file.path(folder, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  al <- set_alignment(set)
  if (is.null(al)) unlink(file.path(folder, "aligned.fasta")) else write_fasta(al$aln, file.path(folder, "aligned.fasta"))
  yaml::write_yaml(set_yaml(set, al), file.path(folder, "set.yaml"))
  set$folder <- normalizePath(folder)
  set$save_in <- dirname(set$folder)
  set$dirty <- FALSE
  set
}

# Open a saved set folder.
set_load <- function(folder) {
  bad <- function(why) stop("The folder ", folder, " is not a FHInumerics set that can be opened (", why, ").", call. = FALSE)
  if (!file.exists(file.path(folder, "set.yaml"))) bad("it has no set.yaml")
  y <- tryCatch(yaml::read_yaml(file.path(folder, "set.yaml")), error = function(e) bad("set.yaml cannot be read"))
  tab <- tryCatch(utils::read.delim(file.path(folder, "metadata.tsv"), colClasses = "character", check.names = FALSE,
                                    na.strings = "", quote = ""),
                  error = function(e) bad("metadata.tsv is missing or cannot be read"))
  seqs <- tryCatch(read_fasta(file.path(folder, "sequences.fasta")), error = function(e) bad("sequences.fasta is missing or damaged"))
  if (!all(SET_TABLE_FIELDS %in% names(tab)) || !setequal(tab$id, names(seqs)) || anyDuplicated(tab$id)) {
    bad("sequences.fasta and metadata.tsv do not list the same sequences")
  }
  short <- unlist(lapply(y$sources, function(s) if (identical(s$type, "database")) stats::setNames(s$short_name, s$database))) %||% character(0)
  db <- tab$source == "database"
  keys <- ifelse(db, paste("db", tab$database, tab$original_id, sep = "\t"), paste0("s", seq_len(nrow(tab))))
  members <- data.frame(key = keys, tab[, setdiff(SET_TABLE_FIELDS, "id"), drop = FALSE], check.names = FALSE)
  members$db_short <- ifelse(db, unname(short[tab$database]), NA_character_)
  lost <- db & is.na(members$db_short)
  members$db_short[lost] <- vapply(tab$database[lost], db_short_name, "")
  members <- cbind(members[, SET_FIELDS], tab[, setdiff(names(tab), SET_TABLE_FIELDS), drop = FALSE])
  set <- set_new(y$name %||% basename(folder))
  set$created <- y$created %||% set$created
  set$history <- unlist(y$history) %||% character(0)
  set$members <- members
  set$seqs <- stats::setNames(unname(seqs[tab$id]), keys)
  set$next_key <- nrow(tab) + 1L
  aln_file <- file.path(folder, "aligned.fasta")
  if (file.exists(aln_file) && !is.null(y$alignment)) {
    aln <- read_alignment(aln_file)
    turned <- unlist(y$alignment$reverse_complemented)
    key <- keys[match(names(aln), tab$id)]
    # Rows that are not a member, or whose letters are not the member's sequence, are not used.
    ok <- !is.na(key) & vapply(seq_along(aln), function(i) {
      s <- gsub("-", "", aln[[i]], fixed = TRUE)
      !is.na(key[i]) && identical(if (names(aln)[i] %in% turned) reverse_complement(s) else s, set$seqs[[key[i]]])
    }, TRUE)
    if (any(ok)) {
      set$alignment <- list(aln = stats::setNames(unname(aln[ok]), key[ok]), turned = key[ok & names(aln) %in% turned],
                            made = y$alignment$made, from = y$alignment$from, method = y$alignment$method,
                            versions = y$alignment$versions)
    }
  }
  set$folder <- normalizePath(folder)
  set$save_in <- dirname(set$folder)
  set
}

# Remove a saved set's files (only SET_FILES), then its folder if nothing else is in it.
# Returns the files left in the folder.
set_delete <- function(folder) {
  if (!file.exists(file.path(folder, "set.yaml"))) stop("This folder is not a saved set: ", folder, call. = FALSE)
  unlink(file.path(folder, SET_FILES))
  left <- list.files(folder, all.files = TRUE, no.. = TRUE)
  if (!length(left)) unlink(folder, recursive = TRUE)
  left
}

# The sets saved in a folder (one subfolder per set): name, folder, number of sequences.
set_list_saved <- function(dir) {
  yamls <- if (dir.exists(dir)) file.path(list.dirs(dir, recursive = FALSE), "set.yaml") else character(0)
  yamls <- yamls[file.exists(yamls)]
  rows <- lapply(yamls, function(f) {
    y <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
    if (is.null(y)) return(NULL)
    data.frame(name = y$name %||% basename(dirname(f)), folder = normalizePath(dirname(f)), n = y$n_sequences %||% NA)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame(name = character(0), folder = character(0), n = integer(0)) else out[order(out$name), ]
}

# ---- Aligning (in a background R process in the window) ----------------------------------

# Six complete genomes, one per subgenotype (data/demo/SOURCE.md): only references for
# OrientNucleotides, never part of the alignment.
orientation_panel <- function() read_fasta(file.path(getOption("fhinumerics.root", "."), "data", "orientation_panel.fasta"))

# Orient (OrientNucleotides against the panel; reverse complement only) and align
# (DECIPHER AlignSeqs, as in runs). Sequences get internal IDs (S1, S2, ...) for the tools.
# Returns the aligned sequences (names as given), the names turned, and how it was made.
set_align <- function(seqs, threads = 1, verbose = FALSE) {
  if (length(seqs) < 2) stop("At least 2 sequences are needed for an alignment.", call. = FALSE)
  step <- function(x) if (verbose) cat("STEP ", x, "\n", sep = "")
  panel <- orientation_panel()
  ids <- paste0("S", seq_along(seqs))
  step("Orienting")
  dna <- Biostrings::DNAStringSet(c(stats::setNames(unname(panel), paste0("R", seq_along(panel))),
                                    stats::setNames(unname(seqs), ids)))
  o <- DECIPHER::OrientNucleotides(dna, reference = seq_along(panel), orientation = "both", type = "both",
                                   verbose = FALSE, processors = threads)
  turned <- names(seqs)[o[[1]][-seq_along(panel)] %in% "rc"]
  aln <- DECIPHER::AlignSeqs(o[[2]][-seq_along(panel)], processors = threads, verbose = verbose)
  step("Saving")
  pkg <- function(p) as.character(utils::packageVersion(p))
  list(aln = stats::setNames(as.character(aln)[match(ids, names(aln))], names(seqs)), turned = turned,
       made = format(Sys.time(), "%Y-%m-%d %H:%M"),
       from = "by FHInumerics with DECIPHER",
       method = paste("DECIPHER OrientNucleotides (reverse complement, against 6 HAV genomes, one per subgenotype),",
                      "then DECIPHER AlignSeqs"),
       versions = list(FHInumerics = app_version(), R = R.version.string, DECIPHER = pkg("DECIPHER"),
                       Biostrings = pkg("Biostrings")))
}

# Start set_align() in a background R process (app/align_worker.R). Returns the job (an
# environment) for align_poll() and align_cancel().
align_start <- function(seqs, threads = 1, dir = tempfile("fhinumerics-align-")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(seqs = seqs, threads = threads, libs = .libPaths()), file.path(dir, "input.rds"))
  root <- getOption("fhinumerics.root", ".")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  job <- new.env()
  job$dir <- dir
  job$n <- length(seqs)
  job$started <- Sys.time()
  job$text <- ""
  job$proc <- processx::process$new(rscript, c(file.path(root, "app", "align_worker.R"), root, dir),
                                    stdout = "|", stderr = file.path(dir, "stderr.txt"),
                                    cleanup = TRUE, cleanup_tree = TRUE, supervise = TRUE)
  job
}

# Where a background alignment is: list(done, value (0-1), text) and, when done, either
# `result` (as set_align()) or `error` (a plain message).
align_poll <- function(job) {
  alive <- job$proc$is_alive()
  out <- tryCatch(job$proc$read_output(), error = function(e) "")
  job$text <- paste0(job$text, out)
  p <- align_progress(job$text)
  if (alive) return(c(list(done = FALSE), p))
  result_file <- file.path(job$dir, "result.rds")
  res <- if (file.exists(result_file)) readRDS(result_file)
  unlink(job$dir, recursive = TRUE)
  if (is.null(res)) {
    return(list(done = TRUE, error = sprintf(paste("The alignment stopped unexpectedly (exit status %s).",
                                                   "Try again; if it happens again, the computer may be short of memory."),
                                             job$proc$get_exit_status())))
  }
  if (!is.null(res$error)) return(list(done = TRUE, error = paste("The alignment could not be made:", res$error)))
  list(done = TRUE, result = res)
}

align_cancel <- function(job) {
  if (job$proc$is_alive()) tryCatch(job$proc$kill_tree(), error = function(e) job$proc$kill())
  job$proc$wait(2000)
  unlink(job$dir, recursive = TRUE)
  invisible(!job$proc$is_alive())
}

# Progress from what the background process printed so far: our STEP lines and DECIPHER
# AlignSeqs' own (a heading line, then a bar of up to 80 "="). Returns list(value, text).
align_progress <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  heads <- grep("^(STEP |Iteration |Determining|Clustering|Reclustering|Aligning|Realigning|Refining)", lines)
  if (!length(heads)) return(list(value = 0, text = "Starting"))
  last <- heads[length(heads)]
  head <- lines[last]
  within <- min(1, nchar(gsub("[^=]", "", paste(lines[-seq_len(last)], collapse = ""))) / 80)
  iter <- grep("^Iteration [0-9]+ of [0-9]+", lines[seq_len(last)], value = TRUE)
  round <- if (length(iter)) as.integer(sub("^Iteration ([0-9]+) of.*", "\\1", iter[length(iter)])) else 0L
  rounds <- if (length(iter)) as.integer(sub("^Iteration [0-9]+ of ([0-9]+).*", "\\1", iter[length(iter)])) else 2L
  span <- function(from, to, what) list(value = from + (to - from) * within, text = what)
  if (startsWith(head, "STEP Orienting")) return(span(0, 0.05, "Finding reverse complements"))
  if (startsWith(head, "STEP Saving")) return(list(value = 0.99, text = "Finishing"))
  if (startsWith(head, "Refining")) return(span(0.95, 0.99, "Final check of the alignment"))
  if (round == 0) {
    if (startsWith(head, "Determining")) return(span(0.05, 0.12, "Comparing the sequences"))
    if (startsWith(head, "Clustering")) return(span(0.12, 0.15, "Grouping similar sequences"))
    return(span(0.15, 0.6, "Aligning"))
  }
  from <- 0.6 + 0.35 * (round - 1) / rounds
  size <- 0.35 / rounds
  part <- if (startsWith(head, "Iteration") || startsWith(head, "Determining")) c(0, 0.1) else
    if (startsWith(head, "Reclustering")) c(0.1, 0.15) else c(0.15, 1)
  span(from + size * part[1], from + size * part[2], sprintf("Improving the alignment (round %d of %d)", round, rounds))
}
