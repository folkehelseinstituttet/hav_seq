# Reference databases: a FASTA file plus a metadata table, as named in config.yaml.
# The search index is built automatically into the app's cache folder, so the
# database files themselves can live on a read-only share (e.g. a network drive).

load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  base <- dirname(normalizePath(path))
  # Relative paths in config.yaml are relative to the config file itself.
  resolve <- function(p) {
    p <- path.expand(gsub("\\\\", "/", p))
    if (is_absolute_path(p)) p else file.path(base, p)
  }
  cfg$output_dir <- resolve(cfg$output_dir %||% "~/FHInumerics-results")
  cfg$cache_dir <- resolve(cfg$cache_dir %||% "cache")
  cfg$python <- find_python(if (!is.null(cfg$python)) resolve(cfg$python))
  if (!is.null(cfg$blast_dir)) cfg$blast_dir <- resolve(cfg$blast_dir)
  cfg$databases <- lapply(cfg$databases, function(db) {
    db$fasta <- resolve(db$fasta)
    db$metadata <- resolve(db$metadata)
    db$id_column <- db$id_column %||% "id"
    db$country_column <- db$country_column %||% "country"
    db
  })
  names(cfg$databases) <- vapply(cfg$databases, `[[`, "", "name")
  cfg
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# "/home/...", "C:/...", "C:\..." and "\\server\share" are absolute.
is_absolute_path <- function(p) grepl("^(/|[A-Za-z]:[\\\\/]|\\\\\\\\|//)", p)

read_metadata <- function(db) {
  sep <- if (grepl("\\.csv$", db$metadata, ignore.case = TRUE)) "," else "\t"
  meta <- utils::read.delim(db$metadata, sep = sep, colClasses = "character",
                            check.names = FALSE, na.strings = c("", "NA"), quote = "\"")
  if (!db$id_column %in% names(meta)) {
    stop("Metadata file has no column '", db$id_column, "': ", db$metadata, call. = FALSE)
  }
  names(meta)[names(meta) == db$id_column] <- "id"
  meta[!duplicated(meta$id), , drop = FALSE]
}

# Columns offered for labels/colouring: the configured ones, else everything but the ID.
display_columns <- function(db, meta) {
  cols <- db$display_columns %||% setdiff(names(meta), "id")
  intersect(cols, names(meta))
}

# Build (or reuse) the search index for a database. Sequences get internal IDs
# (D1, D2, ...) so that unusual characters in real IDs never reach the search,
# alignment or tree. Returns a list with the sequences, ID map and index.
#   search = "decipher": DECIPHER k-mer index (R only), kept in index.rds
#   search = "blast":    BLAST database for rBLAST (needs BLAST+ on this computer)
prepare_db <- function(db, cfg, search = "decipher", log = message) {
  ref <- load_db_seqs(db, cfg)
  dir <- ref$dir

  if (search == "decipher") {
    index_rds <- file.path(dir, "index.rds")
    if (!file.exists(index_rds)) {
      log("Indexing database '", db$name, "' (first use only)...")
      saveRDS(DECIPHER::IndexSeqs(db_stringset(ref), K = 11, maskRepeats = FALSE, maskLCRs = FALSE,
                                  maskNumerous = FALSE, verbose = FALSE), index_rds)
    }
    ref$index <- readRDS(index_rds)
  } else if (search == "blast") {
    setup_blast(cfg$blast_dir)
    ref$prefix <- blast_path(file.path(dir, "db"))
    if (!file.exists(paste0(ref$prefix, ".nsq"))) {
      log("Building BLAST database for '", db$name, "' (first use only)...")
      write_fasta(stats::setNames(ref$seqs, ref$ids$internal), paste0(ref$prefix, ".fasta"))
      status <- rBLAST::makeblastdb(paste0(ref$prefix, ".fasta"), db_name = ref$prefix, dbtype = "nucl",
                                    hash_index = FALSE, verbose = FALSE)
      if (!identical(as.integer(status), 0L)) stop("makeblastdb failed (exit status ", status, ")", call. = FALSE)
    }
  } else {
    stop("Unknown search method '", search, "' (use decipher or blast).", call. = FALSE)
  }
  ref
}

# The database sequences and their internal IDs, cached in cache/<name>_<md5 prefix>/.
load_db_seqs <- function(db, cfg) {
  if (!file.exists(db$fasta)) stop("Database FASTA not found: ", db$fasta, call. = FALSE)
  checksum <- unname(tools::md5sum(db$fasta))
  dir <- file.path(cfg$cache_dir, paste0(gsub("[^A-Za-z0-9]+", "_", db$name), "_", substr(checksum, 1, 12)))
  seqs_rds <- file.path(dir, "seqs.rds")
  if (!file.exists(seqs_rds)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    seqs <- read_fasta(db$fasta)
    ids <- data.frame(internal = paste0("D", seq_along(seqs)), id = names(seqs))
    saveRDS(list(seqs = seqs, ids = ids, checksum = checksum), seqs_rds)
  }
  ref <- readRDS(seqs_rds)
  ref$dir <- dir
  ref
}

db_stringset <- function(ref) {
  Biostrings::DNAStringSet(stats::setNames(ref$seqs, ref$ids$internal))
}

# ---- Optional: BLAST+ through rBLAST ---------------------------------------------
# Only used with `search: blast` in config.yaml, when BLAST+ is installed on this
# computer. rBLAST finds the programs on the PATH; blast_dir is added in front.
setup_blast <- function(blast_dir = NULL) {
  if (!is.null(blast_dir) && dir.exists(blast_dir)) {
    dir <- normalizePath(blast_dir)
    if (.Platform$OS.type == "windows") dir <- utils::shortPathName(dir)  # system2() and spaces
    Sys.setenv(PATH = paste(c(dir, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  }
  if (!requireNamespace("rBLAST", quietly = TRUE) || !rBLAST::has_blast()) {
    stop("search: blast needs the rBLAST package and BLAST+ (blastn) on this computer. ",
         "Set blast_dir in config.yaml, or use search: decipher.", call. = FALSE)
  }
}

# BLAST cannot use paths with spaces; on Windows use the short (8.3) form of the folder.
blast_path <- function(p) {
  p <- file.path(normalizePath(dirname(p), mustWork = FALSE), basename(p))
  if (.Platform$OS.type == "windows" && dir.exists(dirname(p))) {
    p <- file.path(utils::shortPathName(dirname(p)), basename(p))
  }
  if (grepl(" ", p)) {
    stop("BLAST cannot use folders with spaces in the name: ", dirname(p),
         ". Set cache_dir in config.yaml to a folder without spaces.", call. = FALSE)
  }
  p
}
