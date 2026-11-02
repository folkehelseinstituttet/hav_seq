# Database viewer: dates, sample IDs and summary charts of a database's metadata.
# The window (app/app.R) only shows what these functions return.

# "2023-05-14", "2023-05" and "2023" -> Date (a partial date counts as the first day of
# that month or year). Anything else becomes NA.
parse_partial_date <- function(x) {
  x <- gsub("/", "-", trimws(as.character(x)))
  x <- ifelse(grepl("^\\d{4}$", x), paste0(x, "-01-01"),
              ifelse(grepl("^\\d{4}-\\d{1,2}$", x), paste0(x, "-01"), x))
  x[!grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", x)] <- NA
  as.Date(x, format = "%Y-%m-%d")
}

# Sample IDs typed, pasted or read from a file: one per line, or separated by spaces,
# commas or semicolons. Returns the IDs found in the database and those that are not.
match_sample_ids <- function(text, ids) {
  wanted <- unique(unlist(strsplit(paste(text, collapse = "\n"), "[[:space:],;]+")))
  wanted <- sub("^>", "", wanted[nzchar(wanted)])
  list(found = intersect(wanted, ids), missing = setdiff(wanted, ids))
}

# The sequences of the given database samples (in the order given), for FASTA download
# or as queries. Samples without a sequence in the FASTA are left out.
db_sequences <- function(db, cfg, ids) {
  ref <- load_db_seqs(db, cfg)
  ref$seqs[intersect(ids, names(ref$seqs))]
}

# A table as CSV for Excel with Norwegian settings: ";" between columns, decimal comma,
# UTF-8 with a byte order mark (so Excel reads æ, ø, å). Metadata columns are text; one
# where every value is a decimal number (like 12.5) gets decimal commas too. Other
# values are left as they are (leading zeros stay).
write_csv2 <- function(df, file) {
  df[] <- lapply(df, function(x) {
    if (!is.character(x)) return(x)
    v <- x[!is.na(x) & nzchar(x)]
    if (length(v) && all(grepl("^-?[0-9]+(\\.[0-9]+)?$", v)) && any(grepl(".", v, fixed = TRUE))) sub(".", ",", x, fixed = TRUE) else x
  })
  con <- file(file, open = "wb")
  on.exit(close(con))
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  text <- utils::capture.output(utils::write.csv2(df, stdout(), row.names = FALSE, na = ""))
  writeLines(enc2utf8(text), con, sep = "\r\n", useBytes = TRUE)
  invisible(file)
}

# ---- Charts ------------------------------------------------------------------------
# Many different values make an unreadable legend: keep the most common, rest "Other".
lump_values <- function(x, n = 12) {
  x <- ifelse(is.na(x), "(missing)", x)
  top <- names(sort(table(x), decreasing = TRUE))[seq_len(min(n, length(unique(x))))]
  factor(ifelse(x %in% top, x, "Other"), levels = c(setdiff(top, "Other"), "Other"))
}

# Bars of samples per month or year, optionally split by a metadata column.
plot_over_time <- function(meta, dates, unit = "month", colour_by = "none") {
  keep <- !is.na(dates)
  if (!any(keep)) stop("None of these samples has a date.", call. = FALSE)
  d <- data.frame(period = as.Date(cut(dates[keep], unit)))
  d$group <- if (colour_by == "none") "all" else lump_values(meta[[colour_by]][keep])
  p <- ggplot2::ggplot(d, ggplot2::aes(period, fill = group)) +
    ggplot2::geom_bar(width = if (unit == "year") 300 else 25) +
    ggplot2::labs(x = NULL, y = "Number of samples", fill = colour_by) +
    ggplot2::theme_minimal(base_size = 13)
  if (colour_by == "none") p <- p + ggplot2::scale_fill_manual(values = "#3182bd", guide = "none")
  p
}

# Horizontal bars: number of samples per value of a column (the 30 most common values).
plot_counts <- function(meta, column, colour_by = "none") {
  d <- data.frame(value = ifelse(is.na(meta[[column]]), "(missing)", meta[[column]]))
  d$group <- if (colour_by == "none") "all" else lump_values(meta[[colour_by]])
  top <- names(sort(table(d$value), decreasing = TRUE))[seq_len(min(30, length(unique(d$value))))]
  d <- d[d$value %in% top, , drop = FALSE]
  d$value <- factor(d$value, levels = rev(top))
  p <- ggplot2::ggplot(d, ggplot2::aes(y = value, fill = group)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Number of samples", y = column, fill = colour_by) +
    ggplot2::theme_minimal(base_size = 13)
  if (colour_by == "none") p <- p + ggplot2::scale_fill_manual(values = "#3182bd", guide = "none")
  p
}

# Counts of one column against another, with totals.
cross_table <- function(meta, rows, cols) {
  r <- ifelse(is.na(meta[[rows]]), "(missing)", meta[[rows]])
  k <- ifelse(is.na(meta[[cols]]), "(missing)", meta[[cols]])
  tab <- stats::addmargins(table(r, k), FUN = list(Total = sum), quiet = TRUE)
  out <- as.data.frame.matrix(tab)
  cbind(stats::setNames(data.frame(rownames(out), check.names = FALSE), rows), out)
}

# Country names in metadata that differ from the names in the 'maps' package.
COUNTRY_ALIASES <- c(
  "United States" = "USA", "United States of America" = "USA", "US" = "USA",
  "United Kingdom" = "UK", "Great Britain" = "UK", "England" = "UK", "Scotland" = "UK", "Wales" = "UK",
  "Viet Nam" = "Vietnam", "Russian Federation" = "Russia", "Czechia" = "Czech Republic",
  "Korea" = "South Korea", "Republic of Korea" = "South Korea", "Republic of the Congo" = "Republic of Congo",
  "Iran, Islamic Republic of" = "Iran", "Turkiye" = "Turkey", "Hong Kong" = "China", "Cote d'Ivoire" = "Ivory Coast"
)

# Samples per country on a world map (circle size = number of samples).
# Returns the plot, the circles (one row per country: region, long, lat, n, the country
# names as written in the metadata, and the counts per value of `split_by`), and the
# country values that could not be placed on the map.
plot_country_map <- function(meta, column, zoom = TRUE, split_by = "none") {
  world <- ggplot2::map_data("world")
  meta <- meta[!is.na(meta[[column]]), , drop = FALSE]
  if (!nrow(meta)) stop("None of these samples has a country.", call. = FALSE)
  raw <- meta[[column]]
  name <- trimws(sub("\\s*\\(.*\\)\\s*$", "", raw))   # "Egypt (suspected)" -> "Egypt"
  name <- ifelse(name %in% names(COUNTRY_ALIASES), COUNTRY_ALIASES[name], name)
  unknown <- sort(unique(raw[!name %in% world$region]))
  on_map <- name %in% world$region
  if (!any(on_map)) stop("None of the country names could be found on the map.", call. = FALSE)
  meta <- meta[on_map, , drop = FALSE]
  name <- name[on_map]

  circles <- data.frame(region = sort(unique(name)))
  circles$n <- as.vector(table(name)[circles$region])
  circles$names <- vapply(circles$region, function(r) paste(sort(unique(meta[[column]][name == r])), collapse = ", "), "")
  circles$split <- lapply(circles$region, function(r) {
    if (split_by == "none") return(NULL)
    v <- meta[[split_by]][name == r]
    sort(table(ifelse(is.na(v), "(missing)", v)), decreasing = TRUE)
  })

  # A circle sits in the middle of the country's largest piece of land.
  pieces <- split(world[world$region %in% circles$region, ], ~ region + group, drop = TRUE)
  biggest <- do.call(rbind, lapply(split(pieces, vapply(pieces, function(p) p$region[1], "")),
                                   function(ps) ps[[which.max(vapply(ps, nrow, 1L))]]))
  centre <- stats::aggregate(cbind(long, lat) ~ region, biggest, function(v) mean(range(v)))
  circles <- merge(circles, centre, by = "region", sort = FALSE)
  # Big circles first, so small ones are drawn on top and can still be pointed at.
  circles <- circles[order(-circles$n), ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = world, ggplot2::aes(long, lat, group = group),
                          fill = "grey90", colour = "white", linewidth = 0.2) +
    ggplot2::geom_point(data = circles, ggplot2::aes(long, lat, size = n),
                        colour = "#d7301f", alpha = 0.6) +
    ggplot2::scale_size_area(max_size = 18, name = "Samples") +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  if (zoom) {
    pad <- 12
    p <- p + ggplot2::coord_quickmap(xlim = range(circles$long) + c(-pad, pad),
                                     ylim = pmax(pmin(range(circles$lat) + c(-pad, pad), 85), -60))
  } else {
    p <- p + ggplot2::coord_quickmap(ylim = c(-60, 85))
  }
  list(plot = p, circles = circles, unknown = unknown)
}

# What the hover box shows for one map circle: a named character vector.
circle_details <- function(circle, split_by = "none") {
  out <- c(Country = circle$names, Samples = as.character(circle$n))
  counts <- circle$split[[1]]
  if (split_by != "none" && length(counts)) {
    shown <- utils::head(counts, 10)
    out <- c(out, stats::setNames(as.character(shown), paste0(split_by, " ", names(shown))))
    if (length(counts) > 10) out <- c(out, stats::setNames(as.character(sum(utils::tail(counts, -10))), "other values"))
  }
  out
}

# Columns that hold dates: at least 80% of the filled-in values read as a (partial) date.
date_columns <- function(meta) {
  cols <- setdiff(names(meta), "id")
  ok <- vapply(cols, function(col) {
    v <- meta[[col]][!is.na(meta[[col]]) & nzchar(meta[[col]])]
    length(v) > 0 && mean(!is.na(parse_partial_date(v))) >= 0.8
  }, TRUE)
  cols[ok]
}
