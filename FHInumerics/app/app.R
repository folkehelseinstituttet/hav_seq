# FHInumerics: compare new HAV sequences to a reference database, and browse a database.
# Start with FHInumerics.bat (Windows) or ./FHInumerics.sh (Linux/Mac).

library(shiny)
library(bslib)

root <- getOption("fhinumerics.root", normalizePath(".."))   # runApp() starts in app/
options(fhinumerics.root = root)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
cfg <- load_config(Sys.getenv("FHINUMERICS_CONFIG", file.path(root, "config.yaml")))
defaults <- utils::modifyList(DEFAULT_SETTINGS, cfg$settings %||% list())
roots <- vapply(Filter(Negate(is.null), cfg$browse_roots %||% list(Home = "~")), path.expand, "")
roots <- roots[dir.exists(roots)]

# Paths typed into the output box may use Windows backslashes (C:\Users\...).
clean_path <- function(p) path.expand(gsub("\\\\", "/", trimws(p)))

# The menu at the top of the left side: one link per page. The page shown is input$page.
PAGES <- c("Outbreak detection" = "outbreak", "Database viewer" = "viewer")
menu_js <- HTML("
$(function() {
  function setPage(page) {
    $('.fhi-menu .nav-link').removeClass('active');
    $('.fhi-menu .nav-link[data-page=\"' + page + '\"]').addClass('active');
    Shiny.setInputValue('page', page);
  }
  $(document).on('click', '.fhi-menu .nav-link', function(e) { e.preventDefault(); setPage($(this).data('page')); });
  $(document).on('shiny:connected', function() { setPage('outbreak'); });
  Shiny.addCustomMessageHandler('set_page', setPage);
  Shiny.addCustomMessageHandler('clear_file_input', function(id) {
    $('#' + id).closest('.input-group').find('input[type=text]').val('');
  });
});")
menu <- tags$ul(class = "nav nav-pills flex-column fhi-menu mb-2",
  lapply(names(PAGES), function(label) tags$li(class = "nav-item",
    tags$a(class = paste("nav-link", if (PAGES[[label]] == "outbreak") "active"),
           href = "#", `data-page` = PAGES[[label]], label))))

outbreak_inputs <- tagList(
  textAreaInput("fasta_text", "Paste sequence(s) in FASTA format", rows = 8, width = "100%",
                placeholder = ">sample1\nTTCAAGAGGGGTCTCCGG...\n>sample2\n..."),
  uiOutput("sent_note"),
  fileInput("fasta_file", "...or choose a FASTA file", accept = c(".fa", ".fasta", ".fas", ".fna", ".txt")),
  selectInput("db", "Database", choices = names(cfg$databases)),
  sliderInput("n_neighbours", "Closest sequences per query", 5, 50, defaults$n_neighbours, step = 5),
  accordion(open = FALSE, accordion_panel(
    "More settings",
    sliderInput("min_coverage", "Minimum % of the query a match must cover", 50, 100, defaults$min_coverage, step = 5),
    selectInput("aligner", "Alignment", ALIGNERS, defaults$aligner),
    checkboxInput("bootstrap", "Bootstrap support for the tree (100 replicates)", defaults$bootstrap)
  )),
  textInput("out_dir", "Save results in", cfg$output_dir),
  shinyFiles::shinyDirButton("browse_out", "Browse...", "Choose output folder"),
  actionButton("run", "Run comparison", class = "btn-primary btn-lg", width = "100%")
)

viewer_inputs <- tagList(
  selectInput("view_db", "Database", choices = names(cfg$databases)),
  uiOutput("date_ui"),
  h6("Choose samples"),
  p(class = "small text-muted",
    "Click rows in the table to choose samples (hold Ctrl or Shift to choose several).",
    "Or paste sample IDs, or choose a text file with one ID per line, and press",
    em("Choose these IDs."), "Filter the table with the boxes under the column names."),
  textAreaInput("id_text", NULL, rows = 3, width = "100%", placeholder = "Paste sample IDs here"),
  fileInput("id_file", NULL, accept = c(".txt", ".csv", ".tsv"), buttonLabel = "Text file...",
            placeholder = "No file chosen"),
  actionButton("add_ids", "Choose these IDs", width = "100%"),
  div(class = "d-flex gap-2 mt-2",
      actionButton("select_shown", "Choose all shown", class = "btn-sm"),
      actionButton("clear_chosen", "Clear chosen", class = "btn-sm")),
  uiOutput("chosen_count"),
  radioButtons("use", "Charts and downloads use",
               c("All samples shown in the table" = "shown", "Only the chosen samples" = "chosen")),
  div(class = "d-flex gap-2 flex-wrap mb-2",
      downloadButton("dl_csv", "Table (CSV)", class = "btn-sm"),
      downloadButton("dl_fasta", "Sequences (FASTA)", class = "btn-sm")),
  actionButton("send", "Send chosen samples to Outbreak detection", class = "btn-primary", width = "100%")
)

# A chart tab: controls on top, then a note and the chart. With `hover`, pointing at the
# chart shows a box (output <id>_hover_box) next to the pointer.
chart_panel <- function(title, id, ..., hover = FALSE) {
  plot <- plotOutput(paste0(id, "_plot"), height = "auto",
                     hover = if (hover) hoverOpts(paste0(id, "_hover"), delay = 80, delayType = "throttle"))
  nav_panel(title,
    layout_columns(col_widths = 4, ...,
                   downloadButton(paste0("dl_", id), "Download PDF", class = "btn-sm")),
    uiOutput(paste0(id, "_note")),
    if (hover) div(style = "position: relative;", plot, uiOutput(paste0(id, "_hover_box"))) else plot)
}

# The box shown when pointing at a point in a plot: a ring around the point, and a table
# of details next to the pointer, on the side with the most room.
hover_box <- function(hover, hit, details) {
  x <- hover$coords_css$x
  y <- hover$coords_css$y
  width <- hover$range$right / (hover$img_css_ratio$x %||% 1)
  height <- hover$range$bottom / (hover$img_css_ratio$y %||% 1)
  side <- paste(
    if (x > width / 2) sprintf("right: %.0fpx;", width - x + 14) else sprintf("left: %.0fpx;", x + 14),
    if (y > height / 2) sprintf("bottom: %.0fpx;", height - y + 14) else sprintf("top: %.0fpx;", y + 14))
  r <- hit$r + 4
  tagList(
    div(class = "fhi-ring", style = sprintf("left: %.0fpx; top: %.0fpx; width: %.0fpx; height: %.0fpx;",
                                            hit$x - r, hit$y - r, 2 * r, 2 * r)),
    div(class = "fhi-tip", style = side,
        tags$table(lapply(names(details), function(k) tags$tr(tags$td(k), tags$td(details[[k]])))))
  )
}

ui <- page_sidebar(
  title = "FHInumerics - HAV sequence comparison",
  fillable = FALSE,   # tables and plots take their own height; the page scrolls
  tags$head(tags$script(menu_js), tags$style(HTML("
    .fhi-tip { position: absolute; z-index: 1000; pointer-events: none; background: white;
               border: 1px solid #bbb; border-radius: 4px; box-shadow: 0 2px 6px rgba(0,0,0,.2);
               padding: 6px 8px; font-size: 12px; max-width: 380px; }
    .fhi-tip td { padding: 0 6px 0 0; vertical-align: top; }
    .fhi-tip td:first-child { color: #666; white-space: nowrap; }
    .fhi-ring { position: absolute; z-index: 999; pointer-events: none; border-radius: 50%;
                border: 3px solid #ff7f00; box-shadow: 0 0 0 3px rgba(255,127,0,.3); }"))),
  sidebar = sidebar(
    width = 380,
    menu,
    hr(class = "my-1"),
    conditionalPanel("input.page !== 'viewer'", outbreak_inputs),
    conditionalPanel("input.page === 'viewer'", viewer_inputs)
  ),
  navset_hidden(
    id = "main",
    nav_panel_hidden("outbreak", navset_card_tab(
      id = "tabs",
      nav_panel("Tree",
        layout_columns(col_widths = c(3, 6, 3),
          selectInput("colour_by", "Colour by", choices = "none"),
          checkboxGroupInput("label_cols", "Show in labels", choices = NULL, inline = TRUE),
          div(checkboxInput("show_support", "Show bootstrap support (>= 70)", TRUE),
              downloadButton("dl_tree", "Download PDF", class = "btn-sm"))
        ),
        uiOutput("tree_ui")
      ),
      nav_panel("SNP distances", uiOutput("snp_note"), DT::DTOutput("snp_table", fill = FALSE)),
      nav_panel("Closest matches", uiOutput("ties_note"), DT::DTOutput("hits_table", fill = FALSE)),
      nav_panel("Run log", verbatimTextOutput("log"))
    )),
    nav_panel_hidden("viewer", navset_card_tab(
      id = "view_tabs",
      nav_panel("Table", uiOutput("missing_note"), DT::DTOutput("db_table", fill = FALSE)),
      chart_panel("Over time", "time",
        selectInput("time_unit", "Bars per", c("Month" = "month", "Year" = "year"), "year"),
        selectInput("time_colour", "Colour by", "none")),
      chart_panel("Counts", "counts",
        selectInput("counts_col", "Count samples per", NULL),
        selectInput("counts_colour", "Colour by", "none")),
      nav_panel("Cross-table",
        layout_columns(col_widths = 4,
          selectInput("ct_rows", "Rows", NULL),
          selectInput("ct_cols", "Columns", NULL)),
        DT::DTOutput("ct_table", fill = FALSE)),
      chart_panel("Map", "map",
        selectInput("map_split", "In the box, count per", "none"),
        checkboxInput("map_zoom", "Zoom to the samples", TRUE), hover = TRUE)
    ))
  )
)

server <- function(input, output, session) {
  result <- reactiveVal(NULL)
  log_lines <- reactiveVal("No analysis run yet.")
  if (is.null(cfg$python) && defaults$aligner == "mafft") {
    showNotification(paste("Python was not found, so MAFFT is not available.",
                           "Run setup, or choose the DECIPHER alignment under More settings."),
                     type = "warning", duration = NULL)
  }
  observeEvent(input$page, nav_select("main", input$page))

  # ==== Outbreak detection =============================================================
  shinyFiles::shinyDirChoose(input, "browse_out", roots = roots, session = session)
  observeEvent(input$browse_out, {
    path <- shinyFiles::parseDirPath(roots, input$browse_out)
    if (length(path)) updateTextInput(session, "out_dir", value = path)
  })

  # A chosen FASTA file is used instead of the pasted text, until samples are sent here
  # from the Database viewer.
  use_file <- reactiveVal(FALSE)
  observeEvent(input$fasta_file, use_file(TRUE))
  # Text sent from the Database viewer; while it is unchanged, each query's own database
  # entry is left out of its matches.
  sent_text <- reactiveVal(NULL)
  from_viewer <- reactive(!use_file() && !is.null(sent_text()) &&
                            identical(trimws(input$fasta_text), trimws(sent_text())))
  output$sent_note <- renderUI({
    if (!from_viewer()) return(NULL)
    p(class = "small text-muted",
      "From the Database viewer: each sample's own database entry is left out of its closest matches.")
  })

  observeEvent(input$run, {
    lines <- character(0)
    log <- function(...) lines <<- c(lines, paste0(format(Sys.time(), "%H:%M:%S  "), ...))
    res <- tryCatch({
      text <- if (use_file()) {
        paste(readLines(input$fasta_file$datapath, warn = FALSE), collapse = "\n")
      } else input$fasta_text
      queries <- parse_fasta_text(text)
      log("Queries: ", paste(names(queries), collapse = ", "))
      cfg$output_dir <- clean_path(input$out_dir)
      settings <- list(n_neighbours = input$n_neighbours, min_coverage = input$min_coverage,
                       aligner = input$aligner, bootstrap = input$bootstrap,
                       search = defaults$search, threads = defaults$threads,
                       exclude_self = from_viewer())
      withProgress(message = "Running", value = 0, {
        steps <- 7
        run_compare(queries, cfg$databases[[input$db]], cfg, settings, log = log,
                    progress = function(step, detail = NULL) incProgress(1 / steps, step, detail))
      })
    }, error = function(e) {
      log("ERROR: ", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      NULL
    })
    log_lines(lines)
    if (!is.null(res)) {
      result(res)
      db <- cfg$databases[[input$db]]
      updateSelectInput(session, "colour_by", choices = c("none", res$display_columns),
                        selected = db$default_colour %||% "none")
      updateCheckboxGroupInput(session, "label_cols", choices = res$display_columns,
                               selected = intersect(db$default_labels %||% character(0), res$display_columns),
                               inline = TRUE)
      showNotification(paste("Done. Results saved in", res$out_dir), type = "message", duration = 8)
    }
  })

  tree_plot <- reactive({
    req(result()$tree)
    plot_tree(result(), input$colour_by, input$label_cols %||% character(0), input$show_support)
  })
  tree_height <- reactive(max(400, 22 * length(result()$tree$tip.label)))
  output$tree_ui <- renderUI({
    if (is.null(result())) return(p("Paste or choose sequences on the left and press Run comparison."))
    if (is.null(result()$tree)) return(p("Too few sequences for a tree."))
    tagList(
      p(class = "small text-muted", "Point at the dot of a sample in the tree to see its details."),
      div(style = "position: relative;",
          plotOutput("tree", height = tree_height(),
                     hover = hoverOpts("tree_hover", delay = 80, delayType = "throttle")),
          uiOutput("tree_tip"))
    )
  })
  output$tree <- renderPlot(tree_plot())
  output$dl_tree <- downloadHandler(
    filename = function() paste0("tree_", input$colour_by, ".pdf"),
    content = function(file) ggplot2::ggsave(file, tree_plot(), width = 11, height = tree_height() / 72, limitsize = FALSE)
  )

  output$tree_tip <- renderUI({
    hit <- tip_at(tree_plot()$data, input$tree_hover)
    if (is.null(hit)) return(NULL)
    hover_box(input$tree_hover, hit, tip_details(result(), hit$label))
  })

  output$snp_note <- renderUI({
    req(result())
    p(sprintf("SNP differences over %d aligned positions; only positions where both sequences have A/C/G/T are counted. Rows and columns follow the tree order.",
              result()$manifest$compared_columns))
  })
  output$snp_table <- DT::renderDT({
    res <- req(result())
    snp <- res$snp
    if (!is.null(res$tree)) {
      d <- ggtree::fortify(res$tree)
      ord <- d$label[d$isTip][order(d$y[d$isTip], decreasing = TRUE)]
      snp <- snp[ord, ord]
    }
    df <- data.frame(sequence = rownames(snp), snp, check.names = FALSE)
    brks <- c(0, 1, 3, 6, 11, 21)
    DT::formatStyle(
      DT::datatable(df, fillContainer = FALSE, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 100, dom = "t")),
      colnames(snp), backgroundColor = DT::styleInterval(brks - 0.5, c("white", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "white")),
      color = DT::styleInterval(c(-0.5, 2.5), c("black", "white", "black"))
    )
  })

  output$ties_note <- renderUI({
    ties <- req(result())$ties
    many <- ties[ties$n_equally_close > ties$n_shown, ]
    if (!NROW(many)) return(NULL)
    div(class = "alert alert-warning",
        lapply(seq_len(nrow(many)), function(i) div(sprintf(
          "%d database sequences are equally close to %s (%.2f%% identity); %d are in the tree. All are listed below.",
          many$n_equally_close[i], many$query[i], many$best_pident[i], many$n_shown[i]))))
  })
  output$hits_table <- DT::renderDT({
    hits <- req(result())$hits
    DT::datatable(hits, fillContainer = FALSE, rownames = FALSE, filter = "top",
                  options = list(scrollX = TRUE, pageLength = 25))
  })

  output$log <- renderText(paste(log_lines(), collapse = "\n"))

  # ==== Database viewer ================================================================
  view_db <- reactive(cfg$databases[[req(input$view_db)]])
  view_meta <- reactive({
    db <- view_db()
    meta <- tryCatch(read_metadata(db), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      NULL
    })
    req(meta)
    list(meta = meta, date_cols = date_columns(meta))
  })
  # Dates of all samples, from the date column chosen in the sidebar (NULL: no dates).
  view_dates <- reactive({
    col <- input$date_col
    if (is.null(col) || !col %in% view_meta()$date_cols) return(NULL)
    parse_partial_date(view_meta()$meta[[col]])
  })
  NO_DATES <- "No column in this database holds dates (like 2023-05-14, 2023-05 or 2023)."

  # Date column: date_column from config.yaml, else a column called "date", else the
  # first column that holds dates.
  output$date_ui <- renderUI({
    cols <- view_meta()$date_cols
    if (!length(cols)) return(p(class = "small text-muted", NO_DATES, "There is no date filter."))
    first <- c(intersect(c(view_db()$date_column, "date"), cols), cols)[1]
    tagList(selectInput("date_col", "Date column", cols, first), uiOutput("date_range_ui"))
  })
  output$date_range_ui <- renderUI({
    dates <- req(view_dates())
    tagList(
      dateRangeInput("date_range", "Sample date", start = min(dates, na.rm = TRUE), end = max(dates, na.rm = TRUE),
                     min = min(dates, na.rm = TRUE), max = max(dates, na.rm = TRUE), weekstart = 1),
      checkboxInput("include_undated", sprintf("Include samples without a date (%d)", sum(is.na(dates))), TRUE)
    )
  })

  # Row numbers (in the metadata) of the samples in the table, after the date filter.
  table_rows <- reactive({
    dates <- view_dates()
    keep <- rep(TRUE, nrow(view_meta()$meta))
    range <- input$date_range
    if (!is.null(dates) && length(range) == 2 && !anyNA(range)) {
      dated <- !is.na(dates)
      keep <- (dated & dates >= range[1] & dates <= range[2]) |
        (!dated & isTRUE(input$include_undated %||% TRUE))
    }
    which(keep)
  })
  table_ids <- reactive(view_meta()$meta$id[table_rows()])

  # Chosen sample IDs. Clicking rows chooses and un-chooses samples in the table; chosen
  # samples outside the date filter stay chosen.
  chosen <- reactiveVal(character(0))
  not_found <- reactiveVal(character(0))
  observeEvent(input$view_db, { chosen(character(0)); not_found(character(0)) })
  table_proxy <- DT::dataTableProxy("db_table")
  select_in_table <- function(ids) DT::selectRows(table_proxy, which(table_ids() %in% ids))

  output$db_table <- DT::renderDT({
    meta <- view_meta()$meta[table_rows(), , drop = FALSE]
    DT::datatable(meta, fillContainer = FALSE, rownames = FALSE, filter = "top",
                  selection = list(mode = "multiple", selected = which(meta$id %in% isolate(chosen()))),
                  options = list(scrollX = TRUE, pageLength = 25))
  }, server = FALSE)
  observeEvent(input$db_table_rows_selected, ignoreNULL = FALSE, {
    ids <- table_ids()
    chosen(union(setdiff(chosen(), ids), ids[input$db_table_rows_selected]))
  })

  observeEvent(input$add_ids, {
    text <- c(input$id_text, if (!is.null(input$id_file)) readLines(input$id_file$datapath, warn = FALSE))
    ids <- match_sample_ids(text, view_meta()$meta$id)
    if (!length(ids$found) && !length(ids$missing)) {
      showNotification("Paste sample IDs or choose a text file first.", type = "warning")
      return()
    }
    not_found(ids$missing)
    chosen(union(chosen(), ids$found))
    select_in_table(chosen())
    showNotification(sprintf("%d samples chosen.", length(ids$found)), type = "message")
  })
  output$missing_note <- renderUI({
    missing <- not_found()
    if (!length(missing)) return(NULL)
    div(class = "alert alert-warning", sprintf(
      "%d of the IDs you gave are not in this database: %s", length(missing),
      paste(c(utils::head(missing, 20), if (length(missing) > 20) "..."), collapse = ", ")))
  })
  observeEvent(input$select_shown, {
    shown <- table_ids()[input$db_table_rows_all]
    chosen(union(chosen(), shown))
    select_in_table(chosen())
  })
  observeEvent(input$clear_chosen, {
    chosen(character(0))
    not_found(character(0))
    DT::selectRows(table_proxy, NULL)
  })

  output$chosen_count <- renderUI({
    n_dated_out <- length(setdiff(chosen(), table_ids()))
    p(class = "small mt-2", strong(length(chosen()), "chosen"),
      sprintf(" - %d shown in the table, %d in the database", length(input$db_table_rows_all),
              nrow(view_meta()$meta)),
      if (n_dated_out) sprintf(" (%d chosen are outside the date range)", n_dated_out))
  })

  # The samples that charts and downloads use: those shown in the table (after all
  # filters), or only the chosen ones.
  use_rows <- reactive({
    meta <- view_meta()$meta
    if (input$use == "chosen") return(which(meta$id %in% chosen()))
    rows <- table_rows()[input$db_table_rows_all %||% seq_along(table_rows())]
    rows[!is.na(rows)]   # the table may still show the previous filter for a moment
  })
  use_meta <- reactive(view_meta()$meta[use_rows(), , drop = FALSE])

  output$dl_csv <- downloadHandler(
    filename = function() paste0("samples_", format(Sys.Date()), ".csv"),
    content = function(file) utils::write.csv(use_meta(), file, row.names = FALSE, na = "")
  )
  output$dl_fasta <- downloadHandler(
    filename = function() paste0("samples_", format(Sys.Date()), ".fasta"),
    content = function(file) write_fasta(db_sequences(view_db(), cfg, use_meta()$id), file)
  )

  observeEvent(input$send, {
    ids <- chosen()
    if (!length(ids)) {
      showNotification("Choose samples first: click rows in the table, or paste sample IDs.", type = "warning")
      return()
    }
    if (length(ids) > 50) {
      showNotification(sprintf("%d samples are chosen. Outbreak detection takes at most 50; choose fewer.",
                               length(ids)), type = "warning", duration = NULL)
      return()
    }
    seqs <- tryCatch(db_sequences(view_db(), cfg, ids), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      NULL
    })
    if (!length(seqs)) return()
    no_seq <- setdiff(ids, names(seqs))
    text <- paste0(">", names(seqs), "\n", seqs, collapse = "\n")
    updateTextAreaInput(session, "fasta_text", value = text)
    updateSelectInput(session, "db", selected = input$view_db)
    sent_text(text)
    use_file(FALSE)
    session$sendCustomMessage("clear_file_input", "fasta_file")
    session$sendCustomMessage("set_page", "outbreak")
    showNotification(paste0(length(seqs), " samples copied to Outbreak detection. Press Run comparison.",
                            if (length(no_seq)) paste0(" No sequence in the database for: ", paste(no_seq, collapse = ", "))),
                     type = if (length(no_seq)) "warning" else "message", duration = if (length(no_seq)) NULL else 8)
  })

  # ---- Charts ----
  observeEvent(view_meta(), {
    db <- view_db()
    cols <- setdiff(names(view_meta()$meta), "id")
    pick <- function(...) intersect(c(...), cols)[1] %||% cols[1]
    colour <- if (isTRUE(db$default_colour %in% cols)) db$default_colour else "none"
    updateSelectInput(session, "time_colour", choices = c("none", cols), selected = colour)
    updateSelectInput(session, "counts_col", choices = cols, selected = pick(db$country_column, db$default_colour))
    updateSelectInput(session, "counts_colour", choices = c("none", cols), selected = "none")
    updateSelectInput(session, "ct_rows", choices = cols, selected = pick(db$default_colour))
    updateSelectInput(session, "ct_cols", choices = cols, selected = pick(db$country_column, rev(cols)))
    updateSelectInput(session, "map_split", choices = c("none", cols),
                      selected = c(intersect(c("genotype", db$default_colour), cols), "none")[1])
  })
  # Chart errors (e.g. "None of these samples has a date") are shown as plain text.
  chart_or_message <- function(expr) tryCatch(expr, error = function(e) validate(need(FALSE, conditionMessage(e))))
  pdf_download <- function(name, plot) downloadHandler(
    filename = function() paste0(name, "_", format(Sys.Date()), ".pdf"),
    content = function(file) ggplot2::ggsave(file, plot(), width = 11, height = 7)
  )
  no_samples <- function() validate(need(nrow(use_meta()) > 0, "No samples to show: choose samples, or change the filters."))

  time_plot <- reactive({
    no_samples()
    dates <- view_dates()
    validate(need(!is.null(dates), NO_DATES))
    chart_or_message(plot_over_time(use_meta(), dates[use_rows()], input$time_unit, input$time_colour))
  })
  output$time_note <- renderUI({
    dates <- view_dates()
    n <- if (is.null(dates)) 0 else sum(is.na(dates[use_rows()]))
    if (n) p(class = "small text-muted", sprintf("%d samples without a date are not shown.", n))
  })
  output$time_plot <- renderPlot(time_plot(), height = 450)
  output$dl_time <- pdf_download("over_time", time_plot)

  counts_plot <- reactive({
    no_samples()
    req(input$counts_col %in% names(use_meta()))
    chart_or_message(plot_counts(use_meta(), input$counts_col, input$counts_colour))
  })
  output$counts_note <- renderUI({
    req(input$counts_col %in% names(use_meta()))
    n <- length(unique(use_meta()[[input$counts_col]]))
    if (n > 30) p(class = "small text-muted", sprintf("The 30 most common of %d values are shown.", n))
  })
  output$counts_plot <- renderPlot(counts_plot(), height = function() {
    n <- if (isTRUE(input$counts_col %in% names(use_meta()))) length(unique(use_meta()[[input$counts_col]])) else 10
    max(300, 24 * min(n, 30) + 100)
  })
  output$dl_counts <- pdf_download("counts", counts_plot)

  output$ct_table <- DT::renderDT({
    no_samples()
    req(input$ct_rows %in% names(use_meta()), input$ct_cols %in% names(use_meta()))
    DT::datatable(cross_table(use_meta(), input$ct_rows, input$ct_cols), fillContainer = FALSE, rownames = FALSE,
                  options = list(scrollX = TRUE, pageLength = 50))
  })

  map_result <- reactive({
    no_samples()
    column <- view_db()$country_column
    validate(need(column %in% names(use_meta()),
                  paste0("This database has no country column ('", column, "'; set country_column in config.yaml).")))
    chart_or_message(plot_country_map(use_meta(), column, input$map_zoom, input$map_split %||% "none"))
  })
  map_radius <- reactive(point_radius_css(map_result()$plot, 2))
  output$map_note <- renderUI({
    unknown <- map_result()$unknown
    tagList(
      p(class = "small text-muted", "Point at a circle to see the country and its number of samples."),
      if (length(unknown)) div(class = "alert alert-warning", "Not on the map (country name not recognised): ",
                               paste(unknown, collapse = ", ")))
  })
  output$map_plot <- renderPlot(map_result()$plot, height = 550)
  output$dl_map <- pdf_download("map", reactive(map_result()$plot))
  output$map_hover_box <- renderUI({
    circles <- map_result()$circles
    hit <- point_at(circles$long, circles$lat, input$map_hover, radius = map_radius() + 2)
    if (is.null(hit)) return(NULL)
    hit$r <- map_radius()[hit$i]
    hover_box(input$map_hover, hit, circle_details(circles[hit$i, ], input$map_split %||% "none"))
  })

  # When started from the launcher, stop the app once the browser tab is closed.
  session$onSessionEnded(function() {
    if (isTRUE(getOption("fhinumerics.stop_on_close"))) stopApp()
  })
}

shinyApp(ui, server)
