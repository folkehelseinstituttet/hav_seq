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
PAGES <- c("Outbreak detection" = "outbreak", "Database viewer" = "viewer", "Alignment" = "alignment")
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
  // SNP table: pointing at a number with a data-tip shows that text in a box next to
  // the pointer (on the side with the most room).
  $(document).on('mouseenter', 'td[data-tip]', function() { $('#fhi-cell-tip').text($(this).attr('data-tip')).show(); });
  $(document).on('mousemove', 'td[data-tip]', function(e) {
    var tip = $('#fhi-cell-tip'), right = e.clientX > window.innerWidth / 2, below = e.clientY < window.innerHeight / 2;
    tip.css({left: right ? '' : e.clientX + 14, right: right ? window.innerWidth - e.clientX + 14 : '',
             top: below ? e.clientY + 14 : '', bottom: below ? '' : window.innerHeight - e.clientY + 14});
  });
  $(document).on('mouseleave', 'td[data-tip]', function() { $('#fhi-cell-tip').hide(); });
  // SNP table: a click on a name chooses that row for the Alignment page; a click on a
  // number in a query's row or column opens that pair there.
  $(document).on('click', '#snp_table td.fhi-pick', function() {
    $(this).toggleClass('fhi-picked');
    Shiny.setInputValue('snp_rows', $('#snp_table td.fhi-picked').map(function() { return $(this).attr('data-name'); }).get(),
                        {priority: 'event'});
  });
  $(document).on('click', '#snp_table td[data-pair-row]', function() {
    Shiny.setInputValue('snp_pair', {row: $(this).attr('data-pair-row'), col: $(this).attr('data-pair-col')}, {priority: 'event'});
  });
  Shiny.addCustomMessageHandler('snp_clear', function(x) { $('#snp_table td.fhi-picked').removeClass('fhi-picked'); });
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

alignment_inputs <- tagList(
  radioButtons("aln_source", "Show", c("Results of the last run" = "run", "An earlier run folder" = "folder",
                                        "An aligned FASTA file" = "file")),
  conditionalPanel("input.aln_source === 'run'", uiOutput("aln_run_ui")),
  conditionalPanel("input.aln_source === 'folder'",
    shinyFiles::shinyDirButton("aln_browse", "Choose a run folder...", "Choose the folder of one run"),
    uiOutput("aln_folder_ui")),
  conditionalPanel("input.aln_source === 'file'",
    fileInput("aln_file", NULL, accept = c(".fa", ".fasta", ".fas", ".fna", ".aln", ".txt"),
              buttonLabel = "Aligned FASTA...", placeholder = "No file chosen")),
  uiOutput("aln_rows_ui"),
  selectInput("aln_anchor", "Compare with", c("Majority consensus" = "")),
  checkboxInput("aln_variable", "Variable columns only"),
  checkboxInput("aln_hide_gaps", "Hide columns that are a gap in the anchor"),
  accordion(open = FALSE, accordion_panel(
    "Consensus settings",
    radioButtons("aln_ties", "When the most common base is not clear (a tie)",
                 c("Ambiguity code (A/G tie = R)" = "iupac", "N" = "n")),
    checkboxInput("aln_use_threshold", "The most common base must reach a threshold"),
    conditionalPanel("input.aln_use_threshold", sliderInput("aln_threshold", NULL, 50, 100, 75, step = 5, post = "%"))
  )),
  h6("Download", class = "mt-3"),
  div(class = "d-flex gap-2 flex-wrap",
      downloadButton("aln_dl_fasta", "FASTA (whole)", class = "btn-sm"),
      downloadButton("aln_dl_region", "FASTA (columns on screen)", class = "btn-sm"),
      downloadButton("aln_dl_pdf", "PDF (on screen)", class = "btn-sm"),
      actionButton("aln_png", "Picture (PNG)", icon = icon("download"), class = "btn-sm"),
      downloadButton("aln_dl_html", "Interactive page (HTML)", class = "btn-sm"))
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
  tags$head(tags$script(menu_js), tags$script(src = "msa.js"), tags$link(rel = "stylesheet", href = "msa.css"),
            tags$style(HTML("
    .fhi-tip { position: absolute; z-index: 1000; pointer-events: none; background: white;
               border: 1px solid #bbb; border-radius: 4px; box-shadow: 0 2px 6px rgba(0,0,0,.2);
               padding: 6px 8px; font-size: 12px; max-width: 380px; }
    .fhi-tip td { padding: 0 6px 0 0; vertical-align: top; }
    .fhi-tip td:first-child { color: #666; white-space: nowrap; }
    #fhi-cell-tip { position: fixed; display: none; white-space: pre-line; }
    td[data-tip]:hover { box-shadow: inset 0 0 0 3px #ff7f00; }
    td.fhi-flag { position: relative; }
    td.fhi-flag::after { content: ''; position: absolute; top: 0; right: 0; border-style: solid;
                         border-width: 0 9px 9px 0; border-color: transparent #ff7f00 transparent transparent; }
    .fhi-ring { position: absolute; z-index: 999; pointer-events: none; border-radius: 50%;
                border: 3px solid #ff7f00; box-shadow: 0 0 0 3px rgba(255,127,0,.3); }"))),
  div(id = "fhi-cell-tip", class = "fhi-tip"),
  sidebar = sidebar(
    width = 380,
    menu,
    hr(class = "my-1"),
    conditionalPanel("input.page !== 'viewer' && input.page !== 'alignment'", outbreak_inputs),
    conditionalPanel("input.page === 'viewer'", viewer_inputs),
    conditionalPanel("input.page === 'alignment'", alignment_inputs)
  ),
  navset_hidden(
    id = "main",
    nav_panel_hidden("outbreak", uiOutput("view_ui"), navset_card_tab(
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
      nav_panel("SNP distances", uiOutput("snp_note"), uiOutput("snp_pick_ui"), DT::DTOutput("snp_table", fill = FALSE)),
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
    )),
    # Not a card: the viewer's control bar stays at the top while the page scrolls.
    nav_panel_hidden("alignment", div(class = "px-1",
      uiOutput("aln_status"),
      div(id = "msa")
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
        steps <- 3 + 3 * (length(queries) + (length(queries) > 1))   # see run_compare()
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

  # With several queries, there is a result for all queries together and one for each
  # query alone; the "Show results for" list above the tabs picks the one shown.
  output$view_ui <- renderUI({
    views <- req(result())$views
    if (length(views) < 2) return(NULL)
    failed <- vapply(views, function(v) !is.null(v$error), TRUE)
    choices <- stats::setNames(names(views), paste0(names(views), ifelse(failed, " (no result)", "")))
    div(class = "d-flex align-items-center gap-3 mb-2",
        selectInput("view", "Show results for", choices, names(views)[!failed][1], width = "360px"),
        p(class = "small text-muted mb-0",
          "Each query is also compared on its own: its own closest sequences, over its whole length."))
  })
  view <- reactive({
    views <- req(result())$views
    v <- views[[input$view %||% ""]]
    if (is.null(v)) v <- Filter(function(x) is.null(x$error), views)[[1]]   # just after a new run
    v
  })
  view_error <- function() if (!is.null(view()$error)) div(class = "alert alert-warning", view()$error)
  # Search hits and tie notes of the queries in the shown result.
  view_rows <- function(tab) if (identical(view()$label, ALL_QUERIES)) tab else tab[tab$query %in% view()$label, , drop = FALSE]

  tree_plot <- reactive({
    req(view()$tree)
    plot_tree(view(), input$colour_by, input$label_cols %||% character(0), input$show_support)
  })
  tree_height <- reactive(max(400, 22 * length(view()$tree$tip.label)))
  output$tree_ui <- renderUI({
    if (is.null(result())) return(p("Paste or choose sequences on the left and press Run comparison."))
    if (!is.null(view()$error)) return(view_error())
    if (is.null(view()$tree)) return(p("Too few sequences for a tree."))
    tagList(
      p(class = "small text-muted", "Point at the dot of a sample in the tree to see its details.",
        "Drag up or down over the tree to choose samples for the Alignment page."),
      uiOutput("tree_pick_ui"),
      div(style = "position: relative;",
          plotOutput("tree", height = tree_height(),
                     hover = hoverOpts("tree_hover", delay = 80, delayType = "throttle"),
                     brush = brushOpts("tree_brush", direction = "y", resetOnNew = TRUE, fill = "#1a73e8")),
          uiOutput("tree_tip"))
    )
  })
  output$tree <- renderPlot(tree_plot())
  output$dl_tree <- downloadHandler(
    filename = function() paste0("tree_", gsub("[^A-Za-z0-9._-]+", "_", view()$label), "_", input$colour_by, ".pdf"),
    content = function(file) ggplot2::ggsave(file, tree_plot(), width = 11, height = tree_height() / 72, limitsize = FALSE)
  )

  output$tree_tip <- renderUI({
    hit <- tip_at(tree_plot()$data, input$tree_hover)
    if (is.null(hit)) return(NULL)
    hover_box(input$tree_hover, hit, tip_details(view(), hit$label))
  })

  output$snp_note <- renderUI({
    v <- view()
    if (!is.null(v$error)) return(view_error())
    tagList(
      p(sprintf(paste("SNP differences over %d aligned positions. A position counts when the two bases cannot be the same:",
                      "ambiguity codes count only when they do not fit (K = G/T differs from A, M = A/C does not); N and gaps never count.",
                      "Rows and columns follow the tree order. Queries are in red."), v$compared_columns)),
      p(class = "small text-muted",
        "Point at a number in a query's row or column to see which positions differ.",
        "An orange corner marks pairs that also differ by gaps, Ns or ambiguity codes that fit (not counted as SNPs)."))
  })
  output$snp_table <- DT::renderDT({
    v <- view()
    snp <- req(v$snp)
    if (!is.null(v$tree)) {
      d <- ggtree::fortify(v$tree)
      ord <- d$label[d$isTip][order(d$y[d$isTip], decreasing = TRUE)]
      snp <- snp[ord, ord]
    }
    df <- data.frame(sequence = rownames(snp), snp, check.names = FALSE)
    brks <- c(0, 1, 3, 6, 11, 21)
    # Query columns: the header is coloured like the query rows (0-based, after "sequence").
    query_cols <- which(colnames(snp) %in% v$queries)
    header_js <- DT::JS(sprintf(
      "function(thead) { $(thead).find('th').each(function(i) { if ([%s].indexOf(i) >= 0) $(this).css({color: '#d7301f', 'background-color': '#fde0dd'}); }); }",
      paste(query_cols, collapse = ",")))
    # Hover text of the query cells (R/snps.R), attached to each cell as it is drawn;
    # cells where the pair also differs by gaps, Ns or ambiguity codes get a corner mark.
    tips <- snp_tips(v, unlist(result()$manifest$reverse_complemented))
    cell_js <- DT::JS(sprintf(
      "(function() { var tips = %s, flags = %s, cols = %s; return function(td, cellData, rowData, row, col) { var key = rowData[0] + '\\t' + cols[col]; if (col === 0) { td.classList.add('fhi-pick'); td.setAttribute('data-name', rowData[0]); } if (tips[key]) { td.setAttribute('data-tip', tips[key]); td.setAttribute('data-pair-row', rowData[0]); td.setAttribute('data-pair-col', cols[col]); } if (flags[key]) td.classList.add('fhi-flag'); }; })()",
      jsonlite::toJSON(as.list(tips$text), auto_unbox = TRUE),
      jsonlite::toJSON(as.list(tips$flagged[tips$flagged]), auto_unbox = TRUE), jsonlite::toJSON(colnames(df))))
    # All rows on one page (no paging), so every sequence can be reached and chosen.
    table <- DT::datatable(df, fillContainer = FALSE, rownames = FALSE, selection = "none",
                           options = list(scrollX = TRUE, paging = FALSE, dom = "t", headerCallback = header_js,
                                          columnDefs = list(list(targets = "_all", createdCell = cell_js))))
    table <- DT::formatStyle(table, colnames(snp),
      backgroundColor = DT::styleInterval(brks - 0.5, c("white", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "white")),
      color = DT::styleInterval(c(-0.5, 2.5), c("black", "white", "black")))
    DT::formatStyle(table, "sequence", fontWeight = DT::styleEqual(v$queries, "bold", default = "normal"),
                    color = DT::styleEqual(v$queries, "#d7301f", default = "black"),
                    backgroundColor = DT::styleEqual(v$queries, "#fde0dd", default = "white"))
  })

  output$ties_note <- renderUI({
    ties <- view_rows(req(result())$ties)
    many <- ties[ties$n_equally_close > ties$n_shown, ]
    if (!NROW(many)) return(NULL)
    div(class = "alert alert-warning",
        lapply(seq_len(nrow(many)), function(i) div(sprintf(
          "%d database sequences are equally close to %s (%.2f%% identity); %d are in the tree. All are listed below.",
          many$n_equally_close[i], many$query[i], many$best_pident[i], many$n_shown[i]))))
  })
  output$hits_table <- DT::renderDT({
    hits <- view_rows(req(result())$hits)
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

  # ==== Alignment ======================================================================
  # The page shows one alignment (aln_base): a result of the last run, a result in an
  # earlier run folder, or an aligned FASTA file. aln_pick holds what was chosen for it
  # elsewhere (rows, anchor), valid only while the same alignment is shown (same key).
  aln_pick <- reactiveVal(list(key = ""))
  pending_pair <- NULL   # a pair from the SNP table: go to its first difference once
  aln_try <- function(expr) tryCatch(expr, error = function(e) {
    showNotification(conditionMessage(e), type = "error", duration = NULL)
    NULL
  })

  output$aln_run_ui <- renderUI({
    views <- result()$views
    if (is.null(views)) return(p(class = "small text-muted", "No comparison has been run yet. Run one on the",
                                 em("Outbreak detection"), "page, or choose another source."))
    ok <- names(views)[vapply(views, function(v) is.null(v$error), TRUE)]
    selectInput("aln_view", "Result", ok, intersect(isolate(view()$label), ok)[1])
  })
  # The Alignment page follows "Show results for", and the other way round.
  observeEvent(input$view, updateSelectInput(session, "aln_view", selected = input$view))
  observeEvent(input$aln_view, {
    if (!is.null(input$view) && !identical(input$view, input$aln_view)) updateSelectInput(session, "view", selected = input$aln_view)
  })

  aln_roots <- function() {
    out <- c(Results = path.expand(cfg$output_dir), roots)
    out[dir.exists(out) & !duplicated(out)]
  }
  shinyFiles::shinyDirChoose(input, "aln_browse", roots = aln_roots, session = session)
  aln_folder <- reactiveVal(NULL)
  observeEvent(input$aln_browse, {
    path <- shinyFiles::parseDirPath(aln_roots(), input$aln_browse)
    if (!length(path)) return()
    folders <- aln_try(msa_run_results(path))
    if (!is.null(folders)) { attr(folders, "dir") <- path; aln_folder(folders) }
  })
  output$aln_folder_ui <- renderUI({
    f <- req(aln_folder())
    tagList(p(class = "small mt-2 mb-1", basename(attr(f, "dir"))),
            selectInput("aln_folder_result", "Result", names(f)))
  })

  aln_base <- reactive({
    src <- input$aln_source %||% "run"
    if (src == "run") {
      res <- req(result())
      v <- res$views[[req(input$aln_view)]]
      req(v, is.null(v$error))
      list(key = paste("run", res$out_dir, v$label), a = msa_from_result(v, unlist(res$manifest$reverse_complemented)))
    } else if (src == "folder") {
      f <- req(aln_folder())
      i <- req(match(input$aln_folder_result, names(f)))
      a <- aln_try(msa_load_run_folder(f[[i]], paste0(basename(attr(f, "dir")), ": ", names(f)[i]),
                                       attr(f, "manifest"), attr(f, "regions")[[i]], attr(f, "dropped")[[i]]))
      list(key = paste("folder", f[[i]]), a = req(a))
    } else {
      file <- req(input$aln_file)
      a <- aln_try(msa_load(msa_read_fasta(file$datapath), file$name))
      list(key = paste("file", file$name, file$size), a = req(a))
    }
  })
  # The alignment as shown: all rows, or the rows chosen for it.
  aln_shown <- reactive({
    b <- aln_base()
    pick <- aln_pick()
    if (!identical(pick$key, b$key) || is.null(pick$rows)) return(list(key = b$key, a = b$a, pick = NULL))
    list(key = b$key, a = msa_subset(b$a, pick$rows), pick = pick)
  })

  # Anchor: none (the majority consensus) unless chosen; right-click a row, the "Compare
  # with" list, or a pair from the SNP table.
  aln_anchor <- reactiveVal(NULL)
  observeEvent(aln_shown(), priority = 10, {
    s <- aln_shown()
    anchor <- s$pick$anchor
    aln_anchor(anchor)
    updateSelectInput(session, "aln_anchor", choices = c("Majority consensus" = "", rownames(s$a$m)),
                      selected = anchor %||% "")
  })
  observeEvent(input$aln_anchor, aln_anchor(if (nzchar(input$aln_anchor)) input$aln_anchor), ignoreInit = TRUE)
  observeEvent(input$msa_anchor_pick, {
    name <- input$msa_anchor_pick$name
    aln_anchor(if (nzchar(name)) name)
    updateSelectInput(session, "aln_anchor", selected = name)
  })

  aln_view_data <- reactive({
    a <- aln_shown()$a
    threshold <- if (isTRUE(input$aln_use_threshold)) input$aln_threshold
    msa_view(a, aln_anchor(), input$aln_ties %||% "iupac", threshold)
  })
  observeEvent(aln_shown(), priority = 5, session$sendCustomMessage("msa_data", msa_payload(aln_shown()$a)))
  observeEvent(aln_view_data(), {
    a <- aln_shown()$a
    v <- aln_view_data()
    msg <- msa_view_payload(a, v)
    if (!is.null(pending_pair) && identical(v$anchor, pending_pair$query)) {
      msg$goto <- msa_first_difference(a, v, pending_pair$other)
      msg$select <- pending_pair$other
      pending_pair <<- NULL
    }
    session$sendCustomMessage("msa_view", msg)
  })
  observe(session$sendCustomMessage("msa_options", list(variable = isTRUE(input$aln_variable),
                                                        hide_gaps = isTRUE(input$aln_hide_gaps))))

  output$aln_status <- renderUI({
    if (is.null(tryCatch(aln_base(), error = function(e) NULL))) {
      return(p("Run a comparison on the", em("Outbreak detection"), "page, or choose an earlier run folder",
               "or an aligned FASTA file on the left."))
    }
    s <- aln_shown()
    a <- s$a
    v <- aln_view_data()
    n_all <- nrow(aln_base()$a$m)
    region <- if (!is.null(a$compared)) {
      cols <- a$cols[a$compared]
      sprintf(" \u00b7 columns %d\u2013%d compared in the run (the SNP table); the rest is greyed", min(cols), max(cols))
    }
    mode <- if (is.null(v$anchor)) {
      div(class = "msa-mode", "Compared with: ", strong("the majority consensus"), " (no anchor sequence chosen).",
          " To compare with one sequence, right-click its name and choose ", em("Use ... as anchor"),
          ", or pick it under ", em("Compare with"), " on the left.")
    } else {
      div(class = "msa-mode", "Compared with: ", strong(v$anchor), " (the anchor, top row).",
          " Right-click and choose ", em("Compare with the majority consensus"), " to go back.")
    }
    note <- function(...) div(class = "msa-note", ...)
    tagList(
      div(class = "msa-head", h5(a$label), div(class = "msa-sub", sprintf(
        "%s of %d sequences \u00b7 %d columns%s \u00b7 the number after each name: SNPs to the anchor, as in the SNP table",
        nrow(a$m), n_all, ncol(a$m), region %||% ""))),
      mode,
      if (!is.null(s$pick$note)) note(s$pick$note),
      if (length(a$dropped)) note(sprintf(
        "%d sequences (grey names) cover less than half of the compared region and were left out of the SNP table and tree.",
        length(a$dropped))),
      if (!is.null(a$note)) note(a$note)
    )
  })

  # Query and its closest sequences (run results).
  output$aln_rows_ui <- renderUI({
    b <- tryCatch(aln_base(), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    tagList(
      if (!is.null(b$a$snp) && length(b$a$queries)) tagList(
        h6("Show a query and its closest"),
        selectInput("aln_query", "Query", b$a$queries),
        numericInput("aln_n", "Number of closest sequences", 10, min = 1, max = 500),
        actionButton("aln_closest", "Show the query and its closest", class = "btn-sm mb-2")),
      actionButton("aln_all", "Show all sequences", class = "btn-sm mb-3")
    )
  })
  observeEvent(input$aln_closest, {
    b <- aln_base()
    n <- max(1, as.integer(input$aln_n %||% 10))
    cl <- msa_closest(b$a$snp, input$aln_query, n)
    note <- sprintf("Showing %s and its %d closest sequences (by SNPs in the SNP table).", input$aln_query, n)
    if (cl$extra) note <- paste(note, sprintf("%d more are as close as number %d and are shown too.", cl$extra, n))
    aln_pick(list(key = b$key, rows = cl$ids, note = note))
  })
  observeEvent(input$aln_all, aln_pick(list(key = "")))

  # Picks from the run results: rows of the SNP table, a box over the tree, a pair.
  show_in_alignment <- function(ids, anchor = NULL, note = NULL) {
    res <- result()
    v <- view()
    aln_pick(list(key = paste("run", res$out_dir, v$label), rows = ids, anchor = anchor, note = note))
    updateRadioButtons(session, "aln_source", selected = "run")
    updateSelectInput(session, "aln_view", selected = v$label)
    session$sendCustomMessage("set_page", "alignment")
  }
  snp_picked <- reactiveVal(character(0))
  observeEvent(input$snp_rows, snp_picked(unlist(input$snp_rows)), ignoreNULL = FALSE)
  observeEvent(view(), snp_picked(character(0)))
  output$snp_pick_ui <- renderUI({
    ids <- snp_picked()
    if (!length(ids)) {
      return(p(class = "small text-muted", "Click names to choose sequences for the Alignment page.",
               "Click a number in a query's row or column to see that pair in the Alignment page."))
    }
    div(class = "d-flex align-items-center gap-2 mb-2", strong(length(ids), "chosen"),
        actionButton("snp_show", "Show in Alignment", class = "btn-sm btn-primary"),
        actionButton("snp_unpick", "Clear", class = "btn-sm"))
  })
  observeEvent(input$snp_show, show_in_alignment(snp_picked(), note = sprintf(
    "Showing the %d sequences chosen in the SNP table.", length(snp_picked()))))
  observeEvent(input$snp_unpick, { snp_picked(character(0)); session$sendCustomMessage("snp_clear", TRUE) })
  observeEvent(input$snp_pair, {
    r <- input$snp_pair$row
    k <- input$snp_pair$col
    q <- if (r %in% view()$queries) r else k
    other <- setdiff(c(r, k), q)
    pending_pair <<- list(query = q, other = other)
    show_in_alignment(c(q, other), anchor = q, note = sprintf(
      "Showing %s and %s, compared with %s, at their first SNP (if any).", q, other, q))
  })

  tree_brushed <- reactive({
    b <- input$tree_brush
    if (is.null(b)) return(character(0))
    d <- tree_plot()$data
    tips <- d[d$isTip, ]
    tips$label[tips$y >= b$ymin & tips$y <= b$ymax]
  })
  output$tree_pick_ui <- renderUI({
    ids <- tree_brushed()
    if (!length(ids)) return(NULL)
    div(class = "d-flex align-items-center gap-2 mb-2", strong(length(ids), "chosen in the tree"),
        actionButton("tree_show", "Show in Alignment", class = "btn-sm btn-primary"))
  })
  observeEvent(input$tree_show, show_in_alignment(tree_brushed(), note = sprintf(
    "Showing the %d sequences chosen in the tree.", length(tree_brushed()))))

  # Downloads: rows in the order shown (anchor first).
  # The order the viewer shows (its "Sort rows"), else R's order.
  aln_order <- function() {
    a <- aln_shown()$a
    v <- aln_view_data()
    shown <- unlist(input$msa_window$order)
    rows <- if (setequal(shown, v$rows)) shown else v$rows
    c(if (!is.null(v$anchor)) match(v$anchor, rownames(a$m)), rows)
  }
  aln_file_name <- function(ext, cols = NULL) {
    a <- aln_shown()$a
    paste0("alignment_", gsub("[^A-Za-z0-9._-]+", "_", sub("\\.[A-Za-z]+$", "", a$label)),
           if (!is.null(cols)) sprintf("_columns_%d-%d", a$cols[min(cols)], a$cols[max(cols)]), ext)
  }
  window_cols <- function() {
    cols <- unlist(input$msa_window$cols)
    validate(need(length(cols), "Open the Alignment page first."))
    cols
  }
  output$aln_dl_fasta <- downloadHandler(
    filename = function() aln_file_name(".fasta"),
    content = function(file) write_fasta(msa_fasta(aln_shown()$a, aln_order()), file)
  )
  output$aln_dl_region <- downloadHandler(
    filename = function() aln_file_name(".fasta", window_cols()),
    content = function(file) {
      cols <- window_cols()
      write_fasta(msa_fasta(aln_shown()$a, aln_order(), min(cols), max(cols)), file)
    }
  )
  output$aln_dl_pdf <- downloadHandler(
    filename = function() aln_file_name(".pdf", window_cols()),
    content = function(file) msa_pdf(aln_shown()$a, aln_view_data(), unlist(input$msa_window$rows), window_cols(), file)
  )
  observeEvent(input$aln_png, session$sendCustomMessage("msa_png", aln_file_name(".png")))
  output$aln_dl_html <- downloadHandler(
    filename = function() aln_file_name(".html"),
    content = function(file) msa_html(aln_shown()$a, aln_view_data(), file, sort = input$msa_window$sort)
  )

  # When started from the launcher, stop the app once the browser tab is closed.
  session$onSessionEnded(function() {
    if (isTRUE(getOption("fhinumerics.stop_on_close"))) stopApp()
  })
}

shinyApp(ui, server)
