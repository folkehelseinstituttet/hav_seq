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
PAGES <- c("Outbreak detection" = "outbreak", "Database viewer" = "viewer", "Alignment" = "alignment", "Sets" = "sets")
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
  // Sets not saved: the browser asks before the tab is closed (closing it stops FHInumerics).
  window.fhiUnsaved = false;
  Shiny.addCustomMessageHandler('sets_unsaved', function(x) { window.fhiUnsaved = !!x; });
  window.addEventListener('beforeunload', function(e) { if (window.fhiUnsaved) { e.preventDefault(); e.returnValue = ''; } });
  // The background alignment: the box under the menu, on every page.
  Shiny.addCustomMessageHandler('align_job', function(j) {
    var box = $('#fhi-job');
    if (!j || !j.label) { box.hide(); return; }
    box.find('.fhi-job-label').text(j.label);
    box.find('.progress-bar').css('width', Math.round(100 * j.value) + '%');
    box.find('.fhi-job-text').text(j.text);
    box.show();
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
  uiOutput("run_set_ui"),
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
  actionButton("send", "Send chosen samples to Outbreak detection", class = "btn-primary", width = "100%"),
  actionButton("view_add", "Add chosen samples to the set", class = "mt-2", width = "100%"),
  uiOutput("view_set_line")
)

alignment_inputs <- tagList(
  radioButtons("aln_source", "Show", c("Results of the last run" = "run", "An earlier run folder" = "folder",
                                        "A FASTA file (aligned or not)" = "file", "A set" = "set")),
  conditionalPanel("input.aln_source === 'run'", uiOutput("aln_run_ui")),
  conditionalPanel("input.aln_source === 'folder'",
    shinyFiles::shinyDirButton("aln_browse", "Choose a run folder...", "Choose the folder of one run"),
    uiOutput("aln_folder_ui")),
  conditionalPanel("input.aln_source === 'file'",
    fileInput("aln_file", NULL, accept = c(".fa", ".fasta", ".fas", ".fna", ".aln", ".txt"),
              buttonLabel = "FASTA file...", placeholder = "No file chosen"),
    uiOutput("aln_file_ui")),
  conditionalPanel("input.aln_source === 'set'", uiOutput("aln_set_ui")),
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

sets_inputs <- tagList(
  uiOutput("set_list_ui"),
  div(class = "d-flex gap-2 flex-wrap mb-2",
      actionButton("set_new", "New set", class = "btn-sm"),
      shinyFiles::shinyDirButton("set_open", "Open a set folder...", "Choose the folder of a saved set", class = "btn-sm")),
  hr(class = "my-2"),
  textInput("set_name", "Name of the open set", width = "100%"),
  textInput("set_dir", "Save in", width = "100%"),
  shinyFiles::shinyDirButton("set_browse", "Browse...", "Choose the folder to save the set in"),
  div(class = "d-flex gap-2 mt-2",
      actionButton("set_save", "Save", class = "btn-primary"),
      actionButton("set_delete", "Delete")),
  h6("Add sequences", class = "mt-3"),
  fileInput("set_file", NULL, accept = c(".fa", ".fasta", ".fas", ".fna", ".aln", ".txt"),
            buttonLabel = "FASTA file...", placeholder = "aligned or not"),
  p(class = "small text-muted",
    "Or press", em("Add to set"), "on the Database viewer, Outbreak detection or Alignment page."),
  h6("Alignment"),
  uiOutput("set_align_ui"),
  h6("Compare with a database", class = "mt-3"),
  actionButton("set_use", "Use in Outbreak detection", class = "btn-sm"),
  p(class = "small text-muted mt-1", "Uses the rows chosen in the table, or the whole set (at most 50)."),
  h6("Download", class = "mt-3"),
  div(class = "d-flex gap-2 flex-wrap",
      downloadButton("set_dl_fasta", "Sequences (FASTA)", class = "btn-sm"),
      downloadButton("set_dl_csv", "Table (CSV)", class = "btn-sm"))
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
                border: 3px solid #ff7f00; box-shadow: 0 0 0 3px rgba(255,127,0,.3); }
    #fhi-job { display: none; margin: 4px 0 8px; padding: 6px 8px; border: 1px solid #d8dee4; border-radius: 6px;
               background: #f6f7f9; font-size: 13px; }
    #fhi-job .progress { height: 8px; margin: 4px 0; }"))),
  div(id = "fhi-cell-tip", class = "fhi-tip"),
  sidebar = sidebar(
    width = 380,
    menu,
    # The background alignment (sets, unaligned FASTA files), shown on every page.
    div(id = "fhi-job",
        div(strong("Aligning "), span(class = "fhi-job-label")),
        div(class = "progress", div(class = "progress-bar", style = "width: 0%")),
        div(class = "d-flex justify-content-between align-items-center gap-2",
            span(class = "fhi-job-text small text-muted"),
            actionButton("align_cancel", "Cancel", class = "btn-sm btn-outline-danger"))),
    hr(class = "my-1"),
    conditionalPanel("input.page !== 'viewer' && input.page !== 'alignment' && input.page !== 'sets'", outbreak_inputs),
    conditionalPanel("input.page === 'viewer'", viewer_inputs),
    conditionalPanel("input.page === 'alignment'", alignment_inputs),
    conditionalPanel("input.page === 'sets'", sets_inputs)
  ),
  navset_hidden(
    id = "main",
    nav_panel_hidden("outbreak", uiOutput("view_ui"), uiOutput("add_result_ui"), navset_card_tab(
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
    )),
    nav_panel_hidden("sets", div(class = "px-1",
      uiOutput("set_head"),
      uiOutput("set_rows_ui"),
      DT::DTOutput("set_table", fill = FALSE)
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
  observeEvent(input$fasta_file, { use_file(TRUE); use_set_in_run("") })
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
      s <- run_set()
      if (!is.null(s)) {
        q <- set_queries(s, input$db, set_input()$keys)
        queries <- q$seqs
        exclude <- q$exclude_self
        log("Queries from the set '", s$name, "'")
      } else {
        text <- if (use_file()) {
          paste(readLines(input$fasta_file$datapath, warn = FALSE), collapse = "\n")
        } else input$fasta_text
        queries <- parse_fasta_text(text)
        exclude <- from_viewer()
      }
      log("Queries: ", paste(names(queries), collapse = ", "))
      cfg$output_dir <- clean_path(input$out_dir)
      settings <- list(n_neighbours = input$n_neighbours, min_coverage = input$min_coverage,
                       aligner = input$aligner, bootstrap = input$bootstrap,
                       search = defaults$search, threads = defaults$threads,
                       exclude_self = exclude)
      withProgress(message = "Running", value = 0, {
        steps <- 3 + 3 * (length(queries) + (length(queries) > 1))   # see run_compare()
        run_compare(queries, cfg$databases[[input$db]], cfg, settings, log = log,
                    progress = function(step, detail = NULL) incProgress(1 / steps, step, detail),
                    input = if (!is.null(s)) q$input)
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
  # A set as the input, used instead of the pasted text and the file: chosen in the list
  # under the file button, or with "Use in Outbreak detection" on the Sets page (then only
  # the rows chosen there, if any: `keys`). Database samples from the searched database
  # leave out their own entry (set_queries() in R/sets.R).
  set_input <- reactiveVal(list(id = "", keys = NULL))
  use_set_in_run <- function(id, keys = NULL) {
    set_input(list(id = id, keys = keys))
    updateSelectInput(session, "run_set", selected = id)
  }
  run_set <- reactive(if (nzchar(set_input()$id)) sets()[[set_input()$id]])
  output$run_set_ui <- renderUI({
    ch <- set_choices()
    if (!length(ch)) return(NULL)
    tagList(selectInput("run_set", "...or use a set", c("(none)" = "", ch), isolate(set_input()$id)),
            uiOutput("run_set_note"))
  })
  observeEvent(input$run_set, {
    v <- input$run_set
    if (identical(v, set_input()$id)) return()
    if (startsWith(v, "saved:")) {
      open_folder(sub("^saved:", "", v))
      v <- open_id() %||% ""
    }
    use_set_in_run(v)
  })
  output$run_set_note <- renderUI({
    s <- run_set()
    if (is.null(s)) return(NULL)
    q <- tryCatch(set_queries(s, input$db, set_input()$keys), error = function(e) e)
    if (inherits(q, "error")) return(div(class = "alert alert-warning small py-1 px-2", conditionMessage(q)))
    p(class = "small text-muted", sprintf("Using %d sequences of the set '%s'%s, instead of the pasted text or file.",
                                          length(q$seqs), s$name, if (length(set_input()$keys)) " (chosen on the Sets page)" else ""),
      if (q$n_own == 1) "The sample from this database leaves out its own entry from its closest matches."
      else if (q$n_own) sprintf("The %d samples from this database leave out their own entry from their closest matches.", q$n_own))
  })

  view <- reactive({
    views <- req(result())$views
    v <- views[[input$view %||% ""]]
    if (is.null(v)) v <- Filter(function(x) is.null(x$error), views)[[1]]   # just after a new run
    v
  })
  view_error <- function() if (!is.null(view()$error)) div(class = "alert alert-warning", view()$error)

  # The shown result as a viewer alignment, knowing where its rows come from (for sets).
  run_alignment <- function() {
    res <- result()
    a <- msa_from_result(view(), unlist(res$manifest$reverse_complemented))
    a$origin <- msa_run_origin(res$manifest, basename(res$out_dir), cfg$databases)
    a
  }
  add_run_rows <- function(ids) {
    v <- view()
    from <- paste0("Run ", basename(result()$out_dir), if (length(result()$views) > 1) paste0(" (", v$label, ")"))
    add_to_set(set_from_msa(run_alignment(), ids, added_from = from))
  }
  output$add_result_ui <- renderUI({
    v <- view()
    if (!is.null(v$error) || is.null(v$aln_full)) return(NULL)
    div(class = "d-flex align-items-center gap-3 mb-2",
        actionButton("add_result", sprintf("Add the %d sequences of this result to the set", length(v$aln_full)),
                     class = "btn-sm"),
        set_line())
  })
  observeEvent(input$add_result, add_run_rows(names(view()$aln_full)))
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
    content = function(file) write_csv2(use_meta(), file)
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
    use_set_in_run("")
    session$sendCustomMessage("clear_file_input", "fasta_file")
    session$sendCustomMessage("set_page", "outbreak")
    showNotification(paste0(length(seqs), " samples copied to Outbreak detection. Press Run comparison.",
                            if (length(no_seq)) paste0(" No sequence in the database for: ", paste(no_seq, collapse = ", "))),
                     type = if (length(no_seq)) "warning" else "message", duration = if (length(no_seq)) NULL else 8)
  })

  observeEvent(input$view_add, {
    ids <- chosen()
    if (!length(ids)) {
      showNotification("Choose samples first: click rows in the table, or paste sample IDs.", type = "warning")
      return()
    }
    cand <- aln_try(set_from_database(view_db(), cfg, ids, view_meta()$meta))
    if (!is.null(cand)) add_to_set(cand)
  })
  output$view_set_line <- renderUI(set_line())

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

  # A FASTA file that is not aligned, or a set without an alignment, gives no alignment
  # (a = NULL) but a message (`wait`) until it is aligned in the background.
  file_alignments <- reactiveVal(list())   # unaligned FASTA files aligned here, by file key
  aln_base <- reactive({
    src <- input$aln_source %||% "run"
    if (src == "run") {
      res <- req(result())
      v <- res$views[[req(input$aln_view)]]
      req(v, is.null(v$error))
      a <- msa_from_result(v, unlist(res$manifest$reverse_complemented))
      a$origin <- msa_run_origin(res$manifest, basename(res$out_dir), cfg$databases)
      list(key = paste("run", res$out_dir, v$label), a = a)
    } else if (src == "folder") {
      f <- req(aln_folder())
      i <- req(match(input$aln_folder_result, names(f)))
      a <- req(aln_try(msa_load_run_folder(f[[i]], paste0(basename(attr(f, "dir")), ": ", names(f)[i]),
                                           attr(f, "manifest"), attr(f, "regions")[[i]], attr(f, "dropped")[[i]])))
      a$origin <- msa_run_origin(attr(f, "manifest"), basename(attr(f, "dir")), cfg$databases)
      list(key = paste("folder", f[[i]]), a = a)
    } else if (src == "file") {
      file <- req(input$aln_file)
      key <- paste("file", file$name, unname(tools::md5sum(file$datapath)))
      seqs <- req(aln_try(msa_read_fasta(file$datapath)))
      done <- NULL
      if (length(unique(nchar(seqs))) == 1) {
        a <- req(aln_try(msa_load(seqs, file$name)))
      } else {
        done <- file_alignments()[[key]]
        if (is.null(done)) {
          return(list(key = key, a = NULL, seqs = gsub("[-.]", "", seqs), wait = tagList(
            "The sequences in this file are not aligned (they have different lengths). Press ",
            em("Align these sequences"), " on the left: it runs in the background, and the window stays usable.")))
        }
        a <- msa_load(done$aln, file$name, revcomp = done$turned, note = paste0(
          sprintf("Aligned by FHInumerics on %s (DECIPHER).", done$made),
          if (length(done$turned)) paste0(" Turned (reverse complement) for the alignment, so numbered on the reverse complement: ",
                                          paste(done$turned, collapse = ", "), ".")))
      }
      a$origin <- list(kind = "file", file = file$name, made = done$made, method = done$method, versions = done$versions)
      list(key = key, a = a)
    } else {
      req(aln_set_a())
    }
  })
  # A set as the source: rebuilt only when that set's members or alignment change.
  aln_set_a <- reactiveVal(NULL)
  observe({
    req(identical(input$aln_source, "set"))
    id <- input$aln_set
    s <- if (!is.null(id)) sets()[[id]]
    if (is.null(s)) { aln_set_a(NULL); return() }
    sig <- paste(id, s$name, s$alignment$made, paste(set_names(s$members), collapse = "\t"), sep = "\n")
    if (identical(sig, isolate(aln_set_a())$sig)) return()
    a <- msa_from_set(s)
    wait <- if (is.null(a)) {
      if (nrow(s$members) < 2) "This set has fewer than 2 sequences: add more on the Sets page or with Add to set." else
        tagList("This set is not aligned yet. Press ", em("Align the set"),
                " on the left: it runs in the background, and the window stays usable.")
    }
    aln_set_a(list(sig = sig, key = paste("set", sig), a = a, wait = wait))
  })
  # The viewer is hidden while there is no alignment to show.
  observe({
    b <- tryCatch(aln_base(), error = function(e) NULL)
    if (is.null(b$a)) session$sendCustomMessage("msa_data", NULL)
  })
  output$aln_file_ui <- renderUI({
    b <- tryCatch(aln_base(), error = function(e) NULL)
    if (is.null(b$seqs) || !identical(input$aln_source, "file")) return(NULL)
    actionButton("aln_file_align", "Align these sequences", class = "btn-sm btn-primary mb-2")
  })
  observeEvent(input$aln_file_align, {
    b <- aln_base()
    start_alignment("file", b$key, b$seqs, paste("the file", input$aln_file$name))
  })
  output$aln_set_ui <- renderUI({
    ch <- set_choices()
    if (!length(ch)) {
      return(p(class = "small text-muted", "No sets yet. Make one on the", em("Sets"), "page, or press",
               em("Add to set"), "on another page."))
    }
    open <- isolate(open_id())
    if (!isTRUE(open %in% ch)) ch <- c("Choose a set" = "", ch)
    tagList(selectInput("aln_set", "Set", ch, if (isTRUE(open %in% ch)) open else ""), uiOutput("aln_set_status"))
  })
  output$aln_set_status <- renderUI({
    s <- req(sets()[[req(input$aln_set)]])
    tagList(p(class = "small mb-1", aligning_text(input$aln_set) %||% set_alignment_text(s)),
            actionButton("aln_set_align", "Align the set", class = "btn-sm mb-2"))
  })
  observeEvent(input$aln_set, {
    v <- input$aln_set
    if (!nzchar(v)) return()
    if (startsWith(v, "saved:")) open_folder(sub("^saved:", "", v)) else if (!identical(v, open_id())) open_id(v)
  })
  observeEvent(open_id(), if (!is.null(open_id())) updateSelectInput(session, "aln_set", selected = open_id()))
  observeEvent(input$aln_set_align, align_set(input$aln_set))

  # The alignment as shown: all rows, or the rows chosen for it.
  aln_shown <- reactive({
    b <- aln_base()
    req(b$a)
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
    b <- tryCatch(aln_base(), error = function(e) NULL)
    if (is.null(b)) {
      return(p("Run a comparison on the", em("Outbreak detection"), "page, or choose an earlier run folder,",
               "a FASTA file or a set on the left."))
    }
    if (is.null(b$a)) return(div(class = "msa-mode", b$wait))
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
    if (is.null(b$a)) return(NULL)
    tagList(
      if (!is.null(b$a$snp) && length(b$a$queries)) tagList(
        h6("Show a query and its closest"),
        selectInput("aln_query", "Query", b$a$queries),
        numericInput("aln_n", "Number of closest sequences", 10, min = 1, max = 500),
        actionButton("aln_closest", "Show the query and its closest", class = "btn-sm mb-2")),
      div(class = "d-flex gap-2 flex-wrap mb-1",
          actionButton("aln_all", "Show all sequences", class = "btn-sm"),
          actionButton("aln_add", "Add the shown rows to the set", class = "btn-sm")),
      div(class = "mb-3", uiOutput("aln_set_line"))
    )
  })
  output$aln_set_line <- renderUI(set_line())
  observeEvent(input$aln_add, {
    a <- aln_shown()$a
    add_to_set(set_from_msa(a, rownames(a$m), added_from = paste("Alignment page:", a$label)))
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
        actionButton("snp_add", "Add to set", class = "btn-sm"),
        actionButton("snp_unpick", "Clear", class = "btn-sm"))
  })
  observeEvent(input$snp_add, add_run_rows(snp_picked()))
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
        actionButton("tree_show", "Show in Alignment", class = "btn-sm btn-primary"),
        actionButton("tree_add", "Add to set", class = "btn-sm"))
  })
  observeEvent(input$tree_add, add_run_rows(tree_brushed()))
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

  # ==== Sets ===========================================================================
  # Sets in memory (by id), the open one (where "Add to set" adds), and the sets saved in
  # the default place (<results folder>/sets), listed until they are opened.
  sets <- reactiveVal(list())
  open_id <- reactiveVal(NULL)
  n_made <- 0L
  default_sets_dir <- function() unique(c(file.path(cfg$output_dir, "sets"),
                                          file.path(clean_path(isolate(input$out_dir) %||% cfg$output_dir), "sets")))
  list_saved <- function() unique(do.call(rbind, lapply(default_sets_dir(), set_list_saved)))
  saved_sets <- reactiveVal(list_saved())
  open_set <- reactive(if (!is.null(open_id())) sets()[[open_id()]])
  put_set <- function(id, s) { all <- sets(); all[[id]] <- s; sets(all) }
  add_set <- function(s) {
    n_made <<- n_made + 1L
    id <- paste0("set", n_made)
    put_set(id, s)
    open_id(id)
    id
  }
  make_new_set <- function() {
    used <- c(vapply(sets(), function(s) s$name, ""), saved_sets()$name)
    k <- 1
    while (paste("Set", k) %in% used) k <- k + 1
    add_set(set_new(paste("Set", k), save_in = default_sets_dir()[length(default_sets_dir())]))
  }
  open_folder <- function(path) {
    path <- normalizePath(path)
    have <- Filter(function(s) identical(s$folder, path), sets())
    if (length(have)) return(open_id(names(have)[1]))
    s <- aln_try(set_load(path))
    if (!is.null(s)) add_set(s)
  }

  # Add candidates (R/sets.R) to the open set, or to a new one.
  add_to_set <- function(cand) {
    if (is.null(cand)) return()
    if (!nrow(cand$members) && !length(cand$missing)) {
      showNotification("Nothing to add.", type = "warning")
      return()
    }
    id <- open_id() %||% make_new_set()
    r <- set_add(sets()[[id]], cand)
    put_set(id, r$set)
    showNotification(set_add_message(r, cand$missing), type = if (length(r$added)) "message" else "warning", duration = 10)
  }
  # The line under each "Add to set" button.
  set_line <- function() {
    s <- open_set()
    p(class = "small text-muted mb-0", if (is.null(s)) "Adds to a new set (see the Sets page)." else
      sprintf("Adds to the set '%s' (%d sequences), open on the Sets page.", s$name, nrow(s$members)))
  }
  set_choices <- reactive({
    mem <- sets()
    saved <- saved_sets()
    saved <- saved[!saved$folder %in% vapply(mem, function(s) s$folder %||% "", ""), , drop = FALSE]
    stats::setNames(c(names(mem), sprintf("saved:%s", saved$folder)),
                    c(vapply(mem, function(s) s$name, ""), sprintf("%s (saved)", saved$name)))
  })

  output$set_list_ui <- renderUI({
    mem <- sets()
    saved <- saved_sets()
    saved <- saved[!saved$folder %in% vapply(mem, function(s) s$folder %||% "", ""), , drop = FALSE]
    if (!length(mem) && !nrow(saved)) {
      return(p(class = "small text-muted", "No sets yet. Press", em("New set"), "or", em("Add to set"),
               "on another page."))
    }
    labels <- c(vapply(mem, function(s) sprintf("%s (%d%s)", s$name, nrow(s$members),
                                                 if (s$dirty) ", not saved" else ", saved"), ""),
                sprintf("%s (%s, saved)", saved$name, saved$n))
    radioButtons("set_pick", "Sets", stats::setNames(c(names(mem), sprintf("saved:%s", saved$folder)), labels),
                 selected = open_id() %||% character(0))
  })
  observeEvent(input$set_pick, {
    v <- input$set_pick
    if (startsWith(v, "saved:")) open_folder(sub("^saved:", "", v)) else if (!identical(v, open_id())) open_id(v)
  })
  observeEvent(input$set_new, make_new_set())
  shinyFiles::shinyDirChoose(input, "set_open", roots = aln_roots, session = session)
  observeEvent(input$set_open, {
    path <- shinyFiles::parseDirPath(aln_roots(), input$set_open)
    if (length(path)) open_folder(path)
  })

  # Name and folder of the open set. A new name is used at the next Save ("Save as"): the
  # set is then saved in a new folder, and the old one stays.
  observeEvent(open_id(), {
    s <- open_set()
    updateTextInput(session, "set_name", value = s$name %||% "")
    updateTextInput(session, "set_dir", value = s$save_in %||% default_sets_dir()[1])
  }, ignoreNULL = FALSE)
  set_name_typed <- debounce(reactive(input$set_name), 600)
  observeEvent(set_name_typed(), {
    s <- open_set()
    name <- trimws(set_name_typed())
    if (is.null(s) || !nzchar(name) || identical(name, s$name)) return()
    if (name %in% vapply(sets()[names(sets()) != open_id()], function(x) x$name, "")) {
      showNotification(sprintf("Another open set is already called '%s'. Choose another name.", name), type = "warning")
      return()
    }
    s$history <- c(s$history, sprintf("%s renamed from '%s'", format(Sys.time(), "%Y-%m-%d %H:%M"), s$name))
    s$name <- name
    s$dirty <- TRUE
    put_set(open_id(), s)
  })
  observeEvent(input$set_dir, {
    s <- open_set()
    if (!is.null(s) && nzchar(input$set_dir) && !identical(clean_path(input$set_dir), s$save_in)) {
      s$save_in <- clean_path(input$set_dir)
      put_set(open_id(), s)
    }
  })
  shinyFiles::shinyDirChoose(input, "set_browse", roots = aln_roots, session = session)
  observeEvent(input$set_browse, {
    path <- shinyFiles::parseDirPath(aln_roots(), input$set_browse)
    if (length(path)) updateTextInput(session, "set_dir", value = path)
  })

  # Save: into <Save in>/<name>/. A different saved set there is replaced only after asking.
  pending_save <- NULL
  do_save <- function(id, folder) {
    s <- aln_try(set_save(sets()[[id]], folder))
    if (is.null(s)) return()
    all <- sets()
    for (k in setdiff(names(all), id)) {   # another open set that was saved there is now unsaved
      if (identical(all[[k]]$folder, s$folder)) { all[[k]]$folder <- NULL; all[[k]]$dirty <- TRUE }
    }
    all[[id]] <- s
    sets(all)
    saved_sets(list_saved())
    showNotification(paste("Set saved in", s$folder), type = "message", duration = 8)
  }
  observeEvent(input$set_save, {
    s <- open_set()
    if (is.null(s)) return(showNotification("No set is open. Press New set, or Add to set on another page.", type = "warning"))
    if (!nrow(s$members)) return(showNotification("The set is empty: add sequences first.", type = "warning"))
    folder <- file.path(clean_path(input$set_dir), set_folder_name(s$name))
    state <- set_folder_state(folder)
    if (state == "other") {
      return(showNotification(paste("The folder", folder, "already exists and holds other files.",
                                    "Choose another name or another folder."), type = "error", duration = NULL))
    }
    if (state == "set" && !identical(normalizePath(folder), s$folder)) {
      pending_save <<- list(id = open_id(), folder = folder)
      return(showModal(modalDialog(title = "Replace the saved set?",
        sprintf("A saved set is already in %s. Replace it with the set '%s'?", folder, s$name),
        footer = tagList(modalButton("Cancel"), actionButton("set_replace", "Replace", class = "btn-danger")))))
    }
    do_save(open_id(), folder)
  })
  observeEvent(input$set_replace, {
    removeModal()
    if (!is.null(pending_save)) do_save(pending_save$id, pending_save$folder)
    pending_save <<- NULL
  })

  # Delete: a saved set's files (only those FHInumerics writes), or an unsaved set.
  observeEvent(input$set_delete, {
    s <- open_set()
    if (is.null(s)) return()
    saved <- !is.null(s$folder) && file.exists(file.path(s$folder, "set.yaml"))
    showModal(modalDialog(title = "Delete the set?",
      if (saved) sprintf("Delete the saved set '%s'? Its files in %s are removed.", s$name, s$folder)
      else sprintf("Remove the set '%s'? It has not been saved, so it will be gone.", s$name),
      footer = tagList(modalButton("Cancel"), actionButton("set_delete_ok", "Delete", class = "btn-danger"))))
  })
  observeEvent(input$set_delete_ok, {
    removeModal()
    s <- open_set()
    left <- character(0)
    if (!is.null(s$folder) && file.exists(file.path(s$folder, "set.yaml"))) {
      left <- aln_try(set_delete(s$folder))
      if (is.null(left)) return()
    }
    all <- sets()
    all[[open_id()]] <- NULL
    sets(all)
    open_id(if (length(all)) names(all)[1])
    saved_sets(list_saved())
    showNotification(paste0("The set '", s$name, "' was deleted.",
                            if (length(left)) paste0(" Other files were left in its folder: ", paste(left, collapse = ", "))),
                     type = "message", duration = if (length(left)) NULL else 6)
  })

  observeEvent(input$set_file, {
    add_to_set(aln_try(set_from_fasta(input$set_file$datapath, input$set_file$name)))
    session$sendCustomMessage("clear_file_input", "set_file")
  })

  # The open set: heading, table of members (click rows to remove them), downloads.
  output$set_head <- renderUI({
    s <- open_set()
    if (is.null(s)) {
      return(p("No set is open. Press", em("New set"), "on the left, or", em("Add to set"),
               "on the Database viewer, Outbreak detection or Alignment page."))
    }
    m <- s$members
    from <- c(sprintf("%d from %s", as.vector(table(m$database[m$source == "database"])),
                      names(table(m$database[m$source == "database"]))),
              if (any(m$source == "query")) sprintf("%d %s", sum(m$source == "query"), if (sum(m$source == "query") == 1) "query" else "queries"),
              if (any(m$source == "file")) sprintf("%d from files", sum(m$source == "file")))
    tagList(
      div(class = "msa-head", h5(s$name), div(class = "msa-sub", paste(c(
        sprintf("%d sequences", nrow(m)), from,
        if (s$dirty) "not saved" else paste("saved in", s$folder)), collapse = " · "))),
      p(class = "small text-muted mt-1",
        "Click rows to choose them, then press Remove. Add sequences with Add to set on the other pages,",
        "or with a FASTA file on the left. Show the set in the Alignment page with Show in Alignment."),
      div(class = "msa-mode", aligning_text(open_id()) %||% set_alignment_text(s))
    )
  })
  output$set_rows_ui <- renderUI({
    n <- length(input$set_table_rows_selected)
    if (!n || is.null(open_set())) return(NULL)
    div(class = "d-flex align-items-center gap-2 mb-2", strong(n, "chosen"),
        actionButton("set_remove", "Remove from the set", class = "btn-sm btn-danger"),
        actionButton("set_use_rows", "Use the chosen in Outbreak detection", class = "btn-sm"))
  })
  use_open_set <- function() {
    s <- open_set()
    if (is.null(s)) return(showNotification("No set is open.", type = "warning"))
    keys <- s$members$key[input$set_table_rows_selected]
    q <- aln_try(set_queries(s, NULL, keys))
    if (is.null(q)) return()
    use_set_in_run(open_id(), if (length(keys)) keys)
    session$sendCustomMessage("set_page", "outbreak")
    showNotification(sprintf("%d sequences of the set '%s' will be compared. Choose the database and press Run comparison.",
                             length(q$seqs), s$name), type = "message", duration = 8)
  }
  observeEvent(input$set_use, use_open_set())
  observeEvent(input$set_use_rows, use_open_set())
  output$set_table <- DT::renderDT({
    s <- req(open_set())
    tab <- set_table(s)
    DT::datatable(tab, class = "display nowrap", fillContainer = FALSE, rownames = FALSE, filter = "top", selection = "multiple",
                  options = list(scrollX = TRUE, pageLength = 25,
                                 columnDefs = list(list(visible = FALSE, targets = match("database_md5", names(tab)) - 1))))
  }, server = FALSE)
  observeEvent(input$set_remove, {
    s <- open_set()
    keys <- s$members$key[input$set_table_rows_selected]
    put_set(open_id(), set_remove(s, keys))
    showNotification(sprintf("%d removed from the set '%s'.", length(keys), s$name), type = "message")
  })
  set_file_name <- function(ext) paste0("set_", gsub("[^A-Za-z0-9._-]+", "_", open_set()$name %||% "set"), ext)
  output$set_dl_fasta <- downloadHandler(
    filename = function() set_file_name(".fasta"),
    content = function(file) write_fasta(set_fasta(req(open_set())), file)
  )
  output$set_dl_csv <- downloadHandler(
    filename = function() set_file_name(".csv"),
    content = function(file) write_csv2(set_table(req(open_set())), file)
  )
  output$set_align_ui <- renderUI({
    s <- open_set()
    if (is.null(s)) return(NULL)
    tagList(p(class = "small mb-1", aligning_text(open_id()) %||% set_alignment_text(s)),
            div(class = "d-flex gap-2 flex-wrap",
                actionButton("set_align", "Align the set", class = "btn-sm"),
                actionButton("set_show", "Show in Alignment", class = "btn-sm btn-primary")))
  })
  observeEvent(input$set_align, align_set(open_id()))
  observeEvent(input$set_show, {
    updateRadioButtons(session, "aln_source", selected = "set")
    updateSelectInput(session, "aln_set", selected = open_id())
    session$sendCustomMessage("set_page", "alignment")
  })
  # Ask before the tab is closed while a set has unsaved changes.
  observe(session$sendCustomMessage("sets_unsaved", any(vapply(sets(), function(s) isTRUE(s$dirty), TRUE))))

  # ---- Background alignment (one at a time): sets and unaligned FASTA files ----
  job <- reactiveVal(NULL)
  pending_align <- NULL
  start_alignment <- function(kind, target, seqs, label) {
    if (!is.null(job())) {
      return(showNotification(sprintf("An alignment is already running (%s). Wait for it to finish, or press Cancel.",
                                      job()$label), type = "warning"))
    }
    if (length(seqs) < 2) return(showNotification("At least 2 sequences are needed for an alignment.", type = "warning"))
    go <- function() {
      j <- aln_try(align_start(seqs, defaults$threads))
      if (is.null(j)) return()
      job(list(env = j, kind = kind, target = target, label = label))
      session$sendCustomMessage("align_job", list(label = sprintf("%s (%d sequences)", label, length(seqs)),
                                                  value = 0, text = "Starting"))
    }
    if (length(seqs) > SET_ALIGN_WARN) {
      pending_align <<- go
      return(showModal(modalDialog(title = "Align many sequences?",
        sprintf(paste("%s has %d sequences. Aligning more than %d can take several minutes (about 1.5 minutes",
                      "for 1,800 HAV sequences on a fast computer). It runs in the background: the window stays",
                      "usable, and Cancel stops it. Align now?"), label, length(seqs), SET_ALIGN_WARN),
        footer = tagList(modalButton("Not now"), actionButton("align_go", "Align", class = "btn-primary")))))
    }
    go()
  }
  aligning_text <- function(id) {
    j <- job()
    if (!is.null(j) && j$kind == "set" && identical(j$target, id)) "Being aligned now: see the progress box under the menu."
  }
  observeEvent(input$align_go, { removeModal(); if (!is.null(pending_align)) pending_align(); pending_align <<- NULL })
  align_set <- function(id) {
    s <- if (!is.null(id)) sets()[[id]]
    if (is.null(s)) return(showNotification("No set is open.", type = "warning"))
    start_alignment("set", id, stats::setNames(unname(s$seqs[s$members$key]), s$members$key), sprintf("the set '%s'", s$name))
  }
  observe({
    j <- job()
    req(j)
    invalidateLater(500)
    st <- align_poll(j$env)
    if (!st$done) {
      secs <- round(as.numeric(difftime(Sys.time(), j$env$started, units = "secs")))
      return(session$sendCustomMessage("align_job", list(label = sprintf("%s (%d sequences)", j$label, j$env$n),
                                                         value = st$value, text = sprintf("%s · %d s", st$text, secs))))
    }
    job(NULL)
    session$sendCustomMessage("align_job", list())
    if (!is.null(st$error)) return(showNotification(st$error, type = "error", duration = NULL))
    res <- st$result
    secs <- round(as.numeric(difftime(Sys.time(), j$env$started, units = "secs")))
    if (j$kind == "set") {
      s <- isolate(sets())[[j$target]]
      if (is.null(s)) return()   # deleted meanwhile
      s$alignment <- res[c("aln", "turned", "made", "from", "method", "versions")]
      s$history <- c(s$history, sprintf("%s aligned (%d sequences)", res$made, length(res$aln)))
      s$dirty <- TRUE
      put_set(j$target, s)
      turned <- set_alignment_status(s)$turned
    } else {
      done <- file_alignments()
      done[[j$target]] <- res
      file_alignments(done)
      turned <- res$turned
    }
    showNotification(paste0(sprintf("Done: %s is aligned (%d s).", j$label, secs),
                            if (length(turned)) paste0(" Turned (reverse complement): ", paste(turned, collapse = ", "), ".")),
                     type = "message", duration = 10)
  })
  observeEvent(input$align_cancel, {
    j <- job()
    if (is.null(j)) return()
    align_cancel(j$env)
    job(NULL)
    session$sendCustomMessage("align_job", list())
    showNotification("The alignment was cancelled.", type = "message")
  })

  # When started from the launcher, stop the app once the browser tab is closed. A running
  # background alignment is stopped in any case.
  session$onSessionEnded(function() {
    j <- isolate(job())
    if (!is.null(j)) align_cancel(j$env)
    if (isTRUE(getOption("fhinumerics.stop_on_close"))) stopApp()
  })
}

shinyApp(ui, server)
