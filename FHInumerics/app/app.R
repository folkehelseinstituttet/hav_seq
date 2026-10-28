# FHInumerics: compare new HAV sequences to a reference database.
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

ui <- page_sidebar(
  title = "FHInumerics - HAV: find the closest sequences in a database",
  sidebar = sidebar(
    width = 380,
    textAreaInput("fasta_text", "Paste sequence(s) in FASTA format", rows = 8, width = "100%",
                  placeholder = ">sample1\nTTCAAGAGGGGTCTCCGG...\n>sample2\n..."),
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
  ),
  navset_card_tab(
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
    nav_panel("SNP distances", uiOutput("snp_note"), DT::DTOutput("snp_table")),
    nav_panel("Closest matches", uiOutput("ties_note"), DT::DTOutput("hits_table")),
    nav_panel("Run log", verbatimTextOutput("log"))
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

  shinyFiles::shinyDirChoose(input, "browse_out", roots = roots, session = session)
  observeEvent(input$browse_out, {
    path <- shinyFiles::parseDirPath(roots, input$browse_out)
    if (length(path)) updateTextInput(session, "out_dir", value = path)
  })

  observeEvent(input$run, {
    lines <- character(0)
    log <- function(...) lines <<- c(lines, paste0(format(Sys.time(), "%H:%M:%S  "), ...))
    res <- tryCatch({
      text <- if (!is.null(input$fasta_file)) {
        paste(readLines(input$fasta_file$datapath, warn = FALSE), collapse = "\n")
      } else input$fasta_text
      queries <- parse_fasta_text(text)
      log("Queries: ", paste(names(queries), collapse = ", "))
      cfg$output_dir <- clean_path(input$out_dir)
      settings <- list(n_neighbours = input$n_neighbours, min_coverage = input$min_coverage,
                       aligner = input$aligner, bootstrap = input$bootstrap,
                       search = defaults$search, threads = defaults$threads)
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
    plotOutput("tree", height = tree_height())
  })
  output$tree <- renderPlot(tree_plot())
  output$dl_tree <- downloadHandler(
    filename = function() paste0("tree_", input$colour_by, ".pdf"),
    content = function(file) ggplot2::ggsave(file, tree_plot(), width = 11, height = tree_height() / 72, limitsize = FALSE)
  )

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
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 100, dom = "t")),
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
    DT::datatable(hits, rownames = FALSE, filter = "top",
                  options = list(scrollX = TRUE, pageLength = 25))
  })

  output$log <- renderText(paste(log_lines(), collapse = "\n"))

  # When started from the launcher, stop the app once the browser tab is closed.
  session$onSessionEnded(function() {
    if (isTRUE(getOption("fhinumerics.stop_on_close"))) stopApp()
  })
}

shinyApp(ui, server)
