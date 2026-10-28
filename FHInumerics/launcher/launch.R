# Start FHInumerics and open it in the browser. The app only listens on
# 127.0.0.1 (this computer); it stops when the browser tab is closed.
# The first time (or when something is missing) it runs setup.R first, so
# double-clicking FHInumerics.bat is all a user has to do.
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- dirname(dirname(normalizePath(file_arg)))   # this file is in launcher/
line <- strrep("=", 72)

# Any problem ends in one plain message; FHInumerics.bat then keeps the window open.
fail <- function(msg) {
  message("\n", line, "\nFHInumerics could not start.\n\n", paste(strwrap(msg, 70), collapse = "\n"),
          "\n\nIf you cannot fix this, take a screenshot of this window and send it\n",
          "to the person who gave you FHInumerics.\n", line)
  quit(save = "no", status = 1)
}

if (getRversion() < "4.4") {
  fail(paste0("FHInumerics needs R 4.4 or newer, and this computer has R ", getRversion(),
              ". Install a newer R from https://cran.r-project.org (or ask IT)."))
}

tryCatch({
  source(file.path(root, "R", "requirements.R"))
  if (setup_needed(root)) {
    message(line, "\nFHInumerics is installing what it needs. This is only done the first time,\n",
            "takes 5-15 minutes and needs the internet.\n",
            "Please wait, and do not close this window.\n", line)
    source(file.path(root, "launcher", "setup.R"), local = new.env())
    message(line, "\nInstallation finished. Starting FHInumerics...\n", line)
  }
  options(fhinumerics.root = root, fhinumerics.stop_on_close = TRUE)
  open_browser <- function(url) {
    message("FHInumerics is running at ", url, "  (close the browser tab to stop it)")
    try(utils::browseURL(url), silent = TRUE)
  }
  shiny::runApp(file.path(root, "app"), host = "127.0.0.1", port = httpuv::randomPort(), launch.browser = open_browser)
}, error = function(e) fail(conditionMessage(e)))
