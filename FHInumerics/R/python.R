# Running the Python step (MAFFT via pymafft; see python/fhinumerics_tools.py).
# Python is the one installed on this computer: `python` in config.yaml, or found
# automatically. The package is installed by setup.R into tools/python-packages,
# which is put first on PYTHONPATH so the pinned version is always the one used.

# pymafft 0.2.0 (python/requirements.txt) has ready-made packages for Python 3.9-3.13
# only. Newer Pythons would need a Rust compiler to install it, so they are not used.
PYTHON_MINOR <- 9:13
PYTHON_RANGE <- "3.9-3.13"

python_ok <- function(version) {
  minor <- suppressWarnings(as.integer(sub("^3\\.(\\d+)\\..*$", "\\1", version)))
  !is.na(minor) & grepl("^3\\.", version) & minor %in% PYTHON_MINOR
}

# Candidate Python programs on this computer, most specific first.
python_candidates <- function() {
  found <- unname(Sys.which(c("python3", "python")))
  if (.Platform$OS.type == "windows") {
    # The Windows "py" launcher knows the default Python.
    if (nzchar(Sys.which("py"))) {
      py <- tryCatch(system2("py", c("-3", "-c", shQuote("import sys; print(sys.executable)")),
                             stdout = TRUE, stderr = FALSE), error = function(e) character(0))
      found <- c(found, py)
    }
    # python.org installers (per user, all users) and the Python install manager.
    local <- Sys.getenv("LOCALAPPDATA")
    found <- c(found, Sys.glob(file.path(local, "Programs", "Python", "Python3*", "python.exe")),
               Sys.glob(file.path(Sys.getenv("ProgramFiles"), "Python3*", "python.exe")),
               Sys.glob(file.path(local, "Python", "pythoncore-3*", "python.exe")))
    found <- found[!grepl("WindowsApps", found)]   # Microsoft Store placeholders
  }
  unique(found[nzchar(found) & file.exists(found)])
}

python_version <- function(python) {
  v <- tryCatch(system2(python, c("-c", shQuote("import platform; print(platform.python_version())")),
                        stdout = TRUE, stderr = FALSE), error = function(e) character(0))
  if (length(v)) v[1] else NA_character_
}

# Path of the Python to use: config.yaml `python`, else the first suitable one found.
find_python <- function(configured = NULL) {
  if (!is.null(configured)) {
    if (!file.exists(configured)) stop("Python not found at '", configured, "' (python in config.yaml).", call. = FALSE)
    return(configured)   # not normalizePath(): that would resolve a virtualenv's python link
  }
  cands <- python_candidates()
  if (!length(cands)) return(NULL)
  ok <- python_ok(vapply(cands, python_version, ""))
  if (any(ok)) cands[ok][1] else cands[1]
}

python_packages_dir <- function() file.path(getOption("fhinumerics.root", "."), "tools", "python-packages")

# Run python/fhinumerics_tools.py with the given arguments; stop with its output if it fails.
run_python <- function(cfg, args, log = message) {
  if (is.null(cfg$python)) {
    stop("Python was not found. Install Python ", PYTHON_RANGE, " and run setup, or set python in ",
         "config.yaml. Without Python, choose the DECIPHER alignment.", call. = FALSE)
  }
  script <- file.path(getOption("fhinumerics.root", "."), "python", "fhinumerics_tools.py")
  env <- c("current", PYTHONPATH = python_packages_dir(), PYTHONNOUSERSITE = "1", PYTHONIOENCODING = "utf-8")
  res <- processx::run(cfg$python, c(script, args), env = env, error_on_status = FALSE, stderr_to_stdout = TRUE)
  if (res$status != 0) {
    log(res$stdout)
    stop("Python step '", args[1], "' failed (exit status ", res$status, "). See the run log.", call. = FALSE)
  }
  res$stdout
}

python_versions <- function(cfg) {
  if (is.null(cfg$python)) return(list(python = "not found"))
  out <- tryCatch(run_python(cfg, "versions", log = function(...) NULL), error = function(e) NULL)
  if (is.null(out)) return(list(python = cfg$python, error = "could not run python/fhinumerics_tools.py versions"))
  jsonlite::fromJSON(utils::tail(strsplit(trimws(out), "\n")[[1]], 1))
}
