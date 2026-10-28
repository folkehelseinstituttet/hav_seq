# Setup: install the R packages, install the Python package for MAFFT (pymafft)
# into tools/python-packages, and create config.yaml.
# This runs by itself the first time FHInumerics is started (launcher/launch.R, from
# FHInumerics.bat or ./FHInumerics.sh), so it does not have to be run by hand.
# To run it anyway: Rscript launcher/setup.R
# Needs R 4.4 or newer (https://cran.r-project.org) and the internet. Python
# (3.9-3.13, https://www.python.org) is only needed with aligner: mafft in config.yaml.
# No other programs are installed or downloaded.

file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(file_arg)) dirname(dirname(normalizePath(file_arg))) else getwd()   # file is in launcher/
options(timeout = 3600, repos = c(CRAN = "https://cloud.r-project.org"))
if (getRversion() < "4.4") {
  stop("FHInumerics needs R 4.4 or newer, and this computer has R ", getRversion(), ". ",
       "Install a newer R from https://cran.r-project.org (or ask IT).", call. = FALSE)
}
source(file.path(root, "R", "requirements.R"))

# ---- R packages ------------------------------------------------------------------
# Into the personal R library (created if needed), so no admin rights are needed
# when R itself is installed in Program Files.
lib <- path.expand(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep)[[1]][1])
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
# A lock left behind when an earlier install was stopped (window closed) blocks every
# later install. R holds it for seconds per package, so an old one is left over.
locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
unlink(locks[difftime(Sys.time(), file.mtime(locks), units = "mins") > 5], recursive = TRUE)

missing <- missing_packages(CRAN_PACKAGES)
if (length(missing)) install.packages(missing)
missing <- missing_packages(BIOC_PACKAGES)
if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE)
still_missing <- missing_packages()
if (length(still_missing)) {
  stop("Could not install the R packages ", paste(still_missing, collapse = ", "), ". ",
       "Check that this computer is connected to the internet (a VPN or firewall can block ",
       "the download) and start FHInumerics again. If FHInumerics is already installing in ",
       "another window, wait until that has finished.", call. = FALSE)
}

# ---- config.yaml -------------------------------------------------------------------
config_path <- file.path(root, "config.yaml")
if (!file.exists(config_path)) {
  file.copy(file.path(root, "config.example.yaml"), config_path)
  message("Created config.yaml (edit it to add your own databases).")
}
options(fhinumerics.root = root)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
cfg <- load_config(config_path)
aligner <- cfg$settings$aligner %||% DEFAULT_SETTINGS$aligner

# ---- Python package (only for aligner: mafft) -------------------------------------
# Installed with this computer's Python into tools/python-packages (not into Python
# itself), at the version in python/requirements.txt. Only ready-made packages are
# used (--only-binary), so nothing is compiled.
install_pymafft <- function(cfg) {
  if (is.null(cfg$python)) {
    stop("Python was not found. Install Python ", PYTHON_RANGE, " from https://www.python.org ",
         "(or ask IT), or set python in config.yaml, then start FHInumerics again.", call. = FALSE)
  }
  version <- python_version(cfg$python)
  message("Python: ", cfg$python, " (", version, ")")
  if (!python_ok(version)) {
    stop("MAFFT (pymafft) needs Python ", PYTHON_RANGE, "; this is ", version, ". ",
         "Install one of those versions (they can sit next to the one you have) and, if it is ",
         "not found automatically, set python in config.yaml.", call. = FALSE)
  }
  if (system2(cfg$python, c("-m", "pip", "--version"), stdout = FALSE, stderr = FALSE) != 0) {
    stop("This Python has no pip, so packages cannot be installed. Use a Python that has pip ",
         "(the python.org installer includes it) and set python in config.yaml.", call. = FALSE)
  }
  message("Installing pymafft into tools/python-packages...")
  status <- system2(cfg$python, c("-m", "pip", "install", "--disable-pip-version-check", "--upgrade",
                                  "--ignore-installed", "--only-binary", ":all:",
                                  "--target", shQuote(python_packages_dir()),
                                  "-r", shQuote(file.path(root, "python", "requirements.txt"))))
  if (status != 0) stop("pip install failed (see the messages above).", call. = FALSE)
  v <- python_versions(cfg)
  if (is.null(v$pymafft)) stop("The Python package could not be loaded.", call. = FALSE)
  message("MAFFT: pymafft ", v$pymafft)
}

if (aligner == "mafft") {
  tryCatch(install_pymafft(cfg), error = function(e) {
    stop(conditionMessage(e), " To work without Python, set aligner: decipher in config.yaml.",
         call. = FALSE)
  })
} else {
  message("Alignment with DECIPHER (R), so Python is not needed.")
}
message("Setup finished.")
