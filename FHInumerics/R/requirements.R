# What launcher/setup.R installs, and whether it (still) has to run. launcher/launch.R
# checks this on every start and runs setup.R when something is missing, so the first start
# installs everything by itself.

CRAN_PACKAGES <- c("shiny", "bslib", "DT", "yaml", "jsonlite", "ggplot2", "ape", "phangorn",
                   "shinyFiles", "httpuv", "processx", "BiocManager", "maps")
BIOC_PACKAGES <- c("Biostrings", "DECIPHER", "ggtree")

missing_packages <- function(pkgs = c(CRAN_PACKAGES, BIOC_PACKAGES)) {
  pkgs[!nzchar(vapply(pkgs, function(p) system.file(package = p), ""))]
}

# TRUE when R packages or config.yaml are missing, or when config.yaml asks for
# MAFFT (aligner: mafft) and its Python package is not installed yet.
setup_needed <- function(root) {
  config <- file.path(root, "config.yaml")
  if (length(missing_packages()) || !file.exists(config)) return(TRUE)
  aligner <- yaml::read_yaml(config)$settings$aligner %||% "decipher"
  aligner == "mafft" && !dir.exists(file.path(root, "tools", "python-packages", "pymafft"))
}
