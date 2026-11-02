# Background alignment for the window (align_start() in R/sets.R):
#   Rscript app/align_worker.R <FHInumerics folder> <job folder>
# Aligns the sequences in <job folder>/input.rds with set_align() and saves the result
# (or a plain error message) as <job folder>/result.rds. Progress goes to stdout, where
# align_poll() reads it. The window stays usable meanwhile; Cancel kills this process.

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
dir <- args[2]
input <- readRDS(file.path(dir, "input.rds"))
.libPaths(input$libs)
options(fhinumerics.root = root, width = 80)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)

res <- tryCatch(set_align(input$seqs, input$threads, verbose = TRUE),
                error = function(e) list(error = conditionMessage(e)))
saveRDS(res, file.path(dir, "result.tmp"))
file.rename(file.path(dir, "result.tmp"), file.path(dir, "result.rds"))
