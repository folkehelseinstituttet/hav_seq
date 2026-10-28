# Changelog

All notable changes to FHInumerics. The version is in `VERSION` and is written into every
run's `run_manifest.yaml`. Newest first.

## 0.1.0 - 2026-10-28

First version in the `hav_seq` repository (previously developed as "hav-compare").

- One double-click starts everything: `FHInumerics.bat` (Windows) or `FHInumerics.sh`.
  The first start installs the R packages into the personal R library (no admin rights)
  and then opens the window; later starts open it directly. Each R version (4.4 and newer)
  gets its own matching packages.
- Plain-language messages, with the window kept open, when R is missing or too old, the
  download fails, `config.yaml` is broken, or the folder is run from inside a `.zip`.
- Alignment with DECIPHER `AlignSeqs` by default, so only R is needed. MAFFT (pymafft,
  Python 3.9-3.13) remains available with `aligner: mafft`.
- Automated tests: `Rscript tests/run_tests.R`.
- New folder layout (`app/`, `launcher/`, `cli/`, `data/`); manual in `MANUAL.md`.
- Tested on a Windows work PC and on Linux.
