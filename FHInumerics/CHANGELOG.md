# Changelog

All notable changes to FHInumerics. The version is in `VERSION` and is written into every
run's `run_manifest.yaml`. Newest first.

## 0.2.0 - 2026-10-28

- Menu at the top left with two pages: **Outbreak detection** (the comparison, as before)
  and the new **Database viewer**.
- Database viewer: the database as a table with a filter per column and a date range;
  choose samples by clicking rows, pasting IDs or choosing a text file of IDs; charts of
  samples over time, counts per value, a cross-table and a map of samples per country
  (each as PDF); download the table (CSV) or the sequences (FASTA); send up to 50 chosen
  samples to Outbreak detection.
- Samples sent from the Database viewer leave out their own database entry from their
  closest matches (`exclude_self` in `run_manifest.yaml`; `--exclude-self` on the
  command line).
- Tree: point at the dot of a sample to see its metadata (for a query: its closest
  sequences); an orange ring marks the dot.
- Map: point at a circle to see the country, its number of samples and the number per
  genotype (or another column); an orange ring marks the circle.
- New settings per database (both optional): `date_column` (the date column is otherwise
  found automatically and can be changed in the window), `country_column`.
- Tables (also SNP distances and Closest matches) use the full height of the window; the
  page scrolls.
- New R package: maps (downloaded automatically at the next start).

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
