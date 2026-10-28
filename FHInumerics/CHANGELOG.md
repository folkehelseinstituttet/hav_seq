# Changelog

All notable changes to FHInumerics. The version is in `VERSION` and is written into every
run's `run_manifest.yaml`. Newest first.

## Unreleased

- **SNP counting takes ambiguity codes into account** (phangorn `dist.hamming` instead of
  ape `dist.dna`): an ambiguity code counts as a SNP only when it cannot match the other
  base (K = G/T vs A counts; M = A/C vs A does not). Before, ambiguity codes were skipped
  like N. N and gaps are still never counted. SNP numbers can go up for sequences with
  ambiguity codes; others are unchanged. Compared positions now include ambiguity codes.
- SNP distances: point at a number in a query's row or column to see the differing
  positions (base number in the query and the two bases, e.g. `52 - A/K`), and what is
  not counted: gaps (missing ends as a range, gaps inside per position), ambiguity codes
  that fit (e.g. `47 - A/M`) and the number of Ns. An orange corner marks numbers where
  the pair differs by any of these.
- Several queries: besides the result for all queries together, each query gets its own
  result (its own closest sequences, alignment over its whole length, tree and SNP
  matrix). A **Show results for** list above the tabs switches between them; the Closest
  matches tab shows only the chosen query's hits. Saved in `per_query/` in the run folder.
- Queries that do not overlap no longer stop the run: there is no result for all
  queries together, but each query's own result is made.
- SNP distances: queries are shown in red (row names and column headers).
- A query with the same name as a database sequence is renamed `<name>_query`; before,
  the two had the same name in the tree and SNP matrix, and the query got the database
  sample's metadata.
- `run_manifest.yaml`: new `renamed` and `results` (per result: folder, neighbours,
  compared positions); `neighbours` and `compared_columns` moved into `results`.

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
