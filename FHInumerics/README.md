# FHInumerics

Compare new HAV sequences to a reference database, in a point-and-click window.

Paste sequences (FASTA, one or more, like NCBI BLAST) or choose a FASTA file. FHInumerics
finds the closest database sequences, builds a SNP distance matrix and a phylogeny, and
shows the database metadata (date, country, clade, lineage, ...) in the tree.

Everything runs on your own computer. The window opens in your browser but only listens
on `127.0.0.1` (this computer): no server, no open network ports, no data leaves the PC.
Only **R** is needed: no WSL, conda, Python or separate programs.

How the analyses work, every setting and every output file: see **[MANUAL.md](MANUAL.md)**.

## Install and start

1. Install **R** 4.4 or newer (https://cran.r-project.org, or your organisation's software
   portal). RStudio is not needed.
2. Get this `FHInumerics` folder (copy, or clone the `hav_seq` repository), for example to
   your own folder on a network drive or `C:\FHInumerics`. If it came as a `.zip`, unzip
   it first.
3. Put your `config.yaml` (your databases, see below) in the folder. Without one, the
   first start makes one from `config.example.yaml` with the demo database only.
4. Start it:
   - **Windows:** double-click **`FHInumerics.bat`**. To put it on the desktop:
     right-click it > *Send to* > *Desktop (create shortcut)*.
   - **Linux / Mac:** `./FHInumerics.sh`

The **first start** installs the R packages it needs (5-15 minutes, needs the internet),
then opens the window in a browser tab. Later starts open it directly. Closing the browser
tab stops FHInumerics.

If something goes wrong (R not installed, R too old, no internet, a broken `config.yaml`,
started from inside a `.zip`), the black window explains it in plain words and stays open.

What the first start does (`launcher/setup.R`; it runs again by itself if something goes
missing, e.g. after changing R version): it downloads the R packages shiny, bslib, DT,
yaml, ggplot2, maps, ape, phangorn, shinyFiles, processx and, from Bioconductor, Biostrings,
DECIPHER and ggtree into your personal R library (no admin rights needed). Each R version
(4.4, 4.5, 4.6, ...) gets its own copy, matching that version. It also creates
`config.yaml` if there is none.

## Your own databases (`config.yaml`)

`config.yaml` is your local settings file and is never committed to git. Start from
`config.example.yaml` and add an entry under `databases:`:

```yaml
  - name: In-house database
    fasta: X:/path/to/sequences.fasta
    metadata: X:/path/to/metadata.csv        # .tsv or .csv, one row per sequence
    id_column: Key                           # column matching the FASTA IDs
    display_columns: [date, country, genotype]
    date_column: date                        # optional: date column in the Database viewer
    country_column: country                  # optional: country column for the map
```

Use forward slashes in Windows paths. The search index is built automatically in `cache/`
the first time (and again when the FASTA changes), so the database files can stay on a
read-only share. Other settings (results folder, default analysis settings) are explained
in `config.example.yaml` and in [MANUAL.md](MANUAL.md).

## Command line (no window)

    Rscript cli/fhinumerics --query my.fasta --out results/

Options: `--db NAME`, `--neighbours 20`, `--min-coverage 90`, `--no-bootstrap`,
`--aligner decipher|mafft`, `--threads 4`, `--exclude-self` (leave out database entries
with the same name as a query), `--config config.yaml`. Start the window once
first (or run `Rscript launcher/setup.R`) so the packages are installed.

## For developers

- **Tests:** `Rscript tests/run_tests.R` runs the test queries against the demo database
  and checks the known answers (PASS/FAIL). Run it after every change.
- **Changes:** see [CHANGELOG.md](CHANGELOG.md); the version is in `VERSION` and is
  written into every run's `run_manifest.yaml`.

| Path | What |
|---|---|
| `FHInumerics.bat`, `FHInumerics.sh` | start files (Windows, Linux/Mac) |
| `launcher/` | `launch.R` (start, first-start install, error messages), `setup.R` (install), `find-r.bat` (finds R on Windows) |
| `app/app.R` | the window (Shiny) |
| `R/compare.R` | the analysis, shared by window, command line and tests |
| `R/search.R`, `R/db.R`, `R/fasta.R`, `R/plot.R` | database search, config and database loading, FASTA reading, tree plot |
| `R/python.R`, `python/` | optional MAFFT alignment (Python 3.9-3.13 + pymafft) |
| `R/requirements.R` | which packages are needed; whether setup has to run |
| `cli/fhinumerics` | command-line version |
| `data/demo/` | demo database (public sequences, see `data/demo/SOURCE.md`) |
| `tests/` | test queries with known answers and `run_tests.R` |
| `scripts/prepare_demo_db.py` | how the demo database was made (not needed to run) |
| `config.yaml`, `cache/`, `tools/` | local only, created on this computer (not in git) |
