# FHInumerics manual

What FHInumerics does with your sequences, step by step; how to use the window; every
setting; and every file a run writes. For installing and starting, see [README.md](README.md).

Contents: [1. A run in short](#1-a-run-in-short) ·
[2. Using the window](#2-using-the-window) ·
[3. The analysis step by step](#3-the-analysis-step-by-step) ·
[4. Output files](#4-output-files) · [5. Settings](#5-settings-configyaml) ·
[6. Reproducibility](#6-reproducibility) · [7. Known limits](#7-known-limits)

## 1. A run in short

| Step | What happens | Tool (runs in R) |
|---|---|---|
| 1. Read input | Check and clean the pasted or chosen FASTA | FHInumerics |
| 2. Search | Find similar database sequences, both strands | DECIPHER `SearchIndex` + `AlignPairs` |
| 3. Closest neighbours | Keep the closest N per query that cover enough of it | FHInumerics |
| 4. Align | Multiple alignment of queries + neighbours, trimmed to the shared region | DECIPHER `AlignSeqs` (or MAFFT) |
| 5. SNP distances | Count differences at positions where both have A/C/G/T | ape `dist.dna` |
| 6. Tree | Neighbour-joining with bootstrap support, midpoint rooted | ape `bionj`, phangorn `midpoint` |

The window, the command line (`cli/fhinumerics`) and the tests all use the same analysis
code (`R/compare.R`), so they give identical results for the same input and settings.

## 2. Using the window

The menu at the top left switches between four pages:

- **Outbreak detection**: compare new sequences to a database (sections 2.1-2.2).
- **Database viewer**: browse, filter and chart a database, and pick samples from it
  (section 2.3).
- **Alignment**: see the aligned sequences, coloured where they differ (section 2.4).
- **Sets**: collect sequences from several places (databases, runs, files) to view,
  align, save and export them together (section 2.5).

### 2.1 Outbreak detection: input and settings

- **Paste sequence(s) in FASTA format**, or **choose a FASTA file** (`.fa`, `.fasta`,
  `.fas`, `.fna`, `.txt`). A chosen file is used instead of pasted text. At most 50
  sequences per run.
- **...or use a set** (2.5): its sequences are the queries, instead of the pasted text and
  the file (choosing a file again, or sending samples from the Database viewer, switches
  back). The note under it says how many sequences are used. Samples in the set that come
  from the database being searched leave out their own entry from their closest matches
  (as samples sent from the Database viewer do); queries, file sequences and samples from
  other databases do not. The run's `run_manifest.yaml` records the set (`input`: name,
  checksum, number of sequences, folder) and, under `exclude_self`, which entry each
  sample left out.
- **Database**: the databases listed in `config.yaml`.
- **Closest sequences per query** (5-50, default 20): how many database sequences per
  query go into the alignment and tree.
- **More settings** (closed by default):
  - *Minimum % of the query a match must cover* (50-100, default 90).
  - *Alignment*: DECIPHER (default, R only) or MAFFT (only if set up, see section 5).
  - *Bootstrap support for the tree* (on by default; off makes the run faster).
- **Save results in**: the folder where a new subfolder is made for each run. Type a path
  or use *Browse...*.
- **Run comparison** starts the run; a progress bar shows the current step. When it is
  done, a message shows where the results were saved. Errors (for example an invalid
  sequence) are shown as a red message and in the *Run log*.
- Samples sent from the Database viewer (see 2.3) arrive in the paste box, with a note
  under it. As long as that text is not changed, each sample's own database entry is
  left out of its closest matches (otherwise every sample would find itself, 0 SNPs).
  This is recorded as `exclude_self: true` in `run_manifest.yaml`.

### 2.2 Outbreak detection: results

With several queries, a **Show results for** list appears above the tabs: *All queries
together* (all queries and their closest sequences in one tree and matrix) or one query
on its own (its own closest sequences, compared over its whole length). The list sets
what the Tree, SNP distances and Closest matches tabs show. A result that could not be
made (for example, queries from different genome regions cannot be compared together)
is marked "(no result)" and its tabs say why.

- **Tree**: the phylogeny. Queries are marked `(QUERY)`. Point at the dot of a sample to
  see a box with all its metadata (an orange ring marks the dot it belongs to); for a
  query, the box shows the closest sequence(s) in the tree and the number of SNPs.
  - *Colour by*: colour the tips by a metadata column (e.g. clade).
  - *Show in labels*: metadata columns added to each label (e.g. date, country).
    Missing values are shown as `?`.
  - *Show bootstrap support (>= 70)*: support values on the branches.
  - *Download PDF*: the tree as shown.
  - Drag up or down over the tree to choose the samples in that part, then press *Show in
    Alignment*.
- **SNP distances**: the SNP matrix, rows and columns in tree order. Queries are in red
  (row names and column headers). Darker blue = fewer differences: 0, 1-2, 3-5, 6-10,
  11-20; 21 or more is white. Point at a number in a query's row or column (an orange
  frame marks it) to see a box with the number of SNPs, over how many positions they were
  counted, and what differs. Positions are base numbers in the query; bases are given as
  row/column. For example:

  ```
  V123 / S456
  0 SNPs over 440 compared positions
  Positions in V123:
  Gaps (not counted):
  1-19 - missing at the start of S456
  46 - C/-
  Ambiguity codes that fit (not counted):
  47 - A/M
  N (not compared): 3 in S456
  ```

  - *SNPs*: the counted differences (see 3.5), e.g. `52 - A/K`.
  - *Gaps (not counted)*: bases missing at the start or end of a sequence, as one range;
    gaps inside a sequence, one line per run of gaps (`46 - C/-`, or `after 120 - -/A`
    when the other sequence has extra bases after query base 120).
  - *Ambiguity codes that fit (not counted)*: an ambiguity code that can be the other base
    (`47 - A/M`), so it is not a SNP, but the base is uncertain.
  - *N (not compared)*: how many Ns each sequence has.

  Up to 20 lines are listed per part. For a reverse-complemented query the numbers are on
  its reverse complement. A small **orange corner** on a number means the pair also
  differs by gaps, Ns or ambiguity codes that fit: "0 SNPs" there does not mean the two
  sequences are identical letter for letter. The note above the table says
  over how many aligned positions the differences were counted.
  All rows are shown on one page.
  - Click **names** to choose rows (click again to un-choose), then *Show in Alignment*.
  - Click a **number** in a query's row or column to open just that pair in the Alignment
    page, compared with the query, at their first SNP.
  - *Add to set* (next to *Show in Alignment*) adds the chosen rows to the open set (2.5);
    the same button appears when you drag over the tree. *Add the N sequences of this
    result to the set*, above the tabs, adds all sequences of the shown result.
- **Closest matches**: all search hits for every query (also those not in the tree), with
  database metadata, filterable per column. `selected = TRUE` means the hit is in the tree.
  If more database sequences are equally close than fit in the tree, a yellow note says so
  (for example "101 database sequences are equally close ...; 20 are in the tree").
- **Run log**: what happened in each step, with times, warnings and errors.

### 2.3 Database viewer

Left side:

- **Database**: which database to show.
- **Date column**: the column used for the date filter and the chart over time. Offered
  are the columns where most values are dates (`2023-05-14`, `2023-05` or `2023`); the
  default is `date_column` from `config.yaml` (section 5), else a column called `date`.
  Without such a column there is no date filter.
- **Sample date**: show samples from this date range. A date given only as a year or month
  (`2019`, `2019-04`) counts as the first day of that year or month. *Include samples
  without a date* (on by default).
- **Choose samples**, in any combination:
  - click rows in the table (Ctrl or Shift for several; click again to un-choose);
  - paste sample IDs, or choose a text file with one ID per line (spaces, commas and
    semicolons also separate IDs), then press *Choose these IDs*. IDs that are not in the
    database are listed in a yellow note above the table;
  - *Choose all shown*: all samples that pass the current filters; *Clear chosen*.

  The line under the buttons counts the chosen samples. Chosen samples stay chosen when
  the date range changes, even if they are no longer shown.
- **Charts and downloads use**: all samples shown in the table (after all filters), or
  only the chosen samples.
- **Table (CSV)** and **Sequences (FASTA)**: download those samples. The CSV file is made
  for Excel with Norwegian settings: columns separated by `;`, decimal comma (`2,5`), and
  UTF-8 with a byte order mark so that æ, ø and å show correctly. (Earlier versions
  used `,` and a decimal point; scripts that read it must use `;`, e.g. R `read.csv2`.)
- **Send chosen samples to Outbreak detection** (at most 50): copies their sequences into
  the paste box of Outbreak detection, selects the same database and switches page.
  Press *Run comparison* there.
- **Add chosen samples to the set**: adds them to the open set (2.5). The line under the
  button says which set that is.

Right side (tabs):

- **Table**: the database metadata; filter with the boxes under the column names.
- **Over time**: number of samples per month or year, optionally coloured by a column
  (the 12 most common values; the rest is "Other"). Samples without a date are counted
  in a note.
- **Counts**: number of samples per value of a column (the 30 most common), optionally
  coloured by another column.
- **Cross-table**: counts of one column against another, with totals.
- **Map**: samples per country (circle size = number of samples), from the
  `country_column`. Point at a circle to see the country, its number of samples and the
  number per value of the column chosen in *In the box, count per* (default: `genotype`,
  else the database's `default_colour`); an orange ring marks the circle. Country names
  the map does not know are listed above it. *Zoom to the samples* (on by default) or
  show the whole world. Works without internet.

Each chart can be downloaded as a PDF.

### 2.4 Alignment

Shows an alignment with every sequence compared with an **anchor**: by default the
**majority consensus** (no anchor sequence chosen). The blue line above the alignment
always says what is compared with what.

Left side (settings):

- **Show**: *Results of the last run* (with a *Result* list that follows *Show results
  for*), *An earlier run folder* (*Choose a run folder...*, then a result in it), *A FASTA
  file (aligned or not)*, or *A set* (2.5; choosing a set here also opens it on the Sets
  page). A file whose sequences all have the same length is shown as an alignment (gaps
  as `-`). Otherwise the page says it is not aligned: press **Align these sequences**. The
  alignment runs in the background (see "Aligning" in 2.5), and the file is then shown,
  with the sequences that were turned (reverse complement) named in the yellow note.
  A set without an alignment has an **Align the set** button here too.
- **Show a query and its closest** (run results): a query and its N closest sequences by
  SNPs. When more are as close as number N, they are all shown, and a note says so.
  **Show all sequences** goes back to all rows. **Add the shown rows to the set** adds the
  rows shown (for example a query and its closest) to the open set (2.5).
- **Compare with**: the majority consensus or one sequence (the anchor). Also: right-click
  a name or a row and choose *Use ... as anchor* (or *Compare with the majority
  consensus*).
- **Variable columns only**: hide columns where no sequence differs from the anchor
  (SNPs, ambiguity codes, indels; no data does not count). A red tick on the column
  ruler marks where columns were skipped.
- **Hide columns that are a gap in the anchor**: hide columns where the anchor has a gap
  inside it (extra bases in other sequences). A purple mark on the anchor row shows
  where; point at it to see how many columns.
- **Consensus settings** (closed by default): per column, the most common of A/C/G/T
  wins; an ambiguity code splits its vote (R = half A, half G); N, gaps and missing ends
  do not vote. A column where most sequences have a gap is `-`; a column without any
  base is N. *When the most common base is not clear*: an ambiguity code (A/G tie = R,
  default) or N. *The most common base must reach a threshold* (50-100%): otherwise the
  ambiguity code of the fewest bases that together reach it (or N).
- **Download**: *FASTA (whole)*, *FASTA (columns on screen)*, *PDF (on screen)*, *Picture
  (PNG)* and *Interactive page (HTML)*: one file that opens in any browser without
  FHInumerics (zoom, move, sort and point at cells work; it keeps the anchor chosen when
  it was made). It contains the sequences and metadata shown, so share it like the data.
  Rows are in the order shown, the anchor first.

Above the alignment:

- A bar with **Sort rows** (tree order for run results, fewest or most SNPs, name, sample
  date), **Go to position** (a base number in the anchor), **Zoom in**, **Zoom out**,
  **Whole alignment**, **Letters**, **Previous / Next difference** (a SNP or indel in any
  row, or only in the row you clicked), and which columns are on screen.
- The **overview**: how many sequences differ from the anchor in each column, along the
  whole alignment (red where 3 or more differ; the compared region of a run is unshaded).
  The blue box is the part on screen: click or drag it to move.
- The colour key.

The alignment:

- The anchor is the top row. Queries are in red. After each name: the sample date (if the
  data has dates) and the number of SNPs to the anchor, counted like the SNP table (over
  the compared region for run results, so with a sequence as anchor it is the number in
  the SNP table; over the whole alignment otherwise).
- Only what differs from the anchor is drawn: **red-orange** a SNP; **dark blue** an
  ambiguity code that cannot match the anchor (counted as a SNP); **light blue** an
  ambiguity code that fits (not counted); **purple** an indel (a gap where the anchor has
  a base, or extra bases where the anchor has a gap; never counted as a SNP); **light
  grey** no data (N, or the missing start or end of a sequence); **grey with a letter**
  the anchor has no base there, so there is nothing to compare. Zoomed in to letters, `·`
  means the same as the anchor. For run results, the alignment is shown before trimming:
  columns outside the compared region are greyed. Sequences left out by trimming
  (covering less than half of the compared region) have grey names.
- **Mouse wheel** over the alignment zooms at the pointer; **drag** moves;
  **double-click** shows the whole alignment. When zoomed out, a pixel covers several
  columns and shows the most important colour among them (SNP first).
- **Point** at a cell: the sequence, its date and metadata, its totals (SNPs, indel
  events and columns, ambiguity codes, % no data), the alignment column, its base number
  in this sequence and in the anchor (as in the SNP table's box; for a reverse-complemented
  query, on its reverse complement), the two bases, what the colour means, and how many
  sequences differ there. **Click** a name to choose that row for Previous/Next difference.

### 2.5 Sets

A set is a collection of sequences from anywhere in FHInumerics: database samples (from one
or more databases), queries and their closest sequences from runs, rows of the Alignment
page, and FASTA files. Sets are kept in the window until you press *Save*.

Adding sequences: press **Add to set** on another page (Database viewer: chosen samples;
Outbreak detection: names chosen in the SNP table, samples chosen in the tree, or the
whole result; Alignment page: the rows shown), or choose a **FASTA file** (aligned or not)
on the Sets page. They are added to the **open set** - the one chosen in the *Sets* list -
or to a new set if none is open. A message says how many were added, and which were
already in the set or renamed:

- The same database sample (same database, same ID) is kept once. A query or file
  sequence with the same name and the same sequence as one in the set is kept once; with
  another sequence it is added as `<name>_2`.
- A query or file sequence named like a database sample becomes `<name>_query` (as in
  runs, 3.1). When two databases both have a sample with one name, both become
  `<name>_<database>`, using the database's `short_name` from `config.yaml` (section 5),
  e.g. `AY649695.1_demo`.
- Queries keep the orientation they were given; a query that a run reverse-complemented
  has a note in the *note* column.

Left side:

- **Sets**: the sets in the window and the sets saved in the results folder
  (`<results folder>/sets/`); click one to open it. *(not saved)* means it has changes
  that are not saved. **New set** makes an empty set; **Open a set folder...** opens a set
  saved somewhere else.
- **Name of the open set** and **Save in** (the folder; type a path or use *Browse...*;
  default `<results folder>/sets`). **Save** writes the set to `<Save in>/<name>/`. If
  another saved set is already there, you are asked whether to replace it. Changing the
  name and saving again saves a copy under the new name; the old folder stays.
- **Delete** asks first, then removes the set from the window and its saved files (only
  the files FHInumerics writes; any other file you put in the folder stays).
- **Add sequences**: a FASTA file. If it is aligned (all the same length), the set keeps
  its alignment when the set was empty.
- **Alignment**: whether the set is aligned, **Align the set**, and **Show in Alignment**
  (opens the set on the Alignment page).
- **Use in Outbreak detection**: makes the open set the input of Outbreak detection (2.1)
  and switches page. If rows are chosen in the table, only those are used (also *Use the
  chosen in Outbreak detection* above the table). A run takes at most 50 sequences: for a
  larger set, choose up to 50 rows.
- **Download**: *Sequences (FASTA)* (as given, not aligned) and *Table (CSV)*, in the same
  CSV form as the Database viewer's (`;`, decimal comma, UTF-8).

Right side: the set's name, how many sequences come from where, whether it is saved and
aligned, and a table of its sequences: `id` (the name in the set), `source` (database,
query or file), `database`, `original_id` (the ID in its database, or the name as given),
`added_from`, `note`, then the metadata columns of all its databases. Click rows to
choose them and press **Remove from the set**.

**Aligning.** *Align the set* (or *Align these sequences* for a FASTA file) first finds
sequences that are reverse complements: DECIPHER `OrientNucleotides` compares each
sequence with six complete HAV genomes, one per subgenotype (IA, IB, IIA, IIB, IIIA,
IIIB; `data/orientation_panel.fasta`), and turns those that match better reversed and
complemented. Then the sequences are aligned with DECIPHER `AlignSeqs`, as in runs. This
runs in the background: a box under the menu shows the progress on every page, the window
stays usable, and **Cancel** stops it. For more than 500 sequences you are asked first
(the whole demo database, 1,773 sequences, took about 1.5 minutes on a fast computer).
When it is done, a message names the sequences that were turned. In the Alignment page
they are numbered on their reverse complement, and the yellow note says so.

The alignment is kept with the set (and saved with it). Removing sequences keeps it
(their rows are taken out). Sequences added later are not in it: the set is then "out of
date", the Alignment page shows the sequences that are in the alignment and names the
others, and *Align the set* makes a new alignment. A set made from run results, Alignment
page rows or an aligned FASTA file starts with that alignment, so it can be shown at once.

**A saved set** is a folder with:

| File | Content |
|---|---|
| `sequences.fasta` | the sequences by their names in the set, as given (not aligned) |
| `metadata.tsv` | the table shown on the Sets page (tab-separated) |
| `set.yaml` | name, dates, FHInumerics version, number of sequences, the set's checksum, its sources (each database with its file checksum and number of samples; queries; files), what was done when, and the alignment's details |
| `aligned.fasta` | the alignment, if the set has one (turned sequences reverse-complemented, as listed in `set.yaml`) |

The set's checksum is the MD5 of one line `name<TAB>SEQUENCE` per sequence, sorted by
name (byte order), each ending in a newline (sequence in capitals, without gaps). It
changes when a sequence or a name changes, not when the order or the metadata changes.

The browser asks "Leave page?" when you close the tab while a set has unsaved changes
(closing the tab stops FHInumerics, and unsaved sets are lost).

Closing the browser tab stops FHInumerics (when started with the start file).

## 3. The analysis step by step

### 3.1 Reading the input (`R/fasta.R`)

- Works like NCBI BLAST: text without a `>` header is one unnamed sequence (`query_1`);
  the name is the first word after `>`.
- Spaces, line breaks, digits, `.`, `/` and `-` are removed (so GenBank flat-file text can
  be pasted); letters are made upper case; `U` becomes `T`.
- Allowed letters: `ACGT` and the IUPAC ambiguity codes `RYKMSWBDHVN`. Anything else,
  empty entries and duplicate names stop the run with a message naming the sequence.
- A query with the same name as a database sequence gets `_query` added to its name
  (e.g. `V123_query`), so both can be told apart in the tree and SNP matrix. The run log
  and the manifest (`renamed`) list these.

### 3.2 Database search (`R/search.R`)

- The database FASTA gets a DECIPHER k-mer index (k = 11), built once and kept in `cache/`;
  it is rebuilt automatically when the FASTA file changes (checked by MD5 checksum).
- Each query **and its reverse complement** is searched (`SearchIndex`, up to 500 hits per
  query), and every hit is aligned to the query (`AlignPairs`). For each query/database
  pair the strand with the better score is kept.
- The hit table has the same columns as `blastn -outfmt 6`. Percent identity (`pident`)
  = identical positions / alignment length including gaps, as in BLAST. Query coverage
  (`qcovhsp`) = % of the query inside the aligned region.
- Optional: `search: blast` uses BLAST+ (blastn through rBLAST) if it is installed.

### 3.3 Choosing the closest neighbours (`R/compare.R`)

- Only hits that cover at least *minimum coverage* % of the query are used.
- Hits are ranked by percent identity, then alignment score; **ties are broken by
  database order**, so the same input always gives the same neighbours.
- The top N per query are kept (*closest sequences per query*); all queries' neighbours
  together form the sequence set for the alignment.
- With several queries, steps 3.4-3.6 are done for **all queries together** and again
  for **each query alone** (that query and its own N neighbours).
- Ties are counted: how many database sequences have the same identity and length as the
  best hit. If more are equally close than are shown, the run log and the *Closest
  matches* tab say so; `closest_summary.tsv` has the numbers.
- A query whose best hit is on the minus strand is **reverse-complemented** before the
  alignment (noted in the run log and the manifest).

### 3.4 Alignment and trimming

- Queries + neighbours are aligned with **DECIPHER `AlignSeqs`** (default), or with
  **MAFFT `--auto`** (pymafft, Python) when `aligner: mafft` is set.
- The alignment is **trimmed to the region every query covers**: from the latest query
  start to the earliest query end. All-gap columns are removed. Database sequences that
  have A/C/G/T in less than 50% of that region are dropped (named in the run log).
  This way all distances are counted over the same positions.
- A query on its own is trimmed to its own length, so it is often compared over more
  positions than in the result for all queries together.
- If the queries do not overlap at all (e.g. different genome regions), there is no
  result for all queries together; each query's own result is still made.

### 3.5 SNP distances (`R/snps.R`)

- `phangorn::dist.hamming(ratio = FALSE, exclude = "none")`: for each pair, the number of
  positions where the two bases **cannot be the same**. IUPAC ambiguity codes count by
  what they can be:

  | Bases | Counted? | Why |
  |---|---|---|
  | A / T | yes | different bases |
  | A / M (A or C) | no | M can be A |
  | A / K (G or T) | yes | K cannot be A |
  | R (A/G) / Y (C/T) | yes | no base in common |
  | R (A/G) / M (A/C) | no | both can be A |
  | anything / N | no | N can be anything |
  | anything / gap | no | gaps are not counted (see below) |

- For two sequences with only A/C/G/T at the compared positions this gives the same numbers
  as before (and as snp-dists); ambiguity codes can only add differences, never remove
  them. Before (up to FHInumerics 0.2.0) ambiguity codes were skipped like N:
  `ape::dist.dna(model = "N", pairwise.deletion = TRUE)`.
- **Gaps are not counted** (an insertion or deletion is not a SNP), and positions where
  one sequence has a gap are not compared. A sequence that is shorter than the query is
  compared only where it has bases.
- `compared_sites.tsv` gives, for each pair, how many positions were compared: positions
  where both have a base other than N or a gap (ambiguity codes included). This is the
  denominator of the SNP count, and the "compared positions" in the hover box.
- The tree is not affected by this: it uses its own distances (3.6), which still skip
  ambiguity codes.
- Checks: `tests/run_tests.R` checks each case in the table, and that the hover box lists
  as many positions as the matrix counts. `tests/validate_snps.R <fasta>` compares the
  count with an independent one and with the old count on a whole database.

### 3.6 Tree

- Neighbour-joining: `ape::bionj` on TN93 distances (pairwise deletion). Negative branch
  lengths are set to 0.
- Bootstrap: 100 replicates (`ape::boot.phylo`, random seed 1, so values are the same
  every run), when bootstrap is on and there are at least 4 sequences. Support is the %
  of replicates.
- Rooted at the midpoint (`phangorn::midpoint`).

## 4. Output files

Each run makes its own folder `<date>-<time>_<first query name>` in the results folder:

| File | Content |
|---|---|
| `query.fasta` | the queries as given (after cleaning) |
| `search_hits.tsv` | all search hits with database metadata; `selected` = in the tree |
| `closest_summary.tsv` | per query: best identity, best mismatches, how many are equally close, how many shown |
| `alignment.fasta` | the trimmed alignment used for SNPs and tree |
| `alignment_untrimmed.fasta` | the whole alignment before trimming, with real names (shown in the Alignment page) |
| `snp_matrix.tsv` | SNP differences between all sequences |
| `compared_sites.tsv` | number of compared positions per pair (both have a base other than N or a gap) |
| `tree.nwk` | the tree (Newick, bootstrap support as node labels) |
| `tree_metadata.tsv` | metadata of every tip in the tree |
| `run_manifest.yaml` | settings, method description, database file and checksum, queries, reverse-complemented queries, FHInumerics version and the version of every package and tool |
| `work/` | intermediate files (internal IDs, untrimmed alignment) |
| `per_query/NN_<name>/` | with several queries: each query's own `alignment.fasta`, `snp_matrix.tsv`, `compared_sites.tsv`, `tree.nwk`, `tree_metadata.tsv` (and `work/`). The files above are for all queries together |

When the queries came from a set (2.1), `run_manifest.yaml` has an `input` block (the
set's name, checksum, number of sequences, folder and whether it was saved; "all" or
"chosen rows"), and `exclude_self` under `settings` lists, for each sample from the
searched database, the entry it left out.

`run_manifest.yaml` lists under `results` each result (all queries together, then each
query) with its folder, number of neighbours and compared positions, the compared region
(`compared_region`, columns of `alignment_untrimmed.fasta`) and the sequences left out by
trimming (`dropped`), or why it failed.

## 5. Settings (`config.yaml`)

`config.yaml` is local to each copy (not in git); `config.example.yaml` is the template
with comments. Relative paths are relative to the folder.

| Setting | Meaning | Default |
|---|---|---|
| `output_dir` | default results folder (`~` = home / Documents) | `~/FHInumerics-results` |
| `cache_dir` | where search indexes are kept | `cache` |
| `browse_roots` | starting folders for the *Browse...* button | home and drives |
| `settings: n_neighbours` | closest sequences per query | 20 |
| `settings: min_coverage` | minimum % of the query a hit must cover | 90 |
| `settings: search` | `decipher`, or `blast` (needs BLAST+ and `blast_dir`) | `decipher` |
| `settings: aligner` | `decipher` (R only), or `mafft` (needs Python 3.9-3.13) | `decipher` |
| `settings: bootstrap` | 100 bootstrap replicates | `true` |
| `settings: threads` | CPU threads for search and alignment | 4 |
| `python` | Python to use for MAFFT, if several are installed | found automatically |
| `databases` | list of databases: `name`, `fasta`, `metadata` (.tsv/.csv), `id_column`, `display_columns`, `default_labels`, `default_colour`, `date_column` (Database viewer: default date column; optional, default `date` or the first column holding dates), `country_column` (map; default `country`), `short_name` (sets: added to a sample's name when two databases in a set have a sample with that name; optional, default the name with spaces and other characters replaced by `_`) | demo database |

**MAFFT (optional):** with `aligner: mafft`, the first start also installs pymafft
(pinned in `python/requirements.txt`) into `tools/python-packages`, using a Python
3.9-3.13 on the computer. Python 3.14 cannot be used (pymafft has no ready-made package for
it). Without a suitable Python, FHInumerics stops with a message explaining this.

## 6. Reproducibility

- The same input, settings, database and package versions give the same neighbours, SNP
  matrix, tree and bootstrap values every time. `run_manifest.yaml` records all of these.
- Different R versions install their own package versions (whatever CRAN/Bioconductor
  provide for that R version); the versions used are in each manifest.

Comparisons on the test queries (`tests/`, 2026-10-28):

- **DECIPHER vs MAFFT alignment:** fragments: identical alignment, SNP matrix and tree.
  Whole genome: same closest match and SNP count (2); distant neighbours (~420 SNPs)
  differ by 2-4 SNPs because gaps are placed differently, so the tree differs slightly.
- **pymafft vs MAFFT program:** byte-identical alignments.
- **ape vs snp-dists:** identical SNP matrices (also with N, ambiguity codes and gaps).
- **DECIPHER search vs blastn:** fragments: same neighbours. Whole genome: same closest
  match and SNP counts; 4 of the 20 more distant neighbours (~94% identity) differed,
  because the two compute percent identity from slightly different local alignments and
  many hits are within 0.1% of each other at that distance.
- **Tree vs the earlier maximum-likelihood version (IQ-TREE):** not the same method, so
  the trees differ; the neighbour-joining tree is fully reproducible.

## 7. Known limits

- If many database sequences are identical, more can be equally close than fit in the
  tree; all are listed under *Closest matches* and ties are broken by database order.
- Mixing whole genomes and fragments in one run compares only the overlapping region.
- No maximum-likelihood trees: IQ-TREE as a Python package (piqtree) did not give
  identical trees on repeated runs in our tests, and FastTree is not available as an R or
  Python package.
- At most 50 query sequences per run.
- The map places a country by its name. Common variants (United Kingdom, Viet Nam,
  United States, ...) and notes in brackets ("Egypt (suspected)") are handled; other
  unknown names are listed above the map.
