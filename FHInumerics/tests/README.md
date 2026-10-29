# Test queries (made from the demo database, 2026-10-22)

Run all checks (uses only the demo database, writes to a temporary folder):

    Rscript tests/run_tests.R                  # default settings (DECIPHER alignment)
    Rscript tests/run_tests.R --aligner mafft  # needs Python 3.9-3.13 and pymafft set up

It prints PASS/FAIL per check and "All checks passed" (exit status 0) when everything is
as expected. Run it after every change.

| File | Query | Expected |
|---|---|---|
| `fragments.fasta` | `test_KY782331_3snps`: KY782331.1 (IA) with 3 substitutions (0-based positions 50, 150, 250) | 3 SNPs to KY782331.1; 101 database sequences tie at this distance |
| `fragments.fasta` | `test_AY649695_1_copy`: exact copy of AY649695.1 (IIIA) | 0 SNPs to AY649695.1 |
| `genome_revcomp.fasta` | LC128713.1 (IB, 7412 bp) with 2 substitutions (positions 3000, 3100), reverse-complemented | detected as reverse complement; 2 SNPs to LC128713.1 |

Results on 2026-10-22 (earlier conda version, default settings): all as expected. Fragments run
64 s, genome 33 s, both files together 240 s (trimmed to the shared 432 columns).

R + Python version on 2026-10-28 (command line, then-default settings: DECIPHER search,
MAFFT via pymafft, ape SNPs, ape neighbour-joining tree): all as expected. Fragments 20 s
(including building the search index the first time), genome 8 s. Same alignments and SNP
matrices as with the IQ-TREE build; the tree is identical on repeated runs.

    Rscript cli/fhinumerics --query tests/fragments.fasta --out test_fragments

Ambiguity-aware SNP counting (R/snps.R), 2026-10-28: `run_tests.R` checks each case
(A/T, K/A, Y/A and R/Y differ; M/A, R/A, N and gaps do not) and that the hover box lists
as many positions as the SNP matrix counts (whole-genome query: positions 3001 and 3101),
and the hover text of a small alignment with a SNP, a missing start, a gap in each
sequence, an ambiguity code that fits and an N.
On a whole database:

    Rscript tests/validate_snps.R data/demo/sequences.fasta 100

Demo database: 26 of 1773 sequences have ambiguity codes. 126 aligned (all 26 + 100
others), 7875 pairs: phangorn = independent count (Biostrings IUPAC table) for every pair;
phangorn = old A/C/G/T-only count + ambiguity-code differences for every pair; 76 pairs
(1.0%) went up (+1: 70, +2: 5, +3: 1), all others unchanged.

FHInumerics 0.1.0 on 2026-10-28 (`tests/run_tests.R`, default settings: DECIPHER search and
DECIPHER alignment): all checks passed. For the MAFFT comparison see MANUAL.md section 6.

Alignment viewer (R/msa.R, decision 0011), 2026-10-29: `run_tests.R` checks each cell
class on a small alignment, the consensus rules, variable and anchor-gap columns, base
numbers, subsets and "closest N" ties, and on the real runs, with every sequence as
anchor, that the differences equal the SNP matrix and that SNP and strong-ambiguity cells
are exactly `snp_columns()`. The reverse-complemented genome's SNPs are at its positions
3001 and 3101, as in the SNP hover box.
