# Demo database: where it comes from

`sequences.fasta` (1773 sequences) and `metadata.tsv` were made on 2026-10-22 with

    python3 scripts/prepare_demo_db.py <hav_seq>/data/nextclade_datasets/hav-vp1-2b-lineages data/demo

from https://github.com/folkehelseinstituttet/hav_seq at commit abb2b92 (2026-10-07):

- sequences: `community_sequences.fasta`
- metadata: tip attributes in `tree.json` (date, geo_loc_name, clade_membership,
  lineage_phylo, author). `geo_loc_name` "Country: place" is split into `country` and
  `location`; partial dates such as `2017-05-XX` become `2017-05`; `20XX-XX-XX` becomes empty.

That dataset is FHI's copy of the NextClade community dataset
`community/masphl-bioinformatics/hav/vp1-2b-junction` (GenBank, accessed 2025-02-03)
with FHI lineage labels added. Public sequences only. Most are VP1-2B fragments
(349-531 bp); about 10% are complete genomes (~7.4 kb).
