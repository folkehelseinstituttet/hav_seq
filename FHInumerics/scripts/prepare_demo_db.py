#!/usr/bin/env python3
"""Build the bundled demo database from a NextClade dataset folder.

Reads community_sequences.fasta (sequences) and tree.json (tip metadata) and
writes <out_dir>/sequences.fasta and <out_dir>/metadata.tsv in the format the
app expects: one FASTA record per ID, and a TSV whose first column is `id`.

Usage: python3 prepare_demo_db.py <nextclade_dataset_dir> <out_dir>
"""
import csv
import json
import re
import sys
from pathlib import Path


def tips(node):
    if "children" not in node:
        yield node
    for child in node.get("children", []):
        yield from tips(child)


def attr(node, key):
    value = node.get("node_attrs", {}).get(key, {}).get("value")
    return "" if value is None else str(value).strip()


def tidy_date(raw):
    """'2017-05-XX' -> '2017-05', '1985-XX-XX' -> '1985', '20XX-XX-XX' -> ''."""
    parts = []
    for part in raw.split("-"):
        if not part.isdigit():
            break
        parts.append(part)
    return "-".join(parts)


def short_author(raw):
    """'Cleary,N.G., Bryant,P.W. and ...; ...' -> 'Cleary et al.'"""
    first_list = raw.split(";")[0].strip()
    if not first_list:
        return ""
    surname = first_list.split(",")[0].strip()
    several = " and " in first_list or first_list.count(",") > 1
    return surname + (" et al." if several else "")


def main(dataset_dir, out_dir):
    dataset_dir, out_dir = Path(dataset_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tree = json.load(open(dataset_dir / "tree.json"))

    rows = {}
    for tip in tips(tree["tree"]):
        geo = attr(tip, "geo_loc_name")
        country, _, location = geo.partition(":")
        date = tidy_date(attr(tip, "date"))
        rows[tip["name"]] = {
            "id": tip["name"],
            "date": date,
            "year": date[:4],
            "country": country.strip(),
            "location": location.strip(),
            "clade": attr(tip, "clade_membership"),
            "lineage": attr(tip, "lineage_phylo"),
            "author": short_author(attr(tip, "author")),
            "study": tip.get("node_attrs", {}).get("author", {}).get("title", "").strip(),
        }

    seqs, name = {}, None
    for line in open(dataset_dir / "community_sequences.fasta"):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line:
            seqs[name].append(line)

    missing = [n for n in seqs if n not in rows]
    if missing:
        sys.exit(f"{len(missing)} FASTA IDs have no metadata, e.g. {missing[:3]}")

    with open(out_dir / "sequences.fasta", "w") as fh:
        for n, chunks in seqs.items():
            seq = "".join(chunks)
            fh.write(f">{n}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    fields = list(next(iter(rows.values())).keys())
    with open(out_dir / "metadata.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for n in seqs:
            writer.writerow(rows[n])
    print(f"Wrote {len(seqs)} sequences and metadata rows to {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
