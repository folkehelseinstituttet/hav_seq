"""Python step of FHInumerics: MAFFT alignment (pymafft).

pymafft runs MAFFT inside Python, so no separate program is needed.
Called from R (R/python.R); can also be run by hand to check a result:

  python fhinumerics_tools.py align IN.fasta OUT.fasta [--strategy auto] [--threads 1]
  python fhinumerics_tools.py versions
"""

import argparse
import json
import platform
import sys
from importlib import metadata


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif line:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(records, path, width=70):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def align(args):
    """Same as: mafft --<strategy> --thread <threads> IN.fasta"""
    import pymafft

    seqs = read_fasta(args.input)
    result = pymafft.align(list(seqs.items()), strategy=args.strategy, seq_type="nuc",
                           threads=args.threads)
    write_fasta([(s.name, s.sequence.upper()) for s in result], args.output)


def versions(args):
    out = {"python": platform.python_version(), "python_path": sys.executable}
    try:
        out["pymafft"] = metadata.version("pymafft")
    except metadata.PackageNotFoundError:
        out["pymafft"] = None
    print(json.dumps(out))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("align", help="multiple alignment with MAFFT (pymafft)")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--strategy", default="auto")
    p.add_argument("--threads", type=int, default=1)
    p.set_defaults(func=align)

    p = sub.add_parser("versions", help="print package versions as JSON")
    p.set_defaults(func=versions)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
