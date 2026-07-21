#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a TRASH CSV table to FASTA plus positional metadata."
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("output_metadata", type=Path)
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix used to create stable sequence IDs.",
    )
    return parser.parse_args()


args = parse_args()

args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
args.output_metadata.parent.mkdir(parents=True, exist_ok=True)

n_sequences = 0

with (
    args.input_csv.open() as inp,
    args.output_fasta.open("w") as fasta,
    args.output_metadata.open("w") as meta,
):
    reader = csv.DictReader(inp)

    required = {
        "Chr",
        "realstart",
        "realend",
        "seq",
        "strand",
        "width",
    }

    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise SystemExit(
            f"Unexpected input columns: {reader.fieldnames}\n"
            f"Required columns: {sorted(required)}"
        )

    print(
        "seq_id",
        "chrom",
        "start0",
        "end",
        "strand",
        "width",
        sep="\t",
        file=meta,
    )

    for i, row in enumerate(reader, 1):
        seq_id = f"{args.prefix}_{i:09d}"

        start0 = int(float(row["realstart"])) - 1
        end = int(float(row["realend"]))
        width = int(float(row["width"]))
        sequence = row["seq"].strip().upper()

        if not sequence:
            continue

        print(f">{seq_id}", file=fasta)
        print(sequence, file=fasta)

        print(
            seq_id,
            row["Chr"],
            start0,
            end,
            row["strand"],
            width,
            sep="\t",
            file=meta,
        )

        n_sequences += 1

print(f"Wrote {n_sequences:,} sequences to {args.output_fasta}")
print(f"Wrote metadata to {args.output_metadata}")
