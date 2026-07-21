#!/usr/bin/env python3
import sys
import csv
from pathlib import Path

members_tsv = Path(sys.argv[1])
meta_tsv = Path(sys.argv[2])

meta = {}

with open(meta_tsv) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        meta[row["seq_id"]] = row

print(
    "chrom",
    "start0",
    "end",
    "cluster_id",
    "seq_id",
    "representative",
    "width",
    "strand",
    sep="\t",
)

with open(members_tsv) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        seq_id = row["seq_id"]

        if seq_id not in meta:
            raise ValueError(f"No metadata found for {seq_id}")

        m = meta[seq_id]

        print(
            m["chrom"],
            m["start0"],
            m["end"],
            row["cluster_id"],
            seq_id,
            row["representative"],
            m["width"],
            m["strand"],
            sep="\t",
        )
