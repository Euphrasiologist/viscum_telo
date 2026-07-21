#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a CD-HIT .clstr file to a membership table."
    )
    parser.add_argument("cluster_file", type=Path)
    return parser.parse_args()


args = parse_args()

cluster = None

print("cluster_id\tseq_id\trepresentative")

with args.cluster_file.open() as handle:
    for line in handle:
        line = line.strip()

        if line.startswith(">Cluster"):
            cluster_number = line.split()[-1]
            cluster = f"cluster_{cluster_number}"
            continue

        match = re.search(r">(.+?)\.\.\.", line)

        if not match:
            continue

        if cluster is None:
            raise RuntimeError("Cluster member encountered before cluster header")

        seq_id = match.group(1)
        representative = "yes" if line.endswith("*") else "no"

        print(cluster, seq_id, representative, sep="\t")
