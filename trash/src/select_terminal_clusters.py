#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select candidate terminal repeat clusters."
    )
    parser.add_argument("cluster_stats", type=Path)
    parser.add_argument(
        "--min-terminal-hits",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--min-terminal-fraction",
        type=float,
        default=0.40,
    )
    parser.add_argument(
        "--min-enrichment",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--min-total-hits",
        type=int,
        default=20,
    )
    return parser.parse_args()


args = parse_args()

with args.cluster_stats.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        total_hits = int(row["total_hits"])
        terminal_hits = int(row["terminal_hits"])
        terminal_fraction = float(row["terminal_fraction"])
        enrichment = float(row["terminal_enrichment"])

        if (
            total_hits >= args.min_total_hits
            and terminal_hits >= args.min_terminal_hits
            and terminal_fraction >= args.min_terminal_fraction
            and enrichment >= args.min_enrichment
        ):
            print(row["cluster_id"])
