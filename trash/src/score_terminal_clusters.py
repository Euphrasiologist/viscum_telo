#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rank CD-HIT clusters by enrichment near the minimum and maximum "
            "observed TRASH coordinates."
        )
    )
    parser.add_argument("members_bed", type=Path)
    parser.add_argument(
        "--terminal-window",
        type=int,
        default=2_000_000,
        help="Distance from either observed boundary in bp [2000000].",
    )
    return parser.parse_args()


args = parse_args()

rows = []

with args.members_bed.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required = {
        "chrom",
        "start0",
        "end",
        "cluster_id",
        "seq_id",
        "representative",
        "width",
        "strand",
    }

    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise SystemExit(
            "Unexpected header.\n"
            f"Observed: {reader.fieldnames}\n"
            f"Expected: {sorted(required)}"
        )

    for row in reader:
        row["start0"] = int(float(row["start0"]))
        row["end"] = int(float(row["end"]))
        row["width"] = int(float(row["width"]))
        rows.append(row)

if not rows:
    raise SystemExit("No rows found")

# These are observed TRASH-call boundaries, not guaranteed assembly boundaries.
bounds = {}

for row in rows:
    chrom = row["chrom"]
    start = row["start0"]
    end = row["end"]

    if chrom not in bounds:
        bounds[chrom] = [start, end]
    else:
        bounds[chrom][0] = min(bounds[chrom][0], start)
        bounds[chrom][1] = max(bounds[chrom][1], end)

print(
    "Observed coordinate boundaries:",
    file=sys.stderr,
)

for chrom, (left, right) in sorted(bounds.items()):
    print(
        f"  {chrom}: {left:,}-{right:,}",
        file=sys.stderr,
    )

cluster_stats = defaultdict(
    lambda: {
        "total": 0,
        "terminal": 0,
        "left": 0,
        "right": 0,
        "widths": [],
        "chromosomes": set(),
    }
)

expected_terminal_fraction = {}

for chrom, (left_bound, right_bound) in bounds.items():
    span = right_bound - left_bound

    if span <= 0:
        expected_terminal_fraction[chrom] = 1.0
    else:
        terminal_span = min(span, 2 * args.terminal_window)
        expected_terminal_fraction[chrom] = terminal_span / span

for row in rows:
    chrom = row["chrom"]
    cluster = row["cluster_id"]
    start = row["start0"]
    end = row["end"]

    left_bound, right_bound = bounds[chrom]

    left_distance = max(0, start - left_bound)
    right_distance = max(0, right_bound - end)

    is_left = left_distance < args.terminal_window
    is_right = right_distance < args.terminal_window
    is_terminal = is_left or is_right

    stat = cluster_stats[cluster]

    stat["total"] += 1
    stat["terminal"] += int(is_terminal)
    stat["left"] += int(is_left)
    stat["right"] += int(is_right)
    stat["widths"].append(row["width"])
    stat["chromosomes"].add(chrom)

output = []

for cluster, stat in cluster_stats.items():
    total = stat["total"]
    terminal = stat["terminal"]

    terminal_fraction = terminal / total

    expected_values = [
        expected_terminal_fraction[chrom]
        for chrom in stat["chromosomes"]
    ]

    expected_fraction = statistics.mean(expected_values)

    enrichment = (
        terminal_fraction / expected_fraction
        if expected_fraction > 0
        else float("nan")
    )

    # Rewards both terminal abundance and terminal enrichment.
    score = terminal * math.log2(enrichment + 1)

    widths = stat["widths"]

    output.append(
        {
            "cluster_id": cluster,
            "total_hits": total,
            "terminal_hits": terminal,
            "left_hits": stat["left"],
            "right_hits": stat["right"],
            "terminal_fraction": terminal_fraction,
            "expected_terminal_fraction": expected_fraction,
            "terminal_enrichment": enrichment,
            "median_width": statistics.median(widths),
            "mean_width": statistics.mean(widths),
            "min_width": min(widths),
            "max_width": max(widths),
            "n_chromosomes": len(stat["chromosomes"]),
            "score": score,
        }
    )

output.sort(
    key=lambda row: (
        row["score"],
        row["terminal_hits"],
        row["terminal_fraction"],
    ),
    reverse=True,
)

columns = [
    "cluster_id",
    "total_hits",
    "terminal_hits",
    "left_hits",
    "right_hits",
    "terminal_fraction",
    "expected_terminal_fraction",
    "terminal_enrichment",
    "median_width",
    "mean_width",
    "min_width",
    "max_width",
    "n_chromosomes",
    "score",
]

writer = csv.DictWriter(
    sys.stdout,
    fieldnames=columns,
    delimiter="\t",
    lineterminator="\n",
)

writer.writeheader()
writer.writerows(output)
