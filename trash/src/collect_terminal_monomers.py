#!/usr/bin/env python3
import csv
from pathlib import Path

terminal_window = 2_000_000
min_width, max_width = 170, 205

out = open("results/all_terminal_170_205.members.tsv", "w")
print(
    "chromosome", "chrom", "start0", "end", "cluster_id",
    "seq_id", "width", "strand",
    sep="\t", file=out
)

for members in sorted(Path("results").glob("SUPER_*/tables/*.members.bed")):
    chromosome = members.parents[1].name
    candidates_file = (
        Path("results") / chromosome / "candidates" /
        f"{chromosome}.terminal_candidate_clusters.txt"
    )

    if not candidates_file.exists():
        continue

    candidates = {
        x.strip() for x in candidates_file.read_text().splitlines()
        if x.strip()
    }

    with members.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        continue

    left = min(int(float(r["start0"])) for r in rows)
    right = max(int(float(r["end"])) for r in rows)

    for row in rows:
        start = int(float(row["start0"]))
        end = int(float(row["end"]))
        width = int(float(row["width"]))

        terminal = (
            start - left < terminal_window or
            right - end < terminal_window
        )

        if (
            terminal
            and row["cluster_id"] in candidates
            and min_width <= width <= max_width
        ):
            print(
                chromosome,
                row["chrom"],
                start,
                end,
                row["cluster_id"],
                row["seq_id"],
                width,
                row["strand"],
                sep="\t",
                file=out,
            )

out.close()
