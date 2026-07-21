#!/usr/bin/env python3

import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path


COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract representatives and sampled members from selected "
            "CD-HIT clusters."
        )
    )
    parser.add_argument("members_bed", type=Path)
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("candidate_clusters", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument(
        "--max-members",
        type=int,
        default=1000,
        help="Maximum sampled members per cluster [1000].",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
    )
    parser.add_argument(
        "--orient-by-strand",
        action="store_true",
        help="Reverse-complement sequences whose TRASH strand is '-'.",
    )
    return parser.parse_args()


args = parse_args()
random.seed(args.seed)

args.output_directory.mkdir(parents=True, exist_ok=True)
sample_directory = args.output_directory / "cluster_samples"
sample_directory.mkdir(parents=True, exist_ok=True)

with args.candidate_clusters.open() as handle:
    candidates = {
        line.strip()
        for line in handle
        if line.strip()
    }

if not candidates:
    raise SystemExit("No candidate clusters were supplied")

# Store metadata for representatives and reservoir-sampled members.
representatives = {}
samples = defaultdict(list)
cluster_counts = defaultdict(int)
selected_metadata = {}

with args.members_bed.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        cluster = row["cluster_id"]

        if cluster not in candidates:
            continue

        seq_id = row["seq_id"]
        cluster_counts[cluster] += 1

        metadata = {
            "cluster": cluster,
            "chrom": row["chrom"],
            "start0": int(float(row["start0"])),
            "end": int(float(row["end"])),
            "strand": row["strand"],
            "width": int(float(row["width"])),
            "representative": row["representative"],
        }

        if row["representative"] == "yes":
            representatives[cluster] = seq_id
            selected_metadata[seq_id] = metadata

        # Reservoir sampling gives each cluster member equal probability.
        current = samples[cluster]
        seen = cluster_counts[cluster]

        if len(current) < args.max_members:
            current.append(seq_id)
            selected_metadata[seq_id] = metadata
        else:
            replacement = random.randrange(seen)

            if replacement < args.max_members:
                old_seq_id = current[replacement]
                current[replacement] = seq_id
                selected_metadata.pop(old_seq_id, None)
                selected_metadata[seq_id] = metadata

# Always retain representatives, even if reservoir replacement removed them.
for cluster, seq_id in representatives.items():
    selected_metadata.setdefault(
        seq_id,
        {
            "cluster": cluster,
            "chrom": "unknown",
            "start0": -1,
            "end": -1,
            "strand": ".",
            "width": -1,
            "representative": "yes",
        },
    )

wanted = set(selected_metadata)

representative_handle = (
    args.output_directory / "candidate_representatives.fa"
).open("w")

cluster_handles = {
    cluster: (
        sample_directory / f"{cluster}.sample.fa"
    ).open("w")
    for cluster in candidates
}


def write_sequence(seq_id, sequence):
    if seq_id not in wanted:
        return

    metadata = selected_metadata[seq_id]

    if args.orient_by_strand and metadata["strand"] == "-":
        sequence = reverse_complement(sequence)
        orientation = "oriented"
    else:
        orientation = "raw"

    cluster = metadata["cluster"]

    header = (
        f"{cluster}|{seq_id}|"
        f"{metadata['chrom']}:{metadata['start0'] + 1}-"
        f"{metadata['end']}|strand={metadata['strand']}|"
        f"width={metadata['width']}|{orientation}"
    )

    if representatives.get(cluster) == seq_id:
        print(f">{header}", file=representative_handle)
        print(sequence, file=representative_handle)

    if seq_id in samples[cluster]:
        print(f">{header}", file=cluster_handles[cluster])
        print(sequence, file=cluster_handles[cluster])


name = None
sequence_parts = []

with args.input_fasta.open() as handle:
    for line in handle:
        line = line.rstrip()

        if line.startswith(">"):
            if name is not None:
                write_sequence(name, "".join(sequence_parts))

            name = line[1:].split()[0]
            sequence_parts = []
        else:
            sequence_parts.append(line)

if name is not None:
    write_sequence(name, "".join(sequence_parts))

representative_handle.close()

for handle in cluster_handles.values():
    handle.close()

summary_path = args.output_directory / "candidate_cluster_extraction.tsv"

with summary_path.open("w") as out:
    print(
        "cluster_id",
        "total_members",
        "sampled_members",
        "representative_seq_id",
        sep="\t",
        file=out,
    )

    for cluster in sorted(candidates):
        print(
            cluster,
            cluster_counts[cluster],
            len(samples[cluster]),
            representatives.get(cluster, "NA"),
            sep="\t",
            file=out,
        )
