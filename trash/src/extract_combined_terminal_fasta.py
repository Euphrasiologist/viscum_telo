#!/usr/bin/env python3
import csv
from pathlib import Path

comp = str.maketrans(
    "ACGTRYMKSWBDHVN",
    "TGCAYRKMSWVHDBN"
)

metadata = {}

with open("results/all_terminal_170_205.members.tsv") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        metadata[row["seq_id"]] = row

with open("results/all_terminal_170_205.fa", "w") as out:
    for fasta in sorted(Path("results").glob("SUPER_*/fasta/*.fa")):
        name = None
        sequence = []

        def emit():
            if name not in metadata:
                return

            row = metadata[name]
            seq = "".join(sequence).upper()

            if row["strand"] == "-":
                seq = seq.translate(comp)[::-1]

            print(
                f'>{row["chromosome"]}|{name}|'
                f'{row["chrom"]}:{int(row["start0"]) + 1}-{row["end"]}|'
                f'cluster={row["cluster_id"]}|width={row["width"]}',
                file=out,
            )
            print(seq, file=out)

        with fasta.open() as handle:
            for line in handle:
                line = line.rstrip()

                if line.startswith(">"):
                    if name is not None:
                        emit()
                    name = line[1:].split()[0]
                    sequence = []
                else:
                    sequence.append(line)

            if name is not None:
                emit()
