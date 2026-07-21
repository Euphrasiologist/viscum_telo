#!/usr/bin/env bash
set -euo pipefail

WINDOW_DIR="telo_windows"
OUT_DIR="dotplots/flexidot_k15"
LOG_DIR="dotplots/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

mapfile -t FASTAS < <(
    find "$WINDOW_DIR" -maxdepth 1 -type f -name '*.fa' | sort
)

if [[ -z "${LSB_JOBINDEX:-}" ]]; then
    echo "Run this as an LSF array job."
    exit 1
fi

index=$((LSB_JOBINDEX - 1))
fasta="${FASTAS[$index]}"
name=$(basename "$fasta" .fa)

echo "[$(date)] Processing: $fasta"

flexidot \
    -i "$fasta" \
    -m 0 \
    -k 15 \
    -o "$OUT_DIR/$name"

echo "[$(date)] Finished: $name"
