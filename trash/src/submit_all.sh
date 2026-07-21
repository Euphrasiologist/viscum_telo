#!/usr/bin/env bash

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

mkdir -p "${ROOT}/logs"

for INPUT in "${ROOT}"/input/trash_csv/mistletoe_SUPER_*_TRASH.bed; do
    [[ -e "${INPUT}" ]] || continue

    BASENAME=$(basename "${INPUT}" .bed)

    if [[ ${BASENAME} =~ (SUPER_[0-9]+)_TRASH$ ]]; then
        CHROM=${BASH_REMATCH[1]}
    else
        echo "Skipping unrecognised filename: ${INPUT}" >&2
        continue
    fi

    echo "Submitting ${CHROM}"

    # bsub \
    #     -J "telo_${CHROM}" \
    #     -q normal \
    #     -n 32 \
    #     -R "span[hosts=1] select[mem>4000] rusage[mem=4000]" \
    #     -M 4000 \
    #     -o "${ROOT}/logs/${CHROM}.%J.out" \
    #     -e "${ROOT}/logs/${CHROM}.%J.err" \
    #     "THREADS=32 bash '${ROOT}/src/run_one_chromosome.sh' '${INPUT}'"

    env THREADS=32 bash "${ROOT}/src/run_one_chromosome.sh" "${INPUT}"
done
