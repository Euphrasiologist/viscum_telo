#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 input/trash_csv/mistletoe_SUPER_N_TRASH.bed" >&2
    exit 1
fi

INPUT=$1

CDHIT=/software/team301/cdhit/cd-hit-est
SEQKIT=/software/team301/seqkit
MAFFT=/software/team301/mafft-7.525-with-extensions/core/mafft

THREADS=${THREADS:-32}
CDHIT_IDENTITY=${CDHIT_IDENTITY:-0.90}
TERMINAL_WINDOW=${TERMINAL_WINDOW:-2000000}

MIN_TOTAL_HITS=${MIN_TOTAL_HITS:-20}
MIN_TERMINAL_HITS=${MIN_TERMINAL_HITS:-20}
MIN_TERMINAL_FRACTION=${MIN_TERMINAL_FRACTION:-0.40}
MIN_TERMINAL_ENRICHMENT=${MIN_TERMINAL_ENRICHMENT:-3.0}

MAX_ALIGNMENT_MEMBERS=${MAX_ALIGNMENT_MEMBERS:-1000}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SRC="${ROOT}/src"

INPUT=$(realpath "${INPUT}")
BASENAME=$(basename "${INPUT}" .bed)

if [[ ${BASENAME} =~ (SUPER_[0-9]+)_TRASH$ ]]; then
    CHROM=${BASH_REMATCH[1]}
else
    echo "Could not obtain SUPER chromosome name from: ${BASENAME}" >&2
    exit 1
fi

OUTDIR="${ROOT}/results/${CHROM}"
FASTA_DIR="${OUTDIR}/fasta"
CDHIT_DIR="${OUTDIR}/cdhit"
TABLE_DIR="${OUTDIR}/tables"
CANDIDATE_DIR="${OUTDIR}/candidates"
ALIGNMENT_DIR="${OUTDIR}/alignments"
CLUSTER_ALIGNMENT_DIR="${ALIGNMENT_DIR}/clusters"

mkdir -p \
    "${FASTA_DIR}" \
    "${CDHIT_DIR}" \
    "${TABLE_DIR}" \
    "${CANDIDATE_DIR}" \
    "${CLUSTER_ALIGNMENT_DIR}"

FASTA="${FASTA_DIR}/${BASENAME}.fa"
META="${FASTA_DIR}/${BASENAME}.meta.tsv"

CDHIT_PREFIX="${CDHIT_DIR}/${BASENAME}.cdhit90"
CDHIT_FASTA="${CDHIT_PREFIX}.fa"
CDHIT_CLSTR="${CDHIT_FASTA}.clstr"

MEMBERS_TSV="${TABLE_DIR}/${BASENAME}.cdhit90.members.tsv"
MEMBERS_BED="${TABLE_DIR}/${BASENAME}.cdhit90.members.bed"
TERMINAL_STATS="${TABLE_DIR}/${BASENAME}.cdhit90.terminal_cluster_stats.tsv"

CANDIDATE_CLUSTERS="${CANDIDATE_DIR}/${CHROM}.terminal_candidate_clusters.txt"
CANDIDATE_SEQUENCES="${CANDIDATE_DIR}/sequences"
REPRESENTATIVES="${CANDIDATE_SEQUENCES}/candidate_representatives.fa"
REPRESENTATIVE_ALIGNMENT="${ALIGNMENT_DIR}/${CHROM}.candidate_representatives.mafft.fa"

echo "============================================================"
echo "Chromosome:             ${CHROM}"
echo "Input:                  ${INPUT}"
echo "Result directory:       ${OUTDIR}"
echo "Threads:                ${THREADS}"
echo "CD-HIT identity:        ${CDHIT_IDENTITY}"
echo "Terminal window:        ${TERMINAL_WINDOW}"
echo "============================================================"

if [[ ! -s "${FASTA}" || ! -s "${META}" ]]; then
    echo "[1/7] Converting TRASH table to FASTA and metadata"

    python "${SRC}/trash_csv_to_fasta.py" \
        "${INPUT}" \
        "${FASTA}" \
        "${META}" \
        --prefix "${BASENAME}"
else
    echo "[1/7] FASTA and metadata already exist; skipping"
fi

if [[ ! -s "${CDHIT_FASTA}" || ! -s "${CDHIT_CLSTR}" ]]; then
    echo "[2/7] Running cd-hit-est"

    "${CDHIT}" \
        -i "${FASTA}" \
        -o "${CDHIT_FASTA}" \
        -c "${CDHIT_IDENTITY}" \
        -n 8 \
        -d 0 \
        -M 0 \
        -T "${THREADS}"
else
    echo "[2/7] CD-HIT output already exists; skipping"
fi

if [[ ! -s "${MEMBERS_TSV}" ]]; then
    echo "[3/7] Parsing CD-HIT cluster membership"

    python "${SRC}/parse_cd_clusters.py" \
        "${CDHIT_CLSTR}" \
        > "${MEMBERS_TSV}"
else
    echo "[3/7] Membership table already exists; skipping"
fi

if [[ ! -s "${MEMBERS_BED}" ]]; then
    echo "[4/7] Linking cluster membership to genomic coordinates"

    python "${SRC}/cdhit_members_to_bed.py" \
        "${MEMBERS_TSV}" \
        "${META}" \
        > "${MEMBERS_BED}"
else
    echo "[4/7] Coordinate-linked membership table exists; skipping"
fi

if [[ ! -s "${TERMINAL_STATS}" ]]; then
    echo "[5/7] Scoring terminal enrichment"

    python "${SRC}/score_terminal_clusters.py" \
        "${MEMBERS_BED}" \
        --terminal-window "${TERMINAL_WINDOW}" \
        > "${TERMINAL_STATS}"
else
    echo "[5/7] Terminal cluster statistics exist; skipping"
fi

if [[ ! -s "${CANDIDATE_CLUSTERS}" ]]; then
    echo "[6/7] Selecting terminal candidate clusters"

    python "${SRC}/select_terminal_clusters.py" \
        "${TERMINAL_STATS}" \
        --min-total-hits "${MIN_TOTAL_HITS}" \
        --min-terminal-hits "${MIN_TERMINAL_HITS}" \
        --min-terminal-fraction "${MIN_TERMINAL_FRACTION}" \
        --min-enrichment "${MIN_TERMINAL_ENRICHMENT}" \
        > "${CANDIDATE_CLUSTERS}"
else
    echo "[6/7] Candidate cluster list exists; skipping"
fi

N_CANDIDATES=$(grep -cve '^[[:space:]]*$' "${CANDIDATE_CLUSTERS}" || true)

echo "Selected candidate clusters: ${N_CANDIDATES}"

if [[ "${N_CANDIDATES}" -eq 0 ]]; then
    echo "No candidate clusters passed the current thresholds."
    echo "Inspect: ${TERMINAL_STATS}"
    exit 0
fi

if [[ ! -s "${REPRESENTATIVES}" ]]; then
    echo "[7/7] Extracting and orienting candidate sequences"

    python "${SRC}/extract_candidate_sequences.py" \
        "${MEMBERS_BED}" \
        "${FASTA}" \
        "${CANDIDATE_CLUSTERS}" \
        "${CANDIDATE_SEQUENCES}" \
        --max-members "${MAX_ALIGNMENT_MEMBERS}" \
        --orient-by-strand
else
    echo "[7/7] Candidate sequences already extracted; skipping"
fi

if [[ ! -s "${REPRESENTATIVE_ALIGNMENT}" ]]; then
    echo "[alignment] Aligning candidate cluster representatives"

    "${MAFFT}" \
        --auto \
        --thread "${THREADS}" \
        "${REPRESENTATIVES}" \
        > "${REPRESENTATIVE_ALIGNMENT}"
else
    echo "[alignment] Representative alignment exists; skipping"
fi

for SAMPLE_FASTA in "${CANDIDATE_SEQUENCES}"/cluster_samples/*.sample.fa; do
    [[ -e "${SAMPLE_FASTA}" ]] || continue

    CLUSTER=$(basename "${SAMPLE_FASTA}" .sample.fa)
    OUTPUT_ALIGNMENT="${CLUSTER_ALIGNMENT_DIR}/${CLUSTER}.sample.mafft.fa"

    if [[ -s "${OUTPUT_ALIGNMENT}" ]]; then
        echo "[alignment] ${CLUSTER} alignment exists; skipping"
        continue
    fi

    N_SEQUENCES=$("${SEQKIT}" stats -T "${SAMPLE_FASTA}" \
        | awk 'NR == 2 {print $4}')

    echo "[alignment] ${CLUSTER}: ${N_SEQUENCES} sampled sequences"

    if [[ "${N_SEQUENCES}" -gt 500 ]]; then
        "${MAFFT}" \
            --thread "${THREADS}" \
            --parttree \
            --retree 1 \
            "${SAMPLE_FASTA}" \
            > "${OUTPUT_ALIGNMENT}"
    else
        "${MAFFT}" \
            --thread "${THREADS}" \
            --auto \
            "${SAMPLE_FASTA}" \
            > "${OUTPUT_ALIGNMENT}"
    fi
done

echo "Pipeline complete for ${CHROM}"
echo "Terminal statistics:      ${TERMINAL_STATS}"
echo "Candidate clusters:       ${CANDIDATE_CLUSTERS}"
echo "Representative alignment: ${REPRESENTATIVE_ALIGNMENT}"
echo "Cluster alignments:       ${CLUSTER_ALIGNMENT_DIR}"
