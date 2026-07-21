#!/usr/bin/env bash
set -euo pipefail

# Run from:
# /lustre/scratch122/tol/teams/blaxter/users/mb39/viscum_telo/trf_out/blast_out

BEDTOOLS=${BEDTOOLS:-/software/team301/bedtools2/bin/bedtools}
SEQKIT=${SEQKIT:-/software/team301/seqkit}
MAFFT=${MAFFT:-/software/team301/mafft-7.525-with-extensions/core/mafft}

WINDOW_FA=../telo_window_db/viscum_telo_windows.fa

mkdir -p monomer_extracts

make_bed_extract_align () {
    local prefix=$1
    local blast_tsv=$2
    local min_len=$3
    local max_len=$4

    echo "Processing ${prefix}"

    # BLAST outfmt:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    awk -v min_len="$min_len" -v max_len="$max_len" '
    BEGIN { OFS="\t" }

    NR==FNR {
        map[$1]=$2
        next
    }

    {
        sid = ($2 in map) ? map[$2] : $2

        s=$9
        e=$10
        strand="+"

        if (s > e) {
            tmp=s; s=e; e=tmp
            strand="-"
        }

        start=s-1
        end=e
        len=end-start

        if (len >= min_len && len <= max_len) {
            name=sid ":" start "-" end ":" strand "|qid="$1"|pid="$3"|aln_len="$4"|bits="$12
            print sid,start,end,name,$12,strand
        }
    }
    ' telo_window_numeric_id_map.tsv "$blast_tsv" \
    | sort -k1,1 -k2,2n \
    > monomer_extracts/${prefix}.raw_hits.bed

    # Merge overlapping/near-overlapping hits on same strand.
    # This reduces duplicate BLAST HSPs from the same monomer.
    $BEDTOOLS merge \
        -s \
        -d 30 \
        -i monomer_extracts/${prefix}.raw_hits.bed \
        -c 4,5,6 \
        -o distinct,max,distinct \
        > monomer_extracts/${prefix}.merged_hits.bed

    # Convert merged BED back to BED6 with useful names.
    awk '
    BEGIN { OFS="\t" }
    {
        strand=$6
        name=$1 ":" $2 "-" $3 ":" strand "|bits="$5
        print $1,$2,$3,name,$5,strand
    }
    ' monomer_extracts/${prefix}.merged_hits.bed \
    > monomer_extracts/${prefix}.merged_hits.bed6

    # Extract strand-aware FASTA.
    $BEDTOOLS getfasta \
        -fi "$WINDOW_FA" \
        -bed monomer_extracts/${prefix}.merged_hits.bed6 \
        -s \
        -name \
        -fo monomer_extracts/${prefix}.monomers.fa

    # Filter by final extracted length.
    $SEQKIT seq \
        -m "$min_len" \
        -M "$max_len" \
        monomer_extracts/${prefix}.monomers.fa \
        > monomer_extracts/${prefix}.monomers.filtered.fa

    # Align, allowing MAFFT to reverse-complement anything still misoriented.
    $MAFFT \
        --auto \
        --adjustdirection \
        monomer_extracts/${prefix}.monomers.filtered.fa \
        > monomer_extracts/${prefix}.monomers.filtered.mafft.fa

    # Basic stats.
    $SEQKIT stats \
        monomer_extracts/${prefix}.monomers.fa \
        monomer_extracts/${prefix}.monomers.filtered.fa \
        monomer_extracts/${prefix}.monomers.filtered.mafft.fa \
        > monomer_extracts/${prefix}.stats.tsv

    echo "Done ${prefix}"
}

# Core repeat: expected shorter hits.
make_bed_extract_align \
    core80 \
    core80_vs_telo_windows.tsv \
    50 \
    120

# Full HOR/monomer repeat.
make_bed_extract_align \
    full188 \
    full188_vs_telo_windows.tsv \
    120 \
    230

echo "All done."
