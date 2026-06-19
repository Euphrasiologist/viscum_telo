# run the blast pipeline to see where the putative TR's cluster

mkdir -p blast_queries blast_out telo_window_db

rm -f blast_queries/viscum_80bp_core*.fa blast_queries/viscum_188bp_full*.fa

python3 alignment_consensus.py viscum_trf_188bp_candidates.mafft_gt50_trimal.fa \
  > blast_queries/viscum_80bp_core.raw.fa

python3 alignment_consensus.py viscum_trf_188bp_candidates.mafft.fa \
  > blast_queries/viscum_188bp_full.raw.fa

/software/team301/seqkit seq -g blast_queries/viscum_80bp_core.raw.fa \
  | /software/team301/seqkit replace -p '.*' -r 'viscum_80bp_core' \
  > blast_queries/viscum_80bp_core.fa

/software/team301/seqkit seq -g blast_queries/viscum_188bp_full.raw.fa \
  | /software/team301/seqkit replace -p '.*' -r 'viscum_188bp_full' \
  > blast_queries/viscum_188bp_full.fa

# The blastdb

cat ../telo_windows/*.fa > telo_window_db/viscum_telo_windows.fa

/software/team301/ncbi-blast-2.16.0+/bin/makeblastdb \
  -in telo_window_db/viscum_telo_windows.fa \
  -dbtype nucl \
  -parse_seqids \
  -out telo_window_db/viscum_telo_windows

# BLAST core and full repeats...

/software/team301/ncbi-blast-2.16.0+/bin/blastn \
  -query blast_queries/viscum_80bp_core.fa \
  -db telo_window_db/viscum_telo_windows \
  -task blastn-short \
  -evalue 1e-5 \
  -perc_identity 65 \
  -dust no \
  -soft_masking false \
  -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  > blast_out/core80_vs_telo_windows.tsv

/software/team301/ncbi-blast-2.16.0+/bin/blastn \
  -query blast_queries/viscum_188bp_full.fa \
  -db telo_window_db/viscum_telo_windows \
  -task blastn \
  -evalue 1e-10 \
  -perc_identity 65 \
  -dust no \
  -soft_masking false \
  -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  > blast_out/full188_vs_telo_windows.tsv

  # Summarise terminal enrichment per window

awk '
  BEGIN { OFS="\t"; print "window","n_hits","total_aligned_bp","best_bitscore" }
  {
    n[$2]++
    bp[$2]+=$4
    if ($12 > best[$2]) best[$2]=$12
  }
  END {
    for (w in n) print w,n[w],bp[w],best[w]
  }
  ' blast_out/core80_vs_telo_windows.tsv \
  | sort -k1,1 > blast_out/core80_telo_window_summary.tsv

# And full repeat
# 
awk '
  BEGIN { OFS="\t"; print "window","n_hits","total_aligned_bp","best_bitscore" }
  {
    n[$2]++
    bp[$2]+=$4
    if ($12 > best[$2]) best[$2]=$12
  }
  END {
    for (w in n) print w,n[w],bp[w],best[w]
  }
  ' blast_out/full188_vs_telo_windows.tsv \
  | sort -k1,1 > blast_out/full188_telo_window_summary.tsv

# Binned summaries

awk '
BEGIN { OFS="\t"; print "window","bin_10kb","n_hits" }
{
  pos = ($9 < $10 ? $9 : $10)
  bin = int(pos / 10000) * 10000
  key = $2 OFS bin
  n[key]++
}
END {
  for (k in n) print k,n[k]
}
' blast_out/core80_vs_telo_windows.tsv \
  | sort -k1,1 -k2,2n > blast_out/core80_10kb_bins.tsv

awk '
BEGIN { OFS="\t"; print "window","bin_10kb","n_hits" }
{
  pos = ($9 < $10 ? $9 : $10)
  bin = int(pos / 10000) * 10000
  key = $2 OFS bin
  n[key]++
}
END {
  for (k in n) print k,n[k]
}
' blast_out/full188_vs_telo_windows.tsv \
  | sort -k1,1 -k2,2n > blast_out/full188_10kb_bins.tsv
