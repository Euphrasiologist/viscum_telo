/software/team301/seqkit fx2tab -n telo_window_db/viscum_telo_windows.fa \
  | awk '
BEGIN { OFS="\t"; print "window","chrom_end" }
{
  id=$1
  end=id
  sub(/\.fa.*/, "", end)
  gsub(/_/, ":", end)
  print NR, end
}
' > telo_window_id_map.tsv

awk '
BEGIN { OFS="\t" }
NR==FNR && FNR>1 {
  map[$1]=$2
  next
}
FNR==1 {
  print "window","chrom_end","n_hits","total_aligned_bp","best_bitscore"
  next
}
{
  print $1,map[$1],$2,$3,$4
}
' telo_window_id_map.tsv blast_out/core80_telo_window_summary.tsv \
  > blast_out/core80_telo_window_summary.mapped.tsv

  awk '
BEGIN { OFS="\t" }
NR==FNR && FNR>1 {
  map[$1]=$2
  next
}
FNR==1 {
  print "window","chrom_end","n_hits","total_aligned_bp","best_bitscore"
  next
}
{
  print $1,map[$1],$2,$3,$4
}
' telo_window_id_map.tsv blast_out/full188_telo_window_summary.tsv \
  > blast_out/full188_telo_window_summary.mapped.tsv
