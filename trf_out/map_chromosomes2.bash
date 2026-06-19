awk '
BEGIN { OFS="\t"; print "accession","chr_part","chr","part","side_label" }
{
  acc=$1
  chrpart=$0
  sub(/^.*chromosome: /, "", chrpart)
  sub(/[ \t]+[0-9]+$/, "", chrpart)

  if (chrpart ~ /^[0-9]+_[0-9]+$/) {
    split(chrpart,a,"_")
    chr=a[1]; part=a[2]
    if (part > maxpart[chr]) maxpart[chr]=part
    acc_chr[acc]=chr
    acc_part[acc]=part
  }
}
END {
  for (acc in acc_chr) {
    chr=acc_chr[acc]
    part=acc_part[acc]
    label="chr" chr "_internal"
    if (part == 1) label="chr" chr "_left"
    if (part == maxpart[chr]) label="chr" chr "_right"
    print acc, chr "_" part, chr, part, label
  }
}
' ./viscum_map.tsv | sort -k3,3n -k4,4n > viscum_accession_chr_end_map.tsv

awk '
BEGIN { OFS="\t" }
NR==FNR && FNR>1 {
  label[$1]=$5
  part[$1]=$2
  next
}
FNR==1 {
  print $0,"chr_part","chr_end_label"
  next
}
{
  acc=$2
  sub(/:.*/, "", acc)
  print $0,part[acc],label[acc]
}
' viscum_accession_chr_end_map.tsv \
  blast_out/core80_telo_window_summary.mapped.tsv \
  > blast_out/core80_telo_window_summary.labelled.tsv

awk '
BEGIN { OFS="\t" }
NR==FNR && FNR>1 {
  label[$1]=$5
  part[$1]=$2
  next
}
FNR==1 {
  print $0,"chr_part","chr_end_label"
  next
}
{
  acc=$2
  sub(/:.*/, "", acc)
  print $0,part[acc],label[acc]
}
' viscum_accession_chr_end_map.tsv \
  blast_out/full188_telo_window_summary.mapped.tsv \
  > blast_out/full188_telo_window_summary.labelled.tsv
