FA=/lustre/scratch122/tol/data/f/2/e/8/4/d/Viscum_album/assembly/release/drVisAlbu1.1/insdc/GCA_963277665.1.fasta.gz
FAI=${FA}.fai
# 2Mb
W=2000000

mkdir -p telo_windows trf_out dotplots tidk_out

awk -v W=$W '
{
  len[$1]=$2
}
END {
  # biological chromosome starts
  print "OY728119.1:1-" W
  print "OY728125.1:1-" W
  print "OY728130.1:1-" W
  print "OY728135.1:1-" W
  print "OY728140.1:1-" W
  print "OY728145.1:1-" W
  print "OY728150.1:1-" W
  print "OY728155.1:1-" W
  print "OY728159.1:1-" W
  print "OY728163.1:1-" W
  print "OY728167.1:1-" W Viscum album genome assembly, chromosome: B1	1020192845


  # biological chromosome ends
  print "OY728124.1:" len["OY728124.1"]-W+1 "-" len["OY728124.1"]
  print "OY728129.1:" len["OY728129.1"]-W+1 "-" len["OY728129.1"]
  print "OY728134.1:" len["OY728134.1"]-W+1 "-" len["OY728134.1"]
  print "OY728139.1:" len["OY728139.1"]-W+1 "-" len["OY728139.1"]
  print "OY728144.1:" len["OY728144.1"]-W+1 "-" len["OY728144.1"]
  print "OY728149.1:" len["OY728149.1"]-W+1 "-" len["OY728149.1"]
  print "OY728154.1:" len["OY728154.1"]-W+1 "-" len["OY728154.1"]
  print "OY728158.1:" len["OY728158.1"]-W+1 "-" len["OY728158.1"]
  print "OY728162.1:" len["OY728162.1"]-W+1 "-" len["OY728162.1"]
  print "OY728166.1:" len["OY728166.1"]-W+1 "-" len["OY728166.1"]
  print "OY728167.1:" len["OY728167.1"]-W+1 "-" len["OY728166.1"]
}' "$FAI" > telo_regions.txt

while read r; do
  safe=$(echo "$r" | tr ':-' '__')
  /software/team301/samtools/samtools faidx "$FA" "$r" > telo_windows/${safe}.fa
done < telo_regions.txt
