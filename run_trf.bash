for f in telo_windows/*.fa; do
  /software/team301/repeat-annotation/trf "$f" 2 7 7 80 10 50 220 -d -h
  mv *.dat trf_out/
done
