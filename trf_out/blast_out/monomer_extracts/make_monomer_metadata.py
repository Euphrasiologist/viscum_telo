import re
from pathlib import Path

fa = Path("monomer_extracts/full188.monomers.filtered.fa")
out = Path("monomer_extracts/full188.metadata.tsv")

pat = re.compile(r"^>([^:]+):(\d+)-(\d+):(\d+)-(\d+):([+-])")

with fa.open() as f, out.open("w") as o:
    o.write("id\twindow\taccession\twindow_start\twindow_end\tmonomer_start\tmonomer_end\tstrand\tchr\tarm\n")
    for line in f:
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        m = pat.match(h)
        if not m:
            continue

        acc, wstart, wend, mstart, mend, strand = m.groups()
        window = f"{acc}:{wstart}-{wend}"

        # infer arm from whether window starts at 1
        arm = "left" if int(wstart) == 1 else "right"

        # chr needs a map later; temporary accession-based grouping
        o.write(f"{h}\t{window}\t{acc}\t{wstart}\t{wend}\t{mstart}\t{mend}\t{strand}\tNA\t{arm}\n")
