# After generating the dats.

Generate the fasta:

```bash
awk '
$3 >= 150 && $3 <= 220 && $4 > 100 {
    seq = $14
    print ">" FILENAME "|start=" $1 "|end=" $2 "|period=" $3 "|copies=" $4 "|cons_size=" $5
    print seq
}
' *.dat > viscum_trf_188bp_candidates.fa
```

Align with mafft

```bash
/software/team301/mafft-7.525-with-extensions/core/mafft --auto --adjustdirection viscum_trf_188bp_candidates.fa \
    > viscum_trf_188bp_candidates.mafft.fa
```

Extract core with trimal:

```bash
/software/team301/trimal/source/trimal -in ./viscum_trf_188bp_candidates.mafft_gt50.fa \
  -out viscum_trf_188bp_candidates.mafft_gt50_trimal.fa \
  -automated1
```

Then run blast:

```bash
bash ./run_blast.bash
bash ./map_chromosomes.bash
bash ./map_chromosomes2.bash
```
