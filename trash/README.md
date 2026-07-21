# Notes on running here

TRASH output from Alex B via Robin.

0. Turn BED into fasta

```bash
python trash_csv_to_fasta.py mistletoe_SUPER_1_TRASH.bed
```

1. Cluster sequences, testing on a single chromosme first.

```bash
bsub \
  -J cdhit_SUPER_1_TRASH \
  -q normal \
  -n 32 \
  -R "span[hosts=1] select[mem>4000] rusage[mem=4000]" \
  -M 4000 \
  -o logs/cdhit_SUPER_1_TRASH.%J.out \
  -e logs/cdhit_SUPER_1_TRASH.%J.err \
  "bash cluster.bash"
```

2. Parse results into a table

```bash
python parse_cd_clusters.py mistletoe_SUPER_1_TRASH.cdhit90.fa.clstr > mistletoe_SUPER_1_TRASH.cdhit90.members.tsv
```

3. Convert results to BED

```bash
python cdhit_members_to_bed.py \
  mistletoe_SUPER_1_TRASH.cdhit90.members.tsv \
  mistletoe_SUPER_1_TRASH.meta.tsv \
> mistletoe_SUPER_1_TRASH.cdhit90.members.bed
```

Files
