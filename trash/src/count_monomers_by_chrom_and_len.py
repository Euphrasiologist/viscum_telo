import pandas as pd

d = pd.read_csv(
    "results/all_terminal_170_205.members.tsv",
    sep="\t"
)

counts = (
    d.groupby("chromosome")
     .size()
     .reset_index(name="n_monomers")
)

counts.to_csv(
    "results/terminal_170_205.counts_by_chromosome.tsv",
    sep="\t",
    index=False
)

d["length_class"] = pd.cut(
    d["width"],
    bins=[0, 169, 179, 189, 199, 205, 10**9],
    labels=["<170", "170-179", "180-189",
            "190-199", "200-205", ">205"]
)

(
    d.groupby(["chromosome", "length_class"], observed=True)
     .size()
     .reset_index(name="n_monomers")
     .to_csv(
         "results/terminal_monomers.by_chromosome_length.tsv",
         sep="\t",
         index=False
     )
)
