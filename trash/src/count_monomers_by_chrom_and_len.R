library(ggplot2)
library(readr)

d <- read_tsv(
  "results/terminal_monomers.by_chromosome_length.tsv"
)

plot <- ggplot(d, aes(chromosome, n_monomers, fill = length_class)) +
  geom_col() +
  labs(x = NULL, y = "Number of terminal monomers", fill = "Length") +
  theme_bw()

ggsave(plot, file = "results/terminal_monomers.by_chromosome_length.png")
