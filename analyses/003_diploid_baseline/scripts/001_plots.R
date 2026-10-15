#!/usr/bin/env Rscript
# SNP/kb per pair by class from results/pair_baseline.tsv.
# Run from analyses/003_diploid_baseline:  Rscript scripts/001_plots.R
# Writes results/plots/snp_per_kb.{png,pdf}. Descriptive only; pairs share genomes.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)})
out <- "results/plots"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
INK2 <- "#52514e"; GRID <- "#e1e0d9"; SURF <- "#fcfcfb"
theme_set(theme_minimal(base_size = 11) + theme(
  plot.background = element_rect(fill = SURF, colour = NA), panel.background = element_rect(fill = SURF, colour = NA),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
  panel.grid.major.x = element_blank(), axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
  legend.position = "top", legend.title = element_blank(), strip.text = element_text(hjust = 0, size = 11)))

d <- read_tsv("results/pair_baseline.tsv", show_col_types = FALSE) %>%
  mutate(cls = factor(cls, c("diploid-diploid", "Cbp-diploid", "Cbp-Cbp")),
         cg = ifelse(grepl("CgMosher", pair), "pair includes CgMosher", "other pair"),
         set = paste(set, "set"))
med <- d %>% group_by(set, cls) %>% summarise(m = median(snp_per_kb), .groups = "drop")
p <- ggplot(d, aes(cls, snp_per_kb)) +
  geom_point(aes(colour = cg, shape = cg), position = position_jitter(width = 0.12, height = 0, seed = 1), size = 2.6) +
  geom_errorbar(data = med, mapping = aes(x = cls, ymin = m, ymax = m), width = 0.5, linewidth = 0.7, inherit.aes = FALSE) +
  facet_wrap(~set) + scale_colour_manual(values = c(`other pair` = "#2a78d6", `pair includes CgMosher` = "#eb6834")) +
  scale_shape_manual(values = c(`other pair` = 16, `pair includes CgMosher` = 18)) +
  labs(x = NULL, y = "SNPs per kb of syntenic sequence", caption = "black bar = median; Cbp-diploid = Cbp subgenome vs diploid of the same set")
for (e in c("png", "pdf")) ggsave(file.path(out, paste0("snp_per_kb.", e)), p, width = 9, height = 3.9, dpi = 200)
