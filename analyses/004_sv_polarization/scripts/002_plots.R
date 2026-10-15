#!/usr/bin/env Rscript
# Plots for 004_sv_polarization from results/polarization_summary.tsv.
# Run from analyses/004_sv_polarization:  Rscript scripts/002_plots.R
# Writes results/02_plots/{indel_by_lineage,indel_ratio_per_pair,outcome_mix}.{png,pdf}
# Descriptive only (no tests); pairs share genomes so are not independent.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})

out <- "results/02_plots"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
INS <- "#2a78d6"; DEL <- "#eb6834"; INK2 <- "#52514e"; GRID <- "#e1e0d9"; SURF <- "#fcfcfb"
OUTCOME <- c(both = "#2a78d6", `single-frame` = "#86b6ef", ambiguous = "#c3c2b7", conflict = "#eb6834")
LIN <- c("Cbp-Co", "Co diploids", "Cbp-Cr", "Cr/Cg diploids")

theme_set(theme_minimal(base_size = 11) + theme(
  plot.background = element_rect(fill = SURF, colour = NA), panel.background = element_rect(fill = SURF, colour = NA),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
  axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
  plot.title = element_text(face = "plain", size = 12), legend.position = "top", legend.title = element_blank()))
save <- function(p, name, w, h) for (e in c("png", "pdf")) ggsave(file.path(out, paste0(name, ".", e)), p, width = w, height = h, dpi = 200,
                                                device = if (e == "pdf") cairo_pdf else NULL)

lineage <- function(g) ifelse(grepl("^Cbp", g), paste0("Cbp-", substring(g, nchar(g) - 1)),
                              ifelse(grepl("^Co", g), "Co diploids", "Cr/Cg diploids"))
df <- read_tsv("results/polarization_summary.tsv", show_col_types = FALSE)

ind <- df %>% filter(kind == "indel", grepl("^(ins|del)_[AB]$", final)) %>%
  mutate(type = substr(final, 1, 3), genome = ifelse(grepl("_A$", final), A, B),
         lineage = factor(lineage(genome), LIN),
         conf = ifelse(confidence == "both", "both frames", "single frame"))
dc <- ind %>% filter(class == "dip_cbp")

# README check
tot <- dc %>% group_by(lineage, type) %>% summarise(n = sum(n), .groups = "drop") %>% pivot_wider(names_from = type, values_from = n)
print(tot)

# 1. del vs ins per lineage by confidence
b <- dc %>% group_by(lineage, type, conf) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(grp = paste0(type, ": ", conf), conf = factor(conf, c("single frame", "both frames")))
lab <- tot %>% mutate(top = pmax(ins, del), lab = sprintf("ins:del %.2f", ins / del))
b <- b %>% mutate(x = as.numeric(lineage) + ifelse(type == "del", -0.21, 0.21))
p1 <- ggplot(b, aes(x, n, fill = type, alpha = conf)) +
  geom_col(width = 0.4, position = position_stack(), colour = SURF, linewidth = 0.8) +
  geom_text(data = lab, aes(as.numeric(lineage), top * 1.05, label = lab), inherit.aes = FALSE, colour = INK2, size = 3.3) +
  scale_x_continuous(breaks = seq_along(LIN), labels = LIN) +
  scale_fill_manual(values = c(del = DEL, ins = INS), labels = c(del = "deletion", ins = "insertion")) +
  scale_alpha_manual(values = c(`both frames` = 1, `single frame` = 0.45)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1), order = 1), alpha = guide_legend(override.aes = list(fill = "#52514e"), order = 2)) +
  labs(x = NULL, y = "polarized indels >= 50 bp (summed over pairs)",
       title = "Insertions vs deletions by lineage (diploid \u00d7 Cbp pairs)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + theme(panel.grid.major.x = element_blank())
save(p1, "indel_by_lineage", 7.5, 4.4)

# 2. per-pair ratio
rp <- dc %>% group_by(pair, lineage, type) %>% summarise(n = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = n) %>% filter(!is.na(ins), !is.na(del)) %>%
  mutate(ratio = ins / del, mosher = ifelse(grepl("Mosher", pair), "pair includes a Mosher assembly", "no Mosher genome in pair"))
p2 <- ggplot(rp, aes(ratio, lineage, colour = mosher, shape = mosher)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#898781") +
  geom_point(position = position_jitter(height = 0.12, width = 0, seed = 1), size = 2.6) +
  scale_x_log10(breaks = c(0.3, 0.5, 1, 2), labels = c("0.3", "0.5", "1", "2")) +
  scale_y_discrete(limits = rev(LIN)) +
  scale_colour_manual(values = c(INS, DEL)) + scale_shape_manual(values = c(16, 18)) +
  labs(x = "ins:del ratio per pair (log scale; 1 = balanced)", y = NULL,
       title = "Per-pair ins:del ratio (diploid × Cbp pairs)") +
  theme(panel.grid.major.y = element_blank())
save(p2, "indel_ratio_per_pair", 7.5, 4.2)

# 3. outcome mix
oc <- df %>% filter(kind %in% c("indel", "INV", "TRANS", "INVTR")) %>%
  mutate(group = ifelse(kind == "indel", "indels", "INV/TRANS/INVTR"),
         outcome = case_when(final %in% c("amb", "none") | confidence == "none" ~ "ambiguous",
                             final == "conflict" | confidence == "conflict" ~ "conflict",
                             confidence == "both" ~ "both", TRUE ~ "single-frame")) %>%
  group_by(class, group, outcome) %>% summarise(n = sum(n), .groups = "drop") %>%
  group_by(class, group) %>% mutate(frac = n / sum(n), tot = sum(n)) %>% ungroup() %>%
  mutate(cls = factor(recode(class, dip_dip = "diploid × diploid", dip_cbp = "diploid × Cbp", cbp_cbp = "Cbp × Cbp"),
                      c("diploid × diploid", "diploid × Cbp", "Cbp × Cbp")),
         outcome = factor(outcome, names(OUTCOME)), group = factor(group, c("indels", "INV/TRANS/INVTR")))
nlab <- oc %>% distinct(cls, group, tot)
p3 <- ggplot(oc, aes(frac, cls, fill = outcome)) +
  geom_col(colour = SURF, linewidth = 0.8, width = 0.75, position = position_stack(reverse = TRUE)) +
  geom_text(data = nlab, aes(1.01, cls, label = paste0("n=", trimws(format(tot, big.mark = ",")))), inherit.aes = FALSE,
            hjust = 0, colour = INK2, size = 3) +
  facet_wrap(~group) + scale_fill_manual(values = OUTCOME) +
  scale_x_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  labs(x = "fraction of events", y = NULL) + theme(panel.grid.major.y = element_blank(), strip.text = element_text(hjust = 0, size = 11))
save(p3, "outcome_mix", 9, 3.6)
