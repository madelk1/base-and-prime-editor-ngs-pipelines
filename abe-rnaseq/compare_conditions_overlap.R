#!/usr/bin/env Rscript

# Compare overlap between ABE8e, HiFi-ABE, and Inlaid-ABE off-target sites

library(readr)
library(dplyr)

setwd("/pub/madelk1/HiFiRNASeq/variantCalls")

cat("========================================\n")
cat("Comparing off-target site overlap\n")
cat("ABE8e vs HiFi-ABE vs Inlaid-ABE\n")
cat("========================================\n\n")

# Load the three condition outputs
ABE <- read_tsv("ABE_offtarget_sites_3stage_2of3.tsv", show_col_types = FALSE)
HiFi <- read_tsv("HiFi_offtarget_sites_3stage_2of3.tsv", show_col_types = FALSE)
Inlaid <- read_tsv("Inlaid_offtarget_sites_3stage_2of3.tsv", show_col_types = FALSE)

cat("Loaded off-target site counts:\n")
cat("  ABE8e:", nrow(ABE), "sites\n")
cat("  HiFi-ABE:", nrow(HiFi), "sites\n")
cat("  Inlaid-ABE:", nrow(Inlaid), "sites\n\n")

# Extract position keys (unique identifiers) - use pos_key not poskey
abe_sites <- ABE$pos_key
hifi_sites <- HiFi$pos_key
inlaid_sites <- Inlaid$pos_key

# ========== PAIRWISE OVERLAPS ==========
cat("============================================================\n")
cat("PAIRWISE OVERLAPS\n")
cat("============================================================\n\n")

# ABE vs HiFi
abe_hifi_overlap <- intersect(abe_sites, hifi_sites)
cat("ABE8e ∩ HiFi-ABE:", length(abe_hifi_overlap), "sites\n")
cat("  % of ABE8e:", round(length(abe_hifi_overlap)/length(abe_sites)*100, 1), "%\n")
cat("  % of HiFi-ABE:", round(length(abe_hifi_overlap)/length(hifi_sites)*100, 1), "%\n\n")

# ABE vs Inlaid
abe_inlaid_overlap <- intersect(abe_sites, inlaid_sites)
cat("ABE8e ∩ Inlaid-ABE:", length(abe_inlaid_overlap), "sites\n")
cat("  % of ABE8e:", round(length(abe_inlaid_overlap)/length(abe_sites)*100, 1), "%\n")
cat("  % of Inlaid-ABE:", round(length(abe_inlaid_overlap)/length(inlaid_sites)*100, 1), "%\n\n")

# HiFi vs Inlaid
hifi_inlaid_overlap <- intersect(hifi_sites, inlaid_sites)
cat("HiFi-ABE ∩ Inlaid-ABE:", length(hifi_inlaid_overlap), "sites\n")
cat("  % of HiFi-ABE:", round(length(hifi_inlaid_overlap)/length(hifi_sites)*100, 1), "%\n")
cat("  % of Inlaid-ABE:", round(length(hifi_inlaid_overlap)/length(inlaid_sites)*100, 1), "%\n\n")

# ========== THREE-WAY OVERLAP ==========
cat("============================================================\n")
cat("THREE-WAY OVERLAP\n")
cat("============================================================\n\n")

all_three_overlap <- intersect(intersect(abe_sites, hifi_sites), inlaid_sites)
cat("ABE8e ∩ HiFi-ABE ∩ Inlaid-ABE:", length(all_three_overlap), "sites\n")
cat("  % of ABE8e:", round(length(all_three_overlap)/length(abe_sites)*100, 1), "%\n")
cat("  % of HiFi-ABE:", round(length(all_three_overlap)/length(hifi_sites)*100, 1), "%\n")
cat("  % of Inlaid-ABE:", round(length(all_three_overlap)/length(inlaid_sites)*100, 1), "%\n\n")

# ========== UNIQUE SITES ==========
cat("============================================================\n")
cat("UNIQUE SITES (Condition-specific off-targets)\n")
cat("============================================================\n\n")

abe_unique <- setdiff(abe_sites, union(hifi_sites, inlaid_sites))
hifi_unique <- setdiff(hifi_sites, union(abe_sites, inlaid_sites))
inlaid_unique <- setdiff(inlaid_sites, union(abe_sites, hifi_sites))

cat("ABE8e-only:", length(abe_unique), "sites (", round(length(abe_unique)/length(abe_sites)*100, 1), "% of ABE8e)\n")
cat("HiFi-ABE-only:", length(hifi_unique), "sites (", round(length(hifi_unique)/length(hifi_sites)*100, 1), "% of HiFi-ABE)\n")
cat("Inlaid-ABE-only:", length(inlaid_unique), "sites (", round(length(inlaid_unique)/length(inlaid_sites)*100, 1), "% of Inlaid-ABE)\n\n")

# ========== VENN-LIKE BREAKDOWN ==========
cat("============================================================\n")
cat("COMPLETE BREAKDOWN (no overlap categories)\n")
cat("============================================================\n\n")

abe_only_cat <- setdiff(abe_sites, union(hifi_sites, inlaid_sites))
hifi_only_cat <- setdiff(hifi_sites, union(abe_sites, inlaid_sites))
inlaid_only_cat <- setdiff(inlaid_sites, union(abe_sites, hifi_sites))
abe_hifi_only <- setdiff(intersect(abe_sites, hifi_sites), inlaid_sites)
abe_inlaid_only <- setdiff(intersect(abe_sites, inlaid_sites), hifi_sites)
hifi_inlaid_only <- setdiff(intersect(hifi_sites, inlaid_sites), abe_sites)
all_three <- intersect(intersect(abe_sites, hifi_sites), inlaid_sites)

cat("ABE8e only:", length(abe_only_cat), "\n")
cat("HiFi-ABE only:", length(hifi_only_cat), "\n")
cat("Inlaid-ABE only:", length(inlaid_only_cat), "\n")
cat("ABE8e + HiFi-ABE (not Inlaid):", length(abe_hifi_only), "\n")
cat("ABE8e + Inlaid-ABE (not HiFi):", length(abe_inlaid_only), "\n")
cat("HiFi-ABE + Inlaid-ABE (not ABE8e):", length(hifi_inlaid_only), "\n")
cat("All three conditions:", length(all_three), "\n\n")

total_check <- length(abe_only_cat) + length(hifi_only_cat) + length(inlaid_only_cat) + 
               length(abe_hifi_only) + length(abe_inlaid_only) + length(hifi_inlaid_only) + 
               length(all_three)
cat("Total unique sites:", total_check, "\n\n")

# ========== SAVE OVERLAP FILES ==========
cat("============================================================\n")
cat("Saving overlap results...\n")
cat("============================================================\n\n")

if(length(all_three) > 0) {
  all_three_data <- ABE %>% filter(pos_key %in% all_three) %>% 
    select(pos_key, chr, POS, REF, ALT, snv, max_treated_edit)
  write_tsv(all_three_data, "overlap_all_three_conditions.tsv")
  cat("Wrote: overlap_all_three_conditions.tsv\n")
}

if(length(abe_only_cat) > 0) {
  abe_unique_data <- ABE %>% filter(pos_key %in% abe_only_cat) %>%
    select(pos_key, chr, POS, REF, ALT, snv, max_treated_edit)
  write_tsv(abe_unique_data, "overlap_ABE8e_only.tsv")
  cat("Wrote: overlap_ABE8e_only.tsv\n")
}

if(length(hifi_only_cat) > 0) {
  hifi_unique_data <- HiFi %>% filter(pos_key %in% hifi_only_cat) %>%
    select(pos_key, chr, POS, REF, ALT, snv, max_treated_edit)
  write_tsv(hifi_unique_data, "overlap_HiFi_only.tsv")
  cat("Wrote: overlap_HiFi_only.tsv\n")
}

if(length(inlaid_only_cat) > 0) {
  inlaid_unique_data <- Inlaid %>% filter(pos_key %in% inlaid_only_cat) %>%
    select(pos_key, chr, POS, REF, ALT, snv, max_treated_edit)
  write_tsv(inlaid_unique_data, "overlap_Inlaid_only.tsv")
  cat("Wrote: overlap_Inlaid_only.tsv\n")
}

cat("\n============================================================\n")
cat("Overlap analysis complete!\n")
cat("============================================================\n")
