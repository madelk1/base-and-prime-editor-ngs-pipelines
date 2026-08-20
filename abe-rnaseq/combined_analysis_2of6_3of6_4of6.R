#!/usr/bin/env Rscript

# Combined RNP + mRNA Off-Target Analysis - CORRECTED
# Loads raw Step 3 files with corrected sample labels
# S01-S03 = mRNA ABE (treated), S04-S06 = mRNA CTRL (control)
# ABE1-3 = RNP ABE (treated), with their own CTRL data

library(readr)
library(dplyr)

setwd("/pub/madelk1/rnaseq/ABE8.8mRNA")

cat("========================================\n")
cat("Combined RNP + mRNA Off-Target Analysis\n")
cat("2/6, 3/6, and 4/6 filtering\n")
cat("CORRECTED: S01-S03=mRNA ABE, ABE1-3=RNP ABE\n")
cat("========================================\n\n")

# ========== LOAD ALL DATA ==========
cat("Loading data...\n\n")

# mRNA files (raw Step 3 output)
s01_file <- "S01.vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt"
s02_file <- "S02.vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt"
s03_file <- "S03.vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt"

# RNP files (at different location)
rnp_dir <- "/dfs6b/pub/madelk1/rnaseq/fastqs/RNAseq_BE_editing/code/GATK4updatedscripts"
abe1_file <- file.path(rnp_dir, "ABE1.vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt")
abe2_file <- file.path(rnp_dir, "ABE2.vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt")
abe3_file <- file.path(rnp_dir, "ABE3.vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt")

# Load mRNA data (S01-S03 = ABE, S04-S06 = CTRL embedded in files)
s01 <- read_tsv(s01_file, show_col_types = FALSE) %>% mutate(sample = "S01", delivery = "mRNA")
s02 <- read_tsv(s02_file, show_col_types = FALSE) %>% mutate(sample = "S02", delivery = "mRNA")
s03 <- read_tsv(s03_file, show_col_types = FALSE) %>% mutate(sample = "S03", delivery = "mRNA")

# Load RNP data
abe1 <- read_tsv(abe1_file, show_col_types = FALSE) %>% mutate(sample = "ABE1", delivery = "RNP")
abe2 <- read_tsv(abe2_file, show_col_types = FALSE) %>% mutate(sample = "ABE2", delivery = "RNP")
abe3 <- read_tsv(abe3_file, show_col_types = FALSE) %>% mutate(sample = "ABE3", delivery = "RNP")

cat("  S01 (mRNA ABE):", nrow(s01), "sites\n")
cat("  S02 (mRNA ABE):", nrow(s02), "sites\n")
cat("  S03 (mRNA ABE):", nrow(s03), "sites\n")
cat("  ABE1 (RNP):", nrow(abe1), "sites\n")
cat("  ABE2 (RNP):", nrow(abe2), "sites\n")
cat("  ABE3 (RNP):", nrow(abe3), "sites\n\n")

# Combine all data (ABE samples only)
all_data <- bind_rows(s01, s02, s03, abe1, abe2, abe3)

cat("Total combined rows:", nrow(all_data), "\n")
cat("Total unique positions:", n_distinct(paste(all_data$chr, all_data$POS, all_data$REF, all_data$ALT)), "\n\n")

# ========== CREATE POSITION KEY ==========
all_data <- all_data %>%
  mutate(pos_key = paste(chr, POS, REF, ALT, sep = ":"))

all_sites <- unique(all_data$pos_key)
cat("Unique genomic positions:", length(all_sites), "\n\n")

# ========== FILTERING THRESHOLDS ==========
MIN_DEPTH <- 10
MAX_CTRL_EDIT_FREQ <- 0.01
MIN_ABE_EDIT_FREQ <- 0.01

# ========== STAGE 1: ADAR FILTERING ==========
cat("============================================================\n")
cat("STAGE 1: Remove ADAR Sites\n")
cat("============================================================\n\n")

# Remove sites where CTRL shows >1% editing
adar_sites <- all_data %>%
  mutate(
    ctrl_editing = case_when(
      ALT == "A" ~ (A.control / total_reads.control) * 100,
      ALT == "C" ~ (C.control / total_reads.control) * 100,
      ALT == "G" ~ (G.control / total_reads.control) * 100,
      ALT == "T" ~ (T.control / total_reads.control) * 100,
      TRUE ~ 0
    )
  ) %>%
  filter(ctrl_editing > MAX_CTRL_EDIT_FREQ) %>%
  pull(pos_key) %>%
  unique()

cat("ADAR sites (>1% in controls):", length(adar_sites), "\n")

# Remove ADAR sites
all_data_filtered <- all_data %>%
  mutate(pos_key = paste(chr, POS, REF, ALT, sep = ":")) %>%
  filter(!(pos_key %in% adar_sites))

cat("Sites after removing ADAR:", n_distinct(all_data_filtered$pos_key), "\n\n")

# ========== STAGE 2: COVERAGE FILTERING ==========
cat("============================================================\n")
cat("STAGE 2: Coverage Filtering\n")
cat("============================================================\n\n")

# Get coverage info for each site across all 6 ABE samples
coverage_stats <- all_data_filtered %>%
  filter(snv %in% c("A_G", "T_C")) %>%
  group_by(pos_key, chr, POS, REF, ALT, snv) %>%
  summarise(
    n_samples_total = n_distinct(sample),
    n_samples_with_10x = sum(total_reads.treated >= MIN_DEPTH, na.rm = TRUE),
    .groups = "drop"
  )

cat("Sites with A>G or T>C after ADAR filter:", nrow(coverage_stats), "\n\n")

# ========== 2/6 FILTERING ==========
cat("THRESHOLD: 2/6 samples with ≥10x coverage\n")
cat("-------------------------------------------\n")

sites_2of6 <- coverage_stats %>%
  filter(n_samples_with_10x >= 2) %>%
  pull(pos_key)

cat("Sites passing 2/6 threshold:", length(sites_2of6), "\n\n")

# Get editing data for 2/6 sites
data_2of6 <- all_data_filtered %>%
  filter(pos_key %in% sites_2of6, snv %in% c("A_G", "T_C")) %>%
  mutate(
    abe_editing = frac_alt_in_beOverexp * 100
  )

final_sites_2of6 <- data_2of6 %>%
  group_by(pos_key, chr, POS, REF, ALT, snv) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_abe_edit = mean(abe_editing, na.rm = TRUE),
    max_abe_edit = max(abe_editing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_abe_edit >= 1) %>%
  arrange(desc(max_abe_edit))

cat("Sites with ≥1% ABE editing:", nrow(final_sites_2of6), "\n\n")

# ========== 3/6 FILTERING ==========
cat("============================================================\n")
cat("THRESHOLD: 3/6 samples with ≥10x coverage\n")
cat("-------------------------------------------\n")

sites_3of6 <- coverage_stats %>%
  filter(n_samples_with_10x >= 3) %>%
  pull(pos_key)

cat("Sites passing 3/6 threshold:", length(sites_3of6), "\n\n")

# Get editing data for 3/6 sites
data_3of6 <- all_data_filtered %>%
  filter(pos_key %in% sites_3of6, snv %in% c("A_G", "T_C")) %>%
  mutate(
    abe_editing = frac_alt_in_beOverexp * 100
  )

final_sites_3of6 <- data_3of6 %>%
  group_by(pos_key, chr, POS, REF, ALT, snv) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_abe_edit = mean(abe_editing, na.rm = TRUE),
    max_abe_edit = max(abe_editing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_abe_edit >= 1) %>%
  arrange(desc(max_abe_edit))

cat("Sites with ≥1% ABE editing:", nrow(final_sites_3of6), "\n\n")

# ========== 4/6 FILTERING ==========
cat("============================================================\n")
cat("THRESHOLD: 4/6 samples with ≥10x coverage\n")
cat("-------------------------------------------\n")

sites_4of6 <- coverage_stats %>%
  filter(n_samples_with_10x >= 4) %>%
  pull(pos_key)

cat("Sites passing 4/6 threshold:", length(sites_4of6), "\n\n")

# Get editing data for 4/6 sites
data_4of6 <- all_data_filtered %>%
  filter(pos_key %in% sites_4of6, snv %in% c("A_G", "T_C")) %>%
  mutate(
    abe_editing = frac_alt_in_beOverexp * 100
  )

final_sites_4of6 <- data_4of6 %>%
  group_by(pos_key, chr, POS, REF, ALT, snv) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_abe_edit = mean(abe_editing, na.rm = TRUE),
    max_abe_edit = max(abe_editing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_abe_edit >= 1) %>%
  arrange(desc(max_abe_edit))

cat("Sites with ≥1% ABE editing:", nrow(final_sites_4of6), "\n\n")

# ========== WRITE OUTPUT FILES ==========
cat("============================================================\n")
cat("Writing output files...\n")
cat("============================================================\n\n")

write_tsv(final_sites_2of6, "combined_RNP_mRNA_2of6_offtargets.tsv")
cat("✓ combined_RNP_mRNA_2of6_offtargets.tsv (", nrow(final_sites_2of6), "sites)\n")

write_tsv(final_sites_3of6, "combined_RNP_mRNA_3of6_offtargets.tsv")
cat("✓ combined_RNP_mRNA_3of6_offtargets.tsv (", nrow(final_sites_3of6), "sites)\n")

write_tsv(final_sites_4of6, "combined_RNP_mRNA_4of6_offtargets.tsv")
cat("✓ combined_RNP_mRNA_4of6_offtargets.tsv (", nrow(final_sites_4of6), "sites)\n\n")

# ========== SUMMARY STATISTICS ==========
cat("============================================================\n")
cat("SUMMARY STATISTICS\n")
cat("============================================================\n\n")

summary_2of6 <- final_sites_2of6 %>%
  group_by(snv) %>%
  summarise(
    n_sites = n(),
    mean_edit = mean(max_abe_edit),
    median_edit = median(max_abe_edit),
    .groups = "drop"
  ) %>%
  mutate(threshold = "2/6")

summary_3of6 <- final_sites_3of6 %>%
  group_by(snv) %>%
  summarise(
    n_sites = n(),
    mean_edit = mean(max_abe_edit),
    median_edit = median(max_abe_edit),
    .groups = "drop"
  ) %>%
  mutate(threshold = "3/6")

summary_4of6 <- final_sites_4of6 %>%
  group_by(snv) %>%
  summarise(
    n_sites = n(),
    mean_edit = mean(max_abe_edit),
    median_edit = median(max_abe_edit),
    .groups = "drop"
  ) %>%
  mutate(threshold = "4/6")

summary_all <- bind_rows(summary_2of6, summary_3of6, summary_4of6)

write_tsv(summary_all, "combined_RNP_mRNA_summary_statistics.tsv")
cat("✓ combined_RNP_mRNA_summary_statistics.tsv\n\n")

# ========== FINAL SUMMARY ==========
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

print(summary_all)

cat("\nFull breakdown by type and threshold:\n")
cat("2/6:\n")
print(final_sites_2of6 %>% group_by(snv) %>% tally())
cat("\n3/6:\n")
print(final_sites_3of6 %>% group_by(snv) %>% tally())
cat("\n4/6:\n")
print(final_sites_4of6 %>% group_by(snv) %>% tally())

cat("\n============================================================\n")
cat("Output files in: /pub/madelk1/rnaseq/ABE8.8mRNA/\n")
cat("- combined_RNP_mRNA_2of6_offtargets.tsv\n")
cat("- combined_RNP_mRNA_3of6_offtargets.tsv\n")
cat("- combined_RNP_mRNA_4of6_offtargets.tsv\n")
cat("- combined_RNP_mRNA_summary_statistics.tsv\n")
cat("============================================================\n")
