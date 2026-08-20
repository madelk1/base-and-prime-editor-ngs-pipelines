#!/usr/bin/env Rscript

# GATK4 Three-Stage Pooled Filtering - RELAXED (2/3 requirement)
# Loads all three conditions: ABE, HiFi-ABE, Inlaid-ABE

library(readr)
library(dplyr)
library(tidyr)

# Working directory
setwd("/pub/madelk1/HiFiRNASeq/variantCalls")

# Control samples (same for all conditions)
CTRL_SAMPLES <- c("001_CTRL1", "002_CTRL2", "003_CTRL3")

# Filtering thresholds
MIN_DEPTH <- 10                # Minimum coverage
MAX_CTRL_EDIT_FREQ <- 0.01     # Controls must have ≤1% editing
MIN_ABE_EDIT_FREQ <- 0.01      # ABE must have ≥1% to count

# RELAXED: Require only 2/3 of samples
MIN_CTRL_WITH_COVERAGE <- 2    # At least 2 out of 3 controls
MIN_SAMPLES_WITH_COVERAGE <- 4 # At least 4 out of 6 total samples

cat("========================================\n")
cat("GATK4 Three-Stage Pooled Filtering\n")
cat("HiFiRNASeq: ABE, HiFi, Inlaid comparison\n")
cat("RELAXED: 2/3 requirement\n")
cat("========================================\n\n")

# ========== DEFINE CONDITIONS ==========
conditions <- list(
  ABE = list(name = "ABE8e", samples = c("004_ABE1", "005_ABE2", "006_ABE3")),
  HiFi = list(name = "HiFi-ABE", samples = c("007_HiFi1", "008_HiFi2", "009_HiFi3")),
  Inlaid = list(name = "Inlaid-ABE", samples = c("010_Inlaid1", "011_Inlaid2", "012_Inlaid3"))
)

# ========== PROCESS EACH CONDITION ==========
for (cond_name in names(conditions)) {
  
  cat("\n============================================================\n")
  cat("Processing condition:", conditions[[cond_name]]$name, "\n")
  cat("============================================================\n\n")
  
  treated_samples <- conditions[[cond_name]]$samples
  all_samples <- c(CTRL_SAMPLES, treated_samples)
  
  # ========== LOAD DATA ==========
  cat("Loading GATK4 data...\n")
  
  all_data_list <- list()
  
  for (i in 1:3) {
    treated_sample <- treated_samples[i]
    ctrl_sample <- CTRL_SAMPLES[i]
    
    # Extract numeric part for filename construction
    num_treated <- as.numeric(substr(treated_sample, 1, 3))
    filename <- paste0(treated_sample, ".vcf.gz.sorted.vcf_with_coverage.be_and_control.txt.all_snvs.stringent.txt")
    
    if (!file.exists(filename)) {
      cat("ERROR: Cannot find", filename, "\n")
      next
    }
    
    data <- read_tsv(filename, show_col_types = FALSE)
    cat("  Loaded", filename, ":", nrow(data), "sites\n")
    
    # Create two rows per site: one for control, one for treated
    ctrl_data <- data %>%
      mutate(
        sample = ctrl_sample,
        is_control = TRUE,
        depth = total_reads.control,
        alt_frac = case_when(
          ALT == "A" ~ A.control / total_reads.control,
          ALT == "C" ~ C.control / total_reads.control,
          ALT == "G" ~ G.control / total_reads.control,
          ALT == "T" ~ T.control / total_reads.control,
          TRUE ~ 0
        )
      ) %>%
      select(chr, POS, REF, ALT, snv, sample, is_control, depth, alt_frac)
    
    treated_data <- data %>%
      mutate(
        sample = treated_sample,
        is_control = FALSE,
        depth = total_reads.treated,
        alt_frac = frac_alt_in_beOverexp
      ) %>%
      select(chr, POS, REF, ALT, snv, sample, is_control, depth, alt_frac)
    
    all_data_list[[length(all_data_list) + 1]] <- ctrl_data
    all_data_list[[length(all_data_list) + 1]] <- treated_data
  }
  
  all_data <- bind_rows(all_data_list)
  
  cat("\nTotal rows loaded:", nrow(all_data), "\n")
  
  # ========== CREATE POSITION KEY ==========
  all_data <- all_data %>%
    mutate(pos_key = paste(chr, POS, REF, ALT, sep = "_"))
  
  # Get unique sites
  all_sites <- unique(all_data$pos_key)
  cat("Unique genomic positions:", length(all_sites), "\n\n")
  
  # ========== STAGE 1: ≥2/3 CONTROLS NEED ≥10x COVERAGE ==========
  cat("STAGE 1: Require ≥", MIN_CTRL_WITH_COVERAGE, "/3 controls to have ≥", MIN_DEPTH, "x coverage\n")
  
  stage1_sites <- all_data %>%
    filter(is_control) %>%
    group_by(pos_key) %>%
    summarise(
      n_ctrls_with_coverage = sum(depth >= MIN_DEPTH),
      .groups = "drop"
    ) %>%
    filter(n_ctrls_with_coverage >= MIN_CTRL_WITH_COVERAGE) %>%
    pull(pos_key)
  
  cat("Sites with ≥", MIN_CTRL_WITH_COVERAGE, "/3 controls ≥", MIN_DEPTH, "x:", length(stage1_sites), "\n\n")
  
  # ========== STAGE 2: ≥4/6 SAMPLES NEED ≥10x COVERAGE ==========
  cat("STAGE 2: Require ≥", MIN_SAMPLES_WITH_COVERAGE, "/6 samples to have ≥", MIN_DEPTH, "x coverage\n")
  
  stage2_sites <- all_data %>%
    filter(pos_key %in% stage1_sites) %>%
    group_by(pos_key) %>%
    summarise(
      n_samples_with_coverage = sum(depth >= MIN_DEPTH),
      .groups = "drop"
    ) %>%
    filter(n_samples_with_coverage >= MIN_SAMPLES_WITH_COVERAGE) %>%
    pull(pos_key)
  
  cat("Sites with ≥", MIN_SAMPLES_WITH_COVERAGE, "/6 samples ≥", MIN_DEPTH, "x:", length(stage2_sites), "\n\n")
  
  # ========== STAGE 3: REMOVE ADAR SITES ==========
  cat("STAGE 3: Remove ADAR sites (controls with >", MAX_CTRL_EDIT_FREQ*100, "% editing)\n")
  
  adar_sites <- all_data %>%
    filter(pos_key %in% stage2_sites,
           is_control,
           depth >= MIN_DEPTH,
           alt_frac > MAX_CTRL_EDIT_FREQ) %>%
    pull(pos_key) %>%
    unique()
  
  cat("Sites with >", MAX_CTRL_EDIT_FREQ*100, "% editing in controls (ADAR):", length(adar_sites), "\n")
  
  stage3_sites <- setdiff(stage2_sites, adar_sites)
  
  cat("Sites after removing ADAR:", length(stage3_sites), "\n\n")
  
  # ========== IDENTIFY OFF-TARGET EDITING SITES ==========
  cat("Identifying off-target editing sites\n")
  
  offtarget_data <- all_data %>%
    filter(pos_key %in% stage3_sites,
           snv %in% c("A_G", "T_C"))
  
  # Find sites with ≥1% editing in ANY treated sample
  abe_editing <- offtarget_data %>%
    filter(!is_control,
           alt_frac >= MIN_ABE_EDIT_FREQ) %>%
    pull(pos_key) %>%
    unique()
  
  cat("Sites with A>G or T>C editing ≥", MIN_ABE_EDIT_FREQ*100, "% in treated:", length(abe_editing), "\n\n")
  
  # ========== CREATE OUTPUT TABLE ==========
  cat("Creating output table...\n")
  
  output_sites <- offtarget_data %>%
    filter(pos_key %in% abe_editing) %>%
    group_by(pos_key, chr, POS, REF, ALT, snv) %>%
    summarise(
      n_samples_total = n(),
      n_samples_with_coverage = sum(depth >= MIN_DEPTH),
      mean_ctrl_depth = mean(depth[is_control], na.rm = TRUE),
      mean_ctrl_edit = mean(alt_frac[is_control], na.rm = TRUE),
      max_ctrl_edit = max(alt_frac[is_control], na.rm = TRUE),
      mean_treated_depth = mean(depth[!is_control], na.rm = TRUE),
      mean_treated_edit = mean(alt_frac[!is_control], na.rm = TRUE),
      max_treated_edit = max(alt_frac[!is_control], na.rm = TRUE),
      # Individual sample values
      CTRL1_edit = ifelse(any(sample == CTRL_SAMPLES[1]), alt_frac[sample == CTRL_SAMPLES[1]][1], NA),
      CTRL2_edit = ifelse(any(sample == CTRL_SAMPLES[2]), alt_frac[sample == CTRL_SAMPLES[2]][1], NA),
      CTRL3_edit = ifelse(any(sample == CTRL_SAMPLES[3]), alt_frac[sample == CTRL_SAMPLES[3]][1], NA),
      Treated1_edit = ifelse(any(sample == treated_samples[1]), alt_frac[sample == treated_samples[1]][1], NA),
      Treated2_edit = ifelse(any(sample == treated_samples[2]), alt_frac[sample == treated_samples[2]][1], NA),
      Treated3_edit = ifelse(any(sample == treated_samples[3]), alt_frac[sample == treated_samples[3]][1], NA),
      .groups = "drop"
    ) %>%
    arrange(desc(max_treated_edit))
  
  # ========== WRITE OUTPUT FILES ==========
  output_prefix <- paste0(cond_name, "_offtarget_sites_3stage_2of3")
  write_tsv(output_sites, paste0(output_prefix, ".tsv"))
  cat("Wrote:", paste0(output_prefix, ".tsv\n"))
  
  # Summary by edit type
  summary_by_type <- output_sites %>%
    group_by(snv) %>%
    summarise(
      n_sites = n(),
      mean_edit_freq = mean(max_treated_edit),
      median_edit_freq = median(max_treated_edit),
      .groups = "drop"
    )
  
  write_tsv(summary_by_type, paste0(output_prefix, "_summary.tsv"))
  cat("Wrote:", paste0(output_prefix, "_summary.tsv\n"))
  
  # Per-sample summary
  per_sample_summary <- offtarget_data %>%
    filter(pos_key %in% abe_editing) %>%
    group_by(sample) %>%
    summarise(
      condition = ifelse(sample %in% CTRL_SAMPLES, "CTRL", conditions[[cond_name]]$name),
      n_sites_with_editing = sum(alt_frac >= MIN_ABE_EDIT_FREQ, na.rm = TRUE),
      mean_editing = mean(alt_frac, na.rm = TRUE),
      median_editing = median(alt_frac, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_tsv(per_sample_summary, paste0(output_prefix, "_per_sample.tsv"))
  cat("Wrote:", paste0(output_prefix, "_per_sample.tsv\n\n"))
  
  # ========== FINAL SUMMARY ==========
  cat("FILTERING SUMMARY\n")
  cat("Starting unique sites:", length(all_sites), "\n")
  cat("After Stage 1:", length(stage1_sites), "\n")
  cat("After Stage 2:", length(stage2_sites), "\n")
  cat("After Stage 3:", length(stage3_sites), "\n")
  cat("A>G + T>C sites with editing:", length(abe_editing), "\n\n")
  
  cat("Edit type breakdown:\n")
  print(summary_by_type)
  cat("\n")
  
  cat("Per-sample summary:\n")
  print(per_sample_summary)
  cat("\n")
}

cat("============================================================\n")
cat("All conditions processed!\n")
cat("============================================================\n")
EOF
