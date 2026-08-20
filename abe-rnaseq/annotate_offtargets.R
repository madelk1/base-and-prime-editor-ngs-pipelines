#!/usr/bin/env Rscript

# Annotate off-target sites with gene/feature information
# Uses bedtools intersect with GENCODE GTF

library(readr)
library(dplyr)
library(stringr)

cat("========================================\n")
cat("Annotating Off-Target Sites\n")
cat("Using bedtools + GENCODE v44\n")
cat("========================================\n\n")

# Paths
output_dir <- "/pub/madelk1/rnaseq/ABE8.8mRNA"
gtf_file <- "/pub/madelk1/rnaseq/refs/hg38/gencode.v44.annotation.gtf"

# Thresholds to process
thresholds <- c("2of6", "3of6", "4of6")

# ========== PROCESS EACH THRESHOLD ==========
for (thresh in thresholds) {
  cat("============================================================\n")
  cat("Processing", thresh, "threshold\n")
  cat("============================================================\n\n")
  
  # Load sites
  input_file <- file.path(output_dir, paste0("combined_RNP_mRNA_", thresh, "_offtargets.tsv"))
  output_file <- file.path(output_dir, paste0("combined_RNP_mRNA_", thresh, "_offtargets_annotated.tsv"))
  bed_file <- file.path(output_dir, paste0("combined_RNP_mRNA_", thresh, "_temp.bed"))
  intersect_file <- file.path(output_dir, paste0("combined_RNP_mRNA_", thresh, "_intersect.txt"))
  
  if (!file.exists(input_file)) {
    cat("WARNING: Cannot find", input_file, "\n\n")
    next
  }
  
  sites <- read_tsv(input_file, show_col_types = FALSE)
  cat("Loaded:", nrow(sites), "sites\n")
  
  # Convert to BED format (0-based, chr start end)
  bed_data <- sites %>%
    mutate(
      start = POS - 1,  # Convert to 0-based
      end = POS,
      name = paste(chr, POS, REF, ALT, sep = ":")
    ) %>%
    select(chr, start, end, name)
  
  # Write BED file
  write_tsv(bed_data, bed_file, col_names = FALSE)
  cat("Created BED file with", nrow(bed_data), "sites\n")
  
  # Run bedtools intersect
  cat("Running bedtools intersect...\n")
  system(paste0("module load bedtools2/2.30.0 && bedtools intersect -a ", bed_file, " -b ", gtf_file, " -wao > ", intersect_file))
  
  # Load intersect results
  intersect_data <- read_tsv(intersect_file, col_names = FALSE, show_col_types = FALSE)
  colnames(intersect_data) <- c("chr", "start", "end", "name", "gtf_chr", "gtf_source", "gtf_feature", "gtf_start", "gtf_end", "score", "strand", "frame", "gtf_attr", "overlap")
  
  cat("Intersect complete:", nrow(intersect_data), "rows\n")
  
  # Parse GTF attributes to extract gene_name and determine feature type
  parse_feature <- function(feature, attr) {
    # Determine detailed feature type based on GTF feature and attributes
    if (is.na(feature) || feature == ".") {
      return("intergenic")
    }
    
    feature_type <- switch(feature,
      "exon" = "exon",
      "CDS" = "protein_coding_exon",
      "UTR" = "UTR",
      "Selenocysteine" = "Selenocysteine",
      "start_codon" = "start_codon",
      "stop_codon" = "stop_codon",
      "transcript" = "intronic",  # Between exons
      feature  # Default: use as-is
    )
    
    return(feature_type)
  }
  
  extract_gene_name <- function(attr) {
    if (is.na(attr) || attr == ".") {
      return(NA_character_)
    }
    
    # Extract gene_name from GTF attribute string
    match <- str_extract(attr, 'gene_name "([^"]+)"')
    if (is.na(match)) {
      return(NA_character_)
    }
    gene_name <- str_extract(match, '[^"]+$')
    return(gene_name)
  }
  
  # Apply parsing
  intersect_data <- intersect_data %>%
    mutate(
      gene_name = sapply(gtf_attr, extract_gene_name),
      feature_type = mapply(parse_feature, gtf_feature, gtf_attr),
      has_annotation = (overlap > 0)
    )
  
  # For each site, keep the best annotation (prioritize exons > UTR > intronic)
  feature_priority <- c("protein_coding_exon" = 1, "exon" = 2, "CDS" = 3, "UTR" = 4, 
                       "start_codon" = 5, "stop_codon" = 6, "intronic" = 7, "intergenic" = 8)
  
  annotation <- intersect_data %>%
    filter(has_annotation) %>%
    mutate(priority = feature_priority[feature_type]) %>%
    arrange(name, priority) %>%
    group_by(name) %>%
    slice(1) %>%
    ungroup() %>%
    select(name, gene_name, feature_type)
  
  # For sites with no annotation, mark as intergenic
  all_names <- unique(bed_data$name)
  annotated_names <- unique(annotation$name)
  unannotated_names <- setdiff(all_names, annotated_names)
  
  if (length(unannotated_names) > 0) {
    unannotated <- data.frame(
      name = unannotated_names,
      gene_name = NA_character_,
      feature_type = "intergenic"
    )
    annotation <- bind_rows(annotation, unannotated)
  }
  
  # Join back to original sites
  output_data <- sites %>%
    mutate(name = paste(chr, POS, REF, ALT, sep = ":")) %>%
    left_join(annotation, by = "name") %>%
    select(-name) %>%
    arrange(desc(max_abe_edit))
  
  # Write annotated output
  write_tsv(output_data, output_file)
  cat("Wrote:", output_file, "\n")
  cat("  Columns:", paste(colnames(output_data), collapse = ", "), "\n\n")
  
  # Cleanup temp files
  system(paste0("rm -f ", bed_file, " ", intersect_file))
}

# ========== SUMMARY STATS ==========
cat("============================================================\n")
cat("SUMMARY STATISTICS\n")
cat("============================================================\n\n")

for (thresh in thresholds) {
  output_file <- file.path(output_dir, paste0("combined_RNP_mRNA_", thresh, "_offtargets_annotated.tsv"))
  
  if (file.exists(output_file)) {
    data <- read_tsv(output_file, show_col_types = FALSE)
    
    cat(thresh, "threshold:\n")
    cat("  Total sites:", nrow(data), "\n")
    cat("  Feature type breakdown:\n")
    
    feature_counts <- data %>%
      group_by(feature_type) %>%
      summarise(n = n(), mean_edit = mean(max_abe_edit), .groups = "drop") %>%
      arrange(desc(n))
    
    for (i in 1:nrow(feature_counts)) {
      cat("    ", feature_counts$feature_type[i], ":", feature_counts$n[i], 
          "sites (mean edit:", sprintf("%.1f%%", feature_counts$mean_edit[i]*100), ")\n")
    }
    
    cat("  Top genes:\n")
    top_genes <- data %>%
      filter(!is.na(gene_name)) %>%
      group_by(gene_name) %>%
      summarise(n = n(), max_edit = max(max_abe_edit), .groups = "drop") %>%
      arrange(desc(n)) %>%
      head(5)
    
    for (i in 1:nrow(top_genes)) {
      cat("    ", top_genes$gene_name[i], ":", top_genes$n[i], 
          "sites (max edit:", sprintf("%.1f%%", top_genes$max_edit[i]*100), ")\n")
    }
    cat("\n")
  }
}

cat("============================================================\n")
cat("Annotation complete!\n")
cat("============================================================\n")
EOF
