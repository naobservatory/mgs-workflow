#!/usr/bin/env Rscript

library(jsonlite)
library(optparse)
library(tidyverse)

# Set arguments
option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL,
              help="Path to multiqc data directory."),
  make_option(c("-s", "--stage"), type="character", default=NULL,
              help="Stage descriptor."),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Path to output directory.")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Set input and output paths
multiqc_json_path <- file.path(opt$input_dir, "multiqc_data.json")
fastqc_tsv_path <- file.path(opt$input_dir, "multiqc_fastqc.txt")
out_path_basic <- file.path(opt$output_dir, paste0(opt$stage, "_qc_basic_stats.tsv"))
out_path_adapters <- file.path(opt$output_dir, paste0(opt$stage, "_qc_adapter_stats.tsv"))
out_path_quality_base <- file.path(opt$output_dir, paste0(opt$stage, "_qc_quality_base_stats.tsv"))
out_path_quality_sequence <- file.path(opt$output_dir, paste0(opt$stage, "_qc_quality_sequence_stats.tsv"))

#=====================#
# AUXILIARY FUNCTIONS #
#=====================#

process_n_bases <- function(n_bases_vec){
  # Function for extracting approximate base-count information from FASTQC TSV
  val = n_bases_vec %>% str_split(" ") %>% sapply(first) %>% as.numeric
  unit = n_bases_vec %>% str_split(" ") %>% sapply(last)
  # Adjust val based on unit
  val_out = ifelse(unit == "Gbp", val * 10^9, val) # TODO: Add other units as they come up
  val_out = ifelse(unit == "Mbp", val_out * 10^6, val_out)
  return(val_out)
}

basic_info_fastqc <- function(fastqc_tsv, multiqc_json){
  # Read in basic stats from multiqc JSON
  stats_json <- multiqc_json$report_general_stats_data
  tab_json <- lapply(names(stats_json), 
                     function(x) stats_json[[x]] %>% mutate(sample=x)) %>% bind_rows() %>%
    mutate(read_pair = sapply(str_split(sample, "_"), last),
           sample = sapply(str_split(sample, "_"), first)) %>%
    group_by(sample) %>%
    summarize(percent_gc = mean(percent_gc),
              mean_seq_len = mean(avg_sequence_length),
              n_read_pairs = total_sequences[1],
              percent_duplicates = mean(percent_duplicates))
  # Read in basic stats from fastqc TSV
  tab_tsv <- fastqc_tsv %>%
    mutate(read_pair = sapply(str_split(Sample, "_"), last),
           sample = sapply(str_split(Sample, "_"), first),
           n_bases_approx = process_n_bases(`Total Bases`)) %>%
    select(sample, read_pair, n_bases_approx, per_base_sequence_quality:adapter_content) %>%
    group_by(sample) %>% summarize_all(function(x) paste(x, collapse="/")) %>%
    select(-read_pair) %>%
    mutate(n_bases_approx = n_bases_approx %>% str_split("/") %>% sapply(as.numeric) %>% colSums())
  # Combine
  tab <- tab_json %>% inner_join(tab_tsv, by="sample")
}

extract_adapter_data_single <- function(adapter_dataset){
  # Convert a single JSON adapter dataset into a tibble
  data <- lapply(1:length(adapter_dataset$name), function(n)
    adapter_dataset$data[[n]] %>% as.data.frame %>% 
      mutate(sample=adapter_dataset$name[n])) %>%
    bind_rows() %>% as_tibble %>%
    rename(position=V1, pc_adapters=V2) %>%
    separate_wider_delim("sample", " - ", names=c("sample", "adapter")) %>%
    mutate(read_pair = sapply(str_split(sample, "_"), last),
           sample = sapply(str_split(sample, "_"), first))
  return(data)
}


extract_adapter_data <- function(multiqc_json){
  # Extract adapter data from multiqc JSON
  datasets <- multiqc_json$report_plot_data$fastqc_adapter_content_plot$datasets
  data_out <- lapply(datasets, extract_adapter_data_single) %>% bind_rows()
  return(data_out)
}

extract_per_base_quality_single <- function(per_base_quality_dataset){
  # Convert a single JSON per-base-quality dataset into a tibble
  data <- lapply(1:length(per_base_quality_dataset$name), function(n)
    per_base_quality_dataset$data[[n]] %>% as.data.frame %>% 
      mutate(sample=per_base_quality_dataset$name[n])) %>%
    bind_rows() %>% as_tibble %>%
    rename(position=V1, mean_phred_score=V2) %>%
    mutate(read_pair = sapply(str_split(sample, "_"), last),
           sample = sapply(str_split(sample, "_"), first))
  return(data)
}

extract_per_base_quality <- function(multiqc_json){
  # Extract per-base sequence quality data from multiqc JSON
  datasets <- multiqc_json$report_plot_data$fastqc_per_base_sequence_quality_plot$datasets
  data_out <- lapply(datasets, extract_per_base_quality_single) %>% bind_rows()
  return(data_out)
}

extract_per_sequence_quality_single <- function(per_sequence_quality_dataset){
  # Convert a single JSON per-sequence-quality dataset into a tibble
  data <- lapply(1:length(per_sequence_quality_dataset$name), function(n)
    per_sequence_quality_dataset$data[[n]] %>% as.data.frame %>% 
      mutate(sample=per_sequence_quality_dataset$name[n])) %>%
    bind_rows() %>% as_tibble %>%
    rename(mean_phred_score=V1, n_sequences=V2) %>%
    mutate(read_pair = sapply(str_split(sample, "_"), last),
           sample = sapply(str_split(sample, "_"), first))
  return(data)
}

extract_per_sequence_quality <- function(multiqc_json){
  # Extract per-base sequence quality data from multiqc JSON
  datasets <- multiqc_json$report_plot_data$fastqc_per_sequence_quality_scores_plot$datasets
  data_out <- lapply(datasets, extract_per_sequence_quality_single) %>% bind_rows()
  return(data_out)
}

#============#
# RUN SCRIPT #
#============#

# Import data
multiqc_json <- fromJSON(multiqc_json_path)
fastqc_tsv <- readr::read_tsv(fastqc_tsv_path, show_col_types = FALSE)

# Process
basic_info <- basic_info_fastqc(fastqc_tsv, multiqc_json) %>% mutate(stage = opt$stage)
adapters <- extract_adapter_data(multiqc_json) %>% mutate(stage = opt$stage)
per_base_quality <- extract_per_base_quality(multiqc_json) %>% mutate(stage = opt$stage)
per_sequence_quality <- extract_per_sequence_quality(multiqc_json) %>% mutate(stage = opt$stage)

# Write tables
write_tsv(basic_info, out_path_basic)
write_tsv(adapters, out_path_adapters)
write_tsv(per_base_quality, out_path_quality_base)
write_tsv(per_sequence_quality, out_path_quality_sequence)
