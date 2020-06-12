#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-i", "--input_tsv"), type="character", default=NULL,
              help="REQUIRED: TSV of variants to filter [default = %default]"),
  make_option(c("-o", "--output_tsv"), type="character", default=NULL,
              help="REQUIRED: filtered TSV of variants [default = %default]"),
  make_option(c("-p", "--input_blacklist"), type="character", default=NULL,
              help="provide a blacklist TSV to use instead of making a new one", metavar="character"),
  make_option(c("-b", "--output_blacklist"), type="character", default=NULL,
              help="output the newly generated or newly combined blacklist as a TSV, otherwise not outputted", metavar="character"),
  make_option(c("-c", "--combine_blacklists"), action="store_true", default=FALSE,
              help="generate new blacklist and combine with previous blacklist provided (must have same column fields) [default = %default]"),
  make_option(c("-f", "--vaf_max"), type="integer", default=0.1,
              help="variants less frequent than this fraction will not be considered for addition to the new blacklist [default = %default]"),
  make_option(c("-d", "--dp_min"), type="integer", default=5, metavar="number",
              help="variants with depth less than this will not be considered for addition to the new blacklist [default = %default]"),
  make_option(c("-m", "--mq_min"), type="integer", default=30, metavar="number",
              help="variants with mapping quality less than this will not be considered for addition to the new blacklist [default = %default]")
);

opt <- parse_args(OptionParser(option_list=option_list, description="Generates a blacklist of variants then removes those variants"))

# opt <- list("vaf_max" = 0.1, "dp_min" = 10, "mq_min" = 30)

if (is.null(opt$input_tsv) || is.null(opt$output_tsv)) {
  stop("You must specify --input_tsv and --output_tsv. Use -h or --help to get the help message.")
}

suppressPackageStartupMessages({
  library(tidyverse)
})

main <- function() {
  blacklist_columns <- c("CHROM", "POS", "REF", "ALT", "MQ", "DP", "cohort_VAF")
  if (!is.null(opt$input_blacklist) & !opt$combine_blacklists) {
    cat("Input blacklist provided. The program will use the provided blacklist for filtering and not generate a new blacklist\n")
    provided_blacklist_df <- read_tsv(opt$input_blacklist)
    filter_with_blacklist(provided_blacklist_df)
  } else if (!is.null(opt$input_blacklist) & opt$combine_blacklists) {
    cat("Input blacklist provided and instructed to generate new blacklist. The program will combine the new blacklist with the provided blacklist\n")
    provided_blacklist_df <- read_tsv(opt$input_blacklist)
    new_blacklist_df <- generate_blacklist(blacklist_columns)
    combined_blacklist_df <- rbind(provided_blacklist_df[,blacklist_columns], new_blacklist_df[,blacklist_columns])
    combined_blacklist_df <- combined_blacklist_df %>%
      distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)
    filter_with_blacklist(combined_blacklist_df)
    if (!is.null(opt$output_blacklist)) {
      write_tsv(combined_blacklist_df, opt$output_blacklist)
      cat("New blacklist outputted\n")
    }
  } else if (is.null(opt$input_blacklist)) {
    cat("Instructed to generate new blacklist. The program will generate n new blacklist\n")
    new_blacklist_df <- generate_blacklist(blacklist_columns)
    filter_with_blacklist(new_blacklist_df)
    if (!is.null(opt$output_blacklist)) {
      write_tsv(new_blacklist_df, opt$output_blacklist)
      cat("New blacklist outputted\n")
    }
  }
}

generate_blacklist <- function(columns_to_keep) {
  unfiltered_variants_df <- read_tsv(opt$input_tsv)
  total_samples <- length(unique(unfiltered_variants_df$VCF_SAMPLE_ID))
  
  new_blacklist_df <- unfiltered_variants_df %>%
    group_by(CHROM, POS, REF, ALT) %>%
      mutate(cohort_VAF = n()/total_samples) %>%
    ungroup() %>%
    filter(cohort_VAF >= opt$vaf_max & DP_1 >= opt$dp_min & MQ >= opt$mq_min) %>%
    distinct(CHROM, POS, REF, ALT, .keep_all = TRUE) %>%
    select(all_of(columns_to_keep))
  return(new_blacklist_df)
}

filter_with_blacklist <- function(blacklist_df) {
  unfiltered_variants_df <- read_tsv(opt$input_tsv)
  
  filtered_variants_df <- anti_join(unfiltered_variants_df, blacklist_df)
  cat("Variants filtered with blacklist\n")
  
  write_tsv(filtered_variants_df, opt$output_tsv)
  cat("Filtered variants outputted\n")
}

main()
