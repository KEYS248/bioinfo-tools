#!/usr/bin/env Rscript
library(optparse)
option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="REQUIRED: prefix for outputted files", metavar="character"),
  make_option(c("-b", "--kb_bin_width"), type="integer", default=500,
              help="bin width in kilobases [default = %default]", metavar="number"),
  make_option(c("-c", "--controls"), type="character", default=NULL,
              help="REQUIRED: path to directory with control samples", metavar="character"),
	make_option(c("--control_regex"), type="character", default=".*.bam",
              help="regex pattern for control samples [default = %default]", metavar="character"),
  make_option(c("--control_suffix"), type="character", default=".bam",
              help="trim this suffix from control sample names in the plots [default = %default]", metavar="character"),
  make_option(c("-e", "--experimentals"), type="character", default=NULL,
              help="REQUIRED: path to directory with experimental samples", metavar="character"),
	make_option(c("--experimental_regex"), type="character", default=".*.bam",
              help="regex pattern for experimental samples [default = %default]", metavar="character"),
  make_option(c("--experimental_suffix"), type="character", default=".bam",
              help="trim this suffix from experimental sample names in the plots [default = %default]", metavar="character")
); 
 
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$project) || is.null(opt$controls) || is.null(opt$experimentals)){
  stop("You must specify --project, --controls, and --experimentals. Use -h or --help to get the help message.")
}

suppressPackageStartupMessages({
  library(tidyverse) # Used for general purpose data manipulation
  library(SCOPE)
  library(WGSmapp) # Required for SCOPE
  library(BSgenome.Hsapiens.UCSC.hg38) # Required for SCOPE
  library(reshape2) # Used to reformatting data for CNV plot
  library(ggdendro) # Used for plotting cluster tree plot
  library(gridExtra) # Used for combinding tree and CNV plots for export
})

# Variables used in plotting the CNV results
ggplot_theme <- theme(axis.line.y = element_line(size=.1,color = "black"), axis.line.x = element_line(size=.1,color = "black"),
                     axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
                     panel.grid.major = element_line(color = "black"), panel.background = element_rect(fill="white"),
                     panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), 
                     plot.title = element_text(size=15), legend.text=element_text(size=8))
chr_names <- c(paste0("chr", paste(seq(1, 22))), "chrX", "chrY", "chrM")
chr_starts <- c(1, 248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912, 
               1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562, 
               2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417, 3088269832)
chr_middles <- c(124478211, 370053186, 590297730, 784552787, 970429194, 1146601313, 1311677289, 1463919594, 
                 1605686270, 1741782340, 1876224362, 2010405327, 2134225146, 2244929169, 2349446622, 2445611389, 
                 2532409282, 2614224645, 2683720096, 2745250987, 2800828062, 2849592288, 2953021969, 3059656124, 3088278116)


# Functions ------------------------------------------------------
# Modified SCOPE get_bam_bed() function to keep chrX, while the default SCOPE function does not
get_bam_bed_n24chr <- function(bamdir, sampname, hgref = "hg19", resolution = 500) {
  if(!hgref %in% c("hg19", "hg38")){
    stop("Reference genome should be either hg19 or hg38. ")
  }
  if(hgref == "hg19") {
    genome <- BSgenome.Hsapiens.UCSC.hg19
  }else if(hgref == "hg38") {
    genome <- BSgenome.Hsapiens.UCSC.hg38
  }
  if(resolution <= 0){
    stop("Invalid fixed bin length. ")
  }
  bins <- tileGenome(seqinfo(genome), 
                     tilewidth = resolution * 1000, 
                     cut.last.tile.in.chrom = TRUE)
  ref <- bins[which(as.character(seqnames(bins)) %in% paste0("chr", 
                                                             c(seq_len(22), 
                                                               "X")))]
  if (!any(grepl("chr", seqlevels(ref)))) {
    seqlevels(ref) <- paste(c(seq_len(22), "X"), sep = "")
    ref <- sort(ref)
  } else {
    seqlevels(ref) <- paste("chr", c(seq_len(22), "X"), sep = "")
    ref <- sort(ref)
  }
  list(bamdir = bamdir, sampname = sampname, ref = ref)
}

# Function for completing the common analysis and plotting of CNV data for each of the three normalization types
plot_and_export_CNV <- function(Yhat.input = NULL, beta.hat.input = NULL, suffix, alpha.hat.input = NULL) {
  # This plot_EM_fit() function tends to cause problems when using small bin sizes
  # plot_EM_fit(Y_qc = Y_qc, gc_qc = gc_qc, norm_index = control_indices, T = 1:length(samples), ploidyInt = ploidy_data, 
  #             beta0 = beta.hat.input, filename = sprintf("%s.fig_scope_EM_%s.pdf", opt$project, suffix))
  if (is.null(alpha.hat.input)) {
    cat("Performing segment CBS\n")
    start <- Sys.time()
    tryCatch({
      chrs <- unique(as.character(seqnames(ref_qc)))
      segment_cs <- vector('list',length = length(chrs))
      names(segment_cs) <- chrs
      for (chri in chrs) {
        message('\n', chri, '\n')
        # With certain input data, this segment_CBScs() function will cause problems for 2nd pass normalizations (not 1st pass)
        segment_cs[[chri]] <- segment_CBScs(Y = Y_qc, Yhat = Yhat.input, sampname = colnames(Y_qc), 
                                            ref = ref_qc, chr = chri, mode = "integer")
      }
      iCN_mat <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
    }, error = function(err) {
      cat(sprintf("Error occured in %s normalization. Gently ending program.\n", suffix))
      write.table(qc_metric_qc, sprintf("%s.scope_sample_metrics.tsv", opt$project), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      cat("Exported sample data\n")
      cat("Done with SCOPE CNV analysis\n")
      stop()
    })
    cat("Completed segment CBS\n")
    print(Sys.time() - start)
  } else {
    cat("Not performing segment CBS\n")
    colnames(alpha.hat.input) <- samples
    iCN_mat <- alpha.hat.input
  }

  # Round CNV data to whole numbers between 0 and 7, and hierarchically reorder
  iCN_mat <- round(iCN_mat)
  iCN_mat[iCN_mat >= 7] <- 7
  iCN_mat[iCN_mat <= 0] <- 0
  cluster <- hclust(dist(t(iCN_mat)))
  iCN_mat <- iCN_mat[, cluster$order]

  # Save CNV data along with chromosome information
  CNV_df <- cbind(as.data.frame(ref_qc), as.data.frame(iCN_mat))
  write.table(CNV_df, sprintf("%s.scope_CNV_%s.tsv", opt$project, suffix), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat(sprintf("Saved CNVs for %s\n", suffix))

  # Combine chromosomal data with CNV data and reshape for ggplotting
  CNV_df <- CNV_df %>%
    group_by(seqnames, start) %>%
    mutate(abs_start = chr_starts[which(chr_names == seqnames)] + start) %>%
    ungroup() %>%
    melt(id.vars = c("seqnames", "start", "end", "width", "strand", "gc", "mapp", "abs_start"), variable.name = "sample", value.name = "CNV")
  # Plot dendrogram (tree) of hierarchical clustering
  dendro_df <- dendro_data(cluster, type = "rectangle")
  tree_plot <- ggdendrogram(dendro_df, rotate = TRUE, size = 2)
  # Plot CNV heatmap of all samples
  n_chrs <- length(unique(CNV_df$seqnames))
  cnv_colors <- c("#2166AC", "#92C5DE", "#FDFDFD", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
  cnv_plot <- ggplot(CNV_df, aes(abs_start, sample, fill=CNV)) + geom_tile() + xlim(0, chr_starts[(n_chrs+1)]) +
    labs(x = "Chromosome", y = "Sample", title = sprintf("%s CNVs - %skb bins - %s", opt$project, opt$kb_bin_width, suffix)) +
    scale_fill_gradientn(colors = cnv_colors[(min(CNV_df$CNV)+1):(max(CNV_df$CNV)+1)]) +
    scale_x_continuous(labels = chr_names[1:n_chrs], breaks = chr_middles[1:n_chrs], expand = c(0, 0)) +
    geom_vline(xintercept = chr_starts[1:n_chrs]) + ggplot_theme
  # Save two plots together
  pdf(sprintf("%s.fig_scope_CNV_%s.pdf", opt$project, suffix), width = plot_width, height = plot_height)
  grid.arrange(cnv_plot, tree_plot, layout_matrix = matrix(c(rep(1, 9), 2), nrow = 1))
  dev.off()
  cat(sprintf("Plotted CNVs for %s\n", suffix))
  
  # Return hierarchical clustering data
  cluster <- as.data.frame(cluster$order)
  colnames(cluster) <- c(sprintf("%s_hierarchy", suffix))
  return(cluster)
}


# Get sample names and paths ------------------------------------------------------
if (substr(opt$controls, nchar(opt$controls), nchar(opt$controls)) != "/") {
  opt$controls <- paste0(opt$controls, "/")
}
if (substr(opt$experimentals, nchar(opt$experimentals), nchar(opt$experimentals)) != "/") {
  opt$experimentals <- paste0(opt$experimentals, "/")
}
control_samples <- str_remove(list.files(opt$controls, opt$control_regex), opt$control_suffix)
experimental_samples <- str_remove(list.files(opt$experimentals, opt$experimental_regex), opt$experimental_suffix)
samples <- c(control_samples, experimental_samples)
control_indices <- rep(1:length(control_samples))
full_sample_paths <- file.path(c(rep(opt$controls, length(control_samples)), rep(opt$experimentals, length(experimental_samples))),
                               c(list.files(opt$controls, opt$control_regex), list.files(opt$experimentals, opt$experimental_regex)))
cat(sprintf("%s project - %skb bin size\n", opt$project, opt$kb_bin_width))
sample_group_labels <- c(rep("control", length(control_samples)), rep("experimental", length(experimental_samples)))
cat(sprintf("%s samples: %s control, %s experimental\n", length(samples), length(control_samples), length(experimental_samples)))
print(samples)

plot_height <- (length(samples)/10)+7.5
plot_width <- (length(samples)/30)+20


# Get sample bins, mappability, and GC content ------------------------------------------------------
bam_bed_object <- get_bam_bed_n24chr(bamdir = full_sample_paths, sampname = samples, hgref = "hg38", resolution = opt$kb_bin_width)
ref_raw <- bam_bed_object$ref

mapp <- get_mapp(ref_raw, hgref = "hg38")
gc <- get_gc(ref_raw, hgref = "hg38")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
rm(mapp, gc)


# Compute read coverage ------------------------------------------------------
coverage_object <- get_coverage_scDNA(bam_bed_object, mapqthres = 40, seq = 'paired-end', hgref = "hg38")
Y_raw <- coverage_object$Y
qc_metric_raw <- get_samp_QC(bam_bed_object)


# QC filtering and gini coverage inequality index ------------------------------------------------------
qc_object <- perform_qc(Y_raw = Y_raw, sampname_raw = samples, ref_raw = ref_raw, 
                                 QCmetric_raw = qc_metric_raw, mapq20_thresh = 0)
ref_qc <- qc_object$ref
gc_qc <- get_gc(ref_qc, hgref = "hg38")
Y_qc <- qc_object$Y
qc_metric_qc <- as_tibble(qc_object$QCmetric, rownames = "sample")
rm(bam_bed_object, ref_raw, coverage_object, Y_raw, qc_metric_raw)

gini_qc <- get_gini(Y_qc)
gini_df <- data.frame(sample = samples, gini = gini_qc, group = sample_group_labels)
print(gini_df)
qc_metric_qc <- left_join(qc_metric_qc, gini_df, by = c("sample"))
rm(gini_qc, gini_df)


# Multi-pass read depth normalization and plotting ------------------------------------------------------
# First pass fast normalization without latent factors
cat("Performing first pass normalization\n")
start <- Sys.time()
normalized_noK <- normalize_codex2_ns_noK(Y_qc = Y_qc, gc_qc = gc_qc, norm_index = control_indices) # which(gini_qc <= 0.12)
cat("\nCompleted first pass normalization\n")
print(Sys.time() - start)
Yhat.noK <- normalized_noK$Yhat
beta.hat.noK <- normalized_noK$beta.hat
# fGC.hat.noK <- normalized_noK$fGC.hat
N.noK <- normalized_noK$N
qc_metric_qc <- cbind(qc_metric_qc, plot_and_export_CNV(Yhat.input = Yhat.noK, beta.hat.input = beta.hat.noK, suffix = "noK"))

ploidy_data <- initialize_ploidy(Y = Y_qc, Yhat = Yhat.noK, ref = ref_qc)
ploidy_df <- data.frame(sample = samples, ploidy = ploidy_data)
print(ploidy_df)
qc_metric_qc <- left_join(qc_metric_qc, ploidy_df, by = c("sample"))
cat("\nEstimated ploidies for samples\n")
rm(normalized_noK, Yhat.noK, ploidy_df)

# Second pass normalizations let you choose between AIC, BIC, and RSS model selection algorithms for plotting
# Second pass normalization with latent factors and BIC plotting
cat("Performing second pass individual normalization\n")
start <- Sys.time()
normalized_data <- normalize_scope(Y_qc = Y_qc, K = 1, T = 1:7, gc_qc = gc_qc, minCountQC = 20,
                                   beta0 = beta.hat.noK, ploidyInt = ploidy_data, norm_index = control_indices)
cat("\nCompleted second pass individual normalization\n")
print(Sys.time() - start)
Yhat.data <- normalized_data$Yhat[[which.max(normalized_data$BIC)]]
beta.hat.data <- normalized_data$beta.hat[[which.max(normalized_data$BIC)]]
alpha.hat.data <- normalized_data$alpha.hat[[which.max(normalized_data$BIC)]]
qc_metric_qc <- cbind(qc_metric_qc, plot_and_export_CNV(Yhat.input = Yhat.data, beta.hat.input = beta.hat.data, suffix = "BIC_separate"))
rm(normalized_data, Yhat.data, beta.hat.data)

# Second pass normalization with latent factors by group and BIC plotting
cat("Performing second pass group normalization\n")
start <- Sys.time()
normalized_group <- normalize_scope_group(Y_qc = Y_qc, K = 1, T = 1:7, gc_qc = gc_qc, minCountQC = 20,
                                          beta0 = beta.hat.noK, ploidyInt = ploidy_data, norm_index = control_indices, 
                                          groups = sample_group_labels)
cat("\nCompleted second pass group normalization\n")
print(Sys.time() - start)
Yhat.group <- normalized_group$Yhat[[which.max(normalized_group$BIC)]]
beta.hat.group <- normalized_group$beta.hat[[which.max(normalized_group$BIC)]]
alpha.hat.group <- normalized_group$alpha.hat[[which.max(normalized_group$BIC)]]
qc_metric_qc <- cbind(qc_metric_qc, plot_and_export_CNV(Yhat.input = Yhat.group, beta.hat.input = beta.hat.group, suffix = "BIC_group"))

write.table(qc_metric_qc, sprintf("%s.scope_sample_metrics.tsv", opt$project), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("Exported sample data\n")
cat("Done with SCOPE CNV analysis\n")
