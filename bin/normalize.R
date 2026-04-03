#!/usr/bin/env Rscript

library(optparse, quietly=TRUE)

option_list = list(
  make_option(c("-g", "--groups"), type="character", help="Name of groups file.", dest="groups"),
  make_option(c("-c", "--counts"), type="character", help="Name of count file", dest="counts"),
  make_option(c("-p", "--input-prefix"), type="character", default="", help="Input prefix", dest="input_prefix"),
  make_option(c("-e", "--min-counts-per-event"), type="integer", default=10, help="Filter genes with less than this count.", dest="min_counts_per_event"),
  make_option(c("-s", "--min-samples-per-event"), type="integer", default=1, help="Filter genes which do not have at least min-counts-per-event reads in this many samples.", dest="min_samples_per_event"),
  make_option(c("-m", "--min-counts-per-sample"), type="integer", default=1000, help="Filter samples which do not have at least min_counts_per_sample reads.", dest="min_counts_per_sample"),
  make_option(c("-x", "--excluded-events"), type="character", help="Space separated list of events to remove.", dest="excluded_events"),
  make_option(c("-t", "--transformation"), type="character", default="Default", help="Transformation to apply prior to PCA (does not apply to data sent into DESeq). Default will use the rlog transformation for up to 50 samples and vst if there are more than 50 samples.", dest='transformation')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(tidyr, quietly = TRUE))
suppressMessages(library(DESeq2, quietly = TRUE))

read_metadata = function(metadata_file) {
  metadata = read.delim(metadata_file, header=TRUE)
  if (length(metadata) < 2) {
    metadata = read.delim(metadata_file, sep=',', header=TRUE)
  }
  return(metadata)
}

read_counts = function(count_file, samples, excluded_events) {
  
  counts = read.delim(count_file, check.names=FALSE)
  if (length(counts) < 2) {
    counts = read.delim(count_file, check.names=FALSE, sep=',')
  }
  
  counts = counts %>%
    select(feature=1, all_of(samples)) %>%
    filter(!(feature %in% excluded_events)) %>%
    tibble::column_to_rownames('feature')
  counts = as.matrix(counts)
  mode(counts) = 'integer'
  return(counts)
}

check_samples = function(all_counts, total_count_min) {
  bad_samples = (data.frame(Counts=colSums(all_counts)) %>% tibble::rownames_to_column('Sample') %>% filter(Counts < total_count_min))$Sample
  return(bad_samples)
}

sample_filter = function(all_counts, bad_samples) {
  x = as.data.frame(all_counts) %>% select(-all_of(bad_samples))
  return(as.matrix(x))
}

filter_counts = function(input, min_counts_per_event, min_samples_per_event) {
  keep = rowSums(input >= min_counts_per_event) >= min_samples_per_event
  return(input[keep,])
}

run_pca = function(input, retx=TRUE, center=TRUE, scale=TRUE, transformation='Default', write_transform=FALSE) {
  if (transformation == 'None') {
    keep = subset(input, apply(input, 1, var, na.rm = TRUE) >  0)
    return(prcomp(t(keep), retx = retx, center = center, scale. = scale))
  } else if (transformation == 'vst' | (transformation == 'Default' & ncol(input) > 50)) {
    transformed = vst(input)
    if(write_transform) {
      write.table(as.data.frame(transformed) %>% tibble::rownames_to_column("feature"), file='outputs/transformed_counts_vst.tsv', quote=FALSE, sep='	', row.names=FALSE)
    }
    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
  } else {
    transformed = rlog(input)
    if(write_transform) {
      write.table(as.data.frame(transformed) %>% tibble::rownames_to_column("feature"), file='outputs/transformed_counts_rlog.tsv', quote=FALSE, sep='	', row.names=FALSE)
    }
    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
  }
}

if (is.null(opt$excluded_events)) {
  excluded_events = c()
} else {
  excluded_events = strsplit(opt$excluded_events, " ", fixed=TRUE)
}

all_metadata = read_metadata(paste0(opt$input_prefix, opt$groups))
all_counts = read_counts(paste0(opt$input_prefix, opt$counts), all_metadata$sample_name, excluded_events)

bad_samples = check_samples(all_counts, opt$min_counts_per_sample)
all_counts = sample_filter(all_counts, bad_samples)
all_counts_filtered = filter_counts(all_counts, opt$min_counts_per_event, opt$min_samples_per_event)

pca = run_pca(all_counts_filtered, transformation = opt$transformation, write_transform=TRUE)
norm_dds = DESeqDataSetFromMatrix(countData = all_counts_filtered, colData = all_metadata, design = ~1)	
norm = assay(normTransform(norm_dds))
write.table(as.data.frame(norm) %>% tibble::rownames_to_column("feature"), file="outputs/normalized_feature_counts.tsv", sep='\t', quote=FALSE, row.names = FALSE)
write.table(as.data.frame(t(scale(t(norm)))) %>% tibble::rownames_to_column("feature"), file="outputs/expression_z_scores.tsv", sep='\t', quote=FALSE, row.names=FALSE)