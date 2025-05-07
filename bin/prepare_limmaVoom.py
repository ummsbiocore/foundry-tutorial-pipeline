#!/usr/bin/env python3

import sys, os, argparse, subprocess
from collections import defaultdict
from textwrap import dedent
from multiprocessing import Pool

def main(args):

	if os.path.basename(args.groups)[:7] != 'NO_FILE' and os.path.basename(args.comparisons)[:7] != 'NO_FILE': 

		validate_files(args)
	
		comparisons = read_comparisons(args.comparisons)
		highlighted = ', '.join(["'%s'" % i for i in args.highlighted_genes])
	
		if args.convert_names == 'TRUE':
			download_orgdb(args.org)
	
		p = Pool(args.threads)
		p.map(process_comparison, ((comparison[0], comparison[1], comparison[2], comparison[3], highlighted, args) for comparison in comparisons))
		p.close()
		p.join()

		write_params(args)

	else:
		print("Skipped")

def validate_files(args):

	warnings = []
	errors = []

	warnings.extend(check_filenames([args.counts, args.groups, args.comparisons]))
	
	comparisons = check_comparison_file(args.comparisons, errors)
	samples, groups = check_group_file(args.groups, errors)

	for comparison in comparisons:
		treat, control, name, column = comparison
		
		if column not in groups.keys():
			errors.append('\tError: %s from comparison file not a column in groups file.' % (column))
		if treat not in groups[column]:
			errors.append('\tError: %s from comparison file not in %s from groups file.' % (treat, column))
		if control not in groups[column]:
			errors.append('\tError: %s from comparison file not in %s from groups file.' % (control, column))

	data_columns = check_count_file(args.counts, args.exclude_column, errors)

	for sample in samples:
		if sample not in data_columns:
			errors.append('\tError: %s from group file not a column in count file.' % (sample))

	if len(warnings) > 0:
		print('%d Warnings Detected:\n' % len(warnings))
		print('\n\n'.join(warnings))

	if len(errors) > 0:
		error_message = '%d Errors Detected in Input Files:\n\n%s' % (len(errors), '\n\n'.join(errors))
		sys.exit(error_message)

def check_filenames(input_files):
	warnings = []

	problem_files = set()
	for input_file in input_files:
		if " " in input_file:
			problem_files.add(input_file)		
	
	if len(problem_files) > 0:
		warnings.append("\tWarning: Although allowed in the DE module, spaces are generally not recommended in filenames. Affected files: %s" % ('\n'.join((['\t\t%s' % input_file for input_file in problem_files]))))
		
	return(warnings)

def check_count_file(count_file, exclude_column, errors):

	header = open(count_file).readline().rstrip().split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		header = open(count_file).readline().rstrip().split(',')
		if len(header) > 1:
			sep = ','
		else:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (count_file)
			sys.exit(error_message)

	data_columns = set(header)

	with open(count_file) as infile:
		header = infile.readline().rstrip().split(sep)
		num_columns = len(header)
		gene_column = 0

		genes_encountered = set()
		duplicate_gene_names = set()
		missing_values = set()
		num_missing_values = 0
		na_values = set()
		num_na_values = 0
		non_integers = set()
		num_non_integers = 0
		wrong_length = set()
		
		for line in infile:
			cur = line.rstrip().split(sep)
			if len(cur) != num_columns:
				wrong_length.add(cur[gene_column])
			if cur[gene_column] in genes_encountered:
				duplicate_gene_names.add(cur[gene_column])
			genes_encountered.add(cur[gene_column])

			for i, element in enumerate(cur):
				if i != gene_column and i != exclude_column:
					if element == '':
						missing_values.add(cur[gene_column])
						num_missing_values += 1
					elif element == 'NA':
						na_values.add(cur[gene_column])
						num_na_values += 1
					else:
						if not element.isdigit():
							non_integers.add(cur[gene_column])
							num_non_integers += 1

		if len(duplicate_gene_names) > 0:
			errors.append("\tError: Duplicate gene names are not allowed in %s. %s duplicate(s) detected.\n\tExample duplicates: %s" % (count_file, len(duplicate_gene_names), ' '.join(list(duplicate_gene_names)[0:10])))

		if len(wrong_length) > 0:
			errors.append("\tError: Number of columns don't match header. %s issues(s) detected.\n\tExample Genes: %s" % (len(wrong_length), ' '.join(list(wrong_length)[0:10])))

		if len(missing_values) > 0:
			errors.append("\tError: Missing values are not allowed in %s. %s missing values(s) detected.\n\tExample Genes: %s" % (count_file, num_missing_values, ' '.join(list(missing_values)[0:10])))
			
		if len(na_values) > 0:
			errors.append("\tError: NA are not allowed in %s. %s NA values(s) detected.\n\tExample Genes: %s" % (count_file, num_na_values, ' '.join(list(na_values)[0:10])))
	
	return data_columns		
		
def check_group_file(group_infile, errors):

	header = open(group_infile).readline().replace('\ufeff', '').rstrip().split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		header = open(group_infile).readline().replace('\ufeff', '').rstrip().split(',')
		if len(header) > 1:
			sep = ','
		else:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (group_infile)
			sys.exit(error_message)

	samples = set()
	groups = defaultdict(set)

	with open(group_infile) as infile:
		header = infile.readline().replace('\ufeff', '').rstrip().split(sep)

		header_check = [i.lower() for i in header]
		if 'sample_name' not in header:
			sys.exit('Error: "sample_name" must be a column in %s\n%s' % (group_infile, '\n'.join(errors)))

		column_key = {}
		for i, column in enumerate(header):
			column_key[i] = column
		num_columns = len(header)

		for line in infile:
			cur = line.rstrip('\n').split(sep)
			samples.add(cur[0])
			for i, column in enumerate(cur):
				groups[column_key[i]].add(column)

	return samples, groups


def check_comparison_file(comparison_infile, errors):

	comparisons = []
	header = open(comparison_infile).readline().replace('\ufeff', '').rstrip().split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		sep = ','

	with open(comparison_infile) as infile:
		header = infile.readline().replace('\ufeff', '').rstrip().split(sep)
		header = [i.lower() for i in header]
		problem = False

		if 'controls' in header:
			controls_column = header.index('controls')
		else:
			errors.append('\tError: Cannot find "controls" column in %s' % comparison_infile)
			problem = True

		if 'treats' in header:
			treats_column = header.index('treats')
		else:
			errors.append('\tError: Cannot find "treats" column in %s' % comparison_infile)
			problem = True

		if 'names' in header:
			names_column = header.index('names')
		else:
			errors.append('\tError: Cannot find "names" column in %s' % comparison_infile)
			problem = True

		if problem:
			sys.exit('\n'.join(errors))

		if 'grouping_column' in header:
			grouping_column = header.index('grouping_column')

		for line in infile:
			cur = line.rstrip('\n').split(sep)

			if len(cur) == 3:
				comparisons.append((cur[treats_column], cur[controls_column], cur[names_column], 'group'))
			elif len(cur) == 4:

				if 'grouping_column' not in header:
					sys.exit('\tError: When using four column approach, "grouping_column" must exist in %s' % (comparison_infile))

				if cur[grouping_column] == '':
					comparisons.append((cur[treats_column], cur[controls_column], cur[names_column], 'group'))
				else:
					comparisons.append((cur[treats_column], cur[controls_column], cur[names_column], cur[grouping_column]))
			else:
				error_message = '\tError: Incorrect number of columns in comparison file. 3 or 4 required. %d columns in this line: \n\t\t%s' % (len(cur), line)
				sys.exit(error_message)
	
	return(comparisons)

def read_comparisons(comparison_infile):

	comparisons = []

	header = open(comparison_infile).readline().replace('\ufeff', '').rstrip('\n').split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		sep = ','

	header = open(comparison_infile).readline().replace('\ufeff', '').rstrip().split(sep)
	header = [i.lower() for i in header]

	control_column = header.index('controls')
	treats_column = header.index('treats')
	names_column = header.index('names')

	with open(comparison_infile) as infile:
		infile.readline()
		for line in infile:
			cur = line.rstrip('\n').split(sep)

			if len(cur) == 3:
				comparisons.append((cur[treats_column], cur[control_column], cur[names_column], 'group'))
			else:
				grouping_column = header.index('grouping_column')
				if cur[grouping_column] == '':
					comparisons.append((cur[treats_column], cur[control_column], cur[names_column], 'group'))
				else:
					comparisons.append((cur[treats_column], cur[control_column], cur[names_column], cur[grouping_column]))
	
	return comparisons

def read_groups(metadata_infile, column):

	groups = defaultdict(list)

	header = open(metadata_infile).readline().rstrip('\n').split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		sep = ','

	header = open(metadata_infile).readline().rstrip('\n').split(sep)

	if column in header:
		index = header.index(column)
	else:
		index = 1

	with open(metadata_infile) as infile:
		infile.readline()
		for line in infile:
			cur = line.rstrip('\n').split(sep)
			groups[cur[index]].append(cur[0])

	return groups

def process_comparison(parameters):

	treatment, control, name, column, highlighted, args = parameters

	groups = read_groups(args.groups, column)

	build_report(groups, treatment, control, name, column, highlighted, args, '')
	
	if args.excluded_events != '':
		build_report(groups, treatment, control, '%s_excludelist' % (name), column, highlighted, args, args.excluded_events)

def build_report(groups, treatment, control, name, column, highlighted, args, excluded_events):
	output = []
	output.append(print_header(name))
	output.append(print_libraries())
	output.append(print_functions())
	output.append(choose_modules(args))
	output.append(read_counts(groups, treatment, control, args, excluded_events))
	output.append(filter_counts(args))
	output.append(batch_correction(args))
	output.append(sample_info())
	output.append(quality_control())
	output.append(count_distribution(column))
	output.append(all2all())
	output.append(run_pca(args))
	output.append(limma(treatment, control, highlighted, args, name, column))
	output.append(volcano_plot())
	output.append(ma_plot())
	output.append(heatmap())
	output.append(write_tables(name))
	output.append(write_warnings(name))
	output.append(session_info())
	write_output(name, output)
	run_markdown(name)

def print_header(name):

	return(dedent(
	'''
	---
	title: {}
	date: "`r Sys.Date()`"
	output:
	  html_document:
	    code_folding: hide
	---
	'''.format(name)).strip())

def print_libraries():
	
	return(dedent(
	'''
	```{r, load_libraries, message=FALSE, include=FALSE}
	# Load Libraries
	library(dplyr)
	library(ggplot2)
	library(tidyr)
	library(DESeq2)
	library(scales)
	library(ggrepel)
	library(clusterProfiler)
	library(sva)
	library(gplots)
	library(yaml)
	library(DT)
	library(htmltools)
	library(edgeR)
	```'''))

def print_functions():

	return(dedent(
	'''
	```{r local_functions, include=FALSE}
	# Load custom functions
	reverselog = function() {
	  trans_new("reverselog", function(x) -log10(x), function(x) 10^(-x), log_breaks(base = 10), domain = c(1e-1000, Inf))
	}
	
	quiet <- function(x) { 
	  sink(tempfile()) 
	  on.exit(sink()) 
	  invisible(force(x)) 
	} 
	
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

	count_distribution = function(counts, samples, min_counts_per_event, group_by) {
	
	  options(scipen=999)
	
	  needed_colors = nrow(samples %>% distinct(!!sym(group_by)))
	  colors = c("#4682b4", "#000000", "#6d2578", "#DC6D00", "#1e6234", "#FFFF33", "#A65628", "#F781BF", "#999999")

	  df = as.data.frame(counts) %>% tibble::rownames_to_column("feature") %>%
	       pivot_longer(!feature, names_to='sample_name', values_to='count') %>%
	       left_join(samples, by='sample_name') %>%
	       group_by(feature, !!sym(group_by)) %>%
	       summarise(ave_count = mean(count), .groups='drop') %>%
	       filter(ave_count > 0)
	
	  ggplot(df, aes(x=ave_count, fill=as.factor(!!sym(group_by)))) +
	    theme_classic() +
	    theme(strip.background = element_blank(),
	          strip.text = element_blank()) +
	    facet_wrap(as.formula(paste("~", group_by)), ncol=1) +
	    scale_x_continuous(trans='log10', breaks=c(.1,1,10,100,1000,10000,100000,1000000), name='Raw Counts') +
	    scale_y_continuous(name='Number of Features', expand=c(0,0)) +
	    {if (needed_colors <= length(colors)) scale_fill_manual(values = colors, name='')} +
	    geom_histogram(bins=100) +
	    geom_vline(xintercept = min_counts_per_event, linetype=2, color='firebrick')
	}

	filter_counts = function(input, min_counts_per_event, min_samples_per_event) {
	  keep = rowSums(input >= min_counts_per_event) >= min_samples_per_event
	  return(input[keep,])
	}

	run_pca = function(input, retx=TRUE, center=TRUE, scale=TRUE, transformation='Default') {
	  if (transformation == 'None') {
	    keep = subset(input, apply(input, 1, var, na.rm = TRUE) >  0)
	    return(prcomp(t(keep), retx = retx, center = center, scale. = scale))
	  } else if (transformation == 'vst' | (transformation == 'Default' & ncol(input) > 50)) {
	    transformed = vst(input)
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  } else {
	    transformed = rlog(input)
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  }
	}

	batch_correction = function(df, batch_norm_method, batch_correction_column, batch_correction_group_column) {
	  set.seed(1)
	  norm_df = getNormalizedMatrix(df, method=batch_norm_method)
	  mode(norm_df) = 'integer'
	  adjusted_df = apply(norm_df, 2, function(x) as.integer(x) + runif(1, 0, 0.01))
	  colnames(adjusted_df) <- colnames(norm_df)
	  rownames(adjusted_df) <- rownames(norm_df)
	  corrected = quiet(ComBat_seq(adjusted_df, batch = batch_correction_column, group=batch_correction_group_column))
	  mode(corrected) = 'integer'
	  return(corrected)
	}

	call_significance = function(df, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, num_labeled) {
	  labels = df %>%
	           filter(padj < padj_significance_cutoff & abs(log2FoldChange) > fc_significance_cutoff) %>%
	           arrange(-abs(log2FoldChange)) %>%
	           mutate(Rank = row_number()) %>%
	           filter(Rank <= num_labeled) %>%
	           select(feature)
	
	  df = df %>% 
	       mutate(Direction = case_when(padj < padj_significance_cutoff & log2FoldChange > fc_significance_cutoff ~ 'Upregulated',
	                                    padj < padj_significance_cutoff & log2FoldChange < -fc_significance_cutoff ~ 'Downregulated',
	                                    TRUE ~ 'No Change'))
	
	  df = df %>%
	       mutate(Direction = factor(Direction, levels=c("Upregulated", "No Change", "Downregulated"))) %>%
	       mutate(Significant = case_when(Direction == 'No Change' ~ 'Non-Significant',
	                                      TRUE ~ 'Significant')) %>%
	       mutate(padj = case_when(padj < padj_floor ~ padj_floor,
	                               TRUE ~ padj)) %>%
	       mutate(log2FoldChange = case_when(log2FoldChange < 0 & abs(log2FoldChange) > fc_ceiling ~ -fc_ceiling,
	                                         log2FoldChange > 0 & abs(log2FoldChange) > fc_ceiling ~ fc_ceiling,
	                                         TRUE ~ log2FoldChange)
	       ) %>%
	       mutate(Label = case_when(((Significant == 'Significant' & feature %in% labels$feature) | Highlighted == TRUE) ~ alias)) %>%
	       mutate(Group = case_when(Highlighted==TRUE & Direction=='Upregulated' ~ 'Upregulated_Highlighted',
	                                Highlighted==TRUE & Direction=='No Change' ~ 'NoChange_Highlighted',
	                                Highlighted==TRUE & Direction=='Downregulated' ~ 'Downregulated_Highlighted',
	                                Highlighted==FALSE & Direction=='Upregulated' ~ 'Upregulated',
	                                Highlighted==FALSE & Direction=='No Change' ~ 'No Change',
	                                Highlighted==FALSE & Direction=='Downregulated' ~ 'Downregulated'
	                      )
	       )
	  return(df)
	}
	
	add_alias = function(df, add_alias=TRUE, fromType='ENSEMBL', toType='SYMBOL', org='org.Hs.eg.db') {
	  if (add_alias) {
	    mapped_names <- data.frame(bitr(df$feature, fromType=fromType, toType=toType, OrgDb='org.Hs.eg.db'))
	  
	    unique_from = mapped_names %>% 
	                  group_by(!!sym(fromType)) %>% 
	                  summarise(count = n()) %>% 
	                  filter(count == 1) %>%
	                  select(-count)
	    colnames(unique_from) = c('alias')
	    
	    unique_to = mapped_names %>% 
	                  group_by(!!sym(toType)) %>% 
	                  summarise(count = n()) %>% 
	                  filter(count == 1) %>%
	                  select(-count)
	    colnames(unique_to) = c('alias')
	        
	    one_to_one = mapped_names %>% filter(!!sym(toType) %in% unique_to$alias & !!sym(fromType) %in% unique_from$alias)
	  
	    return(df %>%
	           left_join(one_to_one, by=c(feature = fromType)) %>%
	           mutate(alias = case_when(!is.na(!!sym(toType)) ~ !!sym(toType),
	                                    TRUE ~ feature)) %>%
	           select(-!!sym(toType))
	    )
	    } else {
	      return(df %>% mutate(alias=feature))
	    }
	}
	
	add_highlights = function(df, to_highlight) {
	  return(
	    df %>% 
	    mutate(Highlighted = 
	      case_when((alias %in% to_highlight | feature %in% to_highlight) ~ TRUE,
	                 TRUE ~ FALSE
	      )
	    )
	  )
	}

	post_processing = function(res, padj_significance_cutoff, fc_significance_cutoff, num_labeled, highlighted, add_alias=TRUE, fromType='ENSEMBL', toType='SYMBOL', org='org.Hs.eg.db') {
	  post = add_alias(res, add_alias=add_alias, fromType=fromType, toType=toType, org='org.Hs.eg.db')
	  post = add_highlights(post, highlighted)
	  post = call_significance(post, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, num_labeled)
	  return(post)
	}

	pca_plot = function(pca, metadata, x='PC1', y='PC2', color_by='', shape_by='', fill_by='', alpha_by='', label_by='', size=4) {
	
	  col_names <- names(metadata)
	  metadata[,col_names] <- lapply(metadata[,col_names] , factor)
	
	  if (color_by=='NA') { color_by='' }
	  if (shape_by=='NA') { shape_by='' }
	  if (fill_by=='NA') { fill_by='' }
	  if (alpha_by=='NA') { alpha_by='' }
	  if (label_by=='NA') { label_by='' }
	  
	  if (color_by != '') {
	    needed_colors = nrow(metadata %>% distinct(!!sym(color_by)))
	  } else {
	    needed_colors = Inf
	  }
	  colors = c("#b22222", "#4682b4", "#6d2578", "#DC6D00", "#1e6234", "#FFFF33", "#A65628", "#F781BF", "#999999")
	
	  PoV = data.frame(PoV=pca$sdev^2/sum(pca$sdev^2)*100) %>% mutate(PC = paste0("PC", row_number()))
	  pca_results = as.data.frame(pca$x) %>%
	                tibble::rownames_to_column('sample_name') %>%
	                left_join(metadata, by='sample_name')
	  
	  PoV_X = round((PoV %>% filter(PC==x))$PoV, 2)
	  PoV_Y = round((PoV %>% filter(PC==y))$PoV, 2)
	
	  ggplot(pca_results, aes(x=!!sym(x), y=!!sym(y), color=!!sym(color_by), shape=!!sym(shape_by), fill=!!sym(fill_by), alpha=!!sym(alpha_by), label=!!sym(label_by))) +
	    theme_classic() +
	    xlab(paste0(x, ': ', PoV_X, '%')) +
	    ylab(paste0(y, ': ', PoV_Y, '%')) +
	    {if (needed_colors <= length(colors)) scale_color_manual(values = colors)} +
	    geom_point(size=size) +
	    {if (label_by != '') geom_label(label.size=NA, fill=NA, color='black')}
	}
	
	scree_plot = function(pca) {
	
	PoV = data.frame(PoV=pca$sdev^2/sum(pca$sdev^2)*100) %>%
	      mutate(PC = row_number()) %>%
	      mutate(Label = case_when(PoV >=.01 ~ as.character(round(PoV, 2)),
	                               TRUE ~ "<.01"))
	
	ggplot(PoV, aes(x=PC, y=PoV, label=Label)) +
	  theme_classic() +
	  scale_x_continuous(breaks=seq(1, min(nrow(PoV), 20), 1), limits=c(0, min(nrow(PoV), 20))) +
	  scale_y_continuous(name="Percent of Variation", limits=c(0,100), expand=c(0,0)) +
	  geom_bar(stat='identity', fill='black') +
	  geom_label(fill=NA, label.size = NA, vjust=-.05)
	}

	ma_plot = function(deseq_res, Y='log2FoldChange', padj_cutoff=.05, fc_cutoff=1, padj_floor=0, fc_ceiling=Inf, num_labeled=Inf,
	                  upregulated_color='firebrick', noChange_color='grey', downregulated_color='steelblue',
	                  upregulated_highlight_color='#1e6234', noChange_highlight_color='black', downregulated_highlight_color='#1e6234',
	                  upregulated_alpha=1, noChange_alpha=.3, downregulated_alpha=1,
	                  upregulated_highlight_alpha=1, noChange_highlight_alpha=1, downregulated_highlight_alpha=1,
	                  upregulated_size=2, noChange_size=1, downregulated_size=2,
	                  upregulated_highlight_size=2, noChange_highlight_size=2, downregulated_highlight_size=2,
	                  fc_markers=TRUE, fc_marker_color='grey', center_marker=TRUE, center_marker_color='black',
	                  display_upregulated_count = TRUE, display_downregulated_count = TRUE, display_noChange_count=FALSE) {
	
	  options(scipen=999)
	
	  colors = c('Upregulated'=upregulated_color, 'No Change'=noChange_color, 'Downregulated'=downregulated_color, 'Upregulated_Highlighted'=upregulated_highlight_color, 'NoChange_Highlighted'=noChange_highlight_color, 'Downregulated_Highlighted'=downregulated_highlight_color)
	  alphas = c('Upregulated'=upregulated_alpha, 'No Change'=noChange_alpha, 'Downregulated'=downregulated_alpha, 'Upregulated_Highlighted'=upregulated_highlight_alpha, 'NoChange_Highlighted'=noChange_highlight_alpha, 'Downregulated_Highlighted'=downregulated_highlight_alpha)
	  sizes = c('Upregulated'=upregulated_size, 'No Change'=noChange_size, 'Downregulated'=downregulated_size, 'Upregulated_Highlighted'=upregulated_highlight_size, 'NoChange_Highlighted'=noChange_highlight_size, 'Downregulated_Highlighted'=downregulated_highlight_size)
	
	  upregulated_count = nrow(deseq_res %>% filter(Direction == 'Upregulated'))
	  noChange_count = nrow(deseq_res %>% filter(Direction == 'No Change'))
	  downregulated_count = nrow(deseq_res %>% filter(Direction == 'Downregulated'))
	
	  return(
	    ggplot(deseq_res, aes(x=baseMean, y=!!sym(Y), color=Group, alpha=Group, size=Group, label=Label)) +
	      theme_classic() +
	      theme(legend.position = 'none') +
	      scale_x_continuous(trans='log10', name="Expression") +
	      scale_y_continuous(name="Fold Change (log2)") +
	      scale_color_manual(values=colors) +
	      scale_alpha_manual(values=alphas) +
	      scale_size_manual(values=sizes) +
	      geom_point() +
	      {if (fc_markers) geom_hline(yintercept = fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (center_marker) geom_hline(yintercept = 0, linetype=2, color=center_marker_color)} +
	      {if (fc_markers) geom_hline(yintercept = -fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (display_upregulated_count) annotate('text', x=Inf, y=Inf, hjust=1, vjust=1, color=upregulated_color, label=paste0("Upregulated: ", upregulated_count))} +
	      {if (display_downregulated_count) annotate('text', x=Inf, y=-Inf, hjust=1, vjust=-1, color=downregulated_color, label=paste0("Downregulated: ", downregulated_count))} +
	      {if (display_noChange_count) annotate('text', x=Inf, y=fc_cutoff, hjust=1, vjust=1.2, color=noChange_color, label=paste0("No Change: ", noChange_count))} +
	      geom_label_repel(label.size=NA, fill=NA, na.rm=TRUE, max.overlaps = 50, max.time = 5)
	  )
	}
	
	volcano_plot = function(deseq_res, X='log2FoldChange', padj_cutoff=.05, fc_cutoff=1, padj_floor=0, fc_ceiling=Inf, num_labeled=Inf,
	                        upregulated_color='firebrick', noChange_color='grey', downregulated_color='steelblue',
	                        upregulated_highlight_color='#1e6234', noChange_highlight_color='black', downregulated_highlight_color='#1e6234',
	                        upregulated_alpha=1, noChange_alpha=.3, downregulated_alpha=1,
	                        upregulated_highlight_alpha=1, noChange_highlight_alpha=1, downregulated_highlight_alpha=1,
	                        upregulated_size=2, noChange_size=1, downregulated_size=2,
	                        upregulated_highlight_size=2, noChange_highlight_size=2, downregulated_highlight_size=2,
	                        fc_markers=TRUE, fc_marker_color='grey', center_marker=TRUE, center_marker_color='black',
	                        padj_marker=TRUE, padj_marker_color='grey',
	                        display_upregulated_count = TRUE, display_downregulated_count = TRUE, display_noChange_count=FALSE) {
	
	  colors = c('Upregulated'=upregulated_color, 'No Change'=noChange_color, 'Downregulated'=downregulated_color, 'Upregulated_Highlighted'=upregulated_highlight_color, 'NoChange_Highlighted'=noChange_highlight_color, 'Downregulated_Highlighted'=downregulated_highlight_color)
	  alphas = c('Upregulated'=upregulated_alpha, 'No Change'=noChange_alpha, 'Downregulated'=downregulated_alpha, 'Upregulated_Highlighted'=upregulated_highlight_alpha, 'NoChange_Highlighted'=noChange_highlight_alpha, 'Downregulated_Highlighted'=downregulated_highlight_alpha)
	  sizes = c('Upregulated'=upregulated_size, 'No Change'=noChange_size, 'Downregulated'=downregulated_size, 'Upregulated_Highlighted'=upregulated_highlight_size, 'NoChange_Highlighted'=noChange_highlight_size, 'Downregulated_Highlighted'=downregulated_highlight_size)
	
	  upregulated_count = nrow(deseq_res %>% filter(Direction == 'Upregulated'))
	  noChange_count = nrow(deseq_res %>% filter(Direction == 'No Change'))
	  downregulated_count = nrow(deseq_res %>% filter(Direction == 'Downregulated'))
	
	  return(
	    ggplot(deseq_res, aes(x=!!sym(X), y=padj, color=Group, alpha=Group, size=Group, label=Label)) +
	      theme_classic() +
	      theme(legend.position = 'none') +
	      scale_x_continuous(name="Fold Change (log2)") +
	      scale_y_continuous(trans=reverselog(), name='Significance', labels=trans_format('log10',math_format(10^.x))) +
	      scale_color_manual(values=colors) +
	      scale_alpha_manual(values=alphas) +
	      scale_size_manual(values=sizes) +
	      geom_point() +
	      {if (fc_markers) geom_vline(xintercept = fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (center_marker) geom_vline(xintercept = 0, linetype=2, color=center_marker_color)} +
	      {if (fc_markers) geom_vline(xintercept = -fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (padj_marker) geom_hline(yintercept = padj_cutoff, linetype=2, color=padj_marker_color)} +
	      {if (display_upregulated_count) annotate('text', x=Inf, y=0, hjust=1, vjust=1, color=upregulated_color, label=paste0("Upregulated: ", upregulated_count))} +
	      {if (display_downregulated_count) annotate('text', x=-Inf, y=0, hjust=-.01, vjust=1, color=downregulated_color, label=paste0("Downregulated: ", downregulated_count))} +
	      {if (display_noChange_count) annotate('text', x=0, y=0, vjust=1, color=noChange_color, label=paste0("No Change: ", noChange_count))} +
	      geom_label_repel(label.size=NA, fill=NA, na.rm=TRUE, max.overlaps = 50, max.time = 5)
	  )
	}
	
	heatmap_plot = function(df) {
	  norm = getNormalizedMatrix(df, method='TMM')
	  ld = log2(norm+.1)
	  cldt = scale(t(ld), center=TRUE, scale=TRUE)
	  cld = t(cldt)
	  
	  if (nrow(cld) > 1) {
	  	heatmap.2(cld, Rowv=TRUE, dendrogram='column', Colv=TRUE, col=bluered(256), 
	              labRow=NA, density.info='none', trace='none', cexCol=.8,
	              hclust=function(x) hclust(x,method="complete"),
	              distfun=function(x) as.dist((1-cor(t(x)))/2)
	    )
	  } else {
	      x = as.data.frame(cld) %>% tibble::rownames_to_column('Feature') %>%
	          pivot_longer(!Feature, names_to = "Sample", values_to = "Value")
	
	      ggplot(x, aes(x=Sample, y=Feature, fill=Value)) +
	             theme_classic() +
	             theme(axis.line = element_blank(),
	                   axis.ticks= element_blank()) +
	             scale_fill_gradient(low="#0000FF" , high='#FF0000') +
	             geom_tile()
	  }
	}

	all2all <- function(data, cex=2) {
	    pcor <- function(x, y, ...) panel.cor(x, y, cex.cor = cex)
	    nr <- nrow(data)
	    if (nr > 1000)
	        nr <- 1000
	    pairs(log10(data[1:nr, ]), cex = 0.25,
	            diag.panel = panel.hist, lower.panel = pcor)
	}

	panel.hist <- function(x, ...) {
	    usr <- par("usr")
	    on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5))
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks
	    nb <- length(breaks)
	    y <- h$counts
	    y <- y / max(y)
	    rect(breaks[-nb], 0, breaks[-1], y, col = "red", ...)
	}
	
	panel.cor <- function(x, y, prefix = "rho=", cex.cor=2, ...){
	    usr <- par("usr")
	    on.exit(par(usr))
	    par(usr = c(0, 1, 0, 1))
	    r <- cor.test(x, y, method = "spearman",
	        na.rm = TRUE, exact = FALSE)$estimate
	    txt <- round(r, digits = 2)
	    txt <- paste0(prefix, txt)
	    text(0.5, 0.5, txt, cex = cex.cor)
	}

	getNormalizedMatrix <- function(M = NULL, method = "TMM") {
	    if (is.null(M) ) return (NULL)
	    M[is.na(M)] <- 0
	    norm <- M
	    if (!(method == "none" || method == "MRN")){
	        norm.factors <- edgeR::calcNormFactors(M, method = method)
	        norm <- edgeR::equalizeLibSizes(edgeR::DGEList(M,
	            norm.factors = norm.factors))$pseudo.counts
	    }else if(method == "MRN"){
	        columns <- colnames(M)
	        conds <- columns
	        coldata <- prepGroup(conds, columns)
	        M[, columns] <- apply(M[, columns], 2,
	            function(x) as.integer(x))
	        dds <- DESeqDataSetFromMatrix(countData = as.matrix(M),
	            colData = coldata, design = ~group)
	        dds <- estimateSizeFactors(dds)
	        norm <- counts(dds, normalized=TRUE)
	    }
	    return(norm)
	}

	limma_results_table = function(post_res, name) {	
	
	  rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(-abs(log2FoldChange))) %>% select(feature, Rank)
	  df = post_res %>% left_join(rank, by='feature') %>% select(feature, baseMean, pvalue, padj, log2FoldChange, t, B, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange)) %>% arrange(Rank, abs_fc)
	  
	  return(
	    datatable(df,
	      rownames=FALSE,
	      extensions = 'Buttons',
	      options=list(
	        columnDefs=list(list(visible=FALSE, targets=c(7,8,9))),
	        dom = 'lftBipr',
	        buttons = list(
	          list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_limma_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	         )
	        )
	    ) %>% 
	    formatRound(columns=c('baseMean', 'log2FoldChange', 't', 'B', 'pvalue', 'padj'), digits=4) %>%
	    formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	    formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	  )
	}

	write_output = function(post_res, feature_type, prefix, name) {

	  df = post_res %>% relocate(alias) %>% select(-feature)

	  if (feature_type == 'gene') {
	    df = df %>% rename(gene=alias)
	  } else {
	    df = df %>% rename(transcript=alias)
	  }

	  all = df %>% select(-Highlighted, -Direction, -Significant, -Label, -Group)
	  sig = df %>% filter(Significant == 'Significant') %>% select(-Highlighted, -Direction, -Significant, -Label, -Group)

	  write.table(all, file=paste0("outputs/", prefix, "_all_limmaVoom_results.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	  write.table(sig, file=paste0("outputs/", prefix, "_sig_limmaVoom_results.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	}
	```'''))

def choose_modules(args):

	return(dedent(
	'''
	```{{r, choose_modules}}
	include_count_distribution = {}
	include_all2all = {}
	include_pca = {}
	include_batch_correction = {}
	include_limma = {}
	convert_names = {}
	include_ma = {}
	include_volcano = {} 
	include_heatmap = {}
	include_volcano_highlighted = {}
	include_ma_highlighted = {}
	```'''.format(args.include_distribution, args.include_all2all, args.include_pca, args.include_batch_correction, args.include_limma, args.convert_names, args.include_ma, args.include_volcano, args.include_heatmap, args.include_volcano_highlighted, args.include_ma_highlighted)))

def read_counts(groups, treatment, control, args, excluded_events):

	treatments = ', '.join(["'%s'" % sample for sample in groups[treatment]])
	controls = ', '.join(["'%s'" % sample for sample in groups[control]])

	return(dedent(
	'''
	```{{r, load_counts}}
	excluded_events = c({})
	feature_type = "{}"

	# Load count data
	all_metadata = read_metadata("{}{}")
	all_counts = read_counts("{}{}", all_metadata$sample_name, excluded_events)

	min_counts_per_sample = {}
	min_counts_per_event = {}
	min_samples_per_event = {}

	bad_samples = check_samples(all_counts, min_counts_per_sample)
	all_counts = sample_filter(all_counts, bad_samples)

	treatments = c({})
	good_treatments = setdiff(treatments, bad_samples)

	controls = c({})
	good_controls = setdiff(controls, bad_samples)
	
	selected_counts = read_counts("{}{}", c(good_treatments, good_controls), excluded_events)
	selected_metadata = data.frame(sample_name = colnames(selected_counts)) %>% left_join(all_metadata, by='sample_name')
	```'''.format(','.join('"%s"' % (i) for i in excluded_events), args.feature_type, args.input_prefix, args.groups, args.input_prefix, args.counts, args.min_counts_per_sample, args.min_counts_per_event, args.min_samples_per_event, treatments, controls, args.input_prefix, args.counts)))

def filter_counts(args):

	return(dedent(
	'''
	```{{r, filter, echo=FALSE}}
	filter_type = '{}'

	all_counts_filtered = filter_counts(all_counts, min_counts_per_event, min_samples_per_event)
	
	if (filter_type=='local') {{
	  selected_counts_filtered = filter_counts(selected_counts, min_counts_per_event, min_samples_per_event) %>% as.matrix()
	}} else {{
	  selected_counts_filtered = as.data.frame(all_counts_filtered) %>% select(all_of(good_treatments), all_of(good_controls)) %>% as.matrix()
	}}
	```'''.format(args.filter_type)))

def batch_correction(args):

	return(dedent(
	'''
	```{{r, batch_correction, eval=include_batch_correction, echo=FALSE}}
	batch_norm_method = '{}'
	batch_correction_column = as.factor(all_metadata${})
	batch_correction_group_column = as.factor(all_metadata${})
	
	all_counts_batch_corrected = batch_correction(all_counts_filtered, batch_norm_method, batch_correction_column, batch_correction_group_column)
	
	selected_counts_batch_corrected = as.data.frame(all_counts_batch_corrected) %>% 
	                                  select(all_of(good_treatments), all_of(good_controls)) %>% 
	                                  tibble::rownames_to_column("Event") %>%
	                                  filter(Event %in% rownames(selected_counts_filtered)) %>%
	                                  tibble::column_to_rownames("Event") %>%
	                                  as.matrix()
	```'''.format(args.batch_correction_normalization_algorithm, args.batch_correction_column, args.batch_correction_group_column)))

def sample_info():

	return(dedent(
	'''
	```{r, sample_info, echo=FALSE, results='asis'}
	cat("## Sample Info")
	datatable(selected_metadata,
	  rownames = FALSE					
	)
	```'''))

def quality_control():

	return(dedent(
	'''
	```{r, QC, eval=include_count_distribution | include_all2all | include_pca, echo=FALSE, results='asis'}
	cat('## Quality Control {.tabset .tabset-pills}')
	```'''))

def count_distribution(column):

	return(dedent(
	'''
	```{{r, count_distribution, eval=include_count_distribution, echo=FALSE, results='asis'}}
	cat('### Raw Count Distribution\\n')
	cat('In Eukaryotes only a subset of all genes are expressed in a given cell. Expression is therefore a bimodal distribution, with non-expressed genes/transcripts having counts that result from experimental and biological noise. It is important to filter out the features that are not expressed before doing differential gene/transcript expression. You can decide which cutoff separates expressed vs non-expressed features by looking at the following histogram. The current count threshold is at least ', min_counts_per_event, ' raw reads (red line) in at least ', min_samples_per_event, ' samples.\\n')
	
	count_distribution(selected_counts, selected_metadata, min_counts_per_event, "{}")
	```'''.format(column)))

def all2all():

	return(dedent(
	'''
	```{r, reproducibility, eval=include_all2all, echo=FALSE, results='asis', warning=FALSE}
	cat('### Reproducibility\\n')
	

	if (ncol(selected_counts) <= 10) {
	  cat('To check the reproducibility of biological replicates, we use all2all plots.\\n')
	  all2all(selected_counts_filtered, cex=1)
	} else {
	  cat('With more than 10 samples, meaningful reproducibility plots cannot be created due to size restrictions.')
	}
	```'''))

def run_pca(args):

	return(dedent(
	'''
	```{{r, pca, eval=include_pca, echo=FALSE, results='asis', warning=FALSE}}
	transformation = '{}'
	
	cat('### PCA {{ .tabset .tabset-pills }}\\n')
	if (include_batch_correction) {{
	  cat('#### Pre batch correction\\n')
	}}
	
	pca = run_pca(all_counts_filtered, transformation = transformation)
	print(pca_plot(pca, all_metadata, color_by='{}', shape_by='{}', fill_by='{}', alpha_by='{}', label_by='{}'))
	print(scree_plot(pca))

	if (include_batch_correction) {{
	  cat('\\n\\n#### Post batch correction\\n')
	  pca = run_pca(all_counts_batch_corrected, transformation = transformation)
	  print(pca_plot(pca, all_metadata, color_by='{}', shape_by='{}', fill_by='{}', alpha_by='{}', label_by='{}'))
	  print(scree_plot(pca))
	}}
	```'''.format(args.transformation, args.pca_color, args.pca_shape, args.pca_fill, args.pca_transparency, args.pca_label, args.pca_color, args.pca_shape, args.pca_fill, args.pca_transparency, args.pca_label)))

def limma(treatment, control, highlighted, args, name, column):
	return(dedent(
	'''
	```{{r, limmma_info, eval=include_limma, echo=FALSE, results='asis'}}
	cat('## Limma Voom Analysis {{ .tabset .tabset-pills }}\\n')
	cat('The goal of Differential gene expression analysis is to find genes or transcripts whose difference in expression, when accounting for the variance within condition, is higher than expected by chance.')
	```

	```{{r, limma, eval=include_limma, include=include_limma, warning=FALSE, message=FALSE}}
	normalization_method = '{}'
	logratioTrim = {}					# Only affects TMM and TMMwsp methods
	sumTrim = {}							# Only affects TMM and TMMwsp methods
	doWeighting = {}					# Only affects TMM and TMMwsp methods
	Acutoff = {}		# Only affects TMM and TMMwsp methods
	p = {}										# Only affects upperquartile method

	padj_significance_cutoff = {}
	fc_significance_cutoff = {}
	padj_floor = {}
	fc_ceiling = {}
	num_labeled = {}
	highlighted = c({})

	use_batch_correction_in_DE = {}

	if (include_batch_correction & use_batch_correction_in_DE) {{
	  limma_input = selected_counts_batch_corrected
	}} else {{
	  limma_input = selected_counts_filtered 
	}}

	d0 = DGEList(limma_input)
	
	if (normalization_method == 'TMM' | normalization_method == 'TMMwsp') {{
		d0 = normLibSizes(d0, method=normalization_method, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
	}} else if (normalization_method == 'upperquartile') {{
		d0 = normLibSizes(d0, method=normalization_method, p=p)
	}} else {{
		d0 = normLibSizes(d0, method=normalization_method)
	}}
	
	keep = rowSums(limma_input >= min_counts_per_event) >= min_samples_per_event
	d = d0[keep,]
	
	group = as.factor(selected_metadata${})
	
	mm = model.matrix(~0+group)
	
	y=voom(d, mm, plot=FALSE)
	
	fit = lmFit(y, mm)
	
	contrasts = makeContrasts(group{} - group{}, levels=colnames(coef(fit)))
	temp = contrasts.fit(fit,contrasts)
	temp = eBayes(temp)
	results = topTable(temp, sort.by='P', n=Inf) %>% tibble::rownames_to_column('feature') %>% select(feature, baseMean=AveExpr, log2FoldChange=logFC, t, B, pvalue=P.Value, padj=adj.P.Val)
	post_res = post_processing(results, padj_significance_cutoff=padj_significance_cutoff, fc_significance_cutoff=fc_significance_cutoff, num_labeled=num_labeled, highlighted=highlighted, add_alias={}, fromType='{}',  toType='{}', org='org.Hs.eg.db')
	
	if (feature_type == 'transcript') {{
	  transcript_to_gene = read.delim("{}{}") %>% select(feature=1, gene=2)
	  post_res = post_res %>% left_join(transcript_to_gene, by='feature')
	}}

	limma_results_table(post_res, '{}')
	```'''.format(args.normalization_method, args.logratio_trim, args.sum_trim, args.do_weighting, args.A_cutoff, args.p, args.padj_significance_cutoff, args.fc_significance_cutoff, args.padj_floor, args.fc_ceiling, args.num_labeled, highlighted, args.use_bc_in_DE, column, treatment, control, args.convert_names, args.count_file_names, args.converted_names, args.input_prefix, args.counts, name)))

def ma_plot():
	return(dedent(
	'''
	```{r, ma_plot, echo=FALSE, eval=include_ma & include_limma, include=include_ma & include_limma, results='asis', warning=FALSE}
	if(include_ma_highlighted) {
	  cat('### MA {.tabset} \\n')
	  cat('#### All Detected \\n')
	} else {
	  cat('### MA\\n')
	}
	
	ma_plot(post_res, padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	
	if(include_ma_highlighted) {
	  cat('\\n\\n#### Highlighted Only \\n')
	  ma_plot(post_res %>% filter(Highlighted==TRUE), padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	}
	```'''))

def volcano_plot():
	return(dedent(
	'''
	```{r, volcano_plot, echo=FALSE, eval=include_volcano & include_limma, include=include_volcano & include_limma, results='asis', warning=FALSE}
	if(include_volcano_highlighted) {
	  cat('### Volcano {.tabset} \\n')
	  cat('#### All Detected \\n')
	} else {
	  cat('### Volcano\\n')
	}
	
	volcano_plot(post_res, padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	
	if(include_volcano_highlighted) {
	  cat('\\n\\n#### Highlighted Only \\n')
	  volcano_plot(post_res %>% filter(Highlighted==TRUE), padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	}
	```'''))

def heatmap():
	return(dedent(
	'''
	```{r, heatmap, echo=FALSE, eval=include_heatmap & include_limma, include=include_heatmap & include_limma, results='asis'}
	cat('### Heatmap\\n')
	
	significant = (post_res %>% filter(Significant == 'Significant'))$feature
	
	heatmap_data = as.data.frame(selected_counts) %>%
	               tibble::rownames_to_column('feature') %>%
	               filter(feature %in% significant) %>%
	               tibble::column_to_rownames('feature')

	if (nrow(heatmap_data) > 0) {
	  
	  bad_samples = check_samples(heatmap_data, 1)
	  
	  if (length(bad_samples) == 0) {
	    heatmap_plot(heatmap_data)
	  } else {
	    cat("Cannot include the following samples because they have zero counts across all significant features which breaks the normalization function:\\n")
	    cat(bad_samples)
	    heatmap_data_filtered = sample_filter(heatmap_data, bad_samples)
	    heatmap_plot(heatmap_data_filtered)
	  }
	} else {
	  cat("No significant features. Skipping heatmap.")
	}
	```'''))

def write_tables(name):

	return(dedent(
	'''
	```{{r, echo=FALSE, write_tables}}	
	if (include_limma) {{
	  write_output(post_res, feature_type, '{}')
	}}
	```'''.format(name)))

def write_warnings(name):

	return(dedent(
	'''
	```{{r, write_warnings, echo=FALSE}}
	
	prefix = '{}'

	if (length(bad_samples) > 0) {{
	  fail_pca = paste("WARNING -", bad_samples, "is not included in PCA because it has less than", min_counts_per_sample, "counts.", sep=' ')
	}} else {{
	  fail_pca = NULL
	}}

	if (length(c(setdiff(treatments, good_treatments), setdiff(controls, good_controls))) > 0) {{
	  fail_DE = paste("WARNING -", c(setdiff(treatments, good_treatments), setdiff(controls, good_controls)), "is not included in DE because it has less than", min_counts_per_sample, "counts.", sep=' ')
	}} else {{
	  fail_DE = NULL
	}}
	
	if (length(bad_samples) > 0 | length(c(setdiff(treatments, good_treatments), setdiff(controls, good_controls))) > 0) {{
	  writeLines(c(fail_pca, fail_DE), paste0("outputs/", prefix, "_warning_log.txt"))
	}}
	```'''.format(name)))

def session_info():

	return(dedent(
	'''
	```{r, session_info, echo=FALSE, results='asis'}
	cat('## Session Info\\n')
	sessionInfo()
	```'''))

def write_output(name, output):

	with open('%s.Rmd' % (name), 'w') as out:
		out.write('\n'.join(output)) 

def write_params(args):
	output = []
	with open('outputs/limmaVoom_run_params.yaml', 'w') as out:
		for arg in vars(args):
			if getattr(args, arg) == '':
				output.append('%s: ~' % arg)
			elif getattr(args, arg) == 'Inf':
				output.append('%s: .inf' % arg)
			else:
				output.append('%s: %s' % (arg, getattr(args, arg)))
		out.write('\n'.join(output))

def download_orgdb(org):
	cmd = "R -e 'BiocManager::install(\"%s\", update=FALSE, force=TRUE)'" % (org)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print(out.decode())
	print(err.decode(), file=sys.stderr)

def run_markdown(name):
	
	cmd = "Rscript -e 'rmarkdown::render(\"%s.Rmd\", \"html_document\")'" % (name)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print('%s:' % name)
	print(out.decode())
	print('%s:' % name, file=sys.stderr)
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="prepare_voom", description='', usage='%(prog)s [options]')

	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-c', '--counts', required=True, help='Name of count file.', metavar='', dest='counts')
	input_args.add_argument('-g', '--groups', required=True, help='Name of groups file.', metavar='', dest='groups')
	input_args.add_argument('-x', '--comparisons', required=True, help='Name of comparisons file.', metavar='', dest='comparisons')
	input_args.add_argument('-p', '--prefix', default='', help='Path to input directory.', metavar='', dest='input_prefix')
	input_args.add_argument('--feature-type', choices=['gene','transcript'], default='gene', help='Feature type to run LimmaVoom on.', metavar='', dest='feature_type')

	qc_modules = parser.add_argument_group('QC Modules')
	qc_modules.add_argument('--include-distribution', choices=['TRUE', "FALSE"], default='TRUE', help='Create count distribution plot', metavar='', dest='include_distribution')
	qc_modules.add_argument('--include-all2all', choices=['TRUE', "FALSE"], default='TRUE', help='Create reproducibility plots', metavar='', dest='include_all2all')
	qc_modules.add_argument('--include-pca', choices=['TRUE', "FALSE"], default='TRUE', help='Create PCA plots', metavar='', dest='include_pca')

	filtering = parser.add_argument_group('Filtering Criteria')
	filtering.add_argument('--filter-type', choices=["global","local"], default='local', help='Apply filter across all samples (global) or just samples within the specific comparison (local).', metavar='', dest='filter_type')
	filtering.add_argument('--min-counts-per-event', type=int, default=10, help='Filter genes with less than this count.', metavar='', dest='min_counts_per_event')
	filtering.add_argument('--min-samples-per-event', type=int, default=1, help='Filter genes which do not have at least min-counts-per-event reads in this many samples.', metavar='', dest='min_samples_per_event')
	filtering.add_argument('--min-counts-per-sample', type=int, default=1000, help='Filter samples which do not have at least min_counts_per_sample reads.', metavar='', dest='min_counts_per_sample')
	filtering.add_argument('--excluded-events', type=str, nargs='+', default='', help='List of events to remove from analysis.', metavar='', dest='excluded_events')

	pca = parser.add_argument_group('PCA Settings')
	pca.add_argument('--transformation', choices=["Default","rlog","vst","None"], default='Default', help='Transformation to apply prior to PCA (does not apply to data sent into DESeq). Default will use the rlog transformation for up to 50 samples and vst if there are more than 50 samples.', metavar='', dest='transformation')
	pca.add_argument('--pca-color', type=str, default='', help='Column from metadata file by which to determine point color.', metavar='', dest='pca_color')
	pca.add_argument('--pca-shape', type=str, default='', help='Column from metadata file by which to determine point shape.', metavar='', dest='pca_shape')
	pca.add_argument('--pca-fill', type=str, default='', help='Column from metadata file by which to determine point fill.', metavar='', dest='pca_fill')
	pca.add_argument('--pca-transparency', type=str, default='', help='Column from metadata file by which to determine point transparency.', metavar='', dest='pca_transparency')
	pca.add_argument('--pca-label', type=str, default='', help='Column from metadata file by which to determine point label.', metavar='', dest='pca_label')

	batch_correction = parser.add_argument_group('Batch Correction Settings')
	batch_correction.add_argument('--include-batch-correction', choices=["TRUE", "FALSE"], default='FALSE', help='Perform batch correction prior to DESeq', metavar='', dest='include_batch_correction')
	batch_correction.add_argument('--batch-correction-column', type=str, default='batch', help='Column used for batch correction.', metavar='', dest='batch_correction_column')
	batch_correction.add_argument('--batch-correction-group-column', type=str, default='group', help='Column containing biological variable grouping.', metavar='', dest='batch_correction_group_column')
	batch_correction.add_argument('--batch-normalization-algorithm', choices=["MRN","TMM","RLE","upperquartile", "none"], default='TMM', help='Normalization algorithm to use prior to batch correction.', metavar='', dest='batch_correction_normalization_algorithm')

	limma = parser.add_argument_group('Limma Voom Settings')
	limma.add_argument('--include-limma', choices=["TRUE", "FALSE"], default='TRUE', help='Run Limma Voom', metavar='', dest='include_limma')
	limma.add_argument('--use-batch-correction-in-DE', choices=["TRUE", "FALSE"], default='TRUE', help='Use batch corrected values in Limma-Voom', metavar='', dest='use_bc_in_DE')
	limma.add_argument('--normalization-method', choices=["TMM", "TMMwsp", "RLE", "upperquartile", "none"], default='TMM', help='Normalization method to use prior to Limma Voom', metavar='', dest='normalization_method')
	limma.add_argument('--logratio-trim', type=float, default=.3, help='The fraction (0 to 0.5) of observations to be trimmed from each tail of the distribution of log-ratios (M-values) before computing the mean for each pair of samples.', metavar='', dest='logratio_trim')
	limma.add_argument('--sum-trim', type=float, default=.05, help='The fraction (0 to 0.5) of observations to be trimmed from each tail of the distribution of A-values before computing the mean for each pair of samples.', metavar='', dest='sum_trim')
	limma.add_argument('--A-cutoff', type=str, default='-1e10', help='Minimum cutoff applied to A-values for each pair of samples. Count pairs with lower A-values are ignored.', metavar='', dest='A_cutoff')
	limma.add_argument('--do-weighting', choices=["TRUE", "FALSE"], default='TRUE', help='Use (asymptotic binomial precision) weights when computing the mean M-values for each pair of samples.', metavar='', dest='do_weighting')
	limma.add_argument('--p', type=float, default=.75, help='Numeric value between 0 and 1 specifying which quantile of the counts should be used.', metavar='', dest='p')
	limma.add_argument('--include-volcano', choices=["TRUE", "FALSE"], default='TRUE', help='Create volcano plots', metavar='', dest='include_volcano')
	limma.add_argument('--include-ma', choices=["TRUE", "FALSE"], default='TRUE', help='Create MA plots', metavar='', dest='include_ma')
	limma.add_argument('--include-heatmap', choices=["TRUE", "FALSE"], default='TRUE', help='Create heatmap', metavar='', dest='include_heatmap')

	significance = parser.add_argument_group('Significance Settings')
	significance.add_argument('--padj-significance-cutoff', type=float, default=.05, help='Adjusted p-value significance threshold required to be called differentially expressed.', metavar='', dest='padj_significance_cutoff')
	significance.add_argument('--fc-significance-cutoff', type=float, default=1, help='log2 Fold Change threshold required to be called differentially expressed.', metavar='', dest='fc_significance_cutoff')
	significance.add_argument('--padj-floor', type=str, default='0', help='Reassign adjusted p-values below this value to this value. Set at 0 to disable.', metavar='', dest='padj_floor')
	significance.add_argument('--fc-ceiling', type=str, default='Inf', help='Reassign fold changes (log2) greater than this value to this value. Set to Inf to disable.', metavar='', dest='fc_ceiling')

	label = parser.add_argument_group('Label Options')
	label.add_argument('--convert-names', choices=['TRUE', 'FALSE'], default='FALSE', help='Convert gene names', metavar='', dest='convert_names')
	label.add_argument('--count-file-names', choices=["ENSEMBL","ENTREZID","SYMBOL"], default='SYMBOL', help='Name used in count file.', metavar='', dest='count_file_names')
	label.add_argument('--converted-names', choices=["ENSEMBL","ENTREZID","SYMBOL"], default='ENSEMBL', help='Name to use in plots.', metavar='', dest='converted_names')
	label.add_argument('--org-db', default='org.Hs.eg.db', help='OrgDb database', metavar='', dest='org')
	label.add_argument('--num-labeled', type=str, default='Inf', help='Number of significant genes to label. Set to Inf to label all. Set to 0 to turn labels off.', metavar='', dest='num_labeled')
	label.add_argument('--highlighted-genes', type=str, nargs='+', default='', help='Space separated list of genes to label regardless of significance.', metavar='', dest='highlighted_genes')
	label.add_argument('--include-volcano-highlighted', choices=["TRUE", "FALSE"], default='FALSE', help='Create volcano plots of highlighted genes only.', metavar='', dest='include_volcano_highlighted')
	label.add_argument('--include-ma-highlighted', choices=["TRUE", "FALSE"], default='FALSE', help='Create MA plots of highlighted genes only.', metavar='', dest='include_ma_highlighted')

	other = parser.add_argument_group('Other Options')
	other.add_argument('-e', '--exclude', type=int, default=1, help='Column to exclude from value checks. Set to -1 to disable.', metavar='', dest='exclude_column')
	other.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for knitting.', metavar='', dest='threads')

	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)