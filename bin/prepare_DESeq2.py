#!/usr/bin/env python3

import sys, os, argparse, subprocess, csv
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

		normalization_report(args)
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

def check_count_file(count_infile, exclude_column, errors):

	with open(count_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)

		if len(header) <= 1:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (count_infile)
			sys.exit(error_message)
		
		header = [i.replace('\ufeff', '') for i in header]

		data_columns = set(header)
		
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
		
		for line in reader:
			if len(line) != num_columns:
				wrong_length.add(line[gene_column])
			if line[gene_column] in genes_encountered:
				duplicate_gene_names.add(line[gene_column])
			genes_encountered.add(line[gene_column])

			for i, element in enumerate(line):
				if i != gene_column and i != exclude_column:
					if element == '':
						missing_values.add(line[gene_column])
						num_missing_values += 1
					elif element == 'NA':
						na_values.add(line[gene_column])
						num_na_values += 1
					else:
						if not element.isdigit():
							non_integers.add(line[gene_column])
							num_non_integers += 1

		if len(duplicate_gene_names) > 0:
			errors.append("\tError: Duplicate gene names are not allowed in %s. %s duplicate(s) detected.\n\tExample duplicates: %s" % (count_infile, len(duplicate_gene_names), ' '.join(list(duplicate_gene_names)[0:10])))

		if len(wrong_length) > 0:
			errors.append("\tError: Number of columns don't match header. %s issues(s) detected.\n\tExample Genes: %s" % (len(wrong_length), ' '.join(list(wrong_length)[0:10])))

		if len(missing_values) > 0:
			errors.append("\tError: Missing values are not allowed in %s. %s missing values(s) detected.\n\tExample Genes: %s" % (count_infile, num_missing_values, ' '.join(list(missing_values)[0:10])))
			
		if len(na_values) > 0:
			errors.append("\tError: NA are not allowed in %s. %s NA values(s) detected.\n\tExample Genes: %s" % (count_infile, num_na_values, ' '.join(list(na_values)[0:10])))
	
	return data_columns		
		
def check_group_file(group_infile, errors):

	with open(group_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)
		if len(header) <= 1:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (group_infile)
			sys.exit(error_message)

		header = [i.replace('\ufeff', '').lower() for i in header]

		if 'sample_name' not in header:
			sys.exit('Error: "sample_name" must be a column in %s\n%s' % (group_infile, '\n'.join(errors)))

		samples = set()
		groups = defaultdict(set)

		column_key = {}
		for i, column in enumerate(header):
			column_key[i] = column
		num_columns = len(header)

		for line in reader:
			samples.add(line[0])
			for i, column in enumerate(line):
				groups[column_key[i]].add(column)

	return samples, groups

def check_comparison_file(comparison_infile, errors):

	comparisons = []
	
	with open(comparison_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)
		
		if len(header) <= 1:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (comparison_infile)
			sys.exit(error_message)

		header = [i.replace('\ufeff', '').lower() for i in header]

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

		for line in reader:

			if len(line) == 3:
				comparisons.append((line[treats_column], line[controls_column], line[names_column], 'group'))
			elif len(line) == 4:

				if 'grouping_column' not in header:
					sys.exit('\tError: When using four column approach, "grouping_column" must exist in %s' % (comparison_infile))

				if line[grouping_column] == '':
					comparisons.append((line[treats_column], line[controls_column], line[names_column], 'group'))
				else:
					comparisons.append((line[treats_column], line[controls_column], line[names_column], line[grouping_column]))
			else:
				error_message = '\tError: Incorrect number of columns in comparison file. 3 or 4 required. %d columns in this line: \n\t\t%s' % (len(line), '\n'.join(line))
				sys.exit(error_message)
	
	return(comparisons)

def read_comparisons(comparison_infile):

	comparisons = []

	with open(comparison_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)
		

		header = [i.replace('\ufeff', '').lower() for i in header]


		control_column = header.index('controls')
		treats_column = header.index('treats')
		names_column = header.index('names')

		for line in reader:

			if len(line) == 3:
				comparisons.append((line[treats_column], line[control_column], line[names_column], 'group'))
			else:
				grouping_column = header.index('grouping_column')
				if line[grouping_column] == '':
					comparisons.append((line[treats_column], line[control_column], line[names_column], 'group'))
				else:
					comparisons.append((line[treats_column], line[control_column], line[names_column], line[grouping_column]))
		
		return comparisons

def read_groups(group_infile, column):

	groups = defaultdict(list)

	with open(group_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)

		header = [i.replace('\ufeff', '').lower() for i in header]

		sample_name_column = header.index('sample_name')

		if column.lower() in header:
			index = header.index(column)
		else:
			index = 1

		for line in reader:
			groups[line[index]].append(line[sample_name_column])

	return groups

def process_comparison(parameters):

	treatment, control, name, column, highlighted, args = parameters

	groups = read_groups(args.groups, column)

	if args.design.replace('~', '') != column:
		design = '~%s' % column
	else:
		design = args.design

	build_report(groups, design, treatment, control, name, column, highlighted, args, '')
	
	if args.excluded_events != '':
		build_report(groups, design, treatment, control, '%s_excludelist' % (name), column, highlighted, args, args.excluded_events)

def build_report(groups, design, treatment, control, name, column, highlighted, args, excluded_events):
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
	output.append(run_deseq(treatment, control, highlighted, args, name, design, column))
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
	```'''))

def print_functions():

	return(dedent(
	'''
	```{r local_functions, include=FALSE}
	# Load custom functions
	
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
	
	run_pca = function(input, retx=TRUE, center=TRUE, scale=TRUE, transformation='Default', write_transform=FALSE) {
	  if (transformation == 'None') {
	    keep = subset(input, apply(input, 1, var, na.rm = TRUE) >  0)
	    return(prcomp(t(keep), retx = retx, center = center, scale. = scale))
	  } else if (transformation == 'vst' | (transformation == 'Default' & ncol(input) > 50)) {
	    transformed = vst(input)
	    if(write_transform) {
	      write.table(as.data.frame(transformed) %>% tibble::rownames_to_column("feature"), file='outputs/transformed_counts_vst.tsv', quote=FALSE, sep='\t', row.names=FALSE)
	    }
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  } else {
	    transformed = rlog(input)
	    if(write_transform) {
	      write.table(as.data.frame(transformed) %>% tibble::rownames_to_column("feature"), file='outputs/transformed_counts_rlog.tsv', quote=FALSE, sep='\t', row.names=FALSE)
	    }
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

	call_significance = function(df, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, num_labeled, apply_shrinkage) {
	  labels = df %>%
	           filter(padj < padj_significance_cutoff & abs(log2FoldChange) > fc_significance_cutoff) %>%
	           arrange(-abs(log2FoldChange)) %>%
	           mutate(Rank = row_number()) %>%
	           filter(Rank <= num_labeled) %>%
	           select(feature)
	
	  if (apply_shrinkage) {
	    df = df %>% 
	         mutate(Direction = case_when(padj < padj_significance_cutoff & log2FoldChange_shrink > fc_significance_cutoff ~ 'Upregulated',
	                                      padj < padj_significance_cutoff & log2FoldChange_shrink < -fc_significance_cutoff ~ 'Downregulated',
	                                      TRUE ~ 'No Change')) %>%
	         mutate(log2FoldChange_shrink = case_when(log2FoldChange_shrink < 0 & abs(log2FoldChange_shrink) > fc_ceiling ~ -fc_ceiling,
	                                                  log2FoldChange_shrink > 0 & abs(log2FoldChange_shrink) > fc_ceiling ~ fc_ceiling,
	                                                  TRUE ~ log2FoldChange_shrink))
	  } else {
	    df = df %>% 
	         mutate(Direction = case_when(padj < padj_significance_cutoff & log2FoldChange > fc_significance_cutoff ~ 'Upregulated',
	                                      padj < padj_significance_cutoff & log2FoldChange < -fc_significance_cutoff ~ 'Downregulated',
	                                      TRUE ~ 'No Change'))
	  }
	
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
	
	post_processing = function(deseq_res, padj_significance_cutoff, fc_significance_cutoff, num_labeled, highlighted, add_alias=TRUE, fromType='ENSEMBL', toType='SYMBOL', org='org.Hs.eg.db', apply_shrinkage=FALSE) {
	  post = add_alias(deseq_res, add_alias=add_alias, fromType=fromType, toType=toType, org='org.Hs.eg.db')
	  post = add_highlights(post, highlighted)
	  post = call_significance(post, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, num_labeled, apply_shrinkage)
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
	      scale_y_continuous(trans=c("log10", "reverse"), name='Significance', labels=trans_format('log10',math_format(10^.x))) +
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

	deseq_results_table = function(post_res, apply_shrinkage, feature_type, name) {	
	
	  if (apply_shrinkage & feature_type == 'gene') {
	  
	    sketch = htmltools::withTags(table(
	      class = 'display',
	      thead(
	        tr(
	          th(rowspan=2, 'gene'), th(rowspan=2, 'baseMean'), th(rowspan=2, 'pvalue'), th(rowspan=2, 'padj'), th(colspan=3, style = "text-align:center; border-left: solid 1px;", "Standard"), th(colspan=2, style = "text-align:center; border-left: solid 1px;", "Shrinked"), th(rowspan=2, 'Direction'), th(rowspan=2, 'Rank'), th(rowspan=2, '-|log2FoldChange|'),
	        ),
	        tr(
	          th("log2FoldChange", style="border-left: solid 1px;"), th("lfcSE"), th("stat"), th("log2FoldChange", style="border-left: solid 1px;"), th("lfcSE")
	        )
	      )
	    ))
	
	    rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(-abs(log2FoldChange_shrink))) %>% select(feature, Rank)
	    df = post_res %>% left_join(rank, by='feature') %>% select(alias, baseMean, pvalue, padj, log2FoldChange, lfcSE, stat, log2FoldChange_shrink, lfcSE_shrink, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange_shrink)) %>% arrange(Rank, abs_fc)
	    return(
	      datatable(df,
	        rownames=FALSE,
	        extensions = 'Buttons',
	        options=list(
	          columnDefs=list(list(visible=FALSE, targets=c(9, 10, 11))),
	          dom = 'lftBipr',
	          buttons = list(
	            list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_deseq2_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	           )
	          ),
	        container=sketch
	      ) %>% 
	      formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'log2FoldChange_shrink', 'lfcSE_shrink'), digits=4) %>%
	      formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	      formatStyle(c(5, 8), `border-left` = "solid 1px") %>%
	      formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	    )
	  } else if (!apply_shrinkage & feature_type == 'gene') {
	
	    rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(-abs(log2FoldChange))) %>% select(feature, Rank)
	    df = post_res %>% left_join(rank, by='feature') %>% select(gene=alias, baseMean, pvalue, padj, log2FoldChange, lfcSE, stat, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange)) %>% arrange(Rank, abs_fc)
	
	    return(
	      datatable(df,
	        rownames=FALSE,
	        extensions = 'Buttons',
	        options=list(
	          columnDefs=list(list(visible=FALSE, targets=c(7,8,9))),
	          dom = 'lftBipr',
	          buttons = list(
	            list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_deseq2_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	           )
	          )
	      ) %>% 
	      formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), digits=4) %>%
	      formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	      formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	    )

	  } else if (apply_shrinkage & feature_type == 'transcript') {

	    sketch = htmltools::withTags(table(
	      class = 'display',
	      thead(
	        tr(
	          th(rowspan=2, 'transcript'), th(rowspan=2, 'gene'), th(rowspan=2, 'baseMean'), th(rowspan=2, 'pvalue'), th(rowspan=2, 'padj'), th(colspan=3, style = "text-align:center; border-left: solid 1px;", "Standard"), th(colspan=2, style = "text-align:center; border-left: solid 1px;", "Shrinked"), th(rowspan=2, 'Direction'), th(rowspan=2, 'Rank'), th(rowspan=2, '-|log2FoldChange|'),
	        ),
	        tr(
	          th("log2FoldChange", style="border-left: solid 1px;"), th("lfcSE"), th("stat"), th("log2FoldChange", style="border-left: solid 1px;"), th("lfcSE")
	        )
	      )
	    ))
	
	    rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(-abs(log2FoldChange_shrink))) %>% select(feature, Rank)
	    df = post_res %>% left_join(rank, by='feature') %>% select(feature, gene, baseMean, pvalue, padj, log2FoldChange, lfcSE, stat, log2FoldChange_shrink, lfcSE_shrink, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange_shrink)) %>% arrange(Rank, abs_fc)
	    return(
	      datatable(df,
	        rownames=FALSE,
	        extensions = 'Buttons',
	        options=list(
	          columnDefs=list(list(visible=FALSE, targets=c(10, 11, 12))),
	          dom = 'lftBipr',
	          buttons = list(
	            list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_deseq2_results"), extension='.tsv', fieldBoundary='', fieldSeparator='	')
	           )
	          ),
	        container=sketch
	      ) %>% 
	      formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'log2FoldChange_shrink', 'lfcSE_shrink'), digits=4) %>%
	      formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	      formatStyle(c(6, 9), `border-left` = "solid 1px") %>%
	      formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	    )

	  } else {

	    rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(-abs(log2FoldChange))) %>% select(feature, Rank)
	    df = post_res %>% left_join(rank, by='feature') %>% select(transcript=feature, gene, baseMean, pvalue, padj, log2FoldChange, lfcSE, stat, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange)) %>% arrange(Rank, abs_fc)
	
	    return(
	      datatable(df,
	        rownames=FALSE,
	        extensions = 'Buttons',
	        options=list(
	          columnDefs=list(list(visible=FALSE, targets=c(8,9,10))),
	          dom = 'lftBipr',
	          buttons = list(
	            list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_deseq2_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	           )
	          )
	      ) %>% 
	      formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), digits=4) %>%
	      formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	      formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	    )

	  }
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

	  write.table(all, file=paste0("outputs/", prefix, "_all_deseq2_results.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	  write.table(sig, file=paste0("outputs/", prefix, "_sig_deseq2_results.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	}
	```
	'''
	))

def choose_modules(args):

	return(dedent(
	'''
	```{{r, choose_modules}}
	include_count_distribution = {}
	include_all2all = {}
	include_pca = {}
	include_batch_correction = {}
	include_deseq2 = {}
	apply_shrinkage = {}
	convert_names = {}
	include_ma = {}
	include_volcano = {} 
	include_heatmap = {}
	include_volcano_highlighted = {}
	include_ma_highlighted = {}
	```'''.format(args.include_distribution, args.include_all2all, args.include_pca, args.include_batch_correction, args.include_DESeq2, args.apply_shrinkage, args.convert_names, args.include_ma, args.include_volcano, args.include_heatmap, args.include_volcano_highlighted, args.include_ma_highlighted)))

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
	cat('In Eukaryotes only a subset of all genes are expressed in a given cell. Expression is therefore a bimodal distribution, with non-expressed genes having counts that result from experimental and biological noise. It is important to filter out the genes/transcripts that are not expressed before doing differential gene/transcript expression. You can decide which cutoff separates expressed vs non-expressed genes/transcripts by looking at the following histogram. The current count threshold is at least ', min_counts_per_event, ' raw reads (red line) in at least ', min_samples_per_event, ' samples.\\n')
	
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
	
	pca = run_pca(all_counts_filtered, transformation = transformation, write_transform=FALSE)
	print(pca_plot(pca, all_metadata, color_by='{}', shape_by='{}', fill_by='{}', alpha_by='{}', label_by='{}'))
	print(scree_plot(pca))

	if (include_batch_correction) {{
	  cat('\\n\\n#### Post batch correction\\n')
	  pca = run_pca(all_counts_batch_corrected, transformation = transformation, write_transform=FALSE)
	  print(pca_plot(pca, all_metadata, color_by='{}', shape_by='{}', fill_by='{}', alpha_by='{}', label_by='{}'))
	  print(scree_plot(pca))
	}}
	```'''.format(args.transformation, args.pca_color, args.pca_shape, args.pca_fill, args.pca_transparency, args.pca_label, args.pca_color, args.pca_shape, args.pca_fill, args.pca_transparency, args.pca_label)))

def run_deseq(treatment, control, highlighted, args, name, design, column):

	return(dedent(
	'''
	```{{r, DESeq_info, eval=include_deseq2, echo=FALSE, results='asis'}}
	cat('## DESeq Analysis {{ .tabset .tabset-pills }}\\n')
	cat('The goal of Differential gene expression analysis is to find genes or transcripts whose difference in expression, when accounting for the variance within condition, is higher than expected by chance.')
	```

	```{{r, run_DESeq, eval=include_deseq2, include=include_deseq2, warning=FALSE, message=FALSE}}
	design = {}
	fitType = "{}"
	padj_significance_cutoff = {}
	fc_significance_cutoff = {}
	padj_floor = {}
	fc_ceiling = {}
	num_labeled = {}
	input_mode = '{}'
	
	treatment_name = '{}'
	control_name = '{}'
	
	highlighted = c({})

	selected_metadata = selected_metadata %>% mutate({} = factor({}, levels=c(control_name, treatment_name)))

	use_batch_correction_in_DE = {}

	if (include_batch_correction & use_batch_correction_in_DE & input_mode == 'All') {{
	  deseq_count_input = all_counts_batch_corrected
	  deseq_metadata_input = all_metadata
	}} else if (include_batch_correction & use_batch_correction_in_DE & input_mode == 'Comparison-only') {{
	  deseq_count_input = selected_counts_batch_corrected
	  deseq_metadata_input = selected_metadata
	}} else if (input_mode == 'All') {{
	  deseq_count_input = all_counts_filtered
	  deseq_metadata_input = all_metadata 
	}} else {{
	  deseq_count_input = selected_counts_filtered
	  deseq_metadata_input = selected_metadata
	}}

	dds = DESeqDataSetFromMatrix(countData = deseq_count_input, colData = deseq_metadata_input, design = design)
	dds${} = relevel(dds${}, ref = control_name)
	res = DESeq(dds, fitType=fitType)
	results = as.data.frame(results(res, contrast=c('{}', treatment_name, control_name))) %>% mutate(padj = replace_na(padj, 1)) %>% tibble::rownames_to_column('feature')
	
	if (apply_shrinkage) {{
	  shrinkage_type = "{}"
	  shrink_results = as.data.frame(lfcShrink(res, coef=paste('{}', treatment_name, "vs", control_name, sep='_'), type=shrinkage_type)) %>%
	                   tibble::rownames_to_column('feature') %>%
	                   select(feature, log2FoldChange_shrink=log2FoldChange, lfcSE_shrink=lfcSE)
	  results = results %>% left_join(shrink_results, by='feature')
	}}
	
	post_res = post_processing(results, padj_significance_cutoff=padj_significance_cutoff, fc_significance_cutoff=fc_significance_cutoff, num_labeled=num_labeled, highlighted=highlighted, add_alias={}, fromType='{}',  toType='{}', org='{}', apply_shrinkage=apply_shrinkage)
	
	if (feature_type == 'transcript') {{
		transcript_to_gene = read.delim("{}{}") %>% select(feature=1, gene=2)
		post_res = post_res %>% left_join(transcript_to_gene, by='feature')
	}}

	deseq_results_table(post_res, apply_shrinkage, feature_type, '{}')
	```'''.format(design, args.fitType, args.padj_significance_cutoff, args.fc_significance_cutoff, args.padj_floor, args.fc_ceiling, args.num_labeled, args.input_mode, treatment, control, highlighted, column, column, args.use_bc_in_DE, column, column, column, args.shrinkage_type, column, args.convert_names, args.count_file_names, args.converted_names, args.org, args.input_prefix, args.counts, name)))

def ma_plot():
	return(dedent(
	'''
	```{r, ma_plot, echo=FALSE, eval=include_ma & include_deseq2, include=include_ma & include_deseq2, results='asis', warning=FALSE}
	if(include_ma_highlighted) {
	  cat('### MA {.tabset} \\n')
	  cat('#### All Detected \\n')
	} else {
	  cat('### MA\\n')
	}
	
	if(apply_shrinkage) {
	  ma_plot(post_res, Y='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	} else {
	  ma_plot(post_res, padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	}

	if(include_ma_highlighted) {
	  cat('\\n\\n#### Highlighted Only \\n')
	  if(apply_shrinkage) {
	    ma_plot(post_res %>% filter(Highlighted==TRUE), Y='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  } else {
	    ma_plot(post_res %>% filter(Highlighted==TRUE), padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  }
	}
	```'''))

def volcano_plot():
	return(dedent(
	'''
	```{r, volcano_plot, echo=FALSE, eval=include_volcano & include_deseq2, include=include_volcano & include_deseq2, results='asis', warning=FALSE}
	
	if(include_volcano_highlighted) {
	  cat('### Volcano {.tabset} \\n')
	  cat('#### All Detected \\n')
	} else {
	  cat('### Volcano\\n')
	}

	if(apply_shrinkage) {
	  volcano_plot(post_res, X='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	} else {
	  volcano_plot(post_res, padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	}

	if(include_volcano_highlighted) {
	  cat('\\n\\n#### Highlighted Only \\n')
	  if(apply_shrinkage) {
	    volcano_plot(post_res %>% filter(Highlighted==TRUE), X='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  } else {
	    volcano_plot(post_res %>% filter(Highlighted==TRUE), padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  }
	}
	```'''))

def heatmap():
	return(dedent(
	'''
	```{r, heatmap, echo=FALSE, eval=include_heatmap & include_deseq2, include=include_heatmap & include_deseq2, results='asis'}
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
	    cat("Cannot include the following samples because they have zero counts across all significant genes which breaks the normalization function:\\n")
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
	if (include_deseq2) {{
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

def normalization_report(args):
	output = []
	output.append(print_header("Normalization Report"))
	output.append(build_normalization_report(args))
	write_output("Normalization_Report", output)
	run_markdown("Normalization_Report")

def build_normalization_report(args):

	return(dedent(
	'''
	```{{r, report, echo=FALSE}}
	library(dplyr)
	library(tidyr)
	library(DESeq2)
	
	read_metadata = function(metadata_file) {{
	  metadata = read.delim(metadata_file, header=TRUE)
	  if (length(metadata) < 2) {{
	    metadata = read.delim(metadata_file, sep=',', header=TRUE)
	  }}
	  return(metadata)
	}}
	
	read_counts = function(count_file, samples, excluded_events) {{
	
	  counts = read.delim(count_file, check.names=FALSE)
	  if (length(counts) < 2) {{
	    counts = read.delim(count_file, check.names=FALSE, sep=',')
	  }}
	
	  counts = counts %>%
	           select(feature=1, all_of(samples)) %>%
	           filter(!(feature %in% excluded_events)) %>%
	           tibble::column_to_rownames('feature')
	  counts = as.matrix(counts)
	  mode(counts) = 'integer'
	  return(counts)
	}}
	
	check_samples = function(all_counts, total_count_min) {{
	  bad_samples = (data.frame(Counts=colSums(all_counts)) %>% tibble::rownames_to_column('Sample') %>% filter(Counts < total_count_min))$Sample
	  return(bad_samples)
	}}
	
	sample_filter = function(all_counts, bad_samples) {{
	  x = as.data.frame(all_counts) %>% select(-all_of(bad_samples))
	  return(as.matrix(x))
	}}
	
	filter_counts = function(input, min_counts_per_event, min_samples_per_event) {{
	  keep = rowSums(input >= min_counts_per_event) >= min_samples_per_event
	  return(input[keep,])
	}}
	
	run_pca = function(input, retx=TRUE, center=TRUE, scale=TRUE, transformation='Default', write_transform=FALSE) {{
	  if (transformation == 'None') {{
	    keep = subset(input, apply(input, 1, var, na.rm = TRUE) >  0)
	    return(prcomp(t(keep), retx = retx, center = center, scale. = scale))
	  }} else if (transformation == 'vst' | (transformation == 'Default' & ncol(input) > 50)) {{
	    transformed = vst(input)
	    if(write_transform) {{
	      write.table(as.data.frame(transformed) %>% tibble::rownames_to_column("feature"), file='outputs/transformed_counts_vst.tsv', quote=FALSE, sep='	', row.names=FALSE)
	    }}
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  }} else {{
	    transformed = rlog(input)
	    if(write_transform) {{
	      write.table(as.data.frame(transformed) %>% tibble::rownames_to_column("feature"), file='outputs/transformed_counts_rlog.tsv', quote=FALSE, sep='	', row.names=FALSE)
	    }}
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  }}
	}}
	
	excluded_events = c({})
	all_metadata = read_metadata("{}{}")
	all_counts = read_counts("{}{}", all_metadata$sample_name, excluded_events)

	min_counts_per_sample = {}
	min_counts_per_event = {}
	min_samples_per_event = {}
	
	transformation = '{}'
	
	bad_samples = check_samples(all_counts, min_counts_per_sample)
	all_counts = sample_filter(all_counts, bad_samples)
	all_counts_filtered = filter_counts(all_counts, min_counts_per_event, min_samples_per_event)
	
	pca = run_pca(all_counts_filtered, transformation = transformation, write_transform=TRUE)
	norm_dds = DESeqDataSetFromMatrix(countData = all_counts_filtered, colData = all_metadata, design = ~group)	
	norm = assay(normTransform(norm_dds))
	write.table(as.data.frame(norm) %>% tibble::rownames_to_column("feature"), file="outputs/normalized_feature_counts.tsv", sep='\\t', quote=FALSE, row.names = FALSE)
	write.table(as.data.frame(t(scale(t(norm)))) %>% tibble::rownames_to_column("feature"), file="outputs/expression_z_scores.tsv", sep='\\t', quote=FALSE, row.names=FALSE)
	sessionInfo()
	```'''.format(','.join('"%s"' % (i) for i in args.excluded_events), args.input_prefix, args.groups, args.input_prefix, args.counts, args.min_counts_per_sample, args.min_counts_per_event, args.min_samples_per_event, args.transformation)))

def write_output(name, output):

	with open('%s_des.Rmd' % (name), 'w') as out:
		out.write('\n'.join(output)) 

def write_params(args):
	output = []
	with open('outputs/deseq2_run_params.yaml', 'w') as out:
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

	cmd = "Rscript -e 'rmarkdown::render(\"%s_des.Rmd\", \"html_document\")'" % (name)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print('%s:' % name)
	print(out.decode())
	print('%s:' % name, file=sys.stderr)
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="prepare_DESeq2", description='', usage='%(prog)s [options]')
	
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-c', '--counts', required=True, help='Name of count file.', metavar='', dest='counts')
	input_args.add_argument('-g', '--groups', required=True, help='Name of groups file.', metavar='', dest='groups')
	input_args.add_argument('-x', '--comparisons', required=True, help='Name of comparisons file.', metavar='', dest='comparisons')
	input_args.add_argument('-p', '--prefix', default='', help='Path to input directory.', metavar='', dest='input_prefix')
	input_args.add_argument('--feature-type', choices=['gene','transcript'], default='gene', help='Feature type to run DESeq2 on.', metavar='', dest='feature_type')

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

	deseq2 = parser.add_argument_group('DESeq2 Settings')
	deseq2.add_argument('--include-DESeq2', choices=["TRUE", "FALSE"], default='TRUE', help='Run DESeq2', metavar='', dest='include_DESeq2')
	deseq2.add_argument('--design', type=str, default='~group', help='Design formula', metavar='', dest='design')
	deseq2.add_argument('--fitType', choices=["parametric","local","mean","glmGamPoi"], default='parametric', help='fitType used during DESeq2 analysis.', metavar='', dest='fitType')
	deseq2.add_argument('--input-mode', choices=["All", "Comparison-only"], default='All', help='Build DESeq object with all samlpes in group file or only samples from current comparison.', metavar='', dest='input_mode')
	deseq2.add_argument('--apply-shrinkage', choices=["TRUE", "FALSE"], default='TRUE', help='Apply shrinkage to fold change after DESeq2', metavar='', dest='apply_shrinkage')
	deseq2.add_argument('--use-batch-correction-in-DE', choices=["TRUE", "FALSE"], default='TRUE', help='Use batch corrected values in DESeq', metavar='', dest='use_bc_in_DE')
	deseq2.add_argument('--shrinkage-type', choices=["apeglm","ashr","normal"], default='apeglm', help='Shinkage algorithm to use.', metavar='', dest='shrinkage_type')
	deseq2.add_argument('--include-volcano', choices=["TRUE", "FALSE"], default='TRUE', help='Create volcano plots', metavar='', dest='include_volcano')
	deseq2.add_argument('--include-ma', choices=["TRUE", "FALSE"], default='TRUE', help='Create MA plots', metavar='', dest='include_ma')
	deseq2.add_argument('--include-heatmap', choices=["TRUE", "FALSE"], default='TRUE', help='Create heatmap', metavar='', dest='include_heatmap')

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