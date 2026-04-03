#!/usr/bin/env python3

import argparse, os, subprocess, sys, shutil
from textwrap import dedent
from multiprocessing import Pool

def main(args):

	os.makedirs('GSEA/gmt/', exist_ok=True)
	os.makedirs('GSEA/outputs/', exist_ok=True) # Remove me when launched apps can see parent directory
	os.makedirs('outputs/', exist_ok=True)

	gmt_key = read_gmt_key(args.gmt_key)
	copy_gmt_files(args.gmts, gmt_key, args.gmt_source)

	p = Pool(args.threads)
	p.map(build_report, ((infile, gmt_key, args) for infile in args.input))
	p.close()
	p.join()

def build_report(arguments):

	infile, gmt_key, args = arguments

	if args.species == 'human':
		database = 'org.Hs.eg.db'
	elif args.species == 'mouse':
		database = 'org.Mm.eg.db'
	else:
		database = 'NA'

	if infile[-23:] == '_all_deseq2_results.tsv':
		name = os.path.basename(infile[:-23])
	elif infile[-26:] == '_all_limmaVoom_results.tsv':
		name = os.path.basename(infile[:-26])
	else:
		name = "GSEA_results"

	shutil.copy(infile, 'outputs/%s' % (os.path.basename(infile)))
	shutil.copy(infile, 'GSEA/outputs/%s' % (os.path.basename(infile))) # Remove me when launched apps can see parent directory

	output = []
	output.append(print_header(name))

	if database != 'NA':
		output.append(print_libraries())
		output.append(print_functions())
		output.append(set_parameters(args, database))
		output.append(read_data(os.path.basename(infile)))
		output.extend(loop_gmts(args.gmts, gmt_key, name, args.postfix))
		output.append(session_info())
	else:
		output.append("Only human and mouse are supported at this time.")
	
	write_output(name, output, args.postfix)
	run_markdown(name, args.postfix)

def copy_gmt_files(gmts, gmt_key, source):
	for gmt in gmts:
		shutil.copy('%s/%s/%s' % (source, gmt_key[gmt][1], gmt), 'GSEA/gmt/%s' % (gmt))

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
	library(dplyr)
	library(DT)
	library(fgsea)
	library(clusterProfiler)
	library(ggplot2)
	```'''))

def print_functions():
	
	return(dedent(
	'''
	```{r, load_functions, include=FALSE}
	run_gsea = function(pathways, geneList, minSize, maxSize, NES_threshold, padj_threshold, seed) {
	  res = fgsea(pathways=pathways, stats=geneList, eps=0.0, minSize=minSize, maxSize=maxSize) %>% 
	        arrange(-abs(NES)) %>%
	        mutate(Direction = case_when(NES < -NES_threshold & padj < padj_threshold ~ "Downregulated",
	                                     NES >  NES_threshold & padj < padj_threshold ~ "Upregulated",
	                                     TRUE ~ "No Change")) %>%
	        mutate(LeadingEdge = '')
	  
	  for (i in 1:nrow(res)) {
	    le = bitr(res[[i, 'leadingEdge']], fromType='ENTREZID', toType='SYMBOL', OrgDb=species)
	    res[[i, 'LeadingEdge']] = paste0(le[['SYMBOL']], collapse=", ")
	  }
	    
	  return(res %>% select(-leadingEdge))
	}
	
	write_results = function(df, name, prefix, postfix) {
	  write.table(df, paste0(file='outputs/', name, postfix, '_', prefix, '_gsea_results.tsv'), row.names = FALSE, quote=FALSE, sep='\t')
	}
	
	results_table = function(res, name, prefix) {
	  return(datatable(res %>% dplyr::select(-LeadingEdge),
	    rownames=FALSE,
	    extensions = 'Buttons',
	    options=list(
	      columnDefs=list(list(visible=FALSE, targets=c(7))),
	      dom = 'lftBipr',
	      buttons = list(
	        list(extend = 'csvHtml5', text='Download', filename = paste(name, prefix, "gsea_results", sep='_'), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	       )
	      )
	  ) %>% 
	  formatRound(columns=c("log2err", "ES", 'NES'), digits=4) %>%
	  formatSignif(columns=c("pval", "padj"), digits=4) %>%
	  formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	  )
	}
	
	enrichment_plots = function(res) {
	
	  sig_res = res %>% filter(Direction != 'No Change')
	
	  for (pathway in sig_res$pathway) {
	    cat(paste0('\\n#### ', pathway, '\\n'))
	    print(plotEnrichment(pathways[[pathway]], geneList) + labs(title=pathway))
	    cat('\\n')
	  }
	}
	```'''))

def set_parameters(args, database):

	return(dedent(
	'''
	```{{r, parameters, message=FALSE, include=FALSE}}
	species = '{}'
	event_column = '{}'
	fold_change_column = '{}'

	set.seed({})
	minSize = {}
	maxSize = {}
	NES_threshold = {}
	padj_threshold = {}
	```'''.format(database, args.event_column, args.fc_column, args.seed, args.minSize, args.maxSize, args.NES, args.pval)))


def read_data(infile):
	
	return(dedent(
	'''
	```{{r, read_data}}
	raw = read.delim('outputs/{}') # Switch to '../outputs/<>' this line when launched apps can read parent directory
	
	if (fold_change_column == 'log2FoldChange_shrink' & !('log2FoldChange_shrink' %in% colnames(raw)) & 'log2FoldChange' %in% colnames(raw)) {{
	    fold_change_column = 'log2FoldChange'
	    cat('log2FoldChange_shrink is not a column in input data. Using log2FoldChange instead')
	}}
	```
	
	```{{r, prepareList, message=FALSE, include=FALSE}}
	
	mapped_names = data.frame(bitr(raw[[event_column]], fromType='SYMBOL', toType='ENTREZID', OrgDb=species))
	df = raw %>% left_join(mapped_names, by=setNames("SYMBOL", event_column)) %>% filter(!is.na(ENTREZID)) %>% dplyr::select(ENTREZID, all_of(fold_change_column))
	geneList = df[,2]
	names(geneList) = as.character(df[,1])
	geneList = sort(geneList, decreasing = TRUE)
	```'''.format(infile)))

def read_gmt_key(input_file):
	gmt_key = {}
	with open(input_file) as infile:
		infile.readline()
		for line in infile:
			cur = line.rstrip().split('\t')
			gmt_key[cur[2]] = ((cur[0], cur[1], cur[3]))
	return(gmt_key)

def loop_gmts(gmts, gmt_key, name, postfix):
	output = []
	output.append("\n# {.tabset .tabset-pills}")
	for gmt in gmts:
		prefix = gmt_key[gmt][0][:gmt_key[gmt][0].index(':')]
		output.append('\n## %s' % gmt_key[gmt][0])
		output.append(run_gsea('gmt/%s' % (gmt), name, prefix, postfix))
	return(output)

def run_gsea(gmt, name, prefix, postfix):
	
	return(dedent(
	'''
	```{{r, message=FALSE, warning=FALSE, results='asis'}}
	pathways = gmtPathways("{}")
	res = run_gsea(pathways, geneList, minSize, maxSize, NES_threshold, padj_threshold, seed)
	write_results(res, "{}", "{}", "{}")
	cat('### Pathway Enrichment')
	results_table(res, "{}", "{}")
	cat('### Enrichment Plots {{.tabset .tabset-dropdown}}')
	enrichment_plots(res)
	```'''.format(gmt, name, prefix, postfix, name, prefix)))

def session_info():

	return(dedent(
	'''
	```{r, session_info, echo=FALSE, results='asis'}
	cat('# Session Info\\n')
	sessionInfo()
	```'''))

def write_output(name, output, postfix):

	with open('GSEA/%s%s.Rmd' % (name, postfix), 'w') as out:
		out.write('\n'.join(output))

def run_markdown(name, postfix):

	cmd = "Rscript -e 'rmarkdown::render(\"GSEA/%s%s.Rmd\", \"html_document\")'" % (name, postfix)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print('%s:' % name)
	print(out.decode())
	print('%s:' % name, file=sys.stderr)
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="prepare_GSEA", description='', usage='%(prog)s [options]')
	
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-i', '--input', nargs='+', required=True, help='Input file.', metavar='', dest='input')
	input_args.add_argument('-s', '--species', default='human', help='Species. Only human and mouse are currently supported.', metavar='', dest='species')
	input_args.add_argument('-e', '--event-column', default='gene', help='Event column', metavar='', dest='event_column')
	input_args.add_argument('-f', '--fold-change-column', default='log2FoldChange_shrink', help='Log2FoldChange column', metavar='', dest='fc_column')
	
	gmt_args = parser.add_argument_group('Gene Set Parameters')
	gmt_args.add_argument('-k', '--GMT-key', default='data/gmt_key.txt', help='Preloaded GMT file key.', metavar='', dest='gmt_key')
	gmt_args.add_argument('--GMT-source', default='/data/', help='Where gmt files are stored if not running from container.', metavar='', dest='gmt_source')
	gmt_args.add_argument('-g', '--GMT-list', nargs='+', default=['h.all.v2023.1.Hs.entrez.gmt'], help='GMTs to run.', metavar='', dest='gmts')
	gmt_args.add_argument('--minSize', type=int, default=1, help='Minimal gene set size.', metavar='', dest='minSize')
	gmt_args.add_argument('--maxSize', type=int, default=10000, help='Maximal gene set size.', metavar='', dest='maxSize')
	
	significance_args = parser.add_argument_group('Significance Parameters')
	significance_args.add_argument('-n', '--NES', type=float, default=1, help='Normalized Enrichment Score significance threshold.', metavar='', dest='NES')
	significance_args.add_argument('-p', '--pvalue', type=float, default=.05, help='Adjusted pvalue significance threshold', metavar='', dest='pval')

	other_args = parser.add_argument_group('Other Parameters')
	other_args.add_argument('--seed', type=int, default=1, help='Seed.', metavar='', dest='seed')
	other_args.add_argument('--postfix', default='', help='Postfix to append to filenames.', metavar='', dest='postfix')
	other_args.add_argument('-t', '--threads', type=int, default=1, help='Available threads.', metavar='', dest='threads')

	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)