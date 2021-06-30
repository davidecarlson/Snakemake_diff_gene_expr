#!/usr/bin/env python

# take ENSEMBL IDs in deseq2 results add add gene names and symbols, writing out a new table

import argparse
import mygene
import pandas
import multiprocessing

#
parser = argparse.ArgumentParser(description="Take DESeq2 results and get a new table that includes the gene name and symbol for each Ensembl ID")
parser.add_argument('--input', required=True, help='DEseq2 results file to be processed', action='store')
parser.add_argument('--output', required=True, help='Output file to be created', action='store')
parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')

args=parser.parse_args()

deseqs_file = args.input
outfile = args.output
threads = args.threads

# initialize mygene functionality
mg = mygene.MyGeneInfo()

# create function for looking up gene symbol and name from ensembl ID
def lookup_gene(id):
	results = mg.getgenes(id, fields=['symbol', 'name'], as_dataframe = True)
	#print(results)
# set results to "n/a" if they can't be found
	if 'notfound' in results.columns or not 'name' in results.columns or not 'symbol' in results.columns:
		print("There were no results found for gene", id)
		gene_name = "n/a"
		gene_symbol = "n/a"
	elif not 'notfound' in results.columns:
		print("Results found for gene", id)
		gene_name = '\t'.join(results['name'])
		gene_symbol = '\t'.join(results['symbol'])
	return [id, gene_name, gene_symbol]


# read in results and redo the column names
deseq = pandas.read_csv(deseqs_file, na_filter = False)
deseq.columns = ["ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# strip the version number form the ensembl IDs
split_ids = deseq["ID"].str.split(".", n = 1, expand = True)
deseq["Ensembl_ID"] = split_ids[0]

deseq.drop(columns =["ID"], inplace = True)

# create muliprocessing pool to look up gene ids in parallel

p = multiprocessing.Pool(processes=threads)

# use list comprehension to get results
results = [p.apply_async(lookup_gene, [gene])
	for gene in deseq["Ensembl_ID"]]

output = [p.get() for p in results]

p.close()
p.join() # Wait for all child processes to close.


# split the results into 3 lists
ids, names, symbols = map(list, zip(*output))

# check that the output IDs and sort order from multiprocessing matches the input IDs
if list(deseq["Ensembl_ID"]) == ids:
	print("Results look good, preparing output")
	deseq['name'] = names
	deseq['symbol'] = symbols

	# reorder the dataframe
	deseq = deseq[['Ensembl_ID', 'name', 'symbol', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']]
	#print(deseq)

	# write df out to new csv file
	deseq.to_csv(outfile, index = False, sep = "\t")
else:
	print("Error: something went wrong!")
