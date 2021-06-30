#!/usr/bin/env python

import glob
import csv
import argparse

parser = argparse.ArgumentParser(description="Write tab-separated list of fastq files to use in STAR")
parser.add_argument('--inPath', required=True, help='Path to directory with fastq files', action='store')
parser.add_argument('--output', required=True, help='TSV file to be created', action='store')

args=parser.parse_args()

indir = args.inPath
output = args.output

fastq_list = glob.glob(indir + '/*.fastq')

print(fastq_list)


results = []
for fastq in fastq_list:
	sample = fastq.split('.')[0].split('/')[-1]
	combined = [fastq, '-', sample]
	results.append(combined)

print(results)

with open(output, 'wt') as out_file:
	tsv_writer = csv.writer(out_file, delimiter='\t')
	tsv_writer.writerows(results)
