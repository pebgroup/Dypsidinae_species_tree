#!/usr/bin/python3

# Script for splitting fasta files that each contain all contigs from a given sample into 
# fasta files that each contain all samples from a given gene. This requires that all 
# sequences in the original files are named "sample_gene". 
# 
# Requires a text file "filelist.txt" with the names of the input fasta files. I generated
# this using ls *trimmed.fasta > filelist.txt. 
#
# Prints some gene-level stats (gene name, median lenght, number of sequences, number of 
# sequences kept, threshold applied) to stdout. 
# 
# Recommended usage: samples2genes.py > outstats.csv
#
# Wolf Eiserhardt 17.4.2020

from Bio import SeqIO
from statistics import median 

# get list of fasta files to split into genes
files2parse = []
with open("filelist.txt") as files:
	for line in files:
		files2parse.append(line.strip())

# gather sequences from those files into gene dictionary_
# key: gene name
# item: list of seq records
genes = {}
for file2parse in files2parse: 
	for record in SeqIO.parse(file2parse, "fasta"):
		sample_gene = record.id.split("_")
		record.id = sample_gene[0]
		if sample_gene[1] in genes.keys():
			genes[sample_gene[1]].append(record)
		else: 
			genes[sample_gene[1]] = [record]

# print header
print("gene;median_length;sequences;sequences_kept;threshold")

# length-filter sequence sets and write to fasta
prc = .2 #minimum percentage of median length
len_thres = 150 #minimum absolute sequence length
for gene, seq_set in genes.items():
	# get median length
	lens = [len(s.seq) for s in seq_set]
	# get final threshold (combination of percentage and min length)
	thres = max(prc*median(lens), len_thres)
	# select sequences over threshold
	seq_set_keep = [seq for length, seq in zip(lens, seq_set) if length>thres]
	# write out if more than 3 sequences
	if len(seq_set_keep) > 3:
		with open('../seq_sets2/'+gene+'.FNA', "w") as outfile:
			SeqIO.write(seq_set_keep, outfile, "fasta")
	print(gene+';'+str(round(median(lens)))+';'+str(len(seq_set))+';'+str(len(seq_set_keep))+';'+str(round(thres)))