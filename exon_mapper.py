#!/usr/bin/python3

import os, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO

# import recovery statistics
df = pd.read_csv("../assembly/test_seq_lengths.txt", sep="\t")
df = df.drop(0)
df = df[~df['Species'].isin(['0188','0189','0190'])]


# for each aligned locus, identify the two samples with the highest recovery stats
for fn in os.listdir():
	locus = fn.split("_")[1]
	if locus in df.columns:
		df_red = df[['Species',locus]].sort_values(locus, ascending=False)
		df_red = df_red.reset_index()
		# first sample
		sp = df_red["Species"][0]
		exon1 = list(SeqIO.parse("../assembly/"+sp+"/"+locus+"/"+sp+"/sequences/FNA/"+locus+".FNA", "fasta"))[0]
		exon1.id = "exon1"
		# second sample
		sp = df_red["Species"][1]
		exon2 = list(SeqIO.parse("../assembly/"+sp+"/"+locus+"/"+sp+"/sequences/FNA/"+locus+".FNA", "fasta"))[0]
		exon2.id = "exon2"
		with open("temp.fasta", "w") as output_handle:
			SeqIO.write([exon1, exon2], output_handle, "fasta")
		subprocess.call("mafft --add temp.fasta "+fn+" > ../alignments_exon/"+fn ,shell=True)
		subprocess.call("rm temp.fasta", shell=True)
	else:
		print(locus+" not there!!!!")
