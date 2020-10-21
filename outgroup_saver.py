#!/usr/bin/python3

import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# Get all subdirectories in the current working directory. these are the loci recovered by hybpiper
loci = next(os.walk("."))[1]

for l in loci: 
	records = []
	for record in SeqIO.parse(l+"/input_shrunk0.05.fasta", "fasta"):
		records.append(record.id)
	if "1013" not in records:
		# recover the ou outgroup from the original alignment
		for record in SeqIO.parse(l+"/input.fasta", "fasta"):
			if record.id == "1013": 
				outgroup = record
		# load the shrunk alignment:
		aln = []
		for record in SeqIO.parse(l+"/input_shrunk0.05.fasta", "fasta"):
			aln.append(record)
		# add outgroup
		aln.append(outgroup)
		# de-gap
		aln_degapped = []
		for record in aln:
			record.seq = Seq(str(record.seq).replace("-",""))
			aln_degapped.append(record)
		# write complete alignment to file
		with open('../seq_sets2_shrunk/'+l+'_shrunk.fasta', "w") as outfile:
			SeqIO.write(aln_degapped, outfile, "fasta")
	else:
		aln_degapped = []
		for record in SeqIO.parse(l+"/input_shrunk0.05.fasta", "fasta"):
			record.seq = Seq(str(record.seq).replace("-",""))
			aln_degapped.append(record)
		with open('../seq_sets2_shrunk/'+l+'_shrunk.fasta', "w") as outfile:
			SeqIO.write(aln_degapped, outfile, "fasta")
		
