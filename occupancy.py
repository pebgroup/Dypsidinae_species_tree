#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
import os

occ = {}
for file in os.listdir("."):
	if file.endswith("lf.fasta"):
		for record in SeqIO.parse(file, "fasta"):
			if record.id in occ.keys():
				occ[record.id] += 1
			else:
				occ[record.id] = 1

#with open("occupancy.txt", "w") as handle:
#	for i,j in occ.items():
#		print(i+";"+str(j), file=handle)

drop = []
for i,j in occ.items():
	if j < 20: 
		drop.append(i)
		
for file in os.listdir("."):
	if file.endswith("lf.fasta"):
		goodseq = []
		for record in SeqIO.parse(file, "fasta"):
			if record.id not in drop:
				goodseq.append(record)
		with open(file.replace("_lf.fasta","_lf_exl.fasta"), "w") as output_handle:
    			SeqIO.write(goodseq, output_handle, "fasta")
