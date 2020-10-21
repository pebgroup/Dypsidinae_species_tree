#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
import os

occ = {}
for file in os.listdir("."):
	if file.endswith("clean_noempty.fasta"):
		for record in SeqIO.parse(file, "fasta"):
			if record.id in occ.keys():
				occ[record.id] += 1
			else:
				occ[record.id] = 1

with open("occupancy.txt", "w") as handle:
	for i,j in occ.items():
		print(i+";"+str(j), file=handle)
