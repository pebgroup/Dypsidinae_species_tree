#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("endswith")
args = parser.parse_args()
ew = str(args.endswith)

occ = {}
for file in os.listdir("."):
	#if file.endswith("clean.fasta"):
	#if file.endswith("clean_noempty.fasta"):
	#if file.endswith("FNA"):
	if file.endswith(ew):
		for record in SeqIO.parse(file, "fasta"):
			if record.id in occ.keys():
				occ[record.id] += 1
			else:
				occ[record.id] = 1

with open("occupancy.txt", "w") as handle:
	for i,j in occ.items():
		print(i+";"+str(j), file=handle)
