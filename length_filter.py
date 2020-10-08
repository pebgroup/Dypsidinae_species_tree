#!/usr/bin/python3

import argparse, re
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("alignment")
args = parser.parse_args()
alignment = str(args.alignment)

drp = []
for record in SeqIO.parse(alignment, "fasta"):
	sqnc = str(record.seq)
	sqnc_nonull = sqnc.replace("-","")
	if len(sqnc_nonull)/len(sqnc) < 0.5:
		drp.append(record.id)

goodseqs = []
for record in SeqIO.parse(alignment.replace("_clean_70",""), "fasta"):
	if record.id not in drp:
		goodseqs.append(record)
		
with open(alignment.replace("_clean_70","_lf"), "w") as output_handle:
    SeqIO.write(goodseqs, output_handle, "fasta")

print(alignment+";"+str(len(drp))+";"+";".join(drp))