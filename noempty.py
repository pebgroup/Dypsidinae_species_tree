#!/usr/bin/python3

import argparse, re
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("alignment")
args = parser.parse_args()
alignment = str(args.alignment)

goodseqs = []
badseqs_count = 0
for record in SeqIO.parse(alignment, "fasta"):
	if(True in [i in record.seq for i in ["a","g","c","t","A","G","C","T"]]):
		goodseqs.append(record)
	else:
		badseqs_count += 1

with open(re.sub('\.fasta$', '', alignment)+"_noempty.fasta", "w") as output_handle:
    SeqIO.write(goodseqs, output_handle, "fasta")		

if badseqs_count > 0:
	print(re.sub('\.fasta$', '', alignment)+"_noempty.fasta has "+str(badseqs_count)+" empty sequences removed")