#!/usr/bin/python3

# Script for counting the percentage of reads in a fasta file that contain >80% AG or CT.

# Wolf Eiserhardt, 4 April 2020

import argparse, re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("gene")
args = parser.parse_args()
gene = str(args.gene)

baduns_1 = 0
i_1 = 0
for record in SeqIO.parse(gene+"/"+gene+"_spades/split_input/"+gene+"_interleaved_1.fasta", "fasta"):
	i_1+=1
	AGcont = (len(re.findall("A",str(record.seq)))+len(re.findall("G",str(record.seq))))/len(record.seq)
	CTcont = (len(re.findall("C",str(record.seq)))+len(re.findall("T",str(record.seq))))/len(record.seq)
	if AGcont >= 0.8 or CTcont >= 0.8:
		baduns_1 += 1

baduns_2 = 0
i_2 = 0
for record in SeqIO.parse(gene+"/"+gene+"_spades/split_input/"+gene+"_interleaved_2.fasta", "fasta"):
	i_2+=1
	AGcont = (len(re.findall("A",str(record.seq)))+len(re.findall("G",str(record.seq))))/len(record.seq)
	CTcont = (len(re.findall("C",str(record.seq)))+len(re.findall("T",str(record.seq))))/len(record.seq)
	if AGcont >= 0.8 or CTcont >= 0.8:
		baduns_2 += 1
	
print(gene+";"+str(round(baduns_1/i_1,3))+";"+str(round(baduns_2/i_2,3)))