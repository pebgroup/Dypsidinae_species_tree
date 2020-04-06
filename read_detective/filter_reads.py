#!/usr/bin/python3

# Script that extracts reads with Phred information from a FASTQ file based on their 
# occurrence (by ID) in another read file without Phred information (FASTA format).

# Called by gather_reads.sh
# Output: <sample_name>_R*.fastq, separately for forward (1) and reverse (2) reads. 

# Wolf Eiserhardt, 3.4.2020

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("sample")
args = parser.parse_args()

# Get sample ID for file naming purposes
sample = args.sample

# extract the IDs of the reads in the FASTA file
tgts = []
for record in SeqIO.parse("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R1.fasta", "fasta"):#
	tgts.append(record.id)
#print("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R2.fasta")

print(len(tgts))
incr = int(len(tgts)/100)

# create an empty sequence list
fltrd = []

i = 0
# loop through the FASTQ files, first paired then unpaired reads, and append any reads 
# to the sequence list that have an ID that also occurred in the FASTA file. 
for record in SeqIO.parse("/data_vol/wolf/Dypsis/trimmed/"+sample+"_clean-READ1.fastq", "fastq"):
	if record.id in tgts:
		if (i % incr) == 0:
			print(str(i/incr)+"% finished")
		i += 1
		fltrd.append(record)
#print("/data_vol/wolf/Dypsis/trimmed/"+sample+"_clean-READ2.fastq")

for record in SeqIO.parse("/data_vol/wolf/Dypsis/trimmed/"+sample+"_clean-READ1-single.fastq", "fastq"):
	if record.id in tgts:
		if (i % incr) == 0:
			print(str(i/incr)+"% finished")
		i += 1
		fltrd.append(record)
#print("/data_vol/wolf/Dypsis/trimmed/"+sample+"_clean-READ2-single.fastq")
		
# write the sequence list to output file. 
with open("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R1.fastq", "w") as outfile:
 	SeqIO.write(fltrd, outfile, "fastq")
#print("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R2.fastq")
 	
print(len(fltrd))

