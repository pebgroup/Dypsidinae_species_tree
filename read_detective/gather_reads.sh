#!/bin/bash

# Script for extracting and combining all reads that HybPiper has mapped to the given
# targets. In the HybPiper output folder, these are in a set of directories (named after
# the target loci) within a folder <locus name>/<locus_name>_spades/split_input/. Here,
# the reads are stored as fasta files, i.e. not containing Phred information. This script
# extracts those reads and first combines them in one fasta file at the root of the
# HybPiper output directory; then it acccesses the original (trimmed) reads in FASTQ
# format, and extract those that are present in the previously generated fasta file. The
# result is saved in a FASTQ file at the root of the HybPiper output directory, named
# <sample_name>_R*.fastq, separately for forward (1) and reverse (2) reads. 

# This script needs to be run within the HybPiper output directory. 

# Wolf Eiserhardt, 3.4.2020

# name of the sample being processed, for file naming purposes
sample=0152

# creating a list of the target loci as represented by directories
ls -d */ > dirnames.txt
sed -i'.old' -e 's./..g' dirnames.txt

# checking if the output fasta files already exist
if test -f "${sample}_R1.fasta"; then
	echo "${sample}_R1.fasta already exists. Please delete before you execute this script."
	exit 1
fi
if test -f "${sample}_R2.fasta"; then
	echo "${sample}_R2.fasta already exists. Please delete before you execute this script."
	exit 1
fi

# for each locus, first check if the relevant read files exist, then add them to the 
# joint read file. 
while read name; do 
	if test -f "$name/${name}_spades/split_input/${name}_interleaved_2.fasta"; then
		cat "$name/${name}_spades/split_input/${name}_interleaved_2.fasta" >> "${sample}_R2.fasta"
	fi
	if test -f "$name/${name}_spades/split_input/${name}_interleaved_1.fasta"; then
		cat "$name/${name}_spades/split_input/${name}_interleaved_1.fasta" >> "${sample}_R1.fasta"
	fi
done < dirnames.txt

# run the Python3 script that extracts the reads with Phred info from the trimmed read 
# files. 
python3 ~/scripts/dypsidinae/filter_reads.py $sample
