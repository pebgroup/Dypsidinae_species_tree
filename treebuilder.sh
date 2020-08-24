#!/bin/bash

~/scripts/dypsidinae/partitioner.py --smoother 10

for f in *_clean.fasta
do
	~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s $f -T AUTO -ntmax 16 -p ${f/clean.fasta}part.txt -B 1000
	mv ${f/clean.fasta}part.txt.treefile genetrees/${f/clean.fasta}part.txt.tre
	mv ${f/clean.fasta}part.txt* genetrees
	mv ${f/_clean.fasta}.fasta done
	rm $f
done		