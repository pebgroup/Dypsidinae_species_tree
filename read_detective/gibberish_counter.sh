#!/bin/bash

# script for evaluating the reads mapped to each gene in a Hybpiper output for overrepresentation
# of AG or CT motifs. 

# Wolf Eiserhardt, 4 April 2020

# creating a list of the target loci as represented by directories
ls -d */ > dirnames.txt
sed -i'.old' -e 's./..g' dirnames.txt

while read name; do 
	if test -f "$name/${name}_spades/split_input/${name}_interleaved_1.fasta"; then
		~/scripts/dypsidinae/read_detective/gibberish_counter.py $name
	fi
done < dirnames.txt