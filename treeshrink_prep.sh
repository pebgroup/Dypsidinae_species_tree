#!/bin/bash

for f in *.treefile
do 
	mkdir ../treeshrink/${f/_aligned_noempty.fasta.treefile}
	cp $f ../treeshrink/${f/_aligned_noempty.fasta.treefile}/input.tre
	cp ${f/.treefile} ../treeshrink/${f/_aligned_noempty.fasta.treefile}/input.fasta
	cd ../treeshrink/${f/_aligned_noempty.fasta.treefile}
	sed -i'.old' -e "s/-${f/_aligned_noempty.fasta.treefile}//g" input.tre
	sed -i'.old' -e "s/-${f/_aligned_noempty.fasta.treefile}//g" input.fasta
	cd ../../optrimal
done
