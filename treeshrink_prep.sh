#!/bin/bash

#for optrimal

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

# for alignments_trimmed

# for f in *.treefile
# do 
# 	mkdir ../treeshrink_gt0.5/${f/_trimmed_noempty.fasta.treefile}
# 	cp $f ../treeshrink_gt0.5/${f/_trimmed_noempty.fasta.treefile}/input.tre
# 	cp ${f/.treefile} ../treeshrink_gt0.5/${f/_trimmed_noempty.fasta.treefile}/input.fasta
# 	cd ../treeshrink_gt0.5/${f/_trimmed_noempty.fasta.treefile}
# 	sed -i'.old' -e "s/-${f/_trimmed_noempty.fasta.treefile}//g" input.tre
# 	sed -i'.old' -e "s/-${f/_trimmed_noempty.fasta.treefile}//g" input.fasta
# 	cd ../../alignments_trimmed
# done
