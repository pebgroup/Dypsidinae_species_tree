#!/bin/bash

~/scripts/dypsidinae/partitioner.py --smoother 10

for f in *_clean.fasta
do
	~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s $f -T AUTO -ntmax 8 -p ${f/clean.fasta}part.txt -B 1000
	mv ${f/clean.fasta}part.txt.treefile genetrees/${f/clean.fasta}part.txt.tre
	mv ${f/clean.fasta}part.txt* genetrees
	mv ${f/_clean.fasta}.fasta done
	rm $f
done		

cd genetrees
rm -f ../../speciestree/genetrees.tre

for f in *.tre
do 
	~/scripts/dypsidinae/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> ../../speciestree/genetrees.tre 
	rm temp.tre
done

cd ../../speciestree

rm -f astral*
java -jar ~/software/Astral/astral.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree.tre astral_tree_renamed.tre

~/software/QuartetScores -o astral_tree_QS.tre -e genetrees.tre -r astral_tree.tre -v
sed astral_tree_QS.tre -i'.old' -e 's/[0-9]\.*[0-9]*\(:[0-9]\.*[0-9]*\)\[qp-ic:-*[0-9]\.[0-9]*;lq-ic:-*[0-9]\.[0-9]*;eqp-ic:\(-*[0-9]\.[0-9]*\)\]/\2\1/g'
sed astral_tree_QS.tre -i'.old' -e 's/\[eqp-ic:-*[0-9]\.*[0-9]*\]//g'
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree_QS.tre astral_tree_QS_renamed.tre --bs 1