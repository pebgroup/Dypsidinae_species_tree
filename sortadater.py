#!/usr/bin/python3

import dendropy, subprocess

taxa = []
with open("gg", "r") as infile: 
	i = 1
	for line in infile:
		if i > 1:
			line = line.strip()
			LINE = line.split(" ")
			#cmd = "cp genetrees/"+LINE[0]+" trees2date"
			cmd = "cp ../length_filter/"+"_".join(LINE[0].split("_")[:-1])+".fasta alignments"
			subprocess.call(cmd, shell=True)
			cmd = "cp ../length_filter/iqtree/"+"_".join(LINE[0].split("_")[:-1])+"_part.txt alignments"
			subprocess.call(cmd, shell=True)
			tree = dendropy.Tree.get(path="genetrees/"+LINE[0], schema="newick")
			tips = [t.label for t in tree.taxon_namespace]
			taxa = list(set(taxa + tips))
			if len(taxa) == 155:
				break
		i += 1