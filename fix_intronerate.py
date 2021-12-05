#!/usr/bin/python3

import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

for file in os.listdir("."):
	tgt = "../final_tree_nofilter_fix_intronerate/iqtree/reduced_"+file.split(".")[0]+"_aligned_noempty.fasta" #-out.fasta"
	#print(tgt)
	if os.path.exists(tgt):
		# generate list of existing sequences in the edited alignment
		tgt_ids = []
		for record in SeqIO.parse(tgt, "fasta"):
			tgt_ids.append(record.id)
		#print(tgt_ids)
		# check which of the new ones are not in there yet, and create list of them
		toadd = []
		for record in SeqIO.parse(file, "fasta"):
			if record.id not in tgt_ids:
				#print(record.id)
				toadd.append(record)
		#print(toadd)
		if len(toadd) > 0:
			with open('toadd/'+file, "w") as outfile:
				SeqIO.write(toadd, outfile, "fasta")
			cmd = "mafft --add toadd/"+file+" --keeplength "+tgt+" > added/reduced_"+file.split(".")[0]+"_aligned_noempty.fasta"#"-out.fasta"
			subprocess.call(cmd, shell=True)

	
	
	