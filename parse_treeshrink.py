#!/usr/bin/python3

import os, subprocess

genes = [x[0] for x in os.walk(".")][1:]

genes = [x[2:] for x in genes]

with open("overview.txt","w") as outfile:
	for g in genes: 
		with open(g+"/input_shrunk_RS_0.05.txt", "r") as hdl:
			#if len(hdl.readline()) > 0:
			ln = hdl.readline()
			if len(ln) > 0:
				print("edit")
				if g == "32e":
					cmd = "cp /data_vol/wolf/Dypsis/alignments_edited/done/reduced_"+g+"_aligned.fasta /data_vol/wolf/Dypsis/alignments_for_editing2"			
				else:
					cmd = "cp /data_vol/wolf/Dypsis/alignments_edited/done/reduced_"+g+"_aligned_noempty.fasta /data_vol/wolf/Dypsis/alignments_for_editing2"
				subprocess.call(cmd, shell=True)
				print(g+"\t"+ln, file=outfile)
			else:
				print("don't edit")
				if g == "32e":
					cmd = "cp /data_vol/wolf/Dypsis/alignments_edited/done/reduced_"+g+"_aligned.fasta /data_vol/wolf/Dypsis/alignments_edited2"			
				else:
					cmd = "cp /data_vol/wolf/Dypsis/alignments_edited/done/reduced_"+g+"_aligned_noempty.fasta /data_vol/wolf/Dypsis/alignments_edited2"
				subprocess.call(cmd, shell=True)