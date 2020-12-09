#!/usr/bin/python3

import os, subprocess


newgenes = [i.split(".")[0] for i in os.listdir("seq_sets_newspp")]


for fn in os.listdir("alignments_edited2/done"): 
	if fn.split("_")[1] in newgenes:
		cmd = "mafft --add seq_sets_newspp/"+fn.split("_")[1]+".FNA --keeplength alignments_edited2/done/"+fn+" > alignments_edited2_newspp/"+fn 
		subprocess.call(cmd, shell=True)