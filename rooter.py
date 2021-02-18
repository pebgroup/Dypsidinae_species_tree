#!/usr/bin/python3

import argparse, dendropy, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("treefile")
args = parser.parse_args()
treefile = str(args.treefile)

tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

if '1011' in tips and '1012' in tips:
	cmd = "pxrr -t "+treefile+" -g 1011,1012 -o temp.tre"
elif '1011' in tips:
	cmd = "pxrr -t "+treefile+" -g 1011 -o temp.tre"
elif '1012' in tips:
	cmd = "pxrr -t "+treefile+" -g 1012 -o temp.tre"
#if '1013' in tips:
#	cmd = "pxrr -t "+treefile+" -g 1013 -o temp.tre"
else:
	raise ValueError('Outgroup not in tree')

subprocess.call(cmd, shell=True)