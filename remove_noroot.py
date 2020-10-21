#!/usr/bin/python3

import argparse, dendropy, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("treefile")
args = parser.parse_args()
treefile = str(args.treefile)

tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

if '1013' not in tips:
	cmd = "mv "+treefile+" noroot"
	subprocess.call(cmd, shell=True)