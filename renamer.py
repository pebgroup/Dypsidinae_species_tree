#!/usr/bin/python3

import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("mapping")
parser.add_argument("infile")
parser.add_argument("outfile")
parser.add_argument('--bs', default=0, help='1 = tree with bootstrap values')

args = parser.parse_args()

with open(args.infile) as f:
	tree = f.readline().strip()
	
with open(args.mapping, mode='r', encoding='utf-8-sig') as f: 
	for line in f:
		line = line.strip()
		LINE = line.split(";")
		searchterm = LINE[0]
		if args.bs == 0:
			tree = re.sub(rf'\({searchterm},','('+LINE[1]+',',tree)
			tree = re.sub(rf',{searchterm}\)',','+LINE[1]+')',tree)
		else: 
			tree = re.sub(rf'\({searchterm}:','('+LINE[1]+':',tree)
			tree = re.sub(rf',{searchterm}:',','+LINE[1]+':',tree)

with open(args.outfile, "w") as f:
	print(tree,file=f)