#!/usr/bin/python3

from amas import AMAS

with open("filelist.txt") as f:
	infiles = [line.strip() for line in f]

multi_meta_aln = AMAS.MetaAlignment(in_files=infiles, data_type="dna", in_format="fasta", cores=2)

summaries = multi_meta_aln.get_summaries()

with open("alnstats.csv", "w") as outfile:
	print(";".join(summaries[0]), file=outfile)
	for line in summaries[1]:
		print(";".join(line), file=outfile)
	