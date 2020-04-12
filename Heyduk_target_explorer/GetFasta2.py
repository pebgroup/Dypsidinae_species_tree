###########################################################################################################
# Script to get only some sequences from a fasta depending on a list of IDs
# Version 1 
# 1 September 2019
# Originally written by Sidonie BELLOT (s.bellot@kew.org)
# Use and modify as you wish, but please don't hesitate to give feedback!
###########################################################################################################

import sys
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter

# Run the script as such:
# python GetFasta2.py All_genes.fasta name_list.txt output.fasta

# design input and output

queries = sys.argv[1] # fasta file with all the sequences
names = sys.argv[2] # list of names of the sequences we want to keep
outfile = sys.argv[3] # output file containing these sequences



# make a list of the names

GOOD = []

handle = open(names, "r")
lines=handle.readlines()
z = len(lines)
print z

for l in lines:
   name="_" + str(l.split("\n")[0]) + "_"
#   name="_" + str(l.split("\r\n")[0]) + "_"
   GOOD.append(name)
handle.close()
print GOOD

# filter the sequences based on that list (customize depending on your input!)
x = 0
handle = open(queries)  
for seq_record in SeqIO.parse(handle, "fasta"):
   ID = seq_record.id
   print ID
#   if "LOC" in ID:
   if "_" in ID:
#      ID2 = "_" + str(ID.split("LOC")[1].split("\n")[0]) + "_"
      ID2 = "_" + ID + "_"
      print ID2
      if ID2 in GOOD:
         x = x + 1
         with open(outfile, 'a') as fo:
            SeqIO.write(seq_record, fo, "fasta")

handle.close()

print x


