#!/usr/bin/python3

# Script to separate out reads with spurious AG-content

# Wolf Eiserhardt, 4 April 2020

from Bio import SeqIO
import re 

sample = "0152"
n = 478182
incr = int(n/100)

gooduns = []
badunsag = []
badunsct = []

i = 0
for record in SeqIO.parse("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R1.fastq", "fastq"):
	if (i % incr) == 0:
		print(str(i/incr)+"% finished")
	i += 1
	AGcont = (len(re.findall("A",str(record.seq)))+len(re.findall("G",str(record.seq))))/len(record.seq)
	CTcont = (len(re.findall("C",str(record.seq)))+len(re.findall("T",str(record.seq))))/len(record.seq)
#	if "AAAAAAAAA" in record.seq or "GGGGGGGGGG" in record.seq or "AGAGAGAGAG" in record.seq:
	if AGcont >= 0.8:
		badunsag.append(record)
	elif CTcont >= 0.8: 
		badunsct.append(record)
	else:
		gooduns.append(record)
		
with open("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R1_gooduns.fastq", "w") as outfile:
 	SeqIO.write(gooduns, outfile, "fastq")

with open("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R1_badunsag.fastq", "w") as outfile:
 	SeqIO.write(badunsag, outfile, "fastq")

with open("/data_vol/wolf/Dypsis/assembly/"+sample+"/"+sample+"_R1_badunsct.fastq", "w") as outfile:
 	SeqIO.write(badunsct, outfile, "fastq")
