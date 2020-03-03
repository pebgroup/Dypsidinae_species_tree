import numpy as np
import pandas as pd

sampling = pd.read_excel('../sampling.xlsx', sheet_name='sampled_spp', converters={'SECAPR No.':str})

sampling = sampling.loc[sampling['SECAPR No.'].notna()]

f = open("rename4secapr.sh", "w")
print("#!/bin/bash", file=f)

for index, row in sampling.iterrows(): 
	#print(row['Sequence file name'], row['SECAPR No.'])
	print("cp original_data/"+row['Sequence file name']+"_R1_001.fastq for_trimming/"+row['SECAPR No.']+"_R1.fastq", file=f)
	print("cp original_data/"+row['Sequence file name']+"_R2_001.fastq for_trimming/"+row['SECAPR No.']+"_R2.fastq", file=f)
	
f.close()
