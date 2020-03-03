# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 21 February 2020

## 0. Workspace

Data folder on GIS07: `/data_vol/wolf/Dypsis`

Scripts on GIS07: `~/scripts`

Analysis folder on Macbook: `~/Documents/WOLF/PROJECTS/65 Dypsis systematics paper/analysis`

## 1. Preparing data for analysis

Rename read files to four-digit names for compatibility with SECAPR using `rename4secapr.py`

```python
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
``` 
This generates a bash script `rename4secapr.sh` to be run in data directory on GIS07.

## 2. Trimming

`mkdir trimmed` in data directory. 

### Assess pre-trimming data quality

`mkdir fastqc_results`

`mkdir fastqc_results/raw`

SECAPR quality check (!has to be run from within secapr_env!)

```bash
secapr quality_check --input for_trimming --output fastqc_results/raw
```

Trimming: `trimming.sh`

```bash
#!/bin/bash

trim () {
        local f=$1
        java -jar /usr/local/bioinf/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $f ${f/1.fastq}2.fastq ${f/_R1.fastq}_clean-READ1.fastq ${f/_R1.fastq}_clean-READ1-single.fastq ${f/1_R1.fastq}_clean-READ2.fastq ${f/_R1.fastq}_clean-READ2-single.fastq ILLUMINACLIP:/usr/local/bioinf/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
}


for f in *1.fastq; do trim "$f"; done
```