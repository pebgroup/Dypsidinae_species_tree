# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 20 March 2020

## 0. Workspace

Data folder on GIS07: `/data_vol/wolf/Dypsis/`
- `original_data`: raw read files with original naming, cf. sampling.xlsx
- `original_data_renamed`: renamed read files for compatibility with SECAPR (see 1 below)
- ``

Repository location on GIS07: `~/scripts/dypsidinae`

Analysis folder on Macbook: `~/Documents/WOLF/PROJECTS/65 Dypsis systematics paper/analysis`

## 1. Preparing data for analysis

Rename read files to four-digit names for compatibility with SECAPR.

1. Run `rename4secapr.py` to generate a bash script `rename4secapr.sh` with file copy commands. Requires `sampling.xls` (adjust path in script!). This is the reason why a bash script is generated rather than using `subprocess`, as the sampling table is on my local computer but the renaming needs to be done on the server. 

2. Run `rename4secapr.sh` from the data folder (see above). This creates a renamed copy of all files in `original_data`in `original_data_renamed`.

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