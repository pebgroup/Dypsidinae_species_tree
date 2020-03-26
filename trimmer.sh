#!/bin/bash

f=$1

java -jar /usr/local/bioinf/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 $f ${f/1.fastq}2.fastq ../trimmed/${f/_R1.fastq}_clean-READ1.fastq ../trimmed/${f/_R1.fastq}_clean-READ1-single.fastq ../trimmed/${f/_R1.fastq}_clean-READ2.fastq ../trimmed/${f/_R1.fastq}_clean-READ2-single.fastq ILLUMINACLIP:/usr/local/bioinf/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36