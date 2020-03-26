#!/bin/bash

f=$1

cat $f > ../trimmed_for_fastqc/$f
cat ${f/.fastq}-single.fastq >> ../trimmed_for_fastqc/$f

cat ${f/1.fastq}2.fastq > ../trimmed_for_fastqc/${f/1.fastq}2.fastq
cat ${f/1.fastq}2-single.fastq >> ../trimmed_for_fastqc/${f/1.fastq}2.fastq

