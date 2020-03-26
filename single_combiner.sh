#!/bin/bash

f=$1

cat $f > ${f/1-single.fastq}12-single.fastq
cat ${f/1-single.fastq}2-single.fastq >> ${f/1-single.fastq}12-single.fastq

