#!/bin/bash

target="Heyduk_baits/sidonie/Heyduk_palms_exons_final_concatenated_corrected.fasta"
#target="PhyloPalms/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta"

while read name
do /usr/local/bioinf/HybPiper/reads_first.py --cpu 16 -r ../trimmed/"$name"_clean-READ1.fastq ../trimmed/"$name"_clean-READ2.fastq --unpaired ../trimmed/"$name"_clean-READ12-single.fastq -b /data_vol/wolf/"$target" --prefix $name --bwa
done < namelist.txt

#/usr/local/bioinf/HybPiper/reads_first.py --cpu 16 -r ../trimmed/0173_clean-READ1.fastq ../trimmed/0173_clean-READ2.fastq --unpaired ../trimmed/0173_clean-READ12-single.fastq -b /data_vol/wolf/Heyduk_baits/sidonie/Heyduk_palms_exons_final_concatenated_corrected.fasta --prefix "0173" --bwa
