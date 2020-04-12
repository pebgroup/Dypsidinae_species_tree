# Exploring the Heyduk targets

## Questions: 

* What are the Heyduk et al. 2015 targets? 
* How much of them is coding? 
* Do the previously identified repeat regions fall into non-coding parts (UTRs)?

## 1. Extract Elaeis targets only

`python2.7 GetFasta2.py Heyduk_palms.exons.final_corrected.fasta Elaeis_exons.txt Elaeis_exons.fasta`

## 2. Download Elaeis genome and annotation

Navigated to EG5 (GCA_000442705.1) on GenBank: https://www.ncbi.nlm.nih.gov/assembly/GCF_000442705.1  and downloaded `.fna` and `.gff` files. 

## 3. Map Elaeis targets to the genome. 

I first tried doing this using BLAST: 

```bash
makeblastdb -in GCF_000442705.1_EG5_genomic.fna -dbtype nucl -title Elaeis -parse_seqids`
blastn -query /Users/au265104/Desktop/To_send/Elaeis_exons.fasta -db GCF_000442705.1_EG5_genomic.fna -out test.xml -outfmt 5
```

The blasting itself worked, but I did not manage to convert the results to a SAM file. `blast2bam`did not work (somehow the identity of the targets was lost) and running blast with `-outfmt 17` did not produce a viable SAM file either. 

I then moved on to BWA: 

```bash
bwa index GCF_000442705.1_EG5_genomic.fna

bwa mem GCF_000442705.1_EG5_genomic.fna /Users/au265104/Desktop/To_send/Elaeis_exons.fasta > bwamap.sam

samtools view -S -b bwamap.sam > bwamap.bam
samtools sort bwamap.bam -o bwamap.sorted.bam
samtools index bwamap.sorted.bam

bedtools bamtobed -i bwamap.sorted.bam > bwamap.bed
```

## 4. Intersect Heyduk targets with genome annotations and filter out CDS

```bash
bedtools intersect -loj -a bwamap.bed -b ../../genome_assemblies_genome_gff/ncbi-genomes-2020-04-11/GCF_000442705.1_EG5_genomic.gff > bwamap.intersect.txt

grep -e "\tCDS\t" bwamap.intersect.txt > bwamap.intersect.CDS.txt
```

## 5. Analyse results 

`./matchparser.py`

This script generates 

1. `exons_vs_NCBIgenes.csv`: A list of Heyduk *exon* numbers with their corresponding NCBI *gene* number(s) derived from the overlapping CDS feature in the Elaeis `.gff`. Note that some exons overlap with multiple CDS features! 

2. `genes_vs_NCBIgenes2.csv`: A list of Heyduk *gene* numbers with their corresponding NCBI *gene* number(s), and the product of this gene as per the annotation of the CDS feature. If Heyduk genes match multiple NCBI genes (either because individual exons have multiple matches, or because different exons map to different genes - both happens!), these multiple matches are represented as multiple rows. Similarly, if the same NCBI gene has multiple product annotations from different CDS features, these are represented as multiple rows. For convenience, there is a column showing for each Heyduk gene how many NCBI genes it matches, and for each NCBI gene how many product annotations have been found. The table comes in following format: 

Heyduk code | No. of NCBI genes per Heyduk code | NCBI gene | No. of products per NCBI gene | product

_Note that `matchparser.py` ignores different isoforms in product annotations!_

## 6. Conclusions

Manual inspection of `pwamap.intersect.CDS.txt` and comparison with read mapping shows that some repeats are indeed in UTRs, but not all of them. Relatively large parts of some of the Heyduk targets (up to >200bp) are non-coding, and removing all non-coding regions would thus lead to a significant loss of target sequence. 