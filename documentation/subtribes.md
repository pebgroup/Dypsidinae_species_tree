# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 8 February 2021

## 0. Workspace

Data folder on Linospadix: `/data_vol/wolf/Dypsis_subtribes/`
- `original_data`: raw read files with original naming, cf. sampling.xlsx
- `original_data_renamed`: renamed read files for compatibility with SECAPR (see 1. below). Contents deleted after trimming to save space on disk.
- `fastqc_results`: results of fastqc check run via SECAPR
    - `raw`: fastqc results for raw reads (as in `original_data_renamed`)
    - `trimmed`: fastqc results after trimming (as in `trimmed`)
- `trimmed`: trimmed reads
- `assembly`: HybPiper results 
- `coverage`: output of coverage trimming step
- `seq_sets2`: sequence sets after coverage trimming and length filtering
- `alignments`: aligned sequence sets
- `alignments_exon`: aligned sequence sets with mapped exon sequences for partitioning 
- `optrimal`: working directory for dynamic alignment trimming with optrimAl
- `iqtree`: initial gene trees (pre TreeShrink)
- `speciestree_unshrunk`: species tree built from preliminary gene trees (for reference only)
- `treeshrink`: TreeShrink analysis
- `iqtree_shrunk`: final gene trees (post TreeShrink)
- `speciestree`: final species tree (post TreeShrink)

Repository location on GIS07: `~/scripts/dypsidinae`

Analysis folder on Macbook: `~/Documents/WOLF/PROJECTS/65 Dypsis systematics paper/analysis`

## 1. Preparing data for analysis

Rename read files to four-digit names for compatibility with SECAPR.
*NB this has now options for ingroup and outgroup - check before running.*

1. Run `rename4secapr.py` to generate a bash script `rename4secapr.sh` with file copy commands. Requires `sampling.xls` (adjust path in script!). This is the reason why a bash script is generated rather than using `subprocess`, as the sampling table is on my local computer but the renaming needs to be done on the server. 

2. Run `rename4secapr.sh` from the data folder (see above). This creates a renamed copy of all files in `original_data`in `original_data_renamed`.

## 2. Trimming

### Assess pre-trimming data quality

SECAPR quality check (!has to be run from within secapr_env!)

```bash
secapr quality_check --input original_data_renamed --output fastqc_results/raw
```

_This takes about 90 minutes on the server._

PDF results stored in repo in `fastqc_results/raw`.

### Trimming: 

Run in `original_data renamed`:

```bash
ls *R1* | parallel -j 4 ~/scripts/dypsidinae/trimmer.sh
```

Trimmomatic settings used: ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36

Trimmomatic v. 0.39

_This takes <90min on the server._

### Assess post-trimming data quality

Combine paired reads and singles again for comparability (created temporary directory `trimmed_for_fastqc` - this is deleted again after this step to save space). Run from within `trimmed`:

```bash
ls *READ1.fastq | parallel ~/scripts/dypsidinae/combine_posttrim_4_fastqc.sh
```

```bash
secapr quality_check --input trimmed_for_fastqc --output fastqc_results/trimmed
```

## 3. Assembly (HybPiper)

### Combine unpaired reads into a single file: 

Run in `trimmed`:

```bash
ls *1-single.fastq | parallel -j 16 ~/scripts/dypsidinae/single_combiner.sh
```

This merges `####_clean-READ1-single.fastq` and `####_clean-READ2-single.fastq` into a single file, `####_clean-READ12-single.fastq`.

### Generate name list:

Run in `trimmed`:

```bash
ls *READ2.* > namelist.txt
sed -i'.old' -e 's/_clean-READ2.fastq//g' namelist.txt
mv namelist.txt ../assembly/
rm namelist.txt.old
```

### Execute HybPiper:

Run `~/scripts/dypsidinae/piper.sh` from within `assembly`. 

_NB:_ Check that the correct target file (PhyloPalms) is selected in `piper.sh`!

### Get assembly stats: 

From within `assembly` run:

```bash
python /usr/local/bioinf/HybPiper/get_seq_lengths.py /data_vol/wolf/PhyloPalms/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist_full.txt dna > test_seq_lengths.txt

python /usr/local/bioinf/HybPiper/hybpiper_stats.py test_seq_lengths.txt namelist.txt > test_stats.txt
```

Check for paralog warnings: 

```bash
while read i
do
echo $i
python /usr/local/bioinf/HybPiper/paralog_investigator.py $i 2>> paralogreport.txt
done < namelist.txt
sed -i'.old' -e's/ paralogs written for /;/g' paralogreport.txt
```

In R:

```R
data <- read.table("paralogreport.txt", sep=";")
table(as.vector(data$V2))
write.table(unique(as.vector(data$V2)), file="paralogs", row.names=FALSE)
```

| gene | No. paralogs |
| ---- | ------------ |
| EGU105059594 | 19 | 
| EGU105042168 | 12 | 
| EGU105043827 | 10 | 
| EGU105049690 | 7 | 
| EGU105044846 | 6 | 
| HEY362       | 6 | 
| EGU105046168 | 5 | 
| EGU105059636 | 5 | 
| EGU105057015 | 3 | 
| EGU105059479 | 3 | 
| EGU105044758 | 2 | 
| EGU105033626 | 1 | 
| EGU105058687 | 1 | 
| HEY125       | 1 |
| HEY728       | 1 | 

```bash
sed -i'.old' -e's/"//g' paralogs
sed -i '1d' paralogs
```

Exclude one sample (1012) with bad recovery: 

```bash
sed -e '/1012/d' namelist.txt >> namelist_reduced.txt
```

Run `intronerate.py`:

```bash
while read name
do
	echo $name >> intronerate_out_dev.txt
	python /usr/local/bioinf/HybPiper/intronerate_dev.py --prefix $name &>> intronerate_out_dev.txt
done < namelist_reduced.txt
```

_NB_: `intronerate_dev.py` is the development version of this script, as the release version causes an error. See [here](https://github.com/mossmatters/HybPiper/issues/41). 

## 4. Coverage trimming and length filtering

Create directory `coverage` for coverage trimming output. 

In `assembly`, run:

```bash
while read name; do ~/scripts/dypsidinae/coverage.py $name; done < namelist.txt
```

This script does the following: 
- Gather all contigs from each sample in one fasta file: `coverage/sample.fasta`
- Map paired and unpaired reads to that fasta using BWA mem
- Deduplicate reads using Picard
- Calculate depth using samtools
- Mask/strip any bases with coverage *<2*
- Generate a new trimmed sample-level fasta: `coverage/sample_trimmed.fasta`

Then, in `coverage`, run: 

```bash
ls *trimmed.fasta > filelist.txt
~/scripts/dypsidinae/samples2genes.py > outstats.csv
```

This script does the following: 
- Split the sample-level fasta files up and sorts their sequences into genes. 
- Remove any sequences shorter than *150bp* or *20%* of the median sequence length of the gene
- Generate new gene fasta files in `seq_sets2`

These are ready for alignment.

## 5. Alignment

Create directory `alignments`.

Run from `seq_sets2`:

Clean up sequence names:

```bash
for f in *.FNA; do (sed -i'.old' -e $'s/-HEY[0-9]\+[p,n,s,e]* [0-9]\+-HEY[0-9]\+[p,n,s,e]*_HEY[0-9]\+[p,n,s,e]* [0-9]\+-HEY[0-9]\+[p,n,s,e]*//g' $f); done
for f in *.FNA; do (sed -i'.old' -e $'s/-EGU[0-9]\+[p,n,s,e]* [0-9]\+-EGU[0-9]\+[p,n,s,e]*_EGU[0-9]\+[p,n,s,e]* [0-9]\+-EGU[0-9]\+[p,n,s,e]*//g' $f); done
rm *.old 
~/scripts/dypsidinae/occupancy_stats.py
```

Align: 

```bash
for f in *.FNA; do (linsi --adjustdirectionaccurately --thread 16 $f > ../alignments/${f/.FNA}_aligned.fasta); done
```

## 6. Mapping exons to alignments

In `alignments`, run: 

```bash
~/scripts/dypsidinae/exon_mapper.py
```

This creates new alignments in `alignments_exon` that contain the original alignments plus the exon sequences of the two species that had the highest recovery success at each locus. 

## 7. Gap trimming

Copy alignments to new directory `optrimal` (this is necessary as the alignments will get deleted):

```bash
mkdir optrimal
cp alignments_exon/*.fasta optrimal
```

In that directory, generate `cutoff_trim.txt` with desired `-gt` values to be tested. 

Then, from `optrimal`: 

Prepare alignments: 
 
```bash
# replace n's with gaps in alignmenets - this will otherwise trip up TrimAl
for f in *.fasta; do (sed -i'.old' -e 's/n/-/g' $f); done
# change back "exo" to "exon"
for f in *.fasta; do (sed -i'.old' -e 's/exo-/exon/g' $f); done
```

Run optrimal: 

```bash
# create summary tables for all thresholds specified
~/scripts/dypsidinae/PASTA_taster.sh
# create summary table for the raw alignments
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py summary -f fasta -d dna -i *.fasta
mv summary.txt summary_0.txt
rm *.fasta
Rscript --vanilla ~/scripts/dypsidinae/optrimAl.R
```

*NB*: `optrimAL.R` was modified as to NOT discard alignments with data loss exceeding 30% (cf. [Shee et al. 2020](https://doi.org/10.3389/fpls.2020.00258)). Excessive "data loss" is probably an artefact of alignment error. 

## 8. Building initial gene trees

Create directory `iqtree` and copy all trimmed alignments from `optrimal` to this directory. 

Remove paralogous loci: 

```bash
while read l
do
	rm ${l}_aligned.fasta
done < ../assembly/paralogs
```

Then run: 

```bash
for f in *.fasta; do(sed -i'.old' -e 's/ [0-9]\+ bp//g' $f); done
rm *.old
~/scripts/dypsidinae/partitioner.py --smoother 10
for f in *_part.txt; do (cp $f ${f/_part.txt}_clean.part); done
ls *clean.fasta | parallel -j 6 ~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s {} -T AUTO -ntmax 4 -p {.}.part -B 1000
```

_NB_: one gene (EGU105046518) had no intron, resulting in an empty intron partition. The tree for this had to be run manually:

```bash
~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s EGU105046518_aligned_clean.fasta -T AUTO -ntmax 4 -B 1000
```

## 9. Build species tree without further filtering

Create directory `speciestree_unshrunk`.

Remove genetrees that cannot be rooted and thus cannot be used downstream (in `iqtree`):

```bash
mkdir noroot
for f in *.treefile; do (~/scripts/dypsidinae/remove_noroot.py $f); done
```

From `iqtree`, run: 

```bash
for f in *.treefile
do  
	~/scripts/dypsidinae/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> ../speciestree_unshrunk/genetrees.tre
	rm temp.tre
done
```

Then, in `speciestree_unshrunk`, run: 

```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree.tre astral_tree_renamed.tre
```

## 10. Detect branch length outliers with TreeShrink v. 1.3.7

Prepare TreeShrink analysis.

Create directory `treeshrink`. Then, from `iqtree`, run:

```
# Address naming issue:
mv EGU105046518_aligned_clean.fasta.treefile EGU105046518_aligned_clean.part.treefile

# Remove empty sequences that trip up treeshrink
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 1011 -d dna -f fasta -i HEY2291_aligned_clean.fasta -u fasta -g red_
mv red_HEY2291_aligned_clean.fasta-out.fas HEY2291_aligned_clean.fasta
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 1015 -d dna -f fasta -i EGU105046735_aligned_clean.fasta -u fasta -g red_
mv red_EGU105046735_aligned_clean.fasta-out.fas EGU105046735_aligned_clean.fasta

# Build treeshrink file structure
for f in *.treefile
do 
	mkdir ../treeshrink/${f/_aligned_clean.part.treefile}
 	cp $f ../treeshrink/${f/_aligned_clean.part.treefile}/input.tre
 	cp ${f/.part.treefile}.fasta ../treeshrink/${f/_aligned_clean.part.treefile}/input.fasta
 	cd ../treeshrink/${f/_aligned_clean.part.treefile}
 	sed -i'.old' -e $'s/ [0-9]\+ bp//g' input.fasta
 	cd ../../iqtree
done

cd ../treeshrink

python3 ~/software/TreeShrink/run_treeshrink.py -i . -t input.tre -a input.fasta -x 1013
```

Copy alignments to new directory, `iqtree_shrunk`, and fetch corresponding partition files (run from `treeshrink`): 

```bash
for d in *; do(cp $d/input.fasta ../iqtree_shrunk/${d}_aligned_clean.fasta); done
cp ../iqtree/*.part ../iqtree_shrunk
```

_NB_: checked that alignments are same length pre and post TreeShrink, and thus the partition files are still applicable. 

## 11. Building final gene trees: 


From `iqtree_shrunk`, run: 

```bash
ls *clean.fasta | parallel -j 6 ~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s {} -T AUTO -ntmax 4 -p {.}.part -B 1000
```

_NB_: one gene (EGU105046518) had no intron, resulting in an empty intron partition. The tree for this had to be run manually:

```bash
~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s EGU105046518_aligned_clean.fasta -T AUTO -ntmax 4 -B 1000
```

## 12. Building final species tree: 

From `iqtree_shrunk`, run: 

```bash
for f in *.treefile
do  
	~/scripts/dypsidinae/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> ../speciestree/genetrees.tre
	rm temp.tre
done
```

Then, in `speciestree`, run: 

```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree.tre astral_tree_renamed.tre
java -jar ~/software/Astral/astral.5.7.3.jar -q astral_tree.tre -i genetrees.tre -o astral_tree_full_annot.tre -t 2 2> annotation.log
```






















