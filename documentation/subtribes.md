# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 21 October 2020

## 0. Workspace

Data folder on GIS07: `/data_vol/wolf/Dypsis/`
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
- `optrimal`: working directory for dynamic alignment trimming with optrimAl
- `iqtree`: initial gene trees (pre TreeShrink)
- `speciestree_unshrunk`: species tree built from preliminary gene trees (for reference only)
- `treeshrink`: TreeShrink analysis
- `seq_sets2_shrunk`: alignments that have been reduced by TreeShrink and degapped
- `alignments_shrunk`: realigned sequences post TreeShrink
- `iqtree_shrunk`: final gene trees (post TreeShrink)
- `speciestree`: final species tree

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

### Get assembly stats: 

From within `assembly` run:

```bash
python /usr/local/bioinf/HybPiper/get_seq_lengths.py /data_vol/wolf/PhyloPalms/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta namelist_full.txt dna > test_seq_lengths.txt

python /usr/local/bioinf/HybPiper/hybpiper_stats.py test_seq_lengths.txt namelist.txt > test_stats.txt
```

## 4. Coverage trimming and length filtering

Create directory `coverage` for coverage trimming output. 

In `assembly`, run:

```bash
while read name; do ~/scripts/dypsidinae/coverage.py $name; done < namelist.txt
```

*NB* Ensure that "exon" is chosen in the script rather than supercontig. This is currently done by (un)commenting two lines of code. 

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

Run from `seq_sets2`:

```bash
for f in *.FNA; do (linsi --thread 16 $f > ../alignments/${f/.FNA}_aligned.fasta); done
```

## 6. Gap trimming

Copy alignments to new directory `optrimal` (this is necessary as the alignments will get deleted):

```bash
mkdir optrimal
cp alignments/*.fasta optrimal
```

In that directory, generate `cutoff_trim.txt` with desired `-gt` values to be tested. 

Then, from `optrimal`: 

Prepare alignments: 
 
```bash
# replace n's with gaps in alignmenets - this will otherwise trip up TrimAl
for f in *.fasta; do (sed -i'.old' -e 's/n/-/g' $f); done
```

Remove one gene (Hey883n) that causes a weird error in optrimAL- this alignment has essentially no information anyway. 

```bash
rm HEY883n*
rm HEY125_*
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

*NB*: Hey883n lost during optrimal due to weird error - this alignment has essentially no information anyway. 

## 7. Building initial gene trees

Create directory `iqtree` and copy all trimmed alignments from `optrimal` to this directory. Then run: 

```bash
for f in *.fasta; do (~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s $f -T AUTO -ntmax 16 -B 1000 >> iqtree.log); done
```

## 8. Build species tree without further filtering

Create directory `speciestree_unshrunk`.

From `iqtree`, run: 

```bash
for f in *.treefile
do  
	pxrr -t $f -g 1013 -o temp.tre
	nw_ed temp.tre 'i & (b<30)' o >> ../speciestree_unshrunk/genetrees.tre
	rm temp.tre
done
```

Then, in `speciestree_unshrunk`, run: 

```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree.tre astral_tree_renamed.tre

java -jar ~/software/Astral/astral.5.7.3.jar -q astral_tree.tre -i genetrees.tre -o astral_tree_full_annot.tre -t 2 2> annotation.log
```

## 9. Detect branch length outliers with TreeShrink

Remove 15 genetrees that cannot be rooted and thus cannot be used downstream (in `iqtree`):

```bash
mkdir noroot
for f in *.treefile; do (~/scripts/dypsidinae/remove_noroot.py $f); done
```

Prepare TreeShrink analysis.

Create directory `treeshrink`. Then, from `iqtree`, run:

```
for f in *.treefile
do 
	mkdir ../treeshrink/${f/_aligned.fasta.treefile}
 	cp $f ../treeshrink/${f/_aligned.fasta.treefile}/input.tre
 	cp ${f/.treefile} ../treeshrink/${f/_aligned.fasta.treefile}/input.fasta
 	cd ../treeshrink/${f/_aligned.fasta.treefile}
 	sed -i'.old' -e $'s/ [0-9]\+ bp//g' input.fasta
 	cd ../../iqtree
done

cd ../treeshrink

python3 ~/software/TreeShrink/run_treeshrink.py -i . -t input.tre -a input.fasta 
```

Degap shrunk alignments and gather them in new directory `seq_sets2_shrunk` (run from `treeshrink`):

```bash
mkdir ../seq_sets2_shrunk
~/scripts/dypsidinae/outgroup_saver.py
```

## 10. Re-align 

Create directory `alignments_shrunk`, then run from `seq_sets2_shrunk`: 

```bash
for f in *; do (linsi --thread 16 $f > ../alignments_shrunk/${f/.fasta}_aligned.fasta); done
```

## 11. Building final gene trees: 

Create directory `iqtree_shrunk`.

```bash
cp alignments_shrunk/* iqtree_shrunk
```

From `iqtree_shrunk`, run: 

```bash
ls *.fasta | parallel -j 7 ~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s {} -T AUTO -ntmax 4 -B 1000
```

## 12. Building species tree: 

From `iqtree_shrunk`, run: 

```bash
for f in *.treefile
do  
	pxrr -t $f -g 1013 -o temp.tre
	nw_ed temp.tre 'i & (b<30)' o >> ../speciestree/genetrees.tre
	rm temp.tre
done
cd ../speciestree
java -jar ~/software/Astral/astral.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree.tre astral_tree_renamed.tre

java -jar ~/software/Astral/astral.5.7.3.jar -q astral_tree.tre -i genetrees.tre -o astral_tree_full_annot.tre -t 2 2> annotation.log
```

## 13. Building species tree excluding alignments with HybPiper paralog warnings: 

From `iqtree_shrunk`, run: 

```bash
for f in *.treefile
do  
	pxrr -t $f -g 1013 -o temp.tre
	nw_ed temp.tre 'i & (b<30)' o > ../speciestree_noparalog/$f
	rm temp.tre
done
cd ../speciestree_noparalog
while read n
do
	rm "${n}_shrunk_aligned.fasta.treefile"
done < ../assembly/paralogous_loci.txt
for f in *.treefile
do
	cat $f >> genetrees.tre
done

```


























