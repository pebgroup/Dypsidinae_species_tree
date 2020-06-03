# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 3 June 2020

## 0. Workspace

Data folder on GIS07: `/data_vol/wolf/Dypsis/`
- `original_data`: raw read files with original naming, cf. sampling.xlsx
- `original_data_renamed`: renamed read files for compatibility with SECAPR (see 1. below). Contents deleted after trimming to save space on disk.
- `fastqc_results`: results of fastqc check run via SECAPR
    - `raw`: fastqc results for raw reads (as in `original_data_renamed`)
    - `trimmed`: fastqc results after trimming (as in `trimmed`)
    - `trimmed2`: fastqc results after trimming (as in `trimmed2`)
- `trimmed`: trimmed reads (see 2. below)
- `trimmed2`: trimmed reads with alternative trimming criteria (see 2. below)
- `trimmed_for_fastqc`: temporary directory with combined trimmed readfiles for FASTQC
- `assembly`: HybPiper results (see 3. below)
- `assembly_excluded`: samples that have been processed up to HybPiper (step 3) but were subsequently excluded (old outgroup, bad samples)
- `coverage`: output of coverage trimming step (see 7. below)
- `seq_sets2`: sequence sets after coverage trimming and length filtering (see 7. below)
- `alignments2`: aligned sequence sets after coverage trimming and length filtering
- `alignments_trimmed`: alignments after "static" TrimAl trimming (-gt 0.5)
- `optrimal`: working directory for dynamic alignment trimming with optrimAl
- `treeshrink`: working directory for TreeShrink


Repository location on GIS07: `~/scripts/dypsidinae`

Analysis folder on Macbook: `~/Documents/WOLF/PROJECTS/65 Dypsis systematics paper/analysis`

## 1. Preparing data for analysis

Rename read files to four-digit names for compatibility with SECAPR.
*NB this has now options for ingroup and outgroup - check before running.*

1. Run `rename4secapr.py` to generate a bash script `rename4secapr.sh` with file copy commands. Requires `sampling.xls` (adjust path in script!). This is the reason why a bash script is generated rather than using `subprocess`, as the sampling table is on my local computer but the renaming needs to be done on the server. 

2. Run `rename4secapr.sh` from the data folder (see above). This creates a renamed copy of all files in `original_data`in `original_data_renamed`.

3. Manually added a sample that has been resequenced as `Dypsis-heterophylla-SBL179-repooled_*.fastq`. Manually added to `original_data_renamed` as `0201_R*.fastq`

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

### Alternative trimming (more stringent settings, 1.4.2020): 

Run in `original_data renamed`:

```bash
ls *R1* | parallel -j 4 ~/scripts/dypsidinae/trimmer2.sh
```

Trimmomatic settings used: ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36


### Assess post-trimming data quality

Combine paired reads and singles again for comparability (created temporary directory `trimmed_for_fastqc` - this is deleted again after this step to save space). Run from within `trimmed`:

```bash
ls *READ1.fastq | parallel ~/scripts/dypsidinae/combine_posttrim_4_fastqc.sh
```

```bash
secapr quality_check --input trimmed_for_fastqc --output fastqc_results/trimmed
```

Or for alternative trimming settings (see above):

```bash
secapr quality_check --input trimmed_for_fastqc --output fastqc_results/trimmed2
```

PDF results stored in repo in `fastqc_results/trimmed`/`fastqc_results/trimmed2`.

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

*NB* This uses the first trimming settings (i.e., data from `trimmed`). 

### Get assembly stats: 

From within `assembly` run:

```bash
python /usr/local/bioinf/HybPiper/get_seq_lengths.py /data_vol/wolf/Heyduk_baits/sidonie/Heyduk_palms_exons_final_concatenated_corrected.fasta namelist.txt dna > test_seq_lengths.txt

python /usr/local/bioinf/HybPiper/hybpiper_stats.py test_seq_lengths.txt namelist.txt > test_stats.txt
```

### Do intronerate:

(generate new namelist if necessary, e.g. after excluding samples)
```bash
ls -d */ > namelist.txt
sed -i'.old' -e 's/\///g' namelist.txt
```

Run `intronerate.py`:

```bash
while read name; do (python /usr/local/bioinf/HybPiper/intronerate.py --prefix $name &>> intronerate_out.txt); done < namelist.txt
```

## 4. Coverage trimming and length filtering

Create directory `coverage` for coverage trimming output. 

In `assembly`, run:

```bash
while read name; do ~/scripts/dypsidinae/coverage.py $name; done < namelist.txt
```

*NB* Ensure that "supercontig" is chosen in the script rather than exon. This is currently done by (un)commenting two lines of code. 

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

These are ready for further processing (e.g. TrimAl) and phylogenetic analysis. 

## 5. Alignment (MAFFT)

Run from `seq_sets2`:

```bash
for f in *; do (linsi --thread 16 $f > ../alignments2/${f/.FNA}_aligned.fasta); done
```

One alignment (757_aligned.fasta) failed with linsi, and was thus redone with:

`mafft --thread 16 757.FNA > ../alignments2/757_aligned.fasta` 

## 6. Gap trimming

Make `alignments trimmed`. 

In `alignments2` run:

```bash
for f in *; do (~/software/trimal -in $f -out ../alignments_trimmed/${f/_aligned.fasta}_trimmed.fasta -gt 0.5); done
```

Some of the trimmed alignments contain empty sequences. To remove these, run (in `alignments_trimmed`): 

```bash
for f in *.fasta;do(~/scripts/dypsidinae/noempty.py $f);done
1171_trimmed_noempty.fasta has 3 empty sequences removed
120_trimmed_noempty.fasta has 4 empty sequences removed
360_trimmed_noempty.fasta has 6 empty sequences removed
732_trimmed_noempty.fasta has 3 empty sequences removed
874_trimmed_noempty.fasta has 7 empty sequences removed
938_trimmed_noempty.fasta has 5 empty sequences removed
```

Alternatively: use [optrimAl](https://github.com/keblat/bioinfo-utils/blob/master/docs/advice/scripts/optrimAl.txt): 

Copy alignments from `alignments2` into `optrimal`. In that directory, generate `cutoff_trim.txt` with desired `-gt` values to be tested. 

Then, from `optrimal`:  
 
```bash
# create summary tables for all thresholds specified
~/scripts/dypsidinae/PASTA_taster.sh
# create summary table for the raw alignments
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py summary -f fasta -d dna -i *.fasta
mv summary.txt summary_0.txt
rm *.fasta
Rscript --vanilla ~/scripts/dypsidinae/optrimAl.R
```

Some of the alignments generated by optrimal contain empty sequences. To remove these, run (in `optrimal`): 

```bash
for f in *.fasta;do(~/scripts/dypsidinae/noempty.py $f);done
120_aligned_noempty.fasta has 4 empty sequences removed
874_aligned_noempty.fasta has 3 empty sequences removed
938_aligned_noempty.fasta has 4 empty sequences removed
```

## 7. iqtree gene tree inference (with UFbootstrap)

_NB: iqtree is installed in "Dypsis" conda environment._

```bash
for f in *_noempty.fasta;do(iqtree -s $f -m GTR+G10 -B 1000 -T 16); done
```

## 8. Diagnose trees with TreeShrink

Rename some sequences that have been reverse-complemented by MAFFT (in gene 31): 

Run from `optrimal`:

```bash 
sed -i'.old'  -e 's/_R_//g' 31_aligned_noempty.fasta
sed -i'.old'  -e 's/_R_//g' 31_aligned_noempty.fasta.treefile 
```

In `optrimal`, run: `~/scripts/dypsidinae/treeshrink_prep.sh`. This creates a file structure appropriate for TreeShrink in `treeshrink`. 

In `treeshrink`, run: 

```bash
python3 ~/software/TreeShrink/run_treeshrink.py -i . -t input.tre -a input.fasta

Species 0178 only exists in 5 gene trees
Species 0166 only exists in 17 gene trees
Species 0151 only exists in 18 gene trees
Species 0169 only exists in 6 gene trees
Species 0203 only exists in 5 gene trees
Species 0157 only exists in 6 gene trees
Species 0154 only exists in 4 gene trees
Species 0147 only exists in 3 gene trees
Species 0164 only exists in 1 gene trees
Species 0159 only exists in 2 gene trees
Species 0168 only exists in 2 gene trees
```

## 9. Build preliminary species tree

Collapse poorly supported nodes (from `treeshrink`):

```bash
ls -d */ | sed -e 's/\///g' | parallel 'nw_ed {}/input_shrunk_0.05.tre "i & b<=0.3" o > {}/input_shrunk_0.05_collapsed.tre'
```

Collate all trees into one file (this does NOT work properly with parallel):

```bash
ls -d */ | sed -e 's/\///g' > locuslist.txt
while read locus; do(cat $locus/input_shrunk_0.05_collapsed.tre >> iqtrees.tre);done < locuslist.txt
```

Run ASTRAL (from `Dypsis`):

```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i treeshrink/iqtrees.tre -o astral_tree_prelim.tre 2> astral_prelim2.log
```

Rename tip labels:

```bash
~/scripts/dypsidinae/renamer.py rename.csv astral_tree_raxml.tre astral_tree_raxml_renamed.tre
```

Fix some renaming errors:
```bash
sed -i'.old' -e's/3584Dypsis_integra_45014_S41_L001/35840038/g' astral_tree_prelim_renamed.tre
sed -i'.old' -e's/9738Dypsis_vonitrandambo_161115-2_S16_L001/97380074/g' astral_tree_prelim_renamed.tre
sed -i'.old' -e's/7694Dypsis-rakotonasoloi-SBL287/76940178/g' astral_tree_prelim_renamed.tre
sed -i'.old' -e's/3716Dypsis-malcomberi_S58_L001/37160102/g' astral_tree_prelim_renamed.tre

```



[UNTIL HERE]


Diagnose alignments: 

```bash
ls *.fasta > filelist.txt
~/scripts/dypsidinae/pic.py
```

This creates an output file `alnstats.csv` with all sorts of statistics (including PICs) for all alignments. 

## 9. RAxML gene tree inference 

Prepare alignments for RAxML (remove annotation from fasta): 
```bash
for f in *.fasta; do (sed -i'.old' -e $'s/ [0-9]\+ bp//g' $f); done
rm *.old
```

Run RAxML-ng: 

```bash
ls *.fasta | parallel -j 16 raxml-ng --all --msa {} --model GTR+G --tree pars{10} --bs-trees 200 --threads 1
```

Move all RAxML output to `raxmltrees`

## 10. Build species tree

Generate updated locuslist (run from `raxmltrees`):

```bash
ls *.support | sed -e 's/_trimmed.fasta.raxml.support//g' > ../locuslist.updated.txt
```

Collapse very short branches (effectively polytomies, <0.00001): 

```bash
for f in *.support; do (Rscript --vanilla ~/scripts/dypsidinae/di2multi.R $f); done
```

Collapse very poorly supported branches: 

```bash
while read locus; do (nw_ed raxmltrees/"$locus"_trimmed.fasta.raxml.support.noshort 'i & b<=0.1' o > raxmltrees_collapsed/"$locus"_raxml.tre); done < locuslist.updated.txt
```

Gather all trees into one file `raxmltres.tre`

```bash
while read locus; do (cat raxmltrees_collapsed/"$locus"_raxml.tre >> raxmltrees.tre); done < locuslist.updated.txt
```

Remove outlier branches using Treeshrink: 

```bash
python3 ~/software/TreeShrink/run_treeshrink.py -t raxmltrees.tre
```

Output saved in `raxmltrees_treeshrink`

Manually created mapping file to deal with two _Loxococcus_ individuals: `mapfile.txt`
Then run ASTRAL:

(without treeshrink)
```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i raxmltrees.tre -o astral_tree_raxml.tre -a mapfile.txt 2> astral.log
```

(with treeshrink)
```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i raxmltrees_treeshrink/raxmltrees_0.05.tre -o astral_tree_raxml_shrunk.tre -a mapfile.txt 2> astral.log
```

(with reduced sampling: only one _Loxococcus_)
```bash
java -jar ~/software/Astral/astral.5.7.3.jar -i raxmltrees_reduced_treeshrink/raxmltrees_reduced_0.05.tre -o astral_tree_raxml_reduced_shrunk.tre
```

Rename taxa in resulting tree: 

```bash
~/scripts/dypsidinae/renamer.py rename.csv astral_tree_raxml.tre astral_tree_raxml_renamed.tre
```

## 11. Phyparts

```bash
java -jar ~/software/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d test/testd -m test/sp.tre -o out
```

## WORK LOG

[14.4.2020]

HybPiper four additional Dypsidinae samples (0202-0205) that I received from Sidonie using the Heyduk targets. Created new folder `/data_vol/wolf/Dypsis_added`. Within this folder, created same data structure as in main `Dypsis` folder. Moved trimmed reads into `trimmed`, generated name list, ran `piper.sh` from `assembly`.

[15.4.2020]

Reintegrated the four additional Dypsidinae species into main `Dypsis` folder:
```bash
cd /data_vol/wolf/Dypsis_added
mv trimmed/* ../Dypsis/trimmed
cat assembly/namelist.txt >> ../Dypsis/assembly/namelist.txt
rm assembly/namelist.txt
mv assembly/* ../Dypsis/assembly/
```

Preliminary analysis of subtribes (Phylopalms loci): `Dypsis_subtribes`
Extracted sequences. Excluded loci with <10 retrieved sequences:

```bash
grep -e 'Found [0-9] sequences' outstats.txt
Found 2 sequences for HEY168.
Found 8 sequences for HEY763.
Found 1 sequences for HEY111.
Found 8 sequences for HEY363.
mkdir rejected_loci
mv HEY168.FNA rejected_loci/
mv HEY763.FNA rejected_loci/
mv HEY111.FNA rejected_loci/
mv HEY363.FNA rejected_loci/
```

Ran all steps from retrieve sequences till preliminary ASTRAL tree above. 

[16-17.4.2020]

Developed `coverage.py` and `samples2genes.py`. Ran on `Dypsis_subtribes` data. 

[27.4.2020]

Copied Loxococci from `Dypsis_subtribes` to `Dypsis`  - this will be the outgroup

`cp trimmed/1011_* ../Dypsis/trimmed`
`cp trimmed/1012_* ../Dypsis/trimmed`

Moved old outgroups and "bad samples" from `Dypsis/assembly` to `Dypsis/assembly_excluded`:

```bash
mv 0016 ../assembly_excluded/
mv 0056 ../assembly_excluded/
mv 0188 ../assembly_excluded/
mv 0189 ../assembly_excluded/
mv 0190 ../assembly_excluded/
```



