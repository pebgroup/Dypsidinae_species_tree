# Subtribes: TreeShrink step 

10 February 2020

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

python3 ~/software/TreeShrink/run_treeshrink.py -i . -t input.tre -a input.fasta -x 1013 > treeshrink.log
```

Copy alignments to new directory, `iqtree_shrunk`, and fetch corresponding partition files (run from `treeshrink`): 

```bash
for d in *; do(cp $d/output.fasta ../iqtree_shrunk/${d}_aligned_clean.fasta); done
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
## Notes

Sed command to reformat TreeShrink screen output: 

```bash
sed '/:/N;y/\n/\t/' treeshrink.log.red | sed -e 's/:[ \t][ \t]*will be cut in /,/g' | sed 's/ trees where its impact is above /,/g' | sed 's/ for quantile 0.05//g' > treeshrink.log.red2
```

# Species-level analysis (old)

17 February 2021

*NEEDS UPDATING FROM HERE*

Copy alignments to `final_tree_nofilter` and `final_tree_filtered` and `final_tree_filtered_length_only`.

For comparison, copy the alignments that nave _not been edited manually_ to `final_tree_unedited`, and remove the alignments that were excluded during the manual editing step.  

## 13. Remove multiple sequences of same individual

In `final_tree_nofilter` and `final_tree_filtered` and `final_tree_filtered_length_only`, run: 

```bash
for f in *.fasta; do(sed -i'.old' -e 's/ [0-9]\+ bp//g' $f); done
rm *.old
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 0075 0076 0157 0197 -d dna -f fasta -i *.fasta -u fasta -g red_
rm reduced*
for f in *.fas; do (mv $f ${f#red_}ta); done
```

## 14. Filtering

This step does the following: 
- From each alignment, exclude all sequences that cover <50% of "well occupied" alignment columns, defined as columns that have data for >= 70% of species (pers. comm. P. Bailey)  (cf. `lenght_filter.py`).
- Across all alignments, remove any species that is represented in <20 alignments (cf. `occupancy.py`).

From `final_tree_filtered` and `final_tree_filtered_length_only`, run: 

Removing short sequences:

```bash
~/scripts/dypsidinae/partitioner.py --smoother 10
rm *_part.txt
for f in *_clean.fasta; do (~/software/trimal -in $f -out ${f/.fasta}_70.fasta -gt 0.7); done #NB may have to substitute / with % for suffix removal
mkdir bad
mv reduced_417* bad #this alignment has no well-occupied columns
for f in *_clean_70.fasta; do (~/scripts/dypsidinae/length_filter.py $f >> lenght_filter.log); done
```

In `final_tree_filtered` ONLY: dropping all taxa that occur in fewer than 20 of the length filtered alignments: 

```bash
~/scripts/dypsidinae/occupancy.py
```

## 15. Tree building 2nd round

In all tree building subdirectories (`final_tree_nofilter`, `final_tree_filtered`,  `final_tree_filtered_length_only`, and `final_tree_unedited`), move or copy alignments to a subdirectory `iqtree`. Then run:


```bash
~/scripts/dypsidinae/partitioner.py --smoother 10
for f in *clean.fasta; do (~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s $f -T AUTO -ntmax 8 -p ${f/clean.fasta}part.txt -B 1000 >> exl_trees.log); done
```


Build species tree: 

```bash
for f in *.treefile
do 
	~/scripts/dypsidinae/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> ../../speciestree_filtered/genetrees.tre 
	rm temp.tre
done
cd ../../speciestree_filtered
java -jar ~/software/Astral/astral.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
~/scripts/dypsidinae/renamer.py ../rename.csv astral_tree.tre astral_tree_renamed.tre
```

## 16. Dating

Gather genetrees. In `sortadate/genetrees`, run: 

```bash
cp ../../length_filter/iqtree/*treefile .
for f in *.treefile
do 
	~/scripts/dypsidinae/rooter.py $f
	mv temp.tre $f
done
```

In `sortadate`, run: 

```bash
python2 ~/software/SortaDate/src/get_var_length.py genetrees/ --flend .treefile --outf var --outg 1011,1012
python2 ~/software/SortaDate/src/get_bp_genetrees.py genetrees/ ../speciestree_filtered/astral_tree.tre --flend .treefile --outf bp
python2 ~/software/SortaDate/src/combine_results.py var bp --outf comb
python2 ~/software/SortaDate/src/get_good_genes.py comb --max 30 --order 3,1,2 --outf gg
```

Get the minimum number of genes required to cover all 155 species, and copy the alignments and partition files from `lengh_filter` (run in `sortadate`):

```bash
mkdir alignments
~/scripts/dypsidinae/sortadater.py
```

Prepare concatenated alignment: 

```bash
for f in *part.txt 
do
	# Remove prefixes from partition files
	sed -i'.old' -e 's/DNA, //g' $f
	# Remove junk from sequence names
	sed -i'.old' -e 's/ [0-9]\+ bp//g' ${f/_part.txt}.fasta
	# Split alignments into intron and exon parts
	python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py split -f fasta -d dna -i ${f/_part.txt}.fasta -l $f -u fasta
done
# Concatenate all exon alignments...
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py concat -f fasta -d aa -i *exon-out.fas -p partitions_exon.txt -t concatenated_exon.fas
# ... and all intron alignments
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py concat -f fasta -d aa -i *intron-out.fas -p partitions_intron.txt -t concatenated_intron.fas
# concatenate exons and introns
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py concat -f fasta -d aa -i concatenated_intron.fas concatenated_exon.fas
# remove exons from alingment
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x exon1 exon2 -d dna -f fasta -i concatenated.out -g cropped_
sed -i'.old' -e 's///g' cropped_concatenated.out-out.fas
mkdir tree
cp cropped_concatenated.out-out.fas ../tree
cp partitions.txt ../tree 
cd ../tree
```

Build phylogram for dating: 

First, manually edit `partitions.txt` to include `DNA, ` prefix and space around `=`. Then run:  

```bash
~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s cropped_concatenated.out-out.fas -T AUTO -ntmax 16 -p partitions.txt -g ../../speciestree_filtered/astral_tree.tre >> phylogram.log
pxrr -t partitions.txt.treefile -g 1011,1012 -o partitions.txt.treefile.rooted
pxrmt -t partitions.txt.treefile.rooted -n 1011,1012 -o partitions.txt.treefile.pruned
```

Run treepl (after manually creating configuration file `config`):

```bash
/home/au265104/.linuxbrew/Cellar/treepl/2018.05.22/bin/treePL config
~/scripts/dypsidinae/renamer.py ../../rename.csv --bs 1 partitions.txt.treefile.dated partitions.txt.treefile.dated.renamed.tre
```

## 16. Diagnose genetrees 

In `length_filter`, run: 

```bash
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py summary -f fasta -d dna -i *lf_exl.fasta -c 12
```

Older:

## 7. iqtree gene tree inference (with UFbootstrap)

_NB: iqtree is installed in "Dypsis" conda environment._

```bash
for f in *_noempty.fasta;do(iqtree -s $f -m GTR+G10 -B 1000 -T 16); done
```

## 8. Diagnose trees with TreeShrink

*NB: This description is for the opTrimal alignments. For the Trimal -gt 0.5 alignments, use `alignments_trimmed` instead of `optrimal`, and `treeshrink_gt0.5` instead of `treeshrink`. Also, `treeshrink_prep.sh` needs to be adjusted accordingly.*

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
~/scripts/dypsidinae/renamer.py rename.csv astral_tree_prelim.tre astral_tree_prelim_renamed.tre
```

Root to _Loxococcus_:
```bash
pxrr -g Loxococcus-rupicola-SBL8-S7,Loxococcus-rupicola-SBL234-S35 -t astral_tree_prelim_renamed.tre -o astral_tree_prelim_renamed_rerooted.tre
```

Fully annotated ASTRAL tree: 
```bash
java -jar ~/software/Astral/astral.5.7.3.jar -t 2 -i treeshrink/iqtrees.tre -o astral_tree_prelim_t2.tre 2> astral_prelim_t2.log
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

# Even older pipeline...

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
- `assembly`: HybPiper results (see 3. below)
- `assembly_excluded`: samples that have been processed up to HybPiper (step 3) but were subsequently excluded (old outgroup, bad samples)
- `seq_sets`: unaligned sequence sets per locus as retrieved by HybPiper `retrieve_sequences.py_
- `alignments`: aligned sequence sets
- `fasttrees`: preliminary gene trees
- `fasttrees_collapsed`: preliminary gene trees with nodes BS<0.1 collapsed
- `coverage`: output of coverage trimming step (see 7. below)
- `seq_sets2`: sequence sets after coverage trimming and length filtering (see 7. below)
- `alignments2`: aligned sequence sets after coverage trimming and length filtering
- `alignments_trimmed`: alignments after removing gappy sites
- `raxmltrees`: output of raxml analysis
- `seq_sets2_reduced`: `seq_sets2` after removing superfluous species
- `alignments2_reduced`: `alignments2` after removing superfluous species
- `alignments_trimmed_reduced`: `alignments_trimmed` after removing superfluous species
- `raxmltrees_reduced`: same as `raxmltrees` after removing superfluous species
- `raxmltrees_collapsed_reduced`: same as `raxmltrees_collapsed` after removing superfluous species
- `raxmltrees_treeshrink`: treeshrink output

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

### Retrieve sequences: 

Run from `assembly`:

`python /usr/local/bioinf/HybPiper/retrieve_sequences.py /data_vol/wolf/Heyduk_baits/sidonie/Heyduk_palms_exons_final_concatenated_corrected.fasta . supercontig &> outstats.txt`

`outstats.txt` can be used for crude exclusion of loci based on number of seqs retrieved. 

Copied sequence sets to `seq_sets`

## 4. Alignment (MAFFT)

Run from `seq_sets`:

```bash
for f in *; do (linsi --thread 16 $f > ../alignments/${f/.FNA}_aligned.fasta); done
```

## 5. Build preliminary gene trees

Make list of loci: 

`ls alignments | sed 's/_aligned.fasta//g' > locuslist.txt`

Run FastTree on alignments, save trees in `fasttrees`:

`while read locus; do (FastTree -gtr -nt alignments/"$locus"_aligned.fasta > fasttrees/"$locus"_fasttree.tre); done < locuslist.txt`

## 6. Build (preliminary) species tree

Collapse splits with BS <0.1 in preliminary gene trees (results in `fasttrees_collapsed`):

`while read locus; do (nw_ed fasttrees/"$locus"_fasttree.tre 'i & b<=0.1' o > fasttrees_collapsed/"$locus"_fasttree.tre); done < locuslist.txt`

Gather all fasttrees in one tree file, `fasttrees.tre`: 

`while read locus; do (cat fasttrees_collapsed/"$locus"_fasttree.tre >> fasttrees.tre); done < locuslist.txt`

Run ASTRAL on that file: 

`java -jar ~/software/Astral/astral.5.7.3.jar -i fasttrees.tre -o astral_tree.tre 2> astral.log`

## 7. Coverage trimming and length filtering

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

These are ready for further processing (e.g. TrimAl) and phylogenetic analysis. 

## 8. Remove superfluous species 

Remove Loxococcus-rupicola-SBL8-S7 from alignments (SECAPR no 1012). Run in `seq_sets2`:

```bash
for f in *.FNA; do (sed -i'.old' -E 's/ [0-9]{4}_.+//g' $f); done
rm *.old
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 1012 -d dna -f fasta -i *.FNA > amaslog.txt
mkdir ../seq_sets2_reduced
for f in *-out.fas; do (mv $f ${f/-out.fas}); done
for f in reduced_*; do (mv $f ../seq_sets2_reduced/${f/#reduced_}); done
```

Results are in `seq_sets2_reduced`. Run all the following based on that data. 

## 9. Alignment (MAFFT)

see (4) but results in `alignments2` or `alignments2_reduced`

## 10. Gap trimming

In `alignments2` run:

```bash
for f in *; do (~/software/trimal -in $f -out ../alignments_trimmed/${f/_aligned.fasta}_trimmed.fasta -gt 0.5); done
```

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

# Subtribe-level analysis (old)



### Retrieve sequences: 

Run from `assembly`:

`python /usr/local/bioinf/HybPiper/retrieve_sequences.py /data_vol/wolf/Heyduk_baits/sidonie/Heyduk_palms_exons_final_concatenated_corrected.fasta . supercontig &> outstats.txt`

`outstats.txt` can be used for crude exclusion of loci based on number of seqs retrieved. 

Copied sequence sets to `seq_sets`

## 4. Alignment (MAFFT)

Run from `seq_sets`:

```bash
for f in *; do (linsi --thread 16 $f > ../alignments/${f/.FNA}_aligned.fasta); done
```

## 5. Build preliminary gene trees

Make list of loci: 

`ls alignments | sed 's/_aligned.fasta//g' > locuslist.txt`

Run FastTree on alignments, save trees in `fasttrees`:

`while read locus; do (FastTree -gtr -nt alignments/"$locus"_aligned.fasta > fasttrees/"$locus"_fasttree.tre); done < locuslist.txt`

## 6. Build (preliminary) species tree

Collapse splits with BS <0.1 in preliminary gene trees (results in `fasttrees_collapsed`):

`while read locus; do (nw_ed fasttrees/"$locus"_fasttree.tre 'i & b<=0.1' o > fasttrees_collapsed/"$locus"_fasttree.tre); done < locuslist.txt`

Gather all fasttrees in one tree file, `fasttrees.tre`: 

`while read locus; do (cat fasttrees_collapsed/"$locus"_fasttree.tre >> fasttrees.tre); done < locuslist.txt`

Run ASTRAL on that file: 

`java -jar ~/software/Astral/astral.5.7.3.jar -i fasttrees.tre -o astral_tree.tre 2> astral.log`

## 7. Coverage trimming and length filtering

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

These are ready for further processing (e.g. TrimAl) and phylogenetic analysis. 

## 8. Remove superfluous species 

Remove Loxococcus-rupicola-SBL8-S7 from alignments (SECAPR no 1012). Run in `seq_sets2`:

```bash
for f in *.FNA; do (sed -i'.old' -E 's/ [0-9]{4}_.+//g' $f); done
rm *.old
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 1012 -d dna -f fasta -i *.FNA > amaslog.txt
mkdir ../seq_sets2_reduced
for f in *-out.fas; do (mv $f ${f/-out.fas}); done
for f in reduced_*; do (mv $f ../seq_sets2_reduced/${f/#reduced_}); done
```

Results are in `seq_sets2_reduced`. Run all the following based on that data. 

## 9. Alignment (MAFFT)

see (4) but results in `alignments2` or `alignments2_reduced`

## 10. Gap trimming

In `alignments2` run:

```bash
for f in *; do (~/software/trimal -in $f -out ../alignments_trimmed/${f/_aligned.fasta}_trimmed.fasta -gt 0.5); done
```

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


