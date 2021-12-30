# Analysis of Dypsidinae target capture data: main analysis

## 0. Workspace

Data folder on GIS07: `/data_vol/wolf/Dypsis/`
- `original_data`: raw read files with original naming, cf. sampling.xlsx
- `original_data_renamed`: renamed read files for compatibility with SECAPR (see 1. below). 
- `fastqc_results`: results of fastqc check run via SECAPR
    - `raw`: fastqc results for raw reads (as in `original_data_renamed`)
    - `trimmed`: fastqc results after trimming (as in `trimmed`)
    - `trimmed2`: fastqc results after trimming (as in `trimmed2`)
- `trimmed`: trimmed reads (see [below](#2-trimming))
- `trimmed2`: trimmed reads with alternative trimming criteria (see [below](#2-trimming))
- `trimmed_for_fastqc`: temporary directory with combined trimmed readfiles for FASTQC
- `assembly`: HybPiper results (see [below](#3-assembly))
- `coverage`: output of coverage trimming step (see [below](#4-coverage-trimming-and-length-filtering))
- `seq_sets2`: sequence sets after coverage trimming and length filtering [below](#4-coverage-trimming-and-length-filtering)
- `alignments2`: aligned sequence sets after coverage trimming and length filtering (see [below](#6-alignment))
- `alignments_exon`: alignments with added exon sequences for partitioning (see [below](#7-mapping-exons-to-alignments))
- `optrimal`: working directory for dynamic alignment trimming with optrimAl (see [below](#8-gap-trimming))
- `alignments_for_editing`: output of optrimal step, will be manually edited and moved to: 
- `alignments_edited`: manually cleaned alignments (see [below](#9-manual-editing)). Contains subfolders `genetrees` for iqtree results and `done` for processed alignments, allowing batch-wise treebuilding. 
- `alignments_bad`: blacklisted alignments, moved directly from `alignments_for_editing`.
- `speciestree`: intermediate ASTRAL input and output (see [below](#10-tree-building)
- `treeshrink`: output of TreeShrink step (see [below](#11-diagnosing-remaining-errors-using-treeshrink))
- `alignments_edited2`: manually cleaned alignments after 2nd cleaning step. Structure as in `alignments_edited` (see [below](#12-manual-editing-2nd-round))
- `alignments_edited2_newspp/done`: contains cleaned alignments after adding last-minute species
- `final_tree_nofilter`: gene and species trees for the inclusive analysis

Repository location on GIS07: `~/scripts/dypsidinae`

Repository location on Macbook: `~/Documents/WOLF/PROJECTS/65 Dypsis systematics paper/~git/Dypsidinae_species_tree`

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

Trimmomatic settings used:
ILLUMINACLIP:/usr/local/bioinf/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:36 AVGQUAL:30

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

_Based on comparing trimming results, it was decided to use the first trimming settings. I.e., all downstream analyses are based on the contents of_ `trimmed`.

## 3. Assembly

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

Retrieve paralog information ([see](https://github.com/mossmatters/HybPiper/wiki/Paralogs)): 

```bash
while read i
do 
echo $i
python /usr/local/bioinf/HybPiper/paralog_investigator.py $i
done < namelist.txt
```

Paralogs found for 362 825 728 1013 168 985. 

```bash
parallel "python /usr/local/bioinf/HybPiper/paralog_retriever.py namelist.txt {} > {}.paralogs.fasta" ::: 362 825 728 1013 168 985
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

These are ready for blacklisting and alignment.

## 5. Blacklisting 

Run from `seq_sets2` to clean up sequence names:

```bash
for f in *.FNA; do (sed -i'.old' -e $'s/-[0-9]\+[p,n,s,e]* [0-9]\+-[0-9]\+[p,n,s,e]*_[0-9]\+[p,n,s,e]* [0-9]\+-[0-9]\+[p,n,s,e]*//g' $f); done
rm *.old 
```

Remove blacklisted taxa from all sequence sets and tidy up file names. 

```bash
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 0016 0056 0094 0192 0147 0200 0165 0096 0186 0064 0127 0093 -d dna -f fasta -i *FNA -u fasta
for f in *-out.fas; do (mv $f ${f/-out.fas}); done
```

_NB_: The black list is currently hard coded in this command. Add further blacklisted species to the argument `-x`. 

## 6. Alignment

Run from `seq_sets2`:

```bash
for f in reduced_*; do (linsi --thread 16 $f > ../alignments2/${f/.FNA}_aligned.fasta); done
```

## 7. Mapping exons to alignments

In `alignments2`, run: 

```bash
~/scripts/dypsidinae/exon_mapper.py
```

This creates new alignments in `alignments_exon` that contain the original alignments plus the exon sequences of the two species that had the highest recovery success at each locus. 

## 8. Gap trimming

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

Some of the alignments generated by optrimal may contain empty sequences. To remove these, run (in `optrimal`): 

```bash
for f in *.fasta;do(~/scripts/dypsidinae/noempty.py $f);done
reduced_1171_aligned_noempty.fasta has 1 empty sequences removed
reduced_120_aligned_noempty.fasta has 2 empty sequences removed
reduced_326_aligned_noempty.fasta has 2 empty sequences removed
reduced_392_aligned_noempty.fasta has 2 empty sequences removed
reduced_874_aligned_noempty.fasta has 8 empty sequences removed
reduced_938_aligned_noempty.fasta has 3 empty sequences removed
reduced_982_aligned_noempty.fasta has 1 empty sequences removed
```

Copy alignments to `alignments_for_editing`.

## 9. Manual editing

At this stage, each alignment needs to be scrutinised and cleaned by hand as follows: 

1. Move (not copy) the alignment you want to edit from `alignments_for_editing` to a place of your choice. 
2. Make all necessary edits, and save the edited version in `alignments_edited`.

If alignments are found to be overall wrong or doubtful (e.g. alignment patterns indicate the presence of paralogs/chimeric sequences), these should be moved to `alignments_bad` and excluded from further analysis. 

## 10. Tree building

Once a reasonable number of alignments has been saved in `alignments_edited`, run 

```bash
~/scripts/dypsidinae/treebuilder.sh
```

from this folder. This script will

- run `~/scripts/dypsidinae/partitioner.py` with a smoothing parameter of 10bp (i.e. ignoring any mini-partitions <10bp long) to generate RAxML-style partition files called `*_part.txt`, and remove the exon sequences from the alignment (new alignment file saved as `*_clean.fasta`)
- run iqtree with model search and 1000 fast bootstrap replicates
- move iqtree outputs and partition files to `alignments_edited/genetrees`, renaming `*.treefile` to `\*.tre` for convenience
- move the original edited alignment to `alignments_edited/done`
- remove the `*_clean.fasta`
- reroot genetrees on _Loxococcus_, collapse all internal nodes with UFBS<30%, and gather the trees ina file `genetrees.txt` in `Dypsis/speciestree`
- run ASTRAL (output: `speciestree/astral.log` and `speciestree/astral_tree.tre`)
- calculate EQP-IC (output: `speciestree/astral_tree_QS.tre`)
- rename taxa in LPP and EQP-IC trees (output: `speciestree/*renamed.tre`)

*Importantly*, this script will overwrite anything that already exists for an alignment in `alignments_edited/genetrees` or `alignments_edited/done`. This is *on purpose*, as it allows an iterative process: If you check a genetree in `alignments_edited/done` and find that it, e.g., still contains conspicuously long branches, and decide to give the alignment another round of editing, all you need to do is move it back into `alignments_edited`, make your changes, and run `treebuilder.sh` again. 

*Also*, it will overwrite the contents of `speciestree`!!

## 11. Diagnosing remaining errors using TreeShrink

Create directory `Dypsis/treeshrink`.

From `alignments_edit`, run: 

```bash
~/scripts/dypsidinae/treeshrink_prep.sh
```

This creates a folder structure suitable for TreeShrink in `treeshrink`.

From `treeshrink`, run: 

```bash
python3 ~/software/TreeShrink/run_treeshrink.py -i . -t input.tre
```

## 12. Manual editing 2nd round

All alignments that yielded gene trees with anomalously long branches (see previous step) were checked again with focus on the species flagged by TreeShrink. However, alignments in which only outgroup species were flagged by TreeShrink were not checked (assuming that these were likely false positives). Corrected alignments, as well as the alignments not flagged by TreeShrink, were moved to `alignments_edited2`.

## 13. Adding last-minute samples

Some additional samples that were sequenced at the last minute were added to the edited alignments, which were then checked again manually.

The new samples were trimmed, assembled, coverage trimmed and compiled in sequence sets separately, using the same code as described here in 1-4.

The sequence sets of the new samples are stored in `seq_sets_newspp`. 

NB sequence names in sequence sets need cleaning prior to alignment as described in 5. 

To create combined alignments, create `alignments_edited2_newspp` and run from `Dypsis`:

```bash
~/scripts/dypsidinae/add_newspp.py
```

## 14. Manual editing 3rd round

In all new alignments, the newly added species were checked for assembly- or alignment artefacts as described above, then moved to `alignments_edited2_newspp/edited`. NB: After a first round of gene tree building, alignments were checked again and edited some more. The final alignments are now in `alignments_edited2_newspp/edited_again`. 

## 15. Fix intronerate error

In seven species (0003, 0008, 0010, 0017, 0061, 0105, 0110), an unnoticed [error](https://github.com/mossmatters/HybPiper/issues/41) in `intronerate.py` resulted in missing supercontigs. For those species, the intronerate step was repeated with the development version of the script. The contigs were then added to the relevant alignments, and their alignment checked for assembly- or alignment errors as described above. These edited alignments are in `alignments_fix_intronerate`. See also [here](fix_intronerate.md). 

| SECAPR | Species | 
| -------| --------|
| 0003 | _Dypsis albofarinosa_ | 
| 0008 | _Dypsis ankirindro_ |
| 0010 | _Dypsis arenarum_ |
| 0017 | _Dypsis betsimisarakae_ |
| 0061 | _Dypsis rabepierrei_ |
| 0105 | _Dypsis mijoroana_ |
| 0110 | _Dypsis nossibensis_ |

## 16. Gather final alignments

From data folder: 

```bash
mkdir alignments_final
cp alignments_edited2_newspp/edited_again/* alignments_final
cd alignments_edited2/done
for f in *.fasta
do
	if [ ! -f "../../alignments_edited2_newspp/edited/$f" ] ; then
    	cp $f ../../alignments_final
	fi
done 
cd ../../alignments_fix_intronerate
for f in *.fasta
do
	rm ../alignments_final/${f/-out.fasta}
	cp $f ../alignments_final/${f/-out.fasta}
done
```

## 17. Prep for final tree building 

Copy alignments to `final_tree_nofilter/iqtree`. Then, run in this folder: 

```bash
for f in *.fasta; do(sed -i'.old' -e 's/ [0-9]\+ bp//g' $f); done
rm *.old
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py remove -x 0075 0076 0157 0197 0159 0164 2013 2016 0119 2017 -d dna -f fasta -i *.fasta -u fasta -g red_
rm reduced*
for f in *.fas; do (mv $f ${f#red_}ta); done
```

to remove multiple sequences of same individual and blacklist some problematic samples, namely:

 | SECAPR no | species | reason |
 | ----------|---------|--------|
 | 0119, 2013 | _D. sahanofensis_ | polyphyletic, one of them wrongly named? |
 | 2016 | "_D. ifanadianae_" | sample mix-up, most likely _Masoala_ sp. nov. |
 | 0159 | _D humblotiana_ | v. poor data, new sample placed more meaningfully |
 | 0164 | _D. lanceolata_ | v. poor data, new sample placed more meaningfully |
 
_NB_: The alignments that were affected by the [intronerate error](#15-fix-intronerate-error) already had those species removed. This needs to be addressed if those species are ever to be put back into the cleaned alignments. 
 
## 18. Tree building 2nd round

In `final_tree_nofilter/iqtree`, run: 


```bash
~/scripts/dypsidinae/partitioner.py --smoother 10
for f in *_part.txt; do (cp $f ${f/_part.txt}_clean.part); done
ls *clean.fasta | parallel -j 6 ~/software/iqtree-2.0.6-Linux/bin/iqtree2 -s {} -T AUTO -ntmax 4 -p {.}.part -B 1000
```

Build species tree: 

```bash
for f in *.treefile
do 
	~/scripts/dypsidinae/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> ../astral/genetrees.tre 
	rm temp.tre[[:digit:][:punct:]E]\+
done
cd ../astral
java -jar ~/software/Astral/astr[[:digit:]\.\-E]\+al.5.7.3.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
# ~/scripts/dypsidinae/renamer.py ../../rename.csv astral_tree.tre astral_tree_renamed.tre
java -jar ~/software/Astral/astral.5.7.3.jar -q astral_tree.tre -i genetrees.tre -o astral_tree_full_annot.tre -t 2 2> annotation.log
```
Replace:

\[q1=[\d,\.,E,-]+;q2=[\d,\.,E,-]+;q3=[\d,\.,E,-]+;f1=[\d,\.,E,-]+;f2=[\d,\.,E,-]+;f3=[\d,\.,E,-]+;pp1=[\d,\.,E,-]+;pp2=[\d,\.,E,-]+;pp3=[\d,\.,E,-]+;QC=[\d,\.,E,-]+;EN=([\d,\.,E,-]+)\]

With:

\1


## 19. Generate main tree figure

_NB_: All R scripts have to be run locally, NOT on linospadix - otherwise data loss can occur. 

Use `main_figure.R` interactively. 

## 20. Create gene tree figures for supplement

Use `genetrees2.R`.

## 21. Alignment statistics: intron vs exon statistics

In `iqtree`, run:

```bash
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py summary -f fasta -d dna -i *_clean.fasta
#mv summary.txt summary_all.txt
mkdir stats # copy alignments and partition files to separate dir for splitting into exons and introns
cp *.part stats
cp *_clean.fasta stats
cd stats
for f in *.part # reformat partition files for AMAS
do
	sed -i'.old' -e's/DNA, //g' $f
done
for f in *clean.fasta # split alignments into intron and exon
do
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py split -f fasta -d dna -i $f -l ${f/_clean.fasta}_clean.part -u fasta
done
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py summary -f fasta -d dna -i *exon-out.fas
mv summary.txt summary_exon.txt
python3 /home/au265104/.local/lib/python3.6/site-packages/amas/AMAS.py summary -f fasta -d dna -i *intron-out.fas
mv summary.txt summary_intron.txt
```

Then use `alignment_stats.R` to generate raw table. 
