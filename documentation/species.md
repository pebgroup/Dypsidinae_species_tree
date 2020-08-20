# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 20 August 2020

## -1. Current tasks: 



## 0. Workspace

Data folder on GIS07: `/data_vol/wolf/Dypsis/`
- `original_data`: raw read files with original naming, cf. sampling.xlsx
- `original_data_renamed`: renamed read files for compatibility with SECAPR (see 1. below). 
- `fastqc_results`: results of fastqc check run via SECAPR
    - `raw`: fastqc results for raw reads (as in `original_data_renamed`)
    - `trimmed`: fastqc results after trimming (as in `trimmed`)
    - `trimmed2`: fastqc results after trimming (as in `trimmed2`)
- `trimmed`: trimmed reads (see [below](#2-trimming))
- `trimmed2`: trimmed reads with alternative trimming criteria (see 2. below)
- `trimmed_for_fastqc`: temporary directory with combined trimmed readfiles for FASTQC
- `assembly`: HybPiper results (see 3. below)
- `assembly_excluded`: samples that have been processed up to HybPiper (step 3) but were subsequently excluded (old outgroup, bad samples)
- `coverage`: output of coverage trimming step (see 7. below)
- `seq_sets2`: sequence sets after coverage trimming and length filtering (see 7. below)
- `alignments2`: aligned sequence sets after coverage trimming and length filtering
- `alignments_exon`: alignments with added exon sequences for partitioning
- `optrimal`: working directory for dynamic alignment trimming with optrimAl


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

## 6. Alignment (MAFFT)

Run from `seq_sets2`:

```bash
for f in reduced_*; do (linsi --thread 16 $f > ../alignments2/${f/.FNA}_aligned.fasta); done
```

## 7. Map exons to alignemnts

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

Some of the alignments generated by optrimal may contain empty sequences. To remove these, run (in `optrimal`): 

```bash
for f in *.fasta;do(~/scripts/dypsidinae/noempty.py $f);done
reduced_1171_aligned_noempty.fasta has 1 empty sequences removed
reduced_120_aligned_noempty.fasta has 2 empty sequences removed
reduced_874_aligned_noempty.fasta has 8 empty sequences removed
reduced_938_aligned_noempty.fasta has 3 empty sequences removed
```











