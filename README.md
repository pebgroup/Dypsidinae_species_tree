# Analysis of Dypsidinae target capture data

Wolf Eiserhardt (wolf.eiserhardt@bios.au.dk), 15 April 2020

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
- `seq_sets`: unaligned sequence sets per locus as retrieved by HybPiper `retrieve_sequences.py_
- `alignments`: aligned sequence sets
- `fasttrees`: preliminary gene trees
- `fasttrees_collapsed`: preliminary gene trees with nodes BS<0.1 collapsed
- `coverage`: output of coverage trimming step (see 7. below)

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

### Retrieve sequences: 

`python /usr/local/bioinf/HybPiper/retrieve_sequences.py /data_vol/wolf/PhyloPalms/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta . dna`

Saved screen output manually to `outstats.txt`. This should be piped instead. Can be used for crude exclusion of loci based on number of seqs retrieved. 

Copied sequence sets to `seq_sets`

## 4. Alignment (MAFFT)

Run from `seq_sets`:

`for f in *; do (linsi --thread 16 $f > ../alignments/${f/.FNA}_aligned.fasta); done`

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