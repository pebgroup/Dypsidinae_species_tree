# Re-run species tree without manually selected rogue taxa

Wolf Eiserhardt (wolf.eiserhardt@bio.au.dk), 23 November 2021

Taxa: _D. rakotonasoloi_, _D. remotiflora_, _D. rivularis_, _D. curtisii_, _D. corniculata_, _D. montana_ . 

Run locally (`Users/au265104/OneDrive - Aarhus Universitet/ANALYSIS/Dypsis/`):

```
cp -rf final_tree_nofilter final_tree_nofilter_norogues2
cd final_tree_nofilter_norogues2/iqtree
rm *out_clean.*
rm *.txt

python3 ../../AMAS/amas/AMAS.py remove -x 0178 0180 0117 0154 0203 0169 2053 -d dna -f fasta -i *.fasta -u fasta -g red_
rm reduced*
for f in *.fas; do (mv $f ${f#red_}ta); done

for f in *.fasta; do (mv $f ${f/-out.fasta}); done

../../Dypsidinae_species_tree/partitioner.py --smoother 10

rm ../iqtree_32e/reduced_*
cp reduced_32e_aligned.fasta-out_clean.fasta ../iqtree_32e/
rm reduced_32e*

for f in *_part.txt; do (cp $f ${f/_part.txt}_clean.part); done

```

Copy `final_tree_nofilter_norogues2` to GenomeDK. There, run: 

```
gwf status
gwf run

```

Also run batch script in iqtree_32e

Transfer folder back and move files from iqtree_32e to main iqtree folder. Delete contents of astral and run:

```
for f in *.treefile
do 
	python3 ../../Dypsidinae_species_tree/rooter.py $f
	nw_ed temp.tre 'i & (b<30)' o >> ../astral/genetrees.tre 
	rm temp.tre
done
cd ../astral

java -jar ~/software/Astral/astral.5.7.8.jar -i genetrees.tre -o astral_tree.tre  2> astral.log
../../Dypsidinae_species_tree/renamer.py ../../rename.csv astral_tree.tre astral_tree_renamed.tre

``


