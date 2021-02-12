# Fix intronerate error in 7 species

Create list of affected samples: 

`assembly/namelist_redo_int`

0110
0017
0105
0061
0010
0008
0003

Re-run intronerate for those samples with development version: 

```bash
while read name
do
	echo $name >> intronerate_out_redo_int.txt
	python /usr/local/bioinf/HybPiper/intronerate_dev.py --prefix $name &>> intronerate_out_redo_int.txt
done < namelist_redo_int
```

Run `coverage.py` and `samples2genes.py` for those samples. Results in `coverage_fix_intronerate` and `seq_sets2_fix_intronerate`.

Create `final_tree_nofilter_fix_intronerate/iqtree`. Put alignments there as described in main protocol. 

From `seq_sets2_fix_intronerate`, run:

```bash
for f in *.FNA; do (sed -i'.old' -e $'s/-[0-9]\+[p,n,s,e]* [0-9]\+-[0-9]\+[p,n,s,e]*_[0-9]\+[p,n,s,e]* [0-9]\+-[0-9]\+[p,n,s,e]*//g' $f); done
rm *.old 
mkdir toadd
mkdir added
~/scripts/dypsidinae/fix_intronerate.py

```