# Read checking scripts

Some scripts for tracking down the reasons for warnings in Fastqc outputs of trimmed and mapped reads. 

Wolf Eiserhardt, April 2020

`gather_reads.sh`(calling `filter_reads.py`) puts back together reads that have been distributed to targets by HybPiper, and extracts their Phred scores from the original (trimmed) read files. This is probably a really convoluted way of doing something that could have been done from the bam file...

`agfilter.py` separates out reads containing poly-A, poly-G or poly-AG, as encountered in some of the read data. 

`gibberish_counter.sh` (calling `gibberish_counter.py`) counts the percentage of reads in a read file that contain unusually high amounts of "AG" or "CT" as encountered in some of the read data. 