#!/bin/bash
set -e

echo "generate bam using blasr"
blasr -nproc 4 reads.fasta ref.fasta -bestn 1 -sam -clipping hard -minMatch 19 | samtools sort - align.sort 

# echo "calc MD tag"
# samtools calmd -b align.sort.bam ref.fasta > align.sort.md.bam 

# echo "faidx reference"
# samtools faidx ref.fasta

echo "convert bam to m5"
python3 ../bam2m5.py align.sort.bam ref.fasta -5,6,0,5 align.sort.bam.m5 

echo "generate m5 using blasr as standard"
blasr -nproc 4 reads.fasta ref.fasta -bestn 1 -m 5 -minMatch 19 > align.m5  

echo "compare..."
/bin/rm -f m5.diff
sort -k 6,6 -k 8,8n -k 1,1 align.m5 > align.sorted.m5 
sort -k 6,6 -k 8,8n -k 1,1 align.sort.bam.m5 > align.bam.sort.m5 
diff align.sorted.m5 align.bam.sort.m5 > m5.diff
if [ -s m5.dff ]; then
    echo "two m5 file differ, please check"
else
    echo "everything is fine"
fi

