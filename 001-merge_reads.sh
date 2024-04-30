#!/bin/bash


module use /tungstenfs/groups/gbioinfo/Appz/easybuild/modules/all
module use /tungstenfs/groups/gbioinfo/Appz/modules

module load PEAR/0.9.11-GCC-10.2.0
module load cutadapt/1.18-foss-2018b-Python-3.6.6



FORWARD='/tungstenfs/groups/gbioinfo/seqdata/2752F1-1_210310_NB501735_0670_AHJ2YFAFX2_ATCACG-NoIndex_L000_R1_001.fastq.gz'
REVERSE='/tungstenfs/groups/gbioinfo/seqdata/2752F1-1_210310_NB501735_0670_AHJ2YFAFX2_ATCACG-NoIndex_L000_R2_001.fastq.gz'
TRIMMED_FORWARD='tmp/trimmed_jun_f.fastq.gz'
TRIMMED_REVERSE='tmp/trimmed_jun_r.fastq.gz'
MERGED='001-merged_reads/merged_jun'

# error rate of 0.1 will allow 1 errors out of 12 bases
cutadapt -g '^CCTAGG' -G 'GGATCCAAGCTT' -e 0.1 -u -1 -U -1 -O 6 --max-n=0 -j 4 --discard-untrimmed -o $TRIMMED_FORWARD -p $TRIMMED_REVERSE $FORWARD $REVERSE

# expected length is 254. Let's give +- 5
pear -f $TRIMMED_FORWARD -r $TRIMMED_REVERSE -m 259 -v 12 -n 249 -u 0 -j 48 -p 0.05 -o $MERGED


