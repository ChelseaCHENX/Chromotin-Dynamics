#!/bin/bash -ue
TASK_TMP_DIR=$(mktemp -d -p null distiller.tmp.XXXXXXXXXX)
touch SRR6493702.lane1.GRCh38.0.bam 
bwa mem -t 1 -v 3 -SP Homo_sapiens.GRCh38.dna.primary_assembly.fa SRR6493702_1.fastq.gz SRR6493702_2.fastq.gz         | tee >(samtools view -bS > SRR6493702.lane1.GRCh38.0.bam)         | pairtools parse                --add-columns mapq             -c hg38.chrom.sizes             | pairtools sort --nproc 1                              -o SRR6493702.lane1.GRCh38.0.pairsam.gz                              --tmpdir $TASK_TMP_DIR             | cat

rm -rf $TASK_TMP_DIR
