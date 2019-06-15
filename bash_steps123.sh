#!/bin/bash

#$ -N Steps1-3
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)
GZ_DIR="S288C_data/samples/"
GZ_1="_1_paired.fastq.gz"
GZ_2="_2_paired.fastq.gz"

source ~/.miniconda3rc
conda activate hisat2


hisat2 -p 1 --dta -x S288C_data/indexes/S288C -1 ${GZ_DIR}${SAMPLE}${GZ_1} -2 ${GZ_DIR}${SAMPLE}${GZ_2} -S ${SAMPLE}_S288C.sam

samtools sort -@ 1 -o ${SAMPLE}_S288C.bam ${SAMPLE}_S288C.sam 

stringtie -p 1 -G S288C_data/genes/S288C.gtf -o ${SAMPLE}_S288C.gtf -l ${SAMPLE} ${SAMPLE}_S288C.bam


conda deactivate
