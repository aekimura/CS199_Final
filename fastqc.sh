#!/bin/bash

#$ -N FastQC
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

source ~/.miniconda3rc
conda activate hisat2


fastqc ${SAMPLE}_1.fastq.gz
fastqc ${SAMPLE}_2.fastq.gz


conda deactivate
