#!/bin/bash

#$ -N Steps1-3
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)


md5sum ${SAMPLE}_1_paired.fastq.gz
md5sum ${SAMPLE}_2_paired.fastq.gz
