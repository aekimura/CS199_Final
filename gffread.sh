#!/bin/bash

#$ -N Gffread
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

source ~/.miniconda3rc
conda activate hisat2


gffread saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o S288C.gtf


conda deactivate
