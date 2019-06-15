#!/bin/bash

source ~/.miniconda3rc
conda activate hisat2


cat S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]/I/g' | awk '/^>/{print ">" ++i; next}{print}' > S288C_reference_sequence_R64-2-1_20150113_new.fsa

gffread saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o S288C.gtf

hisat2-build -f S288C_reference_sequence_R64-2-1_20150113_new.fsa S288C


conda deactivate
