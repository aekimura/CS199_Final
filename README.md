# Background

For this project, our main goal was to compare levels of gene expression among different mutants and growth conditions.  First, Trimmomatic was used to improve the quality of the 8 Illumina S. cerevisiae S288C fastq samples by dropping the low quality reads and identifying only paired reads.  Then HISAT was applied to map the sample fastq reads to the reference genome. Then the alignments were analyzed by StringTie in order to generate transcripts by estimating expression levels for each gene and isoform in the alignment.  The generated transcripts were then merged to a single annotation file by StringTie to create a common set of transcripts for all the samples. Finally, Ballgown determined which genes and transcripts are differentially expressed taking into account the different conditions of the samples.  For our analysis, the samples were divided into two groups by growth temperature (30°C and 37°C), with the first group grown at 30°C contains wild-type and mutant Isw2 samples and the second group grown at 37°C contains mutant Rsc and mutant Ino80 samples.

**Directory Structure**
.
|-- ballgown
|   |-- temp30
|   |   |-- SRR1257637
|   |   |   |-- e2t.ctab
|   |   |   |-- e_data.ctab
|   |   |   |-- i2t.ctab
|   |   |   |-- i_data.ctab
|   |   |   |-- SRR1257637_S288C.gtf
|   |   |   `-- t_data.ctab
|   |   |-- SRR1257640
|   |   |   |-- e2t.ctab
|   |   |   |-- e_data.ctab
|   |   |   |-- i2t.ctab
|   |   |   |-- i_data.ctab
|   |   |   |-- SRR1257640_S288C.gtf
|   |   |   `-- t_data.ctab
|   |   |-- SRR1257793
|   |   |   |-- e2t.ctab
|   |   |   |-- e_data.ctab
|   |   |   |-- i2t.ctab
|   |   |   |-- i_data.ctab
|   |   |   |-- SRR1257793_S288C.gtf
|   |   |   `-- t_data.ctab
|   |   `-- SRR1259267
|   |       |-- e2t.ctab
|   |       |-- e_data.ctab
|   |       |-- i2t.ctab
|   |       |-- i_data.ctab
|   |       |-- SRR1259267_S288C.gtf
|   |       `-- t_data.ctab
|   `-- temp37
|       |-- SRR1514794
|       |   |-- e2t.ctab
|       |   |-- e_data.ctab
|       |   |-- i2t.ctab
|       |   |-- i_data.ctab
|       |   |-- SRR1514794_S288C.gtf
|       |   `-- t_data.ctab
|       |-- SRR1514795
|       |   |-- e2t.ctab
|       |   |-- e_data.ctab
|       |   |-- i2t.ctab
|       |   |-- i_data.ctab
|       |   |-- SRR1514795_S288C.gtf
|       |   `-- t_data.ctab
|       |-- SRR1515155
|       |   |-- e2t.ctab
|       |   |-- e_data.ctab
|       |   |-- i2t.ctab
|       |   |-- i_data.ctab
|       |   |-- SRR1515155_S288C.gtf
|       |   `-- t_data.ctab
|       `-- SRR1515156
|           |-- e2t.ctab
|           |-- e_data.ctab
|           |-- i2t.ctab
|           |-- i_data.ctab
|           |-- SRR1515156_S288C.gtf
|           `-- t_data.ctab
|-- bash_step6.sh
|-- bash_steps123.sh
|-- bash_steps45.sh
|-- download_data.sh
|-- fastqc_results
|   |-- SRR1257637_1_fastqc.html
|   |-- SRR1257637_1_fastqc.zip
|   |-- SRR1257637_2_fastqc.html
|   |-- SRR1257637_2_fastqc.zip
|   |-- SRR1257640_1_fastqc.html
|   |-- SRR1257640_1_fastqc.zip
|   |-- SRR1257640_2_fastqc.html
|   |-- SRR1257640_2_fastqc.zip
|   |-- SRR1257793_1_fastqc.html
|   |-- SRR1257793_1_fastqc.zip
|   |-- SRR1257793_2_fastqc.html
|   |-- SRR1257793_2_fastqc.zip
|   |-- SRR1259267_1_fastqc.html
|   |-- SRR1259267_1_fastqc.zip
|   |-- SRR1259267_2_fastqc.html
|   |-- SRR1259267_2_fastqc.zip
|   |-- SRR1514794_1_fastqc.html
|   |-- SRR1514794_1_fastqc.zip
|   |-- SRR1514794_2_fastqc.html
|   |-- SRR1514794_2_fastqc.zip
|   |-- SRR1514795_1_fastqc.html
|   |-- SRR1514795_1_fastqc.zip
|   |-- SRR1514795_2_fastqc.html
|   |-- SRR1514795_2_fastqc.zip
|   |-- SRR1515155_1_fastqc.html
|   |-- SRR1515155_1_fastqc.zip
|   |-- SRR1515155_2_fastqc.html
|   |-- SRR1515155_2_fastqc.zip
|   |-- SRR1515156_1_fastqc.html
|   |-- SRR1515156_1_fastqc.zip
|   |-- SRR1515156_2_fastqc.html
|   `-- SRR1515156_2_fastqc.zip
|-- fastqc.sh
|-- fastq_download.sh
|-- gffread.sh
|-- hisat2_env.yml
|-- index_builder.sh
|-- merged.annotated.gtf
|-- merged.loci
|-- merged.stats
|-- merged.stringtie_merged.gtf.refmap
|-- merged.stringtie_merged.gtf.tmap
|-- merged.tracking
|-- preprocessing.sh
|-- r_env.yml
|-- Rplots.pdf
|-- R_temp30.sh
|-- R_temp37.sh
|-- S288C_data
|   |-- genes
|   |   `-- S288C.gtf
|   |-- genome
|   |   |-- S288C_reference_sequence_R64-2-1_20150113_new.fsa
|   |   `-- S288C_reference_sequence_R64-2-1_20150113_new.fsa.fai
|   |-- indexes
|   |   |-- S288C.1.ht2
|   |   |-- S288C.2.ht2
|   |   |-- S288C.3.ht2
|   |   |-- S288C.4.ht2
|   |   |-- S288C.5.ht2
|   |   |-- S288C.6.ht2
|   |   |-- S288C.7.ht2
|   |   `-- S288C.8.ht2
|   |-- mergelist.txt
|   |-- phenotype1.csv
|   |-- phenotype2.csv
|   `-- samples
|       |-- md5sum.o6868713.1
|       |-- md5sum.o6868713.2
|       |-- md5sum.o6868713.3
|       |-- md5sum.o6868713.4
|       |-- md5sum.o6868713.5
|       |-- md5sum.o6868713.6
|       |-- md5sum.o6868713.7
|       |-- md5sum.o6868713.8
|       |-- md5sum.sh
|       |-- samples.txt
|       |-- SRR1257637_1_paired.fastq.gz
|       |-- SRR1257637_2_paired.fastq.gz
|       |-- SRR1257640_1_paired.fastq.gz
|       |-- SRR1257640_2_paired.fastq.gz
|       |-- SRR1257793_1_paired.fastq.gz
|       |-- SRR1257793_2_paired.fastq.gz
|       |-- SRR1259267_1_paired.fastq.gz
|       |-- SRR1259267_2_paired.fastq.gz
|       |-- SRR1514794_1_paired.fastq.gz
|       |-- SRR1514794_2_paired.fastq.gz
|       |-- SRR1514795_1_paired.fastq.gz
|       |-- SRR1514795_2_paired.fastq.gz
|       |-- SRR1515155_1_paired.fastq.gz
|       |-- SRR1515155_2_paired.fastq.gz
|       |-- SRR1515156_1_paired.fastq.gz
|       `-- SRR1515156_2_paired.fastq.gz
|-- S288C_reference_sequence_R64-2-1_20150113.fsa
|-- saccharomyces_cerevisiae_R64-2-1_20150113.gff
|-- samples.txt
|-- setup.sh
|-- SigDE.txt
|-- SRR1257637_S288C.bam
|-- SRR1257637_S288C.gtf
|-- SRR1257637_S288C.sam
|-- SRR1257640_S288C.bam
|-- SRR1257640_S288C.gtf
|-- SRR1257640_S288C.sam
|-- SRR1257793_S288C.bam
|-- SRR1257793_S288C.gtf
|-- SRR1257793_S288C.sam
|-- SRR1259267_S288C.bam
|-- SRR1259267_S288C.gtf
|-- SRR1259267_S288C.sam
|-- SRR1514794_S288C.bam
|-- SRR1514794_S288C.gtf
|-- SRR1514794_S288C.sam
|-- SRR1514795_S288C.bam
|-- SRR1514795_S288C.gtf
|-- SRR1514795_S288C.sam
|-- SRR1515155_S288C.bam
|-- SRR1515155_S288C.gtf
|-- SRR1515155_S288C.sam
|-- SRR1515156_S288C.bam
|-- SRR1515156_S288C.gtf
|-- SRR1515156_S288C.sam
|-- stringtie_merged.gtf
|-- temp30_results
|   |-- MSTRG.194.1.bed
|   |-- MSTRG.194.1.fa.out
|   |-- MSTRG.195.3.bed
|   |-- MSTRG.195.3.fa.out
|   |-- MSTRG.318.1.bed
|   |-- MSTRG.318.1.fa.out
|   |-- MSTRG.347.1.bed
|   |-- MSTRG.347.1.fa.out
|   |-- MSTRG.369.3.bed
|   |-- MSTRG.369.3.fa.out
|   |-- MSTRG.64.2.bed
|   |-- MSTRG.64.2.fa.out
|   |-- MSTRG.64.3.bed
|   `-- MSTRG.64.3.fa.out
|-- temp30_significant.csv
|-- temp37_results
|   |-- MSTRG.106.1.bed
|   |-- MSTRG.106.1.fa.out
|   |-- MSTRG.13.1.bed
|   |-- MSTRG.13.1.fa.out
|   |-- MSTRG.144.2.bed
|   |-- MSTRG.144.2.fa.out
|   |-- MSTRG.179.3.bed
|   |-- MSTRG.179.3.fa.out
|   |-- MSTRG.187.1.bed
|   |-- MSTRG.187.1.fa.out
|   |-- MSTRG.276.4.bed
|   |-- MSTRG.276.4.fa.out
|   |-- MSTRG.28.1.bed
|   |-- MSTRG.28.1.fa.out
|   |-- MSTRG.33.1.bed
|   |-- MSTRG.33.1.fa.out
|   |-- MSTRG.43.5.bed
|   |-- MSTRG.43.5.fa.out
|   |-- MSTRG.89.1.bed
|   `-- MSTRG.89.1.fa.out
|-- temp37_significant.csv
|-- trimmomatic_env.yml
|-- trimmomatic.sh
|-- unpaired_fastq
|   |-- SRR1257637_1_unpaired.fastq.gz
|   |-- SRR1257637_2_unpaired.fastq.gz
|   |-- SRR1257640_1_unpaired.fastq.gz
|   |-- SRR1257640_2_unpaired.fastq.gz
|   |-- SRR1257793_1_unpaired.fastq.gz
|   |-- SRR1257793_2_unpaired.fastq.gz
|   |-- SRR1259267_1_unpaired.fastq.gz
|   |-- SRR1259267_2_unpaired.fastq.gz
|   |-- SRR1514794_1_unpaired.fastq.gz
|   |-- SRR1514794_2_unpaired.fastq.gz
|   |-- SRR1514795_1_unpaired.fastq.gz
|   |-- SRR1514795_2_unpaired.fastq.gz
|   |-- SRR1515155_1_unpaired.fastq.gz
|   |-- SRR1515155_2_unpaired.fastq.gz
|   |-- SRR1515156_1_unpaired.fastq.gz
|   `-- SRR1515156_2_unpaired.fastq.gz
`-- untrimmed_fastq
    |-- SRR1257637_1.fastq.gz
    |-- SRR1257637_2.fastq.gz
    |-- SRR1257640_1.fastq.gz
    |-- SRR1257640_2.fastq.gz
    |-- SRR1257793_1.fastq.gz
    |-- SRR1257793_2.fastq.gz
    |-- SRR1259267_1.fastq.gz
    |-- SRR1259267_2.fastq.gz
    |-- SRR1514794_1.fastq.gz
    |-- SRR1514794_2.fastq.gz
    |-- SRR1514795_1.fastq.gz
    |-- SRR1514795_2.fastq.gz
    |-- SRR1515155_1.fastq.gz
    |-- SRR1515155_2.fastq.gz
    |-- SRR1515156_1.fastq.gz
    `-- SRR1515156_2.fastq.gz

# Pipeline

[**1. Environment Setup**](#1-environment-setup)

[**2. Data Download**](#2-data-download)

[**3. Preprocessing**](#3-preprocessing)

[**4. Map Sample Reads, Assemble Transcripts, and Generate Ballgown Tables**](#4-map-reads)

[**5. R Analysis**](#5-r-analysis)

[**6. Analyze Significant Transcripts**](#6-analyze-significant-transcripts)

# Pipeline Steps

### 1. Environment Setup

>a) Install Miniconda: First Miniconda must be installed and then the provided .yml files can be used to create the environments to complete the rest of the pipeline steps.  To download and set-up Miniconda, enter the following commands in order as they are required. Statments in all caps denote actions to complete.  [Online installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) can also be followed to set up Miniconda if desired.
    
 ```ruby
 qrsh
 mkdir -p bin
 wget https://docs.conda.io/en/latest/miniconda.html/Miniconda3-latest-Linux-x86_64.sh
 pwd >> COPY DIRECTORY
 bash Miniconda3-latest-Linux-x86_64.sh
    COPIED DIRECTORY
    yes >> COPY BACKUP LOCATION
 cat BACKUP LOCATION >> COPY FROM THE END TO WHERE IT SAYS "added by miniconda"
 nano ~/.miniconda3rc >> PASTE COPIED TEXT
 mv ~/.bashrc-miniconda3.bak ~/.bashrc
 source ~/.miniconda3rc
 conda
 conda config --add channels defaults
 conda config --add channels bioconda
 conda config --add channels conda-forge
 ```
 
 >b) Create Environments: After Miniconda set-up is completed do not exit conda.  Then create the three environments (Hisat2, Trimmomatic, and R) can be created using the .yml files.
    
 ```ruby
 conda env create -f hisat2_env.yml
 conda env create -f trimmomatic_env.yml
 conda env create -f r_env.yml
 ```

>c) Download R packages: To complete the setup of the R environment, a few packages need to be downloaded into R.  They include ballgown which is used to complete the final data analysis steps, RSkittleBrewer which controls colors for the output, genefilter which is used for quick statistaical calculations, dplyr which sorts and arranges results, and devtools which is used in the instilation of packages.  The following code should be executed line by line to install the packages, folloing any instructions that appear during the download process.

```
conda activate r_env
R
install.packages("devtools")
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown","genefilter","dplyr","devtools"))
```

### 2. Data Download

>a) Reference Genome: For this experiment the genome [Saccharomyces cerevisiae strain S288C release R64-2-1](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/) was used as a reference. The download of the required reference genome files (fasta and gff) is handled by the script "reference_download.sh".

```
wget https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz
tar xvzf S288C_reference_genome_R64-2-1_20150113.tgz
cd S288C_reference_genome_R64-2-1_20150113
mv saccharomyces_cerevisiae_R64-2-1_20150113.gff ../../
mv S288C_reference_sequence_R64-2-1_20150113.fsa ../../
```

>b) Datasets: Illumina datasets were obtained by searching SRA-Explorer for relevent Illumina datasets for this project.  The data used for this experiment had to be Illumina paired-end data from Saccharomyces cerevisiae strain S288C. Illumina data downloading is handled by the script "fastq_download.sh".

```
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/000/SRR1257640/SRR1257640_1.fastq.gz -o SRR1257640_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/000/SRR1257640/SRR1257640_2.fastq.gz -o SRR1257640_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/003/SRR1257793/SRR1257793_1.fastq.gz -o SRR1257793_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/003/SRR1257793/SRR1257793_2.fastq.gz -o SRR1257793_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1259267/SRR1259267_1.fastq.gz -o SRR1259267_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1259267/SRR1259267_2.fastq.gz -o SRR1259267_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1514795/SRR1514795_1.fastq.gz -o SRR1514795_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1514795/SRR1514795_2.fastq.gz -o SRR1514795_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1515155/SRR1515155_1.fastq.gz -o SRR1515155_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1515155/SRR1515155_2.fastq.gz -o SRR1515155_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1257637/SRR1257637_1.fastq.gz -o SRR1257637_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1257637/SRR1257637_2.fastq.gz -o SRR1257637_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/006/SRR1515156/SRR1515156_1.fastq.gz -o SRR1515156_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/006/SRR1515156/SRR1515156_2.fastq.gz -o SRR1515156_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/004/SRR1514794/SRR1514794_1.fastq.gz -o SRR1514794_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/004/SRR1514794/SRR1514794_2.fastq.gz -o SRR1514794_2.fastq.gz
```

### 3. Preprocessing

>a) Fastq Quality Check: Fastqc version 0.11.8 was used to measure the phred quality scores across the bases of the fastq samples.  The script "fastqc.sh" measures the quality for all the sample fastq files downloaded. Because "samples.txt" contains the name of each sample the script will generate the results for all samples if executed correctly (with SGE_TASK_ID).  The results are included in the folder "fastqc_results".  

```
source ~/.miniconda3rc
conda activate hisat2 

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

fastqc ${SAMPLE}_1.fastq.gz
fastqc ${SAMPLE}_2.fastq.gz

conda deactivate
```

![fastqc results example](https://github.com/aekimura/CS199_Final/blob/master/images/fastqc.png)

>b) Clean Reference Genome Fasta: The script "clean_reference.sh" removes unwanted characters from the fasta headers of the reference genome and replaces them with numbers.

```
cat S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]/I/g' | awk '/^>/{print ">" ++i; next}{print}' > S288C_reference_sequence_R64-2-1_20150113_new.fsa
```

>c) Index Building: Hisat2 v. 2.1.0 builds an index for the reference genome which is utilized later. The relevant script is "index_builder.sh".

```
source ~/.miniconda3rc
conda activate hisat2 

hisat2-build -f S288C_reference_sequence_R64-2-1_20150113_new.fsa S288C

conda deactivate
```

>d) GFF to GTF Conversion: Gffread v. 0.11.4 was used to convert the reference gff file to a gtf file for use later. This is handled by the script "gffread.sh"

```
source ~/.miniconda3rc
conda activate hisat2 

gffread saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o S288C.gtf

conda deactivate
```

>e) Fastq Trimming: Trimmomatic v. 0.39 was used to prepare the 8x Illumina paired-end reads obtained from the Sequence Read Archive for alignment.  Trimmomatic was used to remove adapters, to remove leading and trailing N bases, to scan reads in windows 4 bases long and remove any where the average base quality score is below 15, to remove reads shorter than 5 bases, and to remove reads with an average quality of less than 20.  As with the Fastqc script above the following Trimmomatic script will generate results for all the samples named in "samples.txt".  This step was completed by executing "trimmomatic.sh".

```
SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

source ~/.miniconda3rc
conda activate trimmomatic

java -jar trimmomatic-0.36.jar \
PE -phred33 \
${SAMPLE}_1.fastq.gz \
${SAMPLE}_2.fastq.gz \
${SAMPLE}_1_paired.fastq.gz \
${SAMPLE}_1_unpaired.fastq.gz \
${SAMPLE}_2_paired.fastq.gz \
${SAMPLE}_2_unpaired.fastq.gz \
ILLUMINACLIP:data/users/aekimura/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:0 \
TRAILING:0 \
SLIDINGWINDOW:4:15 \
MINLEN:5 \
AVGQUAL:20

conda deactivate
```

>f) Text File Creation: A few text files that are necessary to the running of the analysis steps of the pipline needed to be manually created at this point.  "samples.txt", "mergelist.txt", and "phenodata.csv" were all generated to fit the samples that were downloaded.  "samples.txt" contains a list of the sample names, one per line.  "mergelist.txt" contains the names of the GTF files that are to be generated later by Stringtie.  "phenodata.csv" contains information about the Illumina samples that were downloaded.  All three files are available to download on github.  

>g) Move Files:  In the final step before beginning to analyze the sample data, the files need to be arranged into the correct directories so that the scripts can be executed.  Directory creation and file moving are handled by the script "setup.sh".  

```
mkdir S288C_data

cd S288C_data

mkdir genes
mkdir genome
mkdir indexes
mkdir samples

cd ../

mv *_new.fsa /S288C_data/genome
mv *.ht2 /S288C_data/indexes
mv *_paired.fastq.gz /S288C_data/samples
mv S228C.gtf /S288C_data/genes
mv mergelist.txt S288C_data
mv S228C_phenodata.csv S288C_data
```

### 4. Map Reads

> Parts A, B, and C were combined into a single script, "bash_steps123.sh" which will generate SAM, BAM, and GTF files for all of the samples named in "samples.txt".

>a) Map Sample Reads to Reference Genome: Hisat2 v. 2.1.0 uses the indicies built in part C of the preprocessing stage in order to generate a Sequence Alignment Map (SAM file) aligning the paired fastq reads from each sample to the reference genome.

>b) Sort and Convert the SAM files to BAM files: Samtools v. 1.9 sorts the SAM file generated by Hisat2 and generates a binary version of the SAM file (BAM file) for each of the samples.

>c) Assemble Transcripts For Each Sample: Stringtie v. 1.3.6 uses the BAM file generated by Samtools to assemble a transcript (GTF file) for each sample.

```
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
```

> Parts D and E are combined into a single script named "bash_steps45.sh" which generates a single merged GTF file for all the samples and a few merged comparison files. 

>d) Merge Transcripts From All Samples: Stringtie v. 1.3.6 uses "mergelist.txt" containing the names of all the GTF files generated for each sample and the reference GTF file to create a merged GTF file for all of the samples.

>e) Compare Transcripts to Reference Annotation: Gffcompare v. 0.11.2 compares the merged GTF file generated by Stringtie and the reference GTF file and generates a report of how similar they are.  The output files begin with whatever immediately follows the "-o" which in this case is "merged".  The output files include "merged.annotated.gtf", "merged.loci", "merged.stats", "merged.stringtie_merged.gtf.refmap", "merged.stringtie_merged.gtf.tmap", and "merged.tracking".  This part is optional and is not necessary to complete the rest of the pipeline.

```
source ~/.miniconda3rc
conda activate hisat2

stringtie --merge -p 1 -G S288C_data/genes/S288C.gtf -o stringtie_merged.gtf  S288C_data/mergelist.txt
gffcompare -r S288C_data/genes/S288C.gtf -G -o merged stringtie_merged.gtf

conda deactivate
```

>e) Estimate Transcript Abundances and Create Table Counts for Ballgown: Stringtie v. 1.3.6 uses the merged GTF file generated by Stringtie in part D and for each sample, generates a GTF file and table counts in a new directory called "ballgown".  This step is accomplished by executing the script "bash_step6.sh".  

```
SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

source ~/.miniconda3rc
conda activate hisat2

stringtie -e -B -p 1 -G stringtie_merged.gtf -o ballgown/${SAMPLE}/${SAMPLE}_S288C.gtf ${SAMPLE}_S288C.bam

conda deactivate
```

>f) Move Ballgown Data: For this experiment, in order to analyze the sample data, the samples were split into two groups and thus the data generated by Stringtie must be separated into two directories.  This is accomplished with the script "move_ballgown.sh".

```
cd ballgown
mkdir temp30
mkdir temp37

mv SRR1257637 temp30
mv SRR1257640 temp30
mv SRR1257793 temp30
mv SRR1259267 temp30
mv * temp37
```

### 5. R Analysis

>The following parts A-R need to be executed for each comparitive analysis of the sample data.  For this experiment, the samples were split into two categories, those grown at 30 degrees Celcius and those grown at 37 degrees Celcius.  The script for the 30 degrees samples is "R_temp30.sh" and the one for the 37 degrees samples is "R_temp37.sh".  The scripts can be executed as they are or by inputting the commands into the R environment individually.  The graphs generated by running this script are placed in a PDF file named "Rplots.pdf".  The results from running the script for the first set of samples must be saved before it is overwritten by running the script for the second set of samples.

>a) Load Required R Packages: The packages that are required for analysis that have previously been downloaded into R should be loaded in R so that they can be used.  The packages are as described in the Enviromental setup stage part C.  

```
R
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
```

>b) Read In Phenotype Data:  This step takes the phenotype data CSV file that was manually created in the preprocessing stage and creates a variable called pheno_data.  

```
pheno_data = read.csv("S288C_data/phenotype1.csv")
```

>c) Read In Expression Data: This step takes in the expression data that was calculated and placed into tables by Stringtie in part E of the read mapping section.  It takes the sample data from where it is stored in the ballgown directory, the common part of the sample names, and the phenotype data that was was read in in the previous step.  The data is stored in the variable bg_genome. 

```
bg_genome = ballgown(dataDir = "ballgown/temp30", samplePattern = "SRR", pData = pheno_data) 
```

>d) Filtering: This step filters out any genes that do not have an expression rate greater than one.  This is important because RNA-seq data often contains genes that have counts of one or lower which are not what we are looking for in our analysis.  

```
bg_genome_filt = subset(bg_genome, "rowVars(texpr(bg_genome)) > 1", genomesubset=TRUE)
```

>e) Transcript Identification: This step uses stattest to identify transcripts that show significantly different expression between the samples.  In this case, the stattest is looking for differences between the different populations of the temp30 samples (WT and Isw2).  getFC is set to True so that the confounder-adjusted fold change is also reported.  

```
results_transcripts = stattest(bg_genome_filt, feature="transcript",covariate="population", getFC=TRUE, meas="FPKM")
```

>f) Gene Identification:  This step uses stattest as before but this time to separate out any genes that show significantly different expression between the different populations. 

```
results_genes = stattest(bg_genome_filt, feature="gene", covariate="population", getFC=TRUE, meas="FPKM")
```

>g) Add gene names and gene IDs to the results_transcripts data frame.  In this experiment gene IDs are vitally important to extracting the transcript sequences of the most significant results.    

```
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_genome_filt), geneIDs=ballgown::geneIDs(bg_genome_filt), results_transcripts)
```

>h) Sort the results from the smallest P value to the largest and write the top 10 results to a CSV.  The transcripts that have the smallest P values are the ones whose expression differs the most between the samples.  This means that analysis of these transcripts may relate to genes that are differentially expressed between the populations of the samples.  The file containing these transcripts will be used later to extract the transcripts.  

```
results_transcripts = arrange(results_transcripts,pval) 
top_results = head(results_transcripts, n=10)
write.csv(top_results, "temp30_results/temp30_significant.csv")
```

>i) Plot Distribution of transcript count per gene:  This step will generate a histogram that maps the number of transcripts that belong to each gene.  The expected result is that most genes will have one transcript that belongs to it, but others will have more.  

```
transcript_gene_table = indexes(bg_genome)$t2g
counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
temp <- hist(counts, breaks=50, col="bisque4", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcripts per gene", ylab="Log10 of Frequency", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL) 
```

>j) Plot Distribution of transcript lengths:  This step generates another histogram that maps out the lengths of this transcripts.  

```
full_table <- texpr(bg_genome, 'all')
temp <- hist(full_table$length, breaks=50, col="steelblue", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcript length (bp)", ylab="Log10 of Frequency", main="Distribution of transcript lengths")
```

![](https://github.com/aekimura/CS199_Final/blob/master/images/transcript_count_length.png)

>k) FPKM values.  This step will provide an output of the maximum FPKM values for each of the temp30 samples.  FPKM stands for fragments per kilobase of exon model per million reads mapped.  

```
gene_expression = as.data.frame(gexpr(bg_genome_filt))
colnames(gene_expression) <- c("SRR1257637","SRR1257640","SRR1257793","SRR1259267")
max(gene_expression[,"SRR1257637"])
max(gene_expression[,"SRR1257640"])
max(gene_expression[,"SRR1257793"])
max(gene_expression[,"SRR1259267"])
```

>l) Plot Distribution of FPKM for all 4 libraries: For the four temp30 samples, this step plots a boxplot of the FPKM values. 

```
data_colors=c("tomato1","tomato2","wheat1","wheat2")
min_nonzero = 1
data_columns=c(1:4)
short_names=c("WT-1","WT-2","Isw2-1","Isw2-2")
boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 4 libraries")
```

![](https://github.com/aekimura/CS199_Final/blob/master/images/FPKM.png)

>m) Plot expression values for a pair of wild-type replicates: This step generates a plot of the FPKM values of the replicates against each other.  First the two WT replicates, then the two Isw2 replicates.  This comparison will assess how similar the data from the replicates are to one another.  
```
x = gene_expression[,"SRR1257637"]
y = gene_expression[,"SRR1257640"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (WT, Replicate 1)", ylab="FPKM (WT, Replicate 2)", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

x = gene_expression[,"SRR1257793"]
y = gene_expression[,"SRR1259267"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (Isw2, Replicate 1)", ylab="FPKM (Isw2, Replicate 2)", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")
```

![](https://github.com/aekimura/CS199_Final/blob/master/images/replicate.png)

>n) MDS distance plot: This step generates a plot that compares the correlation between all of the temp30 samples.  It will sum the FPKM values for each sample and then determine the relative distances between the libraries of each of the samples.  The closer the samples are the more similar they are.  

```
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
i = which(gene_expression[,"sum"] > 5)
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)
```

![](https://github.com/aekimura/CS199_Final/blob/master/images/MDS.png)

>o) View the distribution of differential expression values with significance p < 0.1: In this step, a histogram of differential expression values with a p value of less than 0.1 are shown.  

```
bg_table = texpr(bg_genome_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.1)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) WT vs Isw2", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
```

![](https://github.com/aekimura/CS199_Final/blob/master/images/differential_expression_histogram.png)

>p) Make plots of individual transcripts across samples: In this step can be changed to plot for the most significant transcript by changing the gene ID number to those identified in part H. 

```
ballgown::transcriptNames(bg_genome)[582] ## "MSTRG.278.4" 
plot(fpkm[12,] ~ pheno_data$population, border=c(1,2), main=paste(ballgown::transcriptNames(bg_genome)[582]),pch=19, xlab="Population", ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$population)), col=as.numeric(pheno_data$population))
```

![](https://github.com/aekimura/CS199_Final/blob/master/images/individual.png)

>q) Display the grand expression values and mark those that are significantly differentially expressed: This step generates a plot of the log2 FPKM values for the two temp30 populations against each other highlighting those that are significantly differentially expressed.  

```
gene_expression[,"WT"]=apply(gene_expression[,c(1, 2)], 1, mean)
gene_expression[,"Isw2"]=apply(gene_expression[,c(3, 4)], 1, mean)
x=log2(gene_expression[,"WT"]+min_nonzero)
y=log2(gene_expression[,"Isw2"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="WT FPKM (log2)", ylab="Isw2 FPKM (log2)", main="WT vs Isw2 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)
```

>r) View Significant Results: This step displays the top 25 results after sorting results of significant transcripts with a log2 fold-change >= 2, putting them in order by increasing P value, and then breaking ties using fold-change.  

```
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 1)
sig_tn_de = sigp[sigde,]

o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
output[1:25,c(1,4,5)]
```

### 6. Analyze Significant Transcripts

>a) Without exiting the R environment, if the commands have been input manually, the following script can be executed in order to extract the fasta sequences for the transcripts that have been identified as the most differentially expressed.  If the R analysis part was executed as part of a script, then complete parts A-C manually and use the gene IDs from the significant results output file "temp30_significant.csv" to complete the following. The first step outputs the name of the significant transcript, just change the gene ID.  The second step generates a bed file which allows Bedtools to extract the relevant transcript from the reference fasta file.  

```
ballgown::transcriptNames(bg_genome)[GENEID]
grep -a "MSTRG.347.810" stringtie_merged.gtf | head -n 1 > MSTRG.347.810.bed
bedtools getfasta -fi S288C.fa -bed MSTRG.347.810.bed -fo MSTRG.347.810.bed.fa.out
```

>b) BLAST analysis: Take the sequence obtained in part A and use [UniProt](https://www.uniprot.org/blast/) to blast the transcript against the UniProtKb database.  The results will show which genes in the database the transcript shows similarities to, thus showing genes that have been identified to be differentially expressed by the pipeline.

![](https://github.com/aekimura/CS199_Final/blob/master/images/blast.png)

# Conclusion

The first group of samples, those grown at 30°C , contained two replicates each of the wild-type and the mutant Isw2 strains which are briefly described in Table 2.  The Isw2 mutant has an extra itc1::NatMx allele. The gene itc1 is a subunit of ATP-dependent Isw2p-itc1p chromatin remodeling complex.  It is required for the repression of a-specific genes, of early meiotic genes during mitotic growth, and of the gene INO1.  Thus, we expect that upon comparison, because of the extra itc1 copy present in the Isw2 mutant, genes relating to the Isw2p-itc1p chromatin remodeling complex to be higher in the Isw2 mutant while the genes that are repressed by the complex, meiotic genes or INO1, may be more highly expressed in the wild-type strains.  As shown in Table 8, the transcripts that were identified as the most significantly differentially expressed between the two strains include the chromatin-remodeling ATPase INO80 (MSTRG.347.1) and other chromatin structure remodeling complexes (MSTRG.369.3).  Thus, is seems as though analysis of these Illumina fastq samples with the HISAT-StringTie-Ballgown pipeline was able to with some accuracy identify transcripts relating to genes that we expect to be differentially expressed between wild-type and Isw2 mutants.  

The second group of samples, those grown at 37°C,  contained two replicates each of the mutant Rsc and mutant Ino80 strains which are described in Table 2.  The chromatin structure remodeling (Rsc) mutant and the mutant Ino80 differ by their ATP-dependent chromatin remodellers: Rsc and Ino80.  The Rsc complex belongs the the ATP-dependent chromatin remodeller family and  functions in transcriptional regulation and elongation, chromosome stability, establishing sister chromatid cohesion, and telomere maintenance. Ino80 is an ATPase and nucleosome spacing factor and a subunit of actin and actin related proteins that promotes nucleosome shift in the 3’ direction.  We expected to see differential expression between these two samples in proteins relating to chromatin remodeling or something similar.  As shown by the blast analysis of the transcripts identified as the most significantly differentially expressed between the two strains in Table 9, the identified transcripts aligned with a large variation of genes.  However, the genes that were identified do not seem to have much of a correlation to chromatin remodeling as expected.  This discrepancy could be due to data quality of the original RNA-seq data provided as analyzed by Fastqc (analysis results available on Github).  Another possible reason for these results could be that there are additional confounding factors affecting gene expression in the two mutants being compared. 

Overall, we were able to successfully implement the differential expression analysis pipeline as outlined by Pertea, M. et al. in order to compare gene expression across samples of different conditions using RNA-seq data.  Future implementation of the pipeline would include sample data of higher quality and perhaps more downstream comparison of different strains to the wild-type data using population as the differential variable instead of growth temperature to determine whether growth temperature or genetic background has a greater expression of gene expression in S. cerevisiae S288C.


[**Full Final Report**](https://docs.google.com/document/d/1OVK1lC2Tv07apcZXxRsIEHGQw2ZwCAVIk3lZTvoO_bk)
