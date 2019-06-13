# Background

# Pipeline

[**1. Environment Setup**](#1-environment-setup)

[**2. Data Download**](#2-data-download)

[**3. Preprocessing**](#3-preprocessing)

# Pipeline Steps

### 1. Environment Setup

>First Miniconda must be installed and then the provided .yml files can be used to create the environments to complete the rest of the pipeline steps.  To download and set-up Miniconda, enter the following commands in order as they are required. Statments in all caps denote actions to complete.  [Online installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) can also be followed to set up Miniconda if desired.
    
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
 
 >After Miniconda set-up is completed do not exit conda and the three environments (Hisat2, Trimmomatic, and R) can be created using the .yml files
    
 ```ruby
 conda env create -f hisat2_env.yml
 conda env create -f trimmomatic_env.yml
 conda env create -f r_env.yml
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

>a) Clean Reference Genome Fasta: The script "clean_reference.sh" removes unwanted characters from the fasta headers of the reference genome and replaces them with numbers.

```
cat S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]/I/g' | awk '/^>/{print ">" ++i; next}{print}' > S288C_reference_sequence_R64-2-1_20150113_new.fsa
```

>b) Index Building: Hisat2 v. 2.1.0 builds an index for the reference genome which is utilized later. 

```
source ~/.miniconda3rc
conda activate hisat2 

hisat2-build -f S288C_reference_sequence_R64-2-1_20150113_new.fsa S288C

conda deactivate
```

>c) GFF to GTF Conversion: Gffread v. 0.11.4 was used to convert the reference gff file to a gtf file for use later.

```
source ~/.miniconda3rc
conda activate hisat2 

gffread saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o S288C.gtf

conda deactivate
```

>d) Fastq Trimming: Trimmomatic v. 0.39 was used to prepare the 8x Illumina paired-end reads obtained from the Sequence Read Archive for alignment.  Trimmomatic was used to remove adapters, to remove leading and trailing N bases, to scan reads in windows 4 bases long and remove any where the average base quality score is below 15, to remove reads shorter than 5 bases, and to remove reads with an average quality of less than 20.

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

>e) Fastq Quality Check: Fastqc version 0.11.8 was used to measure the phred quality scores across the bases of the fastq samples.  The script "fastqc.sh" measures the quality for all the sample fastq files downloaded. The results are included in the folder "fastqc_results".  

```
source ~/.miniconda3rc
conda activate hisat2 

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

fastqc ${SAMPLE}_1.fastq.gz
fastqc ${SAMPLE}_2.fastq.gz

conda deactivate
```

# Conclusion
