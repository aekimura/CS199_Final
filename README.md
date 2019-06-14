# Background

# Pipeline

[**1. Environment Setup**](#1-environment-setup)

[**2. Data Download**](#2-data-download)

[**3. Preprocessing**](#3-preprocessing)

[**4. Map Sample Reads, Assemble Transcripts, and Generate Ballgown Tables**](#4-map-reads)

[**5. R Analysis**](#5-r-analysis)

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

>b) Clean Reference Genome Fasta: The script "clean_reference.sh" removes unwanted characters from the fasta headers of the reference genome and replaces them with numbers.

```
cat S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]/I/g' | awk '/^>/{print ">" ++i; next}{print}' > S288C_reference_sequence_R64-2-1_20150113_new.fsa
```

>c) Index Building: Hisat2 v. 2.1.0 builds an index for the reference genome which is utilized later. 

```
source ~/.miniconda3rc
conda activate hisat2 

hisat2-build -f S288C_reference_sequence_R64-2-1_20150113_new.fsa S288C

conda deactivate
```

>d) GFF to GTF Conversion: Gffread v. 0.11.4 was used to convert the reference gff file to a gtf file for use later.

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

### 5. R Analysis

>The following parts A-Z are all included in the script "R_temp30.sh" and can be executed by using the script or by inputting the code into the proper R environment individially.  The output

>a) Load Required R Packages: The packages that are required for analysis that have previously been downloaded into R should be loaded in R so that they can be used.  

```
R
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
```

>b) Read In Expression Data:  This step takes the phenotype data CSV file that was manually created in the preprocessing stage and creates a variable called pheno_data

```
```

>c) Filtering:

```
```

>d) Transcript Identification:

```
```

>e) Gene Identification:

```
```

>f) 

```
```

# Conclusion

[Full Final Report](https://docs.google.com/document/d/1OVK1lC2Tv07apcZXxRsIEHGQw2ZwCAVIk3lZTvoO_bk)
