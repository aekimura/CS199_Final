# Background

# Pipeline

[**1. Environment Setup - Miniconda & .yml**](#1-environment-setup)

[**2. Data Download - SRA-Explorer**](#2-data-download)

# Pipeline Steps

#### 1. Environment Setup

>First Miniconda must be installed and then the provided .yml files can be used to create the environments to complete the rest of the pipeline steps.  To download and set-up Miniconda, enter the following commands in order as they are required. Statments in all caps denote actions to complete.
    
 ```ruby
 ### Miniconda3 Installation ###
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
 ### Environment Set-Up ###
 conda env create -f hisat2_env.yml
 conda env create -f trimmomatic_env.yml
 conda env create -f r_env.yml
 ```

#### 2. Data Download

>Illumina datasets were obtained by searching for relevent Illumina datasets for .

# Conclusion
