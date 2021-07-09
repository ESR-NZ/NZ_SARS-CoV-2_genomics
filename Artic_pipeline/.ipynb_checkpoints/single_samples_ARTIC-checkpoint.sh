#!/bin/bash

#This script is to do a quick and dirty ARTIC v3 assemlby of a single sample. It needs a directory of basecalled fastq reads.   
# It (might) uses ARTIC 1.2.1 release. This is a bit ahead of the version we have been using, some outputs may look different.

## Usage:
# /NGS/active/VIR/SARS-CoV2/NZ_SARS-CoV-2_genomics/Artic_pipeline/single_samples_ARTIC.sh.sh </path/to/fastq_directory>


# Directory of fastqs
FASTQ_PASS=$(realpath $1)

echo $FASTQ_PASS
## Trim trailing / from input path
FASTQ_PASS=$(echo $FASTQ_PASS | sed 's:/*$::')

# boot the conda env
ACTIVATE=/opt/bioinf/anaconda3/anaconda3-5.0.0.1/bin/activate
#ARTIC_MEDAKA=/opt/bioinf/anaconda3/anaconda3-5.0.0.1/envs/artic-fieldbioinformatics_1.2.1_test
ARTIC_MEDAKA=/opt/bioinf/anaconda3/anaconda3-5.0.0.1/envs/artic-ncov2019-medaka_1.1.3_a

source $ACTIVATE $ARTIC_MEDAKA

#make a directory to run the analysis in
ANALYSIS_DIR=single_sample_analysis_$(basename $FASTQ_PASS)
echo 'Output directory is' $ANALYSIS_DIR
[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR 

cd $ANALYSIS_DIR
## Gather all the files together and put in the analysis directory 
artic guppyplex --directory $FASTQ_PASS --min-length 380 --max-length 550 --skip-quality-check\
 --output reads.fastq --prefix $(basename $FASTQ_PASS)



artic minion --medaka --medaka-model r941_min_high_g344 --normalise 1500 --threads 32\
 --scheme-directory /NGS/active/VIR/SARS-CoV2/NZ_SARS-CoV-2_genomics/Artic_pipeline/artic-ncov2019/primer_schemes\
 --read-file reads.fastq nCoV-2019/V3 $(basename $FASTQ_PASS)



# Below for use with artic 1.2.1...
## run the assembly, with strict flag.
#artic minion --medaka --medaka-model r941_min_high_g344 --strict --normalise 1500 --threads 32\
# --scheme-directory /NGS/active/VIR/SARS-CoV2/NZ_SARS-CoV-2_genomics/Artic_pipeline/artic-ncov2019/primer_schemes\
# --read-file reads.fastq nCoV-2019/V3 $(basename $FASTQ_PASS)

## run report 
#multiqc .


