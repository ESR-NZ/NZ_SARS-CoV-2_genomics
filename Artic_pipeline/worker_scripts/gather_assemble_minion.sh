#!/bin/bash

#SBATCH -p covid_urgent
#SBATCH --cpus-per-task=24
#SBATCH -o artic_assemble%a.out
#SBATCH -e artic_assemble%a.err


## Jump into the assembly directroy 
cd $ANALYSIS_DIR/${RUN_ID}_assemblies

N=${SLURM_ARRAY_TASK_ID}
BARCODE=barcode$(printf "%02d" ${N})

## filter by length
srun artic guppyplex --min-length 400 --max-length 700\
 --directory $ANALYSIS_DIR/${RUN_ID}_barcodes/$BARCODE\
 --prefix ${RUN_ID}

wait

srun artic minion --medaka\
 --normalise 1000 --threads 24 --scheme-directory $ARTIC_DIR/primer_schemes\
 --read-file ${RUN_ID}_${BARCODE}.fastq nCoV-2019/$PRIMER_SET ${RUN_ID}_${BARCODE}



