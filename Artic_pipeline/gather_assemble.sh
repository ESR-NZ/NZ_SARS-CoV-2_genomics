#!/bin/bash

#SBATCH -p prod
#SBATCH --cpus-per-task=24
#SBATCH -o artic_assemble.out
#SBATCH -e artic_assemble.err

source /opt/bioinf/anaconda3/anaconda3-5.0.0.1/bin/activate /opt/bioinf/anaconda3/anaconda3-5.0.0.1/envs/artic-ncov2019-medaka

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
 --normalise 200 --threads 24 --scheme-directory $ARTIC_DIR/primer_schemes\
 --read-file ${RUN_ID}_${BARCODE}.fastq nCoV-2019/V3 ${RUN_ID}_${BARCODE}



conda deactivate
