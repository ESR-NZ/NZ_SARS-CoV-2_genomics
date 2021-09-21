#!/bin/bash 
#SBATCH -p covid_urgent
#SBATCH --cpus-per-task=64
#SBATCH -o artic_run_deplex.out
#SBATCH -e artic_run_deplex.err



#Set up environment
INPUT_DIR="$ANALYSIS_DIR/${RUN_ID}_called"
OUT_DIR="$ANALYSIS_DIR/${RUN_ID}_barcodes"

module load guppy/3.6.0


srun guppy_barcoder -t 64 --input_path $INPUT_DIR --save_path $OUT_DIR --barcode_kits SQK-RBK004


module unload guppy/3.6.0
