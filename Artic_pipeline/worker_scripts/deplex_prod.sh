#!/bin/bash 
#SBATCH -p covid_urgent
#SBATCH --cpus-per-task=64
#SBATCH -o artic_run_deplex.out
#SBATCH -e artic_run_deplex.err



#Set up environment
INPUT_DIR="$ANALYSIS_DIR/${RUN_ID}_called"
OUT_DIR="$ANALYSIS_DIR/${RUN_ID}_barcodes"

module load guppy/4.0.15

srun guppy_barcoder -t 64 --input_path $INPUT_DIR --save_path $OUT_DIR --require_barcodes_both_ends --barcode_kits EXP-NBD196

module unload guppy/4.0.15

