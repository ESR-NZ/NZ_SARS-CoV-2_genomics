#!/bin/bash 
#SBATCH -p prod
#SBATCH --cpus-per-task=32
#SBATCH -o deplex_grid.%a.out
#SBATCH -e deplex_grid.%a.err


## CALLEDPATHS exported form parent script to here 
readarray -t PATHS_ARRAY < $LIB_PATHS  ## read in file to array

LIB_PATH=${PATHS_ARRAY[$SLURM_ARRAY_TASK_ID]} ##index path array by the slurm job array number 




module load guppy/3.4.4
srun guppy_barcoder -t 32 --input_path ${LIB_PATH}_basecalled_link --save_path ${LIB_PATH}_barcodes --require_barcodes_both_ends --barcode_kits SQK-RBK004
srun guppy_barcoder -t 32 --input_path ${LIB_PATH}_basecalled_link --save_path ${LIB_PATH}_barcodes --require_barcodes_both_ends --barcode_kits SQK-RBK001
module unload guppy/3.4.4

