#!/bin/bash 
#SBATCH -p covid_urgent
#SBATCH --cpus-per-task=32
#SBATCH -o deplex_grid.%a.out
#SBATCH -e deplex_grid.%a.err


## CALLEDPATHS exported form parent script to here 
readarray -t PATHS_ARRAY < $LIB_PATHS  ## read in file to array

LIB_PATH=${PATHS_ARRAY[$SLURM_ARRAY_TASK_ID]} ##index path array by the slurm job array number 





module load guppy/3.6.0
echo "Running demultiplexing for EXP-NBD114"
guppy_barcoder -t $SLURM_CPUS_PER_TASK --input_path ${LIB_PATH}_basecalled_link --save_path ${LIB_PATH}_barcodes --require_barcodes_both_ends --barcode_kits EXP-NBD114
echo "Running demultiplexing for EXP-NBD104"
guppy_barcoder -t $SLURM_CPUS_PER_TASK --input_path ${LIB_PATH}_basecalled_link --save_path ${LIB_PATH}_barcodes --require_barcodes_both_ends --barcode_kits EXP-NBD104
module unload guppy/3.6.0

## changed to 3.6.0 from 3.4.4
## removed srun from before command
## should run both barcode sets at once, this didnt work for some reason and   
