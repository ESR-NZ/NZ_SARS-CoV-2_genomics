#!/bin/bash
#SBATCH --partition gpu
#SBATCH --cpus-per-task=8 
#SBATCH -o artic_basecall.out
#SBATCH -e artic_basecall.err

## SLURM script to do the basecalling on GPUs

INPUT=$DATA_DIR/$RUN_ID
OUTPUT=$ASSEMBLIES_DIR/${RUN_ID}_analysis/${RUN_ID}_called
 
MODULEPATH=/usr/share/Modules/modulefiles:/etc/modulefiles:/opt/dsc/modulefiles
module load guppy-gpu/3.4.4

guppy_basecaller --disable_pings -c dna_r9.4.1_450bps_fast.cfg -i $INPUT -s $OUTPUT -x 'cuda:all' --recursive --num_callers 8 --gpu_runners_per_device 4 --chunks_per_runner 512 --chunk_size 1000

module unload guppy-gpu/3.4.4
