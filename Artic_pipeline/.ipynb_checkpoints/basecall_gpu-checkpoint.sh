#!/bin/bash
#SBATCH --partition gpu
#SBATCH --cpus-per-task=8 
#SBATCH --array=0,1
#SBATCH -o artic_basecall%a.out
#SBATCH -e artic_basecall%a.err

## SLURM script to do the basecalling on GPUs

INPUT=$DATA_DIR/$RUN_ID
OUTPUT=$ASSEMBLIES_DIR/${RUN_ID}_analysis/${RUN_ID}_called

## Make array of first and second halves of fast5 files 
FILES=(${INPUT}/fast5/*.fast5)  

## devide by two 
num_dirs=2
num_files=${#FILES[@]}
## interger division rounds down in bash
div_2=$(($num_files / $num_dirs))
#Half file arrays
FIRST_HALF=("${FILES[@]:0:$div_2}")
SECOND_HALF=("${FILES[@]:$div_2}")
## write file list
printf "%s\n" "${FIRST_HALF[@]}" > half_0.list
printf "%s\n" "${SECOND_HALF[@]}" > half_1.list

## load some stuff
MODULEPATH=/usr/share/Modules/modulefiles:/etc/modulefiles:/opt/dsc/modulefiles

module load guppy-gpu/3.4.4
#module load guppy-gpu/3.6.0

##Run the jobs
srun guppy_basecaller --disable_pings\
 -c dna_r9.4.1_450bps_fast.cfg\
 --input_file_list half_${SLURM_ARRAY_TASK_ID}.list\
 -s $OUTPUT/${SLURM_ARRAY_TASK_ID} --client_id ${SLURM_ARRAY_TASK_ID}\
 -x "cuda:${SLURM_ARRAY_TASK_ID}"\
 --num_callers 8\
 --gpu_runners_per_device 4\
 --chunks_per_runner 512\
 --chunk_size 1000

wait

module unload guppy-gpu/3.4.4
#module unload guppy-gpu/3.6.0

## Clean up list files 
rm half_{0,1}.list

mv $OUTPUT/{0,1}/* $OUTPUT
rm -r $OUTPUT/{0,1}
