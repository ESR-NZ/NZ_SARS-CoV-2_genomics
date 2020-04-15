#!/bin/bash

##Set variables, hardcoded to work on production server with the strict directory structure.
RUN_ID=$1
VX=V1
DATA_DIR=/NGS/active/VIR/SARS-CoV2/run_links #Where the Nanopore raw data will be stored
ASSEMBLIES_DIR=/NGS/active/VIR/SARS-CoV2/analysis #Where the analysed data will live 
ARTIC_DIR=/NGS/active/VIR/SARS-CoV2/Artic_pipeline/artic-ncov2019 #Where the ref genome and primer scheme live

## check for the existance of command line arg
[ -z "$1" ] && echo "No data directory supplied" && exit
## exit if supplied directory is invalid
[ ! -d $DATA_DIR/$RUN_ID ] && echo "$RUN_ID invalid data directory" && exit
## Check if the data has already been processed, caution overwrite
#[ -d $ASSEMBLIES_DIR/${RUN_ID}_analysis ] && read -p "Warning! Analysis already exists for $RUN_ID. Press ENTER to continue and overwrite!"

## Check for a raw data dir and then make run analysis directory tree
#[ -d $DATA_DIR/$RUN_ID ] && mkdir -p $ASSEMBLIES_DIR/${RUN_ID}_analysis/${RUN_ID}_{called,barcodes,assemblies} && echo "Running $RUN_ID"
echo ""
#Set variable to point to base directory for this runs analysis output 
ANALYSIS_DIR=$ASSEMBLIES_DIR/${RUN_ID}_analysis

## boot the conda env
source /opt/bioinf/anaconda3/anaconda3-5.0.0.1/bin/activate /opt/bioinf/anaconda3/anaconda3-5.0.0.1/envs/artic-ncov2019-medaka


## Call the basecalling SLURM script (on GPU server)
echo "Running guppy gpu basecalling for $RUN_ID"
#sbatch --export=ALL,RUN_ID=$RUN_ID,DATA_DIR=$DATA_DIR,ASSEMBLIES_DIR=$ASSEMBLIES_DIR --wait basecall_gpu.sh 
echo""

## Call the deplexing SLURM script (CPUs on Production servers)
echo "Running barcode demultiplexing for $RUN_ID"
#DEPLEXING_PID=$(sbatch --export=ALL,RUN_ID=$RUN_ID,ANALYSIS_DIR=$ANALYSIS_DIR --wait deplex_prod.sh)

## Count the number of barcoded DIRs made by the deplexing 
BC_LIST="$ANALYSIS_DIR/bc.list"
ls -1 $ANALYSIS_DIR/${RUN_ID}_barcodes | awk '/barcode/{print $0}' > $BC_LIST
NUM_BC=$BC_LIST | wc -l
## This is a string of the two digit barcode numbers passed to sbatch array to tell it what barcodes were used  
BC_ARRAY=$(grep -o '[0-9][0-9]' $BC_LIST | tr -s '\n' ',')


## Call array job to gather and length filter reads in each barcode dir
echo "Gathering reads for assembly for $RUN_ID"
sbatch --export=ALL,RUN_ID=$RUN_ID,ANALYSIS_DIR=$ANALYSIS_DIR,ARTIC_DIR=$ARTIC_DIR,VX=$VX --wait --array=$BC_ARRAY gather_assemble_VX.sh



## Do some clean up and reporting on the results.

# Copy the consensus files to a dir
CONSENSUN_DIR=$ANALYSIS_DIR/consensus
mkdir $CONSENSUN_DIR
cp $ANALYSIS_DIR/${RUN_ID}_assemblies/*.consensus.fasta $CONSENSUN_DIR

## Run the reporting script
report.py --consensus_dir $CONSENSUN_DIR


echo ""
echo "Report for $RUN_ID"
cat $CONSENSUN_DIR/consensus_report.txt

