#!/bin/bash

##Set variables
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config_file)
    CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    *)    
    POSITIONAL+=("$1")
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" 

RUN_ID=$1


source $CONFIG
source $ACTIVATE $ARTIC_MEDAKA

echo "Config file = ${CONFIG}"
echo "DATA_DIR = ${DATA_DIR}"
echo "ASSEMBLIES_DIR = ${ASSEMBLIES_DIR}"


[ ! -z $CONFIG] && echo "Please suppply config file path" && exit 1
[ ! -z $RUN_ID] && echo "Please suppply run-id" && exit 1
## exit if supplied directory is invalid
[ ! -d $DATA_DIR/$RUN_ID ] && echo "$RUN_ID invalid data directory" && exit 1
## Check if the data has already been processed, caution overwrite
[ -d $ASSEMBLIES_DIR/${RUN_ID}_analysis ] && read -p "Warning! Analysis already exists for $RUN_ID. Press ENTER to continue and overwrite!"

## Check for a raw data dir and then make run analysis directory tree
[ -d $DATA_DIR/$RUN_ID ] && mkdir -p $ASSEMBLIES_DIR/${RUN_ID}_analysis/${RUN_ID}_{called,barcodes,assemblies} && echo "Running $RUN_ID"
echo ""
#Set variable to point to base directory for this runs analysis output 
ANALYSIS_DIR=$ASSEMBLIES_DIR/${RUN_ID}_analysis


## Call the basecalling SLURM script (on GPU server)
echo "Running guppy gpu basecalling for $RUN_ID"
sbatch --export=ALL,RUN_ID=$RUN_ID,DATA_DIR=$DATA_DIR,ASSEMBLIES_DIR=$ASSEMBLIES_DIR --wait basecall_gpu.sh 
echo""

## Call the deplexing SLURM script (CPUs on Production servers)
echo "Running barcode demultiplexing for $RUN_ID"
DEPLEXING_PID=$(sbatch --export=ALL,RUN_ID=$RUN_ID,ANALYSIS_DIR=$ANALYSIS_DIR --wait --parsable deplex_prod.sh)

## Count the number of barcoded DIRs made by the deplexing to pass to next script
BC_LIST="$ANALYSIS_DIR/bc.list"
ls -1 $ANALYSIS_DIR/${RUN_ID}_barcodes | awk '/barcode/{print $0}' > $BC_LIST

## This is passed to sbatch array to tell it what barcodes were used  
BC_ARRAY=$(grep -o '[0-9][0-9]' $BC_LIST | tr -s '\n' ',')

NUM_BC=$BC_LIST | wc -l

echo "Number of barcodes is $NUM_BC"
echo ""

## Call array job to gather and length filter reads in each barcode dir
echo "Gathering reads for assembly for $RUN_ID"
sbatch --export=ALL,RUN_ID=$RUN_ID,ANALYSIS_DIR=$ANALYSIS_DIR,ARTIC_DIR=$ARTIC_DIR --wait --array=$BC_ARRAY gather_assemble.sh


## Do some clean up and reporting on the results.

# Copy the consensus files to a dir
CONSENSUN_DIR=$ANALYSIS_DIR/consensus
mkdir $CONSENSUN_DIR
cp $ANALYSIS_DIR/${RUN_ID}_assemblies/*.consensus.fasta $CONSENSUN_DIR



## Run the reporting script
#report.py --consensus_dir $CONSENSUN_DIR --bucket $AWS_BUCKET


echo ""
echo "Report for $RUN_ID"
cat $CONSENSUN_DIR/consensus_report.txt

