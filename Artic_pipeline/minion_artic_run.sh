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
    -l|--location) # can be 'K' for KSC or 'M' for MASC or 'C' for CSC
    LOC="$2"
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

## Default to config in run dir to save typing each time
[ -z "$CONFIG" ] && CONFIG=./Artic_pipeline_config
## bring in the config variables and activate the conda env
source $CONFIG
source $ACTIVATE $ARTIC_MEDAKA
## append DATA_DIR and ASSEMBLIES_DIR with gridion paths


# Location specific paths
if [ -z "$LOC" ] ## if not set, default to KSC settings
then
    LOC_DATA='Minion_runs'
    LOC_ASSEMBLY='nanopore'

elif [ "$LOC" = 'K' ] ## set paths to KSC settings
then
    LOC_DATA='Minion_runs'
    LOC_ASSEMBLY='nanopore'

elif [ "$LOC" = 'M' ] ## set paths to MASC
then
    LOC_DATA='MASC'
    LOC_ASSEMBLY='MASC'

elif [ "$LOC" = 'C' ] ## set paths to MASC
then
    LOC_DATA='CSC'
    LOC_ASSEMBLY='CSC'

fi


#Test for in invalid location parameters
valid='KMC'
[[ "$LOC" =~ [^$valid] ]] && [ ! -z "$LOC" ] && echo "Invalid location supplied" && exit 1


## append DATA_DIR and ASSEMBLIES_DIR with location paths for data in and results out to the right place
# Gridion paths
DATA_DIR=${DATA_DIR}/${LOC_DATA}
ASSEMBLIES_DIR=${ASSEMBLIES_DIR}/${LOC_ASSEMBLY}

#Sanity check
echo "Config file = ${CONFIG}"
echo "DATA_DIR = ${DATA_DIR}/${RUN_ID}"
echo "ASSEMBLIES_DIR = ${ASSEMBLIES_DIR}"
echo "artic_dir = $ARTIC_DIR"
echo "Run_id  = $RUN_ID"
echo ""

##Check for run id 
[ -z "$RUN_ID" ] && echo "Please suppply run-id" && exit 1
## exit if supplied directory is invalid
[ ! -d $DATA_DIR/$RUN_ID ] && echo "$RUN_ID invalid data directory" && exit 1
## Check if the data has already been processed, caution overwrite
[ -d $ASSEMBLIES_DIR/${RUN_ID}_analysis ] && read -p "Warning! Analysis already exists for $RUN_ID. Press ENTER to continue and overwrite!"
## Check for a raw data dir and then make run analysis directory tree
[ -d $DATA_DIR/$RUN_ID ] && mkdir -p $ASSEMBLIES_DIR/${RUN_ID}_analysis/${RUN_ID}_{called,barcodes,assemblies} && echo "Running $RUN_ID"
echo ""
#Set variable to point to base directory for this runs analysis output 
ANALYSIS_DIR=$ASSEMBLIES_DIR/${RUN_ID}_analysis

##Script base dir for finding the resources for the pipeline
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR=${DIR}/worker_scripts
#Sanity check

echo "Config file = ${CONFIG}"
echo "DATA_DIR = ${DATA_DIR}"
echo "ASSEMBLIES_DIR = ${ASSEMBLIES_DIR}"
echo "artic_dir = $ARTIC_DIR"
echo "Run_id  = $RUN_ID"


## Call the basecalling SLURM script (on GPU server)
echo "Running guppy gpu basecalling for $RUN_ID"
sbatch --export=ALL,RUN_ID=$RUN_ID,DATA_DIR=$DATA_DIR,ASSEMBLIES_DIR=$ASSEMBLIES_DIR --wait $DIR/basecall_gpu.sh 

## Call the deplexing SLURM script (CPUs on Production servers)
echo "Running barcode demultiplexing for $RUN_ID"
DEPLEXING_PID=$(sbatch --export=ALL,RUN_ID=$RUN_ID,ANALYSIS_DIR=$ANALYSIS_DIR --wait $DIR/deplex_prod.sh)

## Count the number of barcoded DIRs made by the deplexing to pass to next script
BC_LIST="$ANALYSIS_DIR/bc.list"
ls -1 $ANALYSIS_DIR/${RUN_ID}_barcodes | awk '/barcode/{print $0}' > $BC_LIST

## This is passed to sbatch array to tell it what barcodes were used  
BC_ARRAY=$(grep -o '[0-9][0-9]' $BC_LIST | tr -s '\n' ',')

NUM_BC=$(cat $BC_LIST) 

echo "Barcodes used: $NUM_BC"
echo ""

## Call array job to gather, length filter and assemble reads in each barcode dir
echo "Gathering reads for assembly for $RUN_ID"
sbatch --export=ALL,RUN_ID=$RUN_ID,ANALYSIS_DIR=$ANALYSIS_DIR,ARTIC_DIR=$ARTIC_DIR --wait --array=$BC_ARRAY $DIR/gather_assemble.sh

## Do some clean up and reporting on the results.

# Copy the consensus files to a dir
CONSENSUN_DIR=$ANALYSIS_DIR/consensus
mkdir $CONSENSUN_DIR
cp $ANALYSIS_DIR/${RUN_ID}_assemblies/*.consensus.fasta $CONSENSUN_DIR

## Run the reporting script
$DIR/report.py --consensus_dir $CONSENSUN_DIR #--bucket $AWS_BUCKET

echo ""
echo "Report for $RUN_ID"
cat $CONSENSUN_DIR/consensus_report.txt

conda deactivate

## Find vcf files 
