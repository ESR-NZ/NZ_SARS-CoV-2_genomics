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

## Default to config in run dir to save typing each time
[ -z "$CONFIG" ] && CONFIG=./artic_pipeline_config
## bring in the config variables and activate the conda env
source $CONFIG
source $ACTIVATE $ARTIC_MEDAKA
## append DATA_DIR and ASSEMBLIES_DIR with gridion paths

# Gridion paths
DATA_DIR=${DATA_DIR}/GridIon/
ASSEMBLIES_DIR=${ASSEMBLIES_DIR}/GridIon_assemblies

##Check for run id 
[ -z "$RUN_ID" ] && echo "Please suppply run-id" && exit 1
## exit if supplied directory is invalid
[ ! -d $DATA_DIR/$RUN_ID ] && echo "$RUN_ID invalid data directory" && exit 1


## Check if the data has already been processed, caution overwrite
[ -d $ASSEMBLIES_DIR/GridIon_assemblies/${RUN_ID}_analysis ] && read -p "Warning! Analysis directory already exists for $RUN_ID. Press ENTER to continue and overwrite!"

## Check for a raw data dir and then make run analysis directory tree
[ -d $DATA_DIR/$RUN_ID ] && mkdir -p $ASSEMBLIES_DIR/${RUN_ID}_analysis/

## Make dirctory tree for run
## Get the sub_runs for the gridion run
LIBRARIES=$(ls $DATA_DIR/$RUN_ID)

for library in ${LIBRARIES[@]}
do  
    mkdir -p $ASSEMBLIES_DIR/${RUN_ID}_analysis/${library}/${library}_{barcodes,assemblies,consensus} 
    ln -s $DATA_DIR/$RUN_ID/$library/*/fastq_pass $ASSEMBLIES_DIR/${RUN_ID}_analysis/${library}/${library}_basecalled_link 
done

#Set variable to point to base directory for this runs analysis output base dir
ANALYSIS_DIR=$ASSEMBLIES_DIR/${RUN_ID}_analysis

##Script base dir for finding the resources for the pipeline
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## Call the deplexing SLURM script (CPUs on Production servers)

## File of paths to each lib in run to used in sbatch jobs
LIB_PATHS="$ANALYSIS_DIR/libs.paths"

ls -1 $ANALYSIS_DIR | awk -v var="$ANALYSIS_DIR" '!/libs.paths/ {print var"/"$0"/"$0}' > $LIB_PATHS
readarray -t PATHS_ARRAY < $LIB_PATHS

NUM_LIBS_1=$((${#PATHS_ARRAY[@]}-1))

#echo "Running barcode demultiplexing for run $RUN_ID"
sbatch --export=ALL,ANALYSIS_DIR=$ANALYSIS_DIR,LIB_PATHS=$LIB_PATHS --wait --array=0-$NUM_LIBS_1 $DIR/deplex_prod_grid_V1200.sh  

# ## Count the number of barcoded DIRs made by the deplexing to pass to next script

for lib in ${PATHS_ARRAY[@]}; do
  
    BC_ARRAY=$(ls -1 ${lib}_barcodes | awk '/barcode/{print $0}' | grep -o '[0-9][0-9]'| tr -s '\n' ',')  # ## This is passed to sbatch array to tell it what barcodes were used 
    echo "Running assembly for: $lib"
    echo "With these barcodes: $BC_ARRAY"
    
    ## Call array job to gather, length filter and assemble reads in each barcode dir
    sbatch --export=ALL,ANALYSIS_DIR=$ANALYSIS_DIR,ARTIC_DIR=$ARTIC_DIR,LIB=$lib --wait --array=$BC_ARRAY $DIR/gather_assemble_grid_V1200.sh
    
    ##Clean up and run report script
    cp ${lib}_assemblies/*.consensus.fasta ${lib}_consensus
    $DIR/report_grid.py -c ${lib}_consensus

done

conda deactivate

