#!/bin/bash

#SBATCH -p covid_urgent
#SBATCH --cpus-per-task=24
#SBATCH -o rapid_assemble_grid%a.out
#SBATCH -e rapid_assemble_grid%a.err


## work out what max size filter reads to. I think amplicon size plus 300bps for BCs and adptors 
if [ -z "$PRIMER_SET" ] ## if not set, default to V1200 settings
then
    MAX_LEN=1400
    
elif [ "$PRIMER_SET" = 'V3' ] ## set paths to KSC settings
then
    MAX_LEN=700

elif [ "$PRIMER_SET" = 'V2500' ] ## set paths to MASC
then
    MAX_LEN=2800

elif [ "$PRIMER_SET" = 'V1200' ] ## set paths to MASC
then
    MAX_LEN=1500
fi

## Jump into the assembly directroy 
cd ${LIB}_assemblies

N=${SLURM_ARRAY_TASK_ID}
BARCODE=barcode$(printf "%02d" ${N})

echo "BARCODE = $BARCODE"
echo ""

## filter by length
srun artic guppyplex --min-length 250 --max-length $MAX_LEN\
 --directory ${LIB}_basecalled_link/$BARCODE --prefix $(basename $LIB) 


srun artic minion --medaka\
 --normalise 1000 --threads 24 --scheme-directory $ARTIC_DIR/primer_schemes\
 --read-file $(basename $LIB)_${BARCODE}.fastq nCoV-2019/${PRIMER_SET} $(basename $LIB)_${BARCODE}



