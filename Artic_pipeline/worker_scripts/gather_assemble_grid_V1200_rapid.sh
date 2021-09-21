#!/bin/bash

#SBATCH -p covid_urgent
#SBATCH --cpus-per-task=24
#SBATCH -o rapid_assemble_grid%a.out
#SBATCH -e rapid_assemble_grid%a.err


## Jump into the assembly directroy 

cd ${LIB}_assemblies

N=${SLURM_ARRAY_TASK_ID}
BARCODE=barcode$(printf "%02d" ${N})

echo "BARCODE = $BARCODE"
echo ""

## filter by length
srun artic guppyplex --min-length 300 --max-length 1300\
 --directory ${LIB}_basecalled_link/$BARCODE --prefix $(basename $LIB) 


srun artic minion --medaka\
 --normalise 1000 --threads 24 --scheme-directory $ARTIC_DIR/primer_schemes\
 --read-file $(basename $LIB)_${BARCODE}.fastq nCoV-2019/V1200 $(basename $LIB)_${BARCODE}



