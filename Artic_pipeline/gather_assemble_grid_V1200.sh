#!/bin/bash

#SBATCH -p prod
#SBATCH --cpus-per-task=24
#SBATCH -o artic_assemble_grid%a.out
#SBATCH -e artic_assemble_grid%a.err


## Jump into the assembly directroy 

cd ${LIB}_assemblies

N=${SLURM_ARRAY_TASK_ID}
BARCODE=barcode$(printf "%02d" ${N})

echo "BARCODE = $BARCODE"
echo ""

## filter by length
srun artic guppyplex --min-length 400 --max-length 700\
 --directory ${LIB}_barcodes/$BARCODE --prefix $(basename $LIB) 


srun artic minion --medaka\
 --normalise 200 --threads 24 --scheme-directory $ARTIC_DIR/primer_schemes\
 --read-file $(basename $LIB)_${BARCODE}.fastq nCoV-2019/V1200 $(basename $LIB)_${BARCODE}



