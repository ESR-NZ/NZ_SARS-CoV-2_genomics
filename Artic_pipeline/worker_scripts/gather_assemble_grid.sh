#!/bin/bash

#SBATCH -p covid_urgent
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
artic guppyplex --skip-quality-check --min-length 400 --max-length 700\
 --directory ${LIB}_barcodes/$BARCODE --prefix $(basename $LIB) 


artic minion --medaka --medaka-model r941_min_high_g360\
 --normalise 1000 --threads $SLURM_CPUS_PER_TASK --scheme-directory $ARTIC_DIR/primer_schemes\
 --read-file $(basename $LIB)_${BARCODE}.fastq nCoV-2019/V3 $(basename $LIB)_${BARCODE}



