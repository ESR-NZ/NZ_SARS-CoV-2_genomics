# rsync_examples

## General command structure:
`rsync -zarv --include="*/" --include="*.sh" --exclude="*" "$from" "$to"`

## KSC gridIon
export the run_ID, at KSC we use the following convention

`export RUN_ID=yyyymmdd` 

`rsync -zavmr --include="*/" --exclude="*.fast5" --exclude="*/fastq_fail/*" grid@10.1.30.16:/data/${RUN_ID} run_links/GridIon/`


## MASC gridIon:
for example, a BAU artic run at MASC the run_ID will be something like what is shown here: 

`export RUN_ID=20210203_v3ARTIC_MASC`

`rsync -zavmr --include="*/" --exclude="*.fast5" --exclude="*/fastq_fail/*" grid@10.2.31.23:/data/COVID19/ARTIC/${RUN_ID} run_links/MASC/ARTIC`


