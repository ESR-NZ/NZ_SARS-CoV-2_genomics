# Aritc Assembly Pipeline Docs

## Run_scripts:

These are the scripts that are run by the user. They find the raw data, set up the analysis directories, then they call some sbatch scripts to analyses the data. The scripts are all very similar. The script to call depends on what kind of sequencing library was run, on what type of machine (minion or gridion). The sections below try to explain what one to use. 
The run scripts use a config file to store paths to data and environments ect so these are not exposed to public if we publish the scripts.

The general way to call any run script is:

`./run_script.sh -c ./Artic_config_file -l <str> $RUN_ID`

Where run $RUN_ID is what the run was called on the machine. You would know this from when you rsync'ed the data to production from the sequencer. The `-l` flag is for the location of the run i.e., `K` = KSC, `M` = MASC or `C` = CSC. 

### gridion_artic_96_run.sh
    - This is for the assembly of basecalled ARTIC-v3 data from the gridIon (at KSC) for libraries that have used the ligation sequencing kit with any barcode kit 1-24 and 1-96.
    - Barcoding/demultiplexing should be turned off on the GridIon as it needs to be done on production for ARTIC protocol to work properly
    - Currently active 
    - I will add extra comments to this script so it is clear what some of the weird Bash stuff does.
    - Default location is KSC so remember to add the `-l M` for runs done in MASC 

### gridion_artic_run_V1200_rapid.sh
    - This is for assembly of basecalled and barcoded (demultiplexed) ARTIC-v1200 data produced on the GridIon using the rapid library prep kit. 
    - Basecalling and barcoding must be ON for this to work. 
    - Used for the urgent runs 
    - Currently active
    - This also used the location flag

### minion_artic_run.sh
    - This is for assembly of non-basecalled ARTIC-v3 data produced from the minIon (not gridIon) with the ligation sequencing kit and 1-24 barcode kit. 
    - Uses the GPU server Orac with slurm for basecalling and assumes 1-24 barcoding.
    - This was the first run script, all other run scripts were based on this one. 
    - Not really used anymore since we got a GridIon but could be useful if CSC does COVID sequencing on a minion
    - Could/should upgrade to the 96 barcode kit 

### gridion_artic_run_V1200.sh
    - This is for the assembly of basecalled ARTIC-v1200 data from the gridIon (at KSC) for libraries that have used the ligation sequencing kit/ 1-24 barcode kit.
    - Barcoding/demultiplexing should be turned OFF during the GridIon run as it needs to be done on production for ARTIC protocol.
    - This is rarely used as the v1200 primers set is generally reserved for urgent runs (rapid kit). Sometimes when a genome has a large missing section we may try this protocol 
    - should work but not currently active. Should change the barcodes to 1-96
    - location flag not used in this script will just be default KSC directories.  

## worker_scripts:

These are the worker scripts (sbatch) called by the run scripts:
 
`basecall_gpu.sh`

`deplex_prod_grid_96.sh`

`deplex_prod_grid.sh`

`deplex_prod_grid_V1200.sh`

`deplex_prod_rapid.sh`

`deplex_prod.sh`

`gather_assemble_grid.sh`

`gather_assemble_grid_V1200.sh`

`gridion_artic_run_V1200_rapid.sh`

`gather_assemble.sh`

`report_grid.py report.py`



## run_script_dependencies:

Here are lists of the worker scripts called by each run scripts just for reference. Some of these scripts don't exist anymore.

#### `artic_run.sh`
    basecall_gpu.sh
    gather_assemble.sh
    report.py

#### `CSC_artic_run.sh` 
    basecall_gpu.sh
    gather_assemble.sh
    report.py

#### `gridion_artic_96_run.sh`
    deplex_prod_grid_96.sh
    gather_assemble_grid.sh
    report_grid.py

#### `gridion_artic_run.sh`
    deplex_prod_grid.sh
    gather_assemble_grid.sh
    report_grid.py

#### `gridion_artic_run_V1200.sh`
    deplex_prod_grid_V1200.sh
    gather_assemble_grid_V1200.sh

#### `gridion_artic_run_V1200_rapid.sh`
    gather_assemble_grid_V1200.sh

#### `MASC_gridion_artic_run.sh`
    deplex_prod_grid_96.sh
    gather_assemble_grid.sh
    report_grid.py
