LIB=/NGS/active/VIR/SARS-CoV2/analysis/GridIon_assemblies/20201125_analysis/20201125/20201125
ARTIC_DIR=/NGS/active/VIR/SARS-CoV2/NZ_SARS-CoV-2_genomics/Artic_pipeline/artic-ncov2019/
cd ${LIB}_assemblies

#export BARCODE=barcode11
BAR=("01" "02" "03" "04" "05" "11" "12" "13" "14" "15" "16" "17" "25" "26" "27" "28" "29" "37" "38" "39" "40" "41" "49" "50" "51" "52" "61" "62" "64" "73" "74" "75" "76" "85" "86" "87" "88")
echo ${BAR[2]}
for BARCODE in ${BAR[@]}
do
    echo $BARCODE
    artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${LIB}_barcodes/barcode${BARCODE} --prefix $(basename $LIB)


    artic minion --medaka --normalise 200 --threads 12 --scheme-directory $ARTIC_DIR/primer_schemes --read-file $(basename $LIB)_barcode${BARCODE}.fastq nCoV-2019/V3 $(basename $LIB)_barcode${BARCODE}
done
