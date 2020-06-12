#!/bin/bash

#options to define directories where the reads are
miseq_reads_directory=''
designation_directory=''
name=''

while getopts 'r:d:n:h' opt; do
  case $opt in
    r) miseq_reads_directory=$OPTARG ;;
    d) designation_directory=$OPTARG ;;
    n) name=$OPTARG ;;
    h) echo "Usage: bash new_softclip.sh -r <fullpath miseq_reads_directory> -d <fullpath designation_directory> -n <sequence ending ie _R1.fastq.gz etc> -h <help>"
        exit
        ;;
  esac
done


mkdir ${designation_directory}/trimmed
mkdir ${designation_directory}/bwa_align
mkdir ${designation_directory}/vcf
mkdir ${designation_directory}/consensus
mkdir ${designation_directory}/ivar

trimmed_directory=${designation_directory}/trimmed
bwa_directory=${designation_directory}/bwa_align
vcf_directory=${designation_directory}/vcf
consensus_directory=${designation_directory}/consensus
ivar=${designation_directory}/ivar
analysis_files=/NGS/active/VIR/SARS-CoV2/analysis/illumina/nsw/analysis_pipeline

##Define the files to be processed, for KSC_miseq reads no need to change 

files=${miseq_reads_directory}/*$name
for f in $files
do

##Define the files to be processed, for KSC_miseq reads no need to change 

if [[ "$name" == "1.fastq.gz" ]];
 then
     x=$(basename $f _$name)
     echo ${x}
     read1=${x}_$name
     read2=${x}_2.fastq.gz
     echo $read2
 else
    echo next
fi

if [[ "$name" == "R1.fastq.gz" ]];
 then
     x=$(basename $f _$name)
     echo ${x}
     read1=${x}_$name
     read2=${x}_R2.fastq.gz
     echo $read2
 else
   echo next
fi

if [[ "$name" == "R1.fastq" ]];
 then
     x=$(basename $f _$name)
     echo ${x}
     read1=${x}_$name
     read2=${x}_R2.fastq
     echo $read2
 else
   echo next
fi

if [[ "$name" == "R1_001.fastq.gz" ]];
 then
     x=$(basename $f _$name)
     echo ${x}
     read1=${x}_$name
     read2=${x}_R2_001.fastq.gz
     echo $read2
 else
   echo next
fi

#x=$(basename $f _$name)
#echo ${x}

#read1=${x}_$name
#read2=${x}_

#echo $read2

##quality trimming
trimmomatic PE -threads 30 -phred33 ${miseq_reads_directory}/$read1 \
${miseq_reads_directory}/$read2 \
${trimmed_directory}/${x}_R1_P.fastq ${trimmed_directory}/${x}_R1_uP.fastq \
${trimmed_directory}/${x}_R2_P.fastq ${trimmed_directory}/${x}_R2_uP.fastq  \
ILLUMINACLIP:/home/uren/bin/nullarbor_prokka/bin/../conf/trimmomatic.fa:1:30:11 \
HEADCROP:15 CROP:130 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:25 MINLEN:100

##overall coverage and mark bad regions
bwa mem -t 10 ${analysis_files}/nCoV-2019.reference.fasta ${trimmed_directory}/${x}_R1_P.fastq ${trimmed_directory}/${x}_R2_P.fastq | \
samtools sort -o ${bwa_directory}/${x}_all_sort.bam -
java -jar /opt/bioinf/picard/picard-2.10.10/bin/picard.jar MarkDuplicates I=${bwa_directory}/${x}_all_sort.bam O=${bwa_directory}/${x}_all_sortdup.bam M= ${bwa_directory}/${x}_all_dup_metics.txt
samtools index ${bwa_directory}/${x}_all_sortdup.bam
samtools depth -a ${bwa_directory}/${x}_all_sortdup.bam > ${bwa_directory}/${x}_depth.tsv
awk 'BEGIN {FS="\t"; OFS="\t"} $3 < 20 {print $1, $2-1, $2, "lowCoverage", "ncov"}' ${bwa_directory}/${x}_depth.tsv | cat - ${analysis_files}/CoV-2019_start_end.bed > ${bwa_directory}/${x}_badpostion.bed

##trim primers from alignment with ivar
module load iVar/1.2

ivar trim -i ${bwa_directory}/${x}_all_sortdup.bam -b ${analysis_files}/NSW_primer_A.bed -e -p ${x}_trimA
ivar trim -i ${bwa_directory}/${x}_all_sortdup.bam -b ${analysis_files}/NSW_primer_B.bed -e -p ${x}_trimB
mv ${x}_trimA.bam ${bwa_directory}
mv ${x}_trimB.bam ${bwa_directory}

samtools sort ${bwa_directory}/${x}_trimA.bam > ${bwa_directory}/${x}_trimA_sort.bam
samtools sort ${bwa_directory}/${x}_trimB.bam > ${bwa_directory}/${x}_trimB_sort.bam
samtools index ${bwa_directory}/${x}_trimA_sort.bam
samtools index ${bwa_directory}/${x}_trimB_sort.bam

##align to primer masked fastas##
#bwa mem -t 10 nCoV-2019_maskA.fasta ${trimmed_directory}/${x}_R1_P.fastq ${trimmed_directory}/${x}_R2_P.fastq | \
#samtools sort -o ${bwa_directory}/${x}_maskA_sort.bam -
#java -jar /opt/bioinf/picard/picard-2.10.10/bin/picard.jar MarkDuplicates I=${bwa_directory}/${x}_maskA_sort.bam O=${bwa_directory}/${x}_maskA_sortdup.bam M= ${bwa_directory}/${x}_maskA_dup_metics.txt
#samtools index ${bwa_directory}/${x}_maskA_sortdup.bam
#bwa mem -t 10 nCoV-2019_maskB.fasta ${trimmed_directory}/${x}_R1_P.fastq ${trimmed_directory}/${x}_R2_P.fastq | \
#samtools sort -o ${bwa_directory}/${x}_maskB_sort.bam -
#samtools index ${bwa_directory}/${x}_maskB_sort.bam
#java -jar /opt/bioinf/picard/picard-2.10.10/bin/picard.jar MarkDuplicates I=${bwa_directory}/${x}_maskB_sort.bam #O=${bwa_directory}/${x}_maskB_sortdup.bam M= ${bwa_directory}/${x}_maskB_dup_metics.txt
#samtools index ${bwa_directory}/${x}_maskB_sortdup.bam

##varaint calling##
## non primer regions ##

bcftools mpileup --max-depth 1000 -Ou -f ${analysis_files}/nCoV-2019.reference.fasta ${bwa_directory}/${x}_trimA_sort.bam  -R  ${analysis_files}/NSW_primer_A_comp.bed | bcftools call --ploidy 1 -mv > ${vcf_directory}/${x}_maskA.raw.vcf
bcftools mpileup --max-depth 1000 -Ou -f ${analysis_files}/nCoV-2019.reference.fasta ${bwa_directory}/${x}_trimB_sort.bam  -R  ${analysis_files}/NSW_primer_A_comp.bed | bcftools call --ploidy 1 -mv > ${vcf_directory}/${x}_maskB.raw.vcf
vcffilter -f "QUAL > 30" -f "DP > 19" ${vcf_directory}/${x}_maskA.raw.vcf > ${vcf_directory}/${x}_maskA_f.vcf
vcffilter -f "QUAL > 30" -f "DP > 19" ${vcf_directory}/${x}_maskB.raw.vcf > ${vcf_directory}/${x}_maskB_f.vcf

##primers regions, diploid at least not 1, also higher depth (?) ##
bcftools mpileup --max-depth 1000 -Ou -f ${analysis_files}/nCoV-2019.reference.fasta ${bwa_directory}/${x}_all_sortdup.bam  -R ${analysis_files}/NSW_primer.bed | bcftools call --ploidy 1 -mv > ${vcf_directory}/${x}_primer.raw.vcf
vcffilter -f "QUAL > 30" -f "DP > 38" ${vcf_directory}/${x}_primer.raw.vcf > ${vcf_directory}/${x}_primer.f.vcf

##combine all the vcfs
vcfcombine ${vcf_directory}/${x}_maskA_f.vcf ${vcf_directory}/${x}_maskB_f.vcf ${vcf_directory}/${x}_primer.f.vcf > ${vcf_directory}/${x}_comb_f.vcf
bgzip -c ${vcf_directory}/${x}_comb_f.vcf > ${vcf_directory}/${x}_comb_f.vcf.gz
tabix -p vcf ${vcf_directory}/${x}_comb_f.vcf.gz

##ivar analysis to detect minor alleles should could be interesting later. slightly more strigent than default, because fewer mutations? not sure##
samtools mpileup -aa -d 1000000 -Q 0 -A -f ${analysis_files}/nCoV-2019.reference.fasta ${bwa_directory}/${x}_all_sortdup.bam | \
ivar variants -p ${x}_ivar -q 30 -t 0.05 -m 20 -r ${analysis_files}/nCoV-2019.reference.fasta -g ${analysis_files}/annotation/GCF_009858895.2_ASM985889v3_genomic_mod.gff
mv ${x}_ivar.tsv ${ivar}
bedtools intersect -a ${ivar}/${x}_ivar.tsv -b ${analysis_files}/NSW_primer.bed  > ${ivar}/${x}_ivar_primer.tsv

##filter out alternative alleles as N. 
##talked to joep and matt decided:
## > 80% call alternative allele called by samtools  anyways
## between 20% to 79% call N sometimes called by samtools
## below 20% call reference allele called by ivar

grep -v '+[A,G,T,C]\|-[A,G,T,C]' ${ivar}/${x}_ivar.tsv | awk 'BEGIN {FS="\t";OFS="\t"} $11 < 0.8 && $11 > 0.19 {print $1, $2-1, $2, "IVAR", "freq:"$11}' >> ${bwa_directory}/${x}_badpostion.bed

##make consensus##
cat ${analysis_files}/nCoV-2019.reference.fasta | bcftools consensus  ${vcf_directory}/${x}_comb_f.vcf.gz > ${consensus_directory}/${x}_premask.fasta

##mask low depth regions and confusing alleles in consensus and mask dodgy sites
##http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473

bedtools maskfasta -fi ${consensus_directory}/${x}_premask.fasta -bed ${bwa_directory}/${x}_badpostion.bed -fo ${consensus_directory}/${x}_consensus.fasta
bedtools maskfasta -fi ${consensus_directory}/${x}_consensus.fasta -bed ${analysis_files}/dodgy_sites.bed -fo ${consensus_directory}/${x}_consensus_mask.fasta
sed -i "s/MN908947.3/${x}-NSW-tiling-illumina-consensus/" ${consensus_directory}/${x}_consensus.fasta
sed -i "s/MN908947.3/${x}-NSW-tiling-illumina-consensus_masked/" ${consensus_directory}/${x}_consensus_mask.fasta

##report
grep -v '>' ${consensus_directory}/${x}_consensus.fasta | tr '\n' ' ' | sed 's/ //g' |\
 sed 's/[A,G,C,T]N/\nN/g' | sed 's/N[A,G,C,T]/N\n/g' | grep 'N' > Ns.txt
awk '{print length($1)}' Ns.txt | paste  - Ns.txt > ${consensus_directory}/${x}_gaps.txt
awk '{sum+=$1} END {print FILENAME, sum}' ${consensus_directory}/${x}_gaps.txt >> ${consensus_directory}/gaps_summary.txt

done
