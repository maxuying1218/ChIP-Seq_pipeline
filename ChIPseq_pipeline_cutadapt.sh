#!/usr/bin/env bash

#sh ChIPseq_pipeline.sh input_fq1,input_fq2 out_dir out_label ./genome/bowtie_indexes/mm10 10 0 0 1 genome_name gtf_path

fq=${1//,/ }
out_dir=$2
out_label=$3
bowtie_index=$4
cpu=$5
trim5=$6
trim3=$7
mismatch=$8
genome=$9
gtf=${10}

if [ ! -d ${out_dir} ]
then
    mkdir ${out_dir}
else
    echo ${out_dir}" exist!!!"
fi


out_dir=`cd $out_dir|pwd`/${out_dir}
cut_dir=${out_dir}/cutadapt
bowtie_dir=${out_dir}/bowtie
bed_rm_dir=${out_dir}/bed/rmdup
bw_dir=${out_dir}/bw
macs_dir=${out_dir}/macs
pkan_dir=${out_dir}/pkan

mkdir -p ${cut_dir}/logs ${bowtie_dir}/logs ${bed_norm_dir} ${bed_rm_dir}/logs ${bw_dir} ${pro_tss_dir} ${pro_gb_dir} ${macs_dir} ${pkan_dir}/logs

###cutadapt
#adapt="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapt="GATCGGAAGAGCACA"
echo `date +%F": "%X`,": Cutadapt ${out_label}!"
cut_fq=${cut_dir}/${out_label}.fq
if [[ $fq =~ ".fq.gz" ]]||[[ $fq =~ ".fastq.gz" ]]
then
    gzip -cd $fq| cutadapt -a $adapt -f fastq -m 20 -o ${cut_fq} - > ${cut_dir}/logs/${out_label}.txt
else
    cat $fq| cutadapt -a $adapt -f fastq -m 20 -o ${cut_fq} - > ${cut_dir}/logs/${out_label}.txt
fi


### bowtie_mapping
echo `date +%F": "%X`,": Bowtie mapping ${out_label}!"
bowtie -S $bowtie_index -p $cpu -m 1 -n ${mismatch} -q ${cut_fq} ${bowtie_dir}/${out_label}.sam 2> ${bowtie_dir}/logs/${out_label}.mapping.log
rm ${cut_fq}

echo `date +%F": "%X`,": Sam to Bam ${out_label}!"
sam=${bowtie_dir}/${out_label}.sam 
cpu=$cpu
flag=${sam%.sam}
bam=${flag}.bam
sorted_bam=${flag}.sorted.bam
samtools view -bS -@ $cpu  $sam > $bam
samtools sort -T $flag -@ $cpu -o $sorted_bam $bam
samtools index $sorted_bam
rm $bam
rm ${bowtie_dir}/${out_label}.sam
###

echo `date +%F": "%X`,": Remove duplicates ${out_label}!"
picard MarkDuplicates REMOVE_DUPLICATES=true I=${bowtie_dir}/${out_label}.sorted.bam O=${bowtie_dir}/${out_label}.rmdup.bam > ${bowtie_dir}/logs/${out_label}.rmdup.txt

echo `date +%F": "%X`,": Bam to Bed ${out_label}!"
bamToBed -i ${bowtie_dir}/${out_label}.rmdup.bam|sort -S 5G --parallel=10 -k1,1 -k2,2n > ${bed_norm_dir}/${out_label}.bed

echo `date +%F": "%X`,": Bed to bw ${out_label}!"
cd ${bw_dir}
bedToBigWig -b ${bed_rm_dir}/${out_label}.bed -o ${out_label}.bw -g $genome -e 200

#echo `date +%F": "%X`,": Make geneProfile ${out_label}!"
#cd ${pro_tss_dir}
#geneProfile -rd ../../bw/${out_label}.dup1.200bp.bdg -o ${out_label}.profile -g ${gtf} -n $genome -tss -tu 5000 -td 5000
#cd ${pro_gb_dir}
#geneProfile -rd ../../bw/${out_label}.dup1.200bp.bdg -o ${out_label}.profile -g ${gtf} -n $genome 

echo `date +%F": "%X`,": Macs Call Peak ${out_label}!"
cd ${macs_dir}
macs -t ${bed_rm_dir}/${out_label}.bed --name ${out_label} --tsize 150 --nomodel --nolambda
rm ${out_label}_peaks.xls
sed -i '1d' ${out_label}_peaks.bed
sed -i 's/\t-[0-9]*\t/\t0\t/g' ${out_label}_peaks.bed

echo `date +%F": "%X`,": PeakAnnotate ${out_label}!"
cd ${pkan_dir}
annotatePeaks.pl ${macs_dir}/${out_label}_peaks.bed $gn 1>${out_label}.peak_anno.xlsx 2>logs/${out_label}.pkan.log
cd ${out_dir}

echo `date +%F": "%X`,": Finished ${out_label}!"



