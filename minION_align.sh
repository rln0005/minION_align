#!/bin/bash

####ONLY NEEDED IF USING ASC
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load fastqc/0.10.1
#module load gnu_parallel/201612222
#module load gnu_parallel

#move into data directory
cd ../minION_align/data

#cut file names at underscores and keep unique names stored in list
ls | grep ".fastq" | cut -d "_" -f 1 | sort | uniq > list

#combine all fastq files for a single barcode
while read i
do
	cat "$i"*.fastq > "$i"_all.fastq
done<list

#move the combined fastq file to main directory and cd into that directory
mv *_all.fastq ../
cd ../

#use NanoFilt to filter reads (based on avg length=500, quality=10, remove first & last 50 bp reads)
NanoFilt -l 500 -q 12 --headcrop 100 --tailcrop 100 *_all.fastq > trim.fastq

#use minimap2 to map to reference genome
minimap2 -ax map-ont --MD *_ref.fasta trim.fastq > mapped.sam

#use samtools to convert sam to bam
samtools view -bS mapped.sam > mapped.bam
#use samtools to sort the bam file
samtools sort -o mapped.sorted.bam mapped.bam
#use samtools to create an index file
samtools index mapped.sorted.bam
######idx stats?? check samtools mapping stats
####samtools idxstats
####samtools flagstat
####samtools depth ###huge output
####write function to get depth over a window instead of per site

#use Sniffles to call variants
sniffles -m mapped.sorted.bam -v variants.vcf

#use IGV to visualize
igv
#genomes --> upload the reference fasta file
#file --> upload the vcf file (variants.vcf) AND the sorted bam file (mapped.sorted.bam)



