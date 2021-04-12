#!/bin/bash

################################# minION_align.sh created by Rebecca Nance ########################################


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

#use NanoFilt to filter reads (based on avg length=500, quality=12, remove first 50 bp reads)
NanoFilt -l 500 -q 12 --headcrop 50 < *_all.fastq > trim.fastq
#### fastqc analysis shows these reads have fairly lousy quality scores (this is as much trimming as can be done)


#use medaka to generate consensus seq and generate variants
medaka_haploid_variant -t 2 trim.fastq *_ref.fasta
#annotate medaka's vcf file so it can be uploaded into IGV
medaka tools annotate consensus_to_ref.vcf *_ref.fasta mapped.sorted.bam variants.annotated.vcf

######### IGV LOOKS CRAZY WITH THIS OUTPUT DATA......... ????????? ERROR CORRECTION? FILTERING? RACON?
#use IGV to visualize
#igv
#genomes --> upload the reference fasta file
#file --> upload the vcf file (variants.vcf) AND the sorted bam file (mapped.sorted.bam)



###### Alternative Map/Alignment using Minimap/Miniasm
#use minimap2 to align sequencing reads onto themselves
minimap2 -x ava-ont trim.fastq trim.fastq | gzip -1 > minimap.paf.gz
#use the above minimap output and the trimmed reads to assemble untigs using miniasm
miniasm -f trim.fastq minimap.paf.gz > miniasm.gfa

###### Alternative Map/Alignment using Shasta
#generates a folder called 'ShastaRun' that contains a file called Assembly.gfa
shasta-Linux-0.7.0 --input trim.fastq

#use Bandage to visualize different gfa assembly files
Bandage
#upload each of the gfa files into the Bandage gui



########## THIS HAS NO ERROR CORRECTION ##############
##use minimap2 to map to reference genome
#minimap2 -ax map-ont --MD *_ref.fasta trim.fastq > mapped.sam
##use samtools to convert sam to bam
#samtools view -bS mapped.sam > mapped.bam
#use samtools to sort the bam file
#samtools sort -o mapped.sorted.bam mapped.bam
##use samtools to create an index file
#samtools index mapped.sorted.bam





