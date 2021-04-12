#!/bin/sh

################################# minION_align.sh created by Rebecca Nance ########################################
## sequence data (fastq files) must be located in a /minION_align/data
## a reference sequence file (fasta) must include ref.fasta in the file name and be located in minION_align directory



#make output directory
mkdir -p output
#cd to output directory
cd output
#make other subdirectories
mkdir -p consensus_sequence_fasta
mkdir -p draft_assembly_gfa

#cd back into main directory
cd ../

#move into data directory where fastq seq files are contained
cd data

#cut file names at underscores and keep unique names stored in list
ls | grep ".fastq" | cut -d "_" -f 3 | sort | uniq > list

#combine all fastq files for a single barcode
while read i
do
cat *"$i"*.fastq > "$i"_all.fastq
done<list

mv *_all.fastq ../
cd ../

#use NanoFilt to filter reads (based on avg length=500, quality=12, remove first 50 bp reads)
cat *_all.fastq | NanoFilt -l 500 -q 12 --headcrop 50 > trim.fastq



### Generate first assembly (from reference sequence, using minimap2)
#generate draft assembly using reference seq with minimap2
minimap2 -x map-ont *_ref.fasta trim.fastq | gzip -1 > map_minimap_ref.paf.gz
#assemble untigs from the mapped reads
miniasm -f trim.fastq map_minimap_ref.paf.gz > map_minimap_ref.gfa
#generate consensus seq by converting gfa to fasta
awk '/^S/{print ">"$2"\n"$3}' map_minimap_ref.gfa > map_minimap_ref.fasta
#map trimmed reads onto miniasm assembly
minimap2 map_minimap_ref.fasta trim.fastq | gzip -1 > map_minimap_ref.racon.paf.gz

echo "FIRST ASSEMBLY COMPLETE"

### Generate second assembly (from de novo, using minimap2)
#map reads onto themsevles to identify overlaps (aka untigs)
minimap2 -x ava-ont trim.fastq trim.fastq | gzip -1 > overlaps.paf.gz
#assemble untigs from the mapped overlapping reads
miniasm -f trim.fastq overlaps.paf.gz > untigs.gfa
#generate consensus seq by converting gfa to fasta
awk '/^S/{print ">"$2"\n"$3}' untigs.gfa > untigs.fasta
#map trimmed reads onto miniasm untigs assembly
minimap2 untigs.fasta trim.fastq | gzip -1 > untigs.racon.paf.gz

echo "SECOND ASSEMBLY COMPLETE"

### Generate third assembly (from de novo, using shasta)
./shasta-Linux-0.7.0 --input trim.fastq

echo "THIRD ASSEMBLY COMPLETE"

cd ShastaRun
mv Assembly.gfa ../
mv Assembly.fasta ../
cd ../

#map trimmed reads onto shasta assembly
minimap2 Assembly.fasta trim.fastq | gzip -1 > shasta.racon.paf.gz

### Error correction/consensus sequence with Racon
#build consensus using the first assembly (from reference seq, using minimap2)
racon trim.fastq map_minimap_ref.racon.paf.gz map_minimap_ref.fasta > map_minimap_ref.consensus.fasta
#build consensus using the second assembly (from de novo, using minimap2)
racon trim.fastq untigs.racon.paf.gz untigs.fasta > map_minimap_denovo.racon.consensus.fasta
#build consensus (really just error correction) using the third assembly (from de novo, using shasta)
racon trim.fastq shasta.racon.paf.gz Assembly.fasta > map_shasta_denovo.consensus.fasta



#copy gfa files into the draft_assembly_gfa subdirectory and the consensus.fasta files into the consensus_sequences_fasta subdirectory
cp *.gfa output/draft_assembly_gfa
cp *.consensus.fasta output/consensus_sequence_fasta


echo "Script complete. Consensus sequences generated in /minION_align/output/draft_assembly_gfa. Assembly graphs generated in /minION_align/output/consensus_sequences_fasta. Examine these with Bandage."




############################ OPTIONAL ################################
##use medaka to generate consensus seq and generate variants
#medaka_haploid_variant -t 2 trim.fastq *_ref.fasta
##annotate medaka's vcf file so it can be uploaded into IGV
#medaka tools annotate consensus_to_ref.vcf *_ref.fasta mapped.sorted.bam variants.annotated.vcf

##use IGV to visualize
#igv
#genomes --> upload the reference fasta file
#file --> upload the vcf file (variants.vcf) AND the sorted bam file (mapped.sorted.bam)

##use Bandage to visualize different gfa assembly files
#Bandage
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





