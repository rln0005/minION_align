# minION 'Quick & Dirty' Sequence Analysis 
This script is designed to use sequence data from a Nanopore minION (or any long-read sequence data) specifically outputted from the minKNOW program (in .fastq format). The script will combine all .fastq files for an individual barcode into a single .fastq file, filter the reads using NanoFilt, map reads to a reference sequence using minimap2, convert from .sam to .bam using samtools, call variants using Sniffles to generate a .vcf file, and visualize the output using IGV. The script was designed to be used on a circular haploid plasmid genome of approximately 36kb. It was intended to be used as a 'quick & dirty' reference to determine if CRISPR-modified plasmid genomes were accurately constructed in the lab prior to transfection. Visualization with IGV allows regions of dissimilarity to be further investigated.

## De-Multiplex and Concatenate FastQ Files
Use CL tools (such as grep) to extract barcodes for samples so that they can be processed individually (for example, different reference sequences for barcodes 1 and 2). Combine all fastq files for a given barcode into a single fastq file. 

## Filter Reads
Using NanoFilt with standard/fixed filtering conditions

## Align with Reference Sequence of Choice
Using minimap2 to index and align with reference sequence. Then use samtools to convert files from sam to bam.

## Variant Calling 
Using Sniffles 

## Analyze and Visualize
Using IGV
