# minION 'Quick & Dirty' Sequence Analysis 
This script is designed to use sequence data from a Nanopore minION (or any long-read sequence data) specifically outputted from the minKNOW program (in .fastq format). The script will combine all .fastq files for an individual barcode into a single .fastq file, filter the reads using NanoFilt, map reads to a reference sequence using minimap2, convert from .sam to .bam using samtools, call variants using Sniffles to generate a .vcf file, and visualize the output using IGV. The script was designed to be used on a circular haploid plasmid genome of approximately 36kb. It was intended to be used as a 'quick & dirty' reference to determine if CRISPR-modified plasmid genomes were accurately constructed in the lab prior to transfection. Visualization with IGV allows regions of dissimilarity to be further investigated.

## De-Multiplex and Concatenate fastq Files
Use CL tools (such as grep) to extract barcodes for samples so that they can be processed individually (for example, different reference sequences for barcodes 1 and 2). Combine all fastq files for a given barcode into a single fastq file. 
## Filter Reads
Using NanoFilt with standard/fixed filtering conditions. (https://github.com/wdecoster/nanofilt)
## Align with Reference Sequence of Choice
Using minimap2 to index and align with reference sequence. Then use samtools to convert files from sam to bam. (https://github.com/lh3/minimap2)
## Variant Calling 
Using Sniffles (https://github.com/fritzsedlazeck/Sniffles) 
## Analyze and Visualize
Using IGV (https://github.com/igvteam/igv). The script will end with executing 'igv', which will open the IGV application. From there, the user must manually upload the files generated from this script. Upload the reference fasta file using Genomes -->  upload from file. Upload the vcf file using File --> upload from file. Upload the sorted bam file (mapped.sorted.bam) using File --> upload from file. 

## *PARAMETERS / REQUIREMENTS*
- Requires NanoFilt, Samtools, Minimap2, Sniffles, and IGV which can all be downloaded from their github links. 
- Script is written in bash and can be executed on the command line using 'sh minION_align.sh'
- Script should be executed within a directory titled 'minION_align'. 
- Sequence data files (fastq) should be contained in a subdirectory titled 'data'. 
- A reference genome in fasta format should contain 'ref.fasta' in the label and should be in the 'minION_align' main directory.
- All output files will be deposited in the 'minION_align' main directory.
- The fastq file combining all fastq files for an individual barcode will include the label 'all.fastq' and will be deposited in the 'minION_align' main directory.


