# minION Sequence Analysis 
This script is designed to take multiplexed sequence data from a Nanopore minION (.fastq) and align individually with user-inputted reference sequences to identify any regions of dissimilarity.

## De-Multiplex
Use CL tools (such as grep) to extract barcodes for samples so that they can be processed individually (for example, different reference sequences for barcodes 1 and 2)

## Align with Reference Sequence of Choice
Using HiSAT (or maybe BWA...tbd) to index and align with reference sequence. Then use samtools to convert files from SAM to BAM.

## Variant Calling 
Not sure yet how I will go about doing this

## Analyze
Ideally, I would like to generate an easy to read visual of the alignment - where it aligns and where/if it does not. In the regions that it doesn't align, I would like it to show why it doesn't align (for example, have two rows with the top row containing the sequenced data and the bottom row containing the ref seq...regions that don't align are indicated in red font and regions that do align are green)

