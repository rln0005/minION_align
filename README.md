# minION Alignment Sequence Analysis   
*This script was initially intended to be used on sequence data from a circular haploid plasmid genome of approximately 37 kb that was sequenced using Oxford Nanopore's (ONT) minION device. It was intended to be used as a quick reference to determine if CRISPR-modified plasmid genomes were accurately constructed in the lab prior to transfection. Visualization of assmeblies using Bandage and variants using IGV allows regions of dissimilarity to be further investigated.* 
  
This script (minION_align.sh) is designed to use sequence data from a Nanopore minION (or any long-read sequence data) specifically outputted from the minKNOW program (in .fastq format). The script will combine all .fastq files for an individual barcode into a single .fastq file, filter the reads using NanoFilt, map reads in a few different ways with different packages -- including de novo assembly using Shasta, de novo assembly using Minimap2/Miniasm, and assembly using a reference sequence using Minimap2/Miniasm. Racon will be used for error correction and to generate a consensus sequence for each of the assemblies.  
to a reference sequence using minimap2, convert from .sam to .bam using samtools, call variants using Sniffles to generate a .vcf file, and visualize the output using IGV.  
  
Some difficulties encountered that are specific to the aims of this project include: 1) draft assembly/consensus sequence generation for relatively small CIRCULAR plasmids (some reads span the entire plasmid and more) and 2) difficulty in SNP/indel variant calling using long-read HAPLOID sequence data. As new tools are developed and pre-existing tools are improved, this script/repository will be updated.     
  
**There are 5 general steps in this script: 1) Combine & Filter 2) Align/Map 3) Generate Consensus 4) Variant Call 5) Visualize Assemblies/Variants.** The 5th step is not executed by the script itself (see info below for options to visualize results).    

## Concatenate fastq files from sequence data generated from minKNOW ONT software
Use CL tools (such as grep) to extract barcodes for samples so that they can be processed individually (for example, different reference sequences for barcodes 1 and 2). Combine all fastq files for a given barcode into a single fastq file.  
These reads are already labeled 'pass' -- I think they have already been filtered to an extent. In the future, I plan to implement Guppy as a basecaller in this script.  
## Trim/Filter Reads
Using NanoFilt with standard/fixed filtering conditions. These parameters can be adjusted within the script if desired.  
`NanoFilt -l 500 -q 12 --headcrop 50 < *_all.fastq > trim.fastq `  
  
`-l 500` option filters based on average sequence length of 500 bp  
`-q 12` option filters based on average quality score of 12  
`--headcrop 50` filters by removing first 50bp of sequence  

## ASSEMBLY USING REFERENCE with Minimap/Miniasm:   
Align reads with a reference sequence using minimap2 to align with reference sequence (used for ONT genomic reads).  
`minimap2 -ax map-ont --MD *_ref.fasta trim.fastq > map_MM_ref_mapped.sam`  
  
 `-a` option generates CIGAR and output alignment files in sam format (default is paf format)  
 `-x` in combination with `map-ont` specifies the seq data is derived from an Oxford Nanopore device to enhance read alignment.   
 `--MD` option ensures that the MD tag will be included in the sam output file. This is necessary for Sniffles.  
   
 Then use samtools to convert files from sam to bam:  `samtools view -bS mapped.sam > mapped.bam`  
 Then use samtools to sort & index the bam file: `samtools sort -o mapped.sorted.bam mapped.bam` `samtools index mapped.sorted.bam`  
  
# DE NOVO ASSEMBLY with Shasta (generates contigs):  
This will generate a new subfolder called "ShastaRun" which will contain many files, including Assembly.gfa and Assembly.fasta.   
`shasta-Linus-0.7.0 --input *_all.fastq`  
Error Correction using Racon can be applied to the Shasta mapped file -- to compare consensus assemblies generated using different alignment programs (Shasta vs. Minimap2/Miniasm)   
  
# DE NOVO ASSEMBLY with Minimap/Miniasm (generates untigs):   
Map Reads Onto Themselves -- Identifies Overlaps: `minimap2 -x ava-ont trim.fastq trim.fastq | gzip -1 > overlaps.paf.gz `  
Assemble Untigs: `miniasm -f trim.fastq overlaps.paf.gz > untigs.gfa`    
Convert gfa to fasta: ``awk `/^S/{print ">"$2"\n"$3}' untigs.gfa > untigs.fasta``    
  
# ERROR CORRECTION/GENERATE CONSENSUS SEQUENCE with Racon:  
Racon was developed to complement minimap2/miniasm pipeline but can be used for any long-read assembly. It is used to correct draft assemblies.  
## Map Trimmed Reads onto Miniasm Assembly  
`minimap2 untigs.fasta trim.fastq > untigs.racon.paf`  
## Build Consensus Using Trimmed Reads and Untigs/Minimap.Racon Assembly (ref seq assembly)    
`racon trim.fastq untigs.racon.paf untigs.fasta > untigs.racon.consensus.fasta `  
## Build Consensus Using Trimmed Reads and Contigs/Shasta Assembly (de novo assembly)  
`racon trim.fastq `




## Variant Calling    
It has been difficult to find a reliable SNP/indel variant caller for long read data on a haploid genome. Medaka has a command `medaka_haploid_variant` that is supposed to be used to generate variants on haploid sequences, but there are some bugs in their code. See below.  
`medaka_haploid_variant -t 2 trim.fastq *_ref.fasta `  
****BUG IN MEDAKA CODE -- the output of medaka_haploid_variant includes a vcf file (consensus_to_ref.vcf) but this file is generated with blanks (whitespace) in the "INFO" column -- which is required to be populated to be visualized in IGV. This has already been reported to medaka team (see https://github.com/nanoporetech/medaka/issues/286). I attempted to overcome this with the following: `medaka tools annotate consensus_to_ref.vcf *_ref.fasta mapped.sorted.bam variants.annotated.vcf` This will output a file variants.annotated.vcf which can be uploaded into IGV.****  
I believe `medaka_haploid_variant` includes the `medaka_consensus` step within it, but this is not directly clear from the manual.   
  
  
## Analyze and Visualize Variants Using IGV
Type `igv` on the command line to open the IGV program -- upload files directly.  
Upload the reference fasta file using Genomes -->  upload from file.  
Upload the vcf file using File --> upload from file. Upload the sorted bam file (mapped.sorted.bam) using File --> upload from file.   
  
## Analyze and Visualize Assemblies Using Bandage  
Type `Bandage` on the command line to open the Bandage GUI -- assemblies (gtf files) can be uploaded & visualized directly   
   
  
# OTHER EXPLORED/POTENTIAL AVENUES  
## Structural Variant Calling Using Sniffles  
Sniffles is used for *structural* variants, which is not the intention of this current project. But here's some info on it: 
Using Sniffles (https://github.com/fritzsedlazeck/Sniffles)  
`sniffles -m mapped.sorted.bam -v variants.vcf`  
   
`-m` option is used to indicate the mapped bam file name  
`-v` option is used to indicate the output vcf file name    

## De Novo Assembly Using Canu  
Canu is another popular assembly package. This was attempted using this dataset, but the program ran for approximately 12+ hours (on a standard laptop, not HPC) and was still running. It was purposefully killed (and otherwise might have taken days to complete the run). This option should be explored using HPC.  
`canu -p output genomeSize=37.8k -corMemory=3 -corThreads=3 -redMemory=3 -redThreads=2 -oeaMemory=3 -oeaThreads=2 -nanopore-raw *_all.fastq`  

##  Visualize Differences in Consensus vs. Reference with DNAdiff/MUMmer  
Since we are not expecting major variations between the sequenced data and the reference sequence, MUMmer (https://github.com/mummer4/mummer) was installed and dnadiff was attempted using `dnadiff <reference.fasta> <consensus.fasta>` but repeatedly got "ERROR: Failed to run show-snps, aborting." Consensus quality can also be explored using dnadiff. This tool should be further explored for the purposes of this script.  

## Assembly Using Reference Sequence with Graphmap   
Graphmap (https://github.com/isovic/graphmap) tool to map/align long reads to a reference sequence. Unfortunately I was unable to get this package successfully installed on my machine, but this could be another potential mapping tool to explore.  



# *REQUIRED PACKAGES* 
 - NanoFilt (https://github.com/wdecoster/nanofilt)  
 - Minimap2 (https://github.com/lh3/minimap2)  
 - IGV (https://github.com/igvteam/igv) 
 - Bandage (https://github.com/rrwick/Bandage)  
 - Racon (https://github.com/isovic/racon)
 - Shasta (https://github.com/chanzuckerberg/shasta) 
 - Medaka (https://github.com/nanoporetech/medaka)

# *PARAMETERS / REQUIREMENTS*
- Script is written in bash and can be executed on the command line using 'sh minION_align.sh'
- Script should be executed within a directory titled 'minION_align'. 
- Sequence data files (fastq) should be contained in a subdirectory titled 'data'. (minION_align/data) 
- A reference genome in fasta format should contain 'ref.fasta' in the file name and should be contained in the 'minION_align' main directory. 
- All output files will be deposited in the 'minION_align' main directory unless otherwise specified. 
- The fastq file combining all fastq files for an individual barcode will include the file suffix 'all.fastq' and will be deposited in the 'minION_align' main directory.


