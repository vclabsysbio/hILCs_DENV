# RNAseq_DV_ILCs Project
# Part 1 : raw data processing
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This script show command and parameters setting for raw data processing
# We start with fastq files from illumina sequencer (deposited under ancession GSEXXXXXXX)
# First we trim Nextera adapter using trimmomatic
# Read mapping was done using HISAT2
# TPM was obtained by StingTie. Read counts were obtain by HTseq-count
# Processed data of all samples were deposited on GSEXXXXXXX

# Nextera adapter reads were removed using Trimmomatic v0.36
java -jar <trimmomatic-0.36.jar> PE -threads 8 -phred33 <R1.fastq> <R2.fastq> <R1_U_output.fastq> <R2_U_output.fastq> <R1_P_output.fastq> <R2_P_output.fastq> ILLUMINACLIP:<Trimmomatic-0.36/adapters/NexteraPE-PE.fa>:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Sequence was mapped and aligned using HISAT2 (REF) with GRCh38 reference .
hisat2 -p 24 -k 1 --dta --no-mixed -x indexes/human_NCBI_GRCh38_index -1 data/"$filename"_R1_P.fastq.gz -2 data/"$filename"_R2_P.fastq.gz -S tmp/"$filename"_hisat2_dta_nomixed.sam 2> "$filename"_hisat.txt

# Gene abundance normalized in TPM was obtained by StringTie
stringtie -p 24 -G <ref_ann.gff> -e -B -A <gene_abund.tab> -o <out.gtf> <aligned_reads.bam>

# Raw reads counts were obtained with the HTseq-count
htseq-count -r name -s no -f bam -t exon -i gene_id <alignment_files> <gff_file>
