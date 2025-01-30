###########################################################################################################################################################################
#Programmer: Wenqiao He
#Last update: December 2024
#Purpose: This is a bash script for Nanopore data analysis to investigate mosquito species and blood meal 
###########################################################################################################################################################################

#!/bin/bash
#SBATCH -n 16
#SBATCH -t 3-00:00:00
#SBATCH --mem 16Gb

## 1. Copy all fast5 file (Pass and Fail fast5 files) to a new folder ######################################################################################################
cp raw_read_folder/fast5_pass/*/* new_folder
cp raw_read_folder/fast5_fail/*/* new_folder

## 2. Use Guppy for basecalling and demultiplex ###########################################################################################################################
guppy_basecaller=/ont-guppy/bin/guppy_basecaller
guppy_barcoder=/ont-guppy-cpu/bin/guppy_barcoder
mkdir basecalling
mkdir demultiplex

${guppy_basecaller}     -i new_folder \
                        -s basecalling \
                        -c dna_r10.4.1_e8.2_400bps_5khz_hac_qs6.cfg \ ##(change phred quality score to 6 in the .cfg file)
                        --device auto


${guppy_barcoder} -i basecalling/pass \
                  -s demultiplex \
                  -r \
                  --enable_trim_barcodes \
                  --detect_mid_strand_adapter \
                  --trim_adapters \
                  --barcode_kits SQK-NBD114-96 \
                  -c configuration.cfg

## 3. Blast all pass reads ##############################################################################################################################################
module load blast
module load seqkit

cd demultiplex

for i in */*.fastq; do seqkit fq2fa ${i} -o ${i%.fastq}.fasta; done

for i in */*.fasta; do blastn -query ${i} -db nt -evalue 1e-6 -max_target_seqs 1 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore staxids sscinames scomnames'  -out ${i%.fastq}_blastn_1.txt -num_threads 10; done
