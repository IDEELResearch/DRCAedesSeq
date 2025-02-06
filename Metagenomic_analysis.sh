###########################################################################################################################################################################
#Programmer: Wenqiao He
#Last update: December 2024
#Purpose: Aedes mosquito viral metagenomic analysis
#This is a bash script for metagenomic data analysis
###########################################################################################################################################################################


#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --time=140:00:00
#SBATCH --mem=100G

# 1. Merge paired-end reads   ###########################################################################################################################################

## Set working directory 
cd /path_with_raw_reads/ 

## Merge paired-end reads using bbmerge (version 38.96) with default parameters
module load bbmap
mkdir merged_reads/
mkdir unmerged_reads/
for i in *_R1_001.fastq.gz; do bbmerge.sh in1=${i} in2=${i%_R1_001.fastq.gz}_R2_001.fastq.gz out=/merged_reads/${i%_R1_001.fastq.gz}.fq outu1=/unmerged_reads/${i%_R1_001.fastq.gz}_R1_unmerged.fq outu2=/unmerged_reads/${i%_R1_001.fastq.gz}_R2_unmerged.fq ihist=/output_folder/${i%_R1_001.fastq.gz}_ihist.txt;done

# 2. Quality trimming and adapter removal using trimmomatic  (version 0.36)  #############################################################################################
module load trimmomatic
mkdir trimmomatic
for i in *.fq; do trimmomatic SE -threads 16 -phred33 ${i} /trimmomatic/${i%.fq}_trimmed.fq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;done

# 3. Filtering reads mapping to mosquito/human genomes   ##################################################################################################################
## Index reference genomes (use human genome (GRCh38.p14), Ae. albopictus genome (Aalbo_primary.1), Ae. aegypti (AaegL5.0) genome, and Ae. simpsoni CO1 gene (KT881399.1) from NCBI database as references)
module load bwa-mem2 # Version 2.2.1
cd /reference_genome_folder/
cat Aedes_albo.fna GCA_002204515.1_Aedes_aegypti_L5.0_genomic.fna Aedes_simpsoni_COI.fasta GRCh38_latest_genomic.fna > al_ae_simCOI_human.fasta
bwa-mem2 index al_ae_simCOI_human.fasta

## Align reads to references
cd /merged_reads/trimmomatic/
mkdir alignment/
for i in *_trimmed.fq; do bwa-mem2 mem -t 10 /reference_genome_folder/al_ae_simCOI_human.fasta ${i} -o /alignment/${i%.fq}_pairalignment_al_ae_simCOI_human.bam; done

## Get the mapped/unmapped reads
module load samtools # Version 1.21
cd alignment/
mkdir mapped_reads/
mkdir unmapped_reads/
for i in *.bam; do samtools view -@ 16 -b -F 4 ${i} > /mapped_reads/${i%.bam}_mapped.bam; done
for i in *.bam; do samtools view -@ 16 -b -f 4 ${i} > /unmapped_reads/${i%.bam}_unmapped.bam; done	

## Convert the reads that are not mapped to the mosquito/human genomes from BAM files to FASTQ files
cd unmapped_reads/
module load bedtools # Version 2.30
mkdir bam_to_fastq/
for i in *.bam; do bedtools bamtofastq -i ${i} -fq bam_to_fastq/${i%.bam}.fastq; done

# 4. De novo assemble contigs from lab mosquito and water control pools (These contigs are considered as contaminating contigs) ###########################################
module load python/3.9.6
module load spades # Version 4.0.0
cd bam_to_fastq/
metaspades.py --12 Lab_Mos_unmapped.fastq -o /contaminating_contigs/Lab_mosquito/
metaspades.py --12 Mos_Control_H2O_unmapped.fastq -o /contaminating_contigs/Water_control/

# 5. Align reads in field mosquito pools to contaminating contigs and remove the mapped reads #############################################################################
## Concatenate and index the contaminating contigs from lab mosquito and water control pools as references
cat /contaminating_contigs/Lab_mosquito/contigs.fasta /contaminating_contigs/Water_control/contigs.fasta > cat_Lab_water_control_contigs.fasta
bwa-mem2 index cat_Lab_water_control_contigs.fasta

## Align reads in field mosquito pools to contaminating contigs
mkdir align_contaminating_contigs/
for i in *.fastq; do bwa-mem2 mem -t 10 cat_Lab_water_control_contigs.fasta ${i} -o /align_contaminating_contigs/${i%.fastq}_pairalignment_contaminating_contigs.bam; done
cd align_contaminating_contigs/

## Get the mapped/unmapped reads 
mkdir mapped_reads/
mkdir unmapped_reads/
for i in *.bam; do samtools view -@ 16 -b -F 4 ${i} > /mapped_reads/${i%.bam}_mapped.bam; done
for i in *.bam; do samtools view -@ 16 -b -f 4 ${i} > /unmapped_reads/${i%.bam}_unmapped.bam; done	

## Convert the reads not mapped to contaminating contigs from BAM files to FASTQ files
cd unmapped_reads/
mkdir bam_to_fastq/
for i in *.bam; do bedtools bamtofastq -i ${i} -fq bam_to_fastq/${i%.bam}.fastq; done

# 6. Taxonomic classification of reads remaining in field mosquito pools using KrakenUniq (version 1.0.4) with its default mode and default nt and viral databases ######################

module load perl 
module load jellyfish (version 1.1.12)
# download and build krakenuniq default viral database
krakenuniq-download --db viral refseq/viral/Any viral-neighbors
krakenuniq-build --db viral --kmer-len 31 --threads 10 --taxids-for-genomes --taxids-for-sequences

# Taxonomic classification using viral database
mkdir align_krakenuniq_viral/
mkdir align_krakenuniq_viral/classified_sequences/

for i in *.fastq; do krakenuniq --db viral --classified-out align_krakenuniq_viral/classified_sequences/${i%.fastq}classfied_sequences.fastq --threads 16 --report-file align_krakenuniq_viral/${i%.fastq}_not_mpa.report  --output align_krakenuniq_viral/${i%.fastq}_not_mpa.txt ${i} ; done 

# download and build krakenuniq default nt database
krakenuniq-download --db nt --threads 10 nt
krakenuniq-build --db nt --kmer-len 31 --threads 10 --taxids-for-genomes --taxids-for-sequences

# Taxonomic classification using nt database
mkdir align_krakenuniq_nt/
mkdir align_krakenuniq_nt/classified_sequences/

for i in *.fastq; do krakenuniq --db nt --classified-out align_krakenuniq_nt/classified_sequences/${i%.fastq}classfied_sequences.fastq --threads 16 --report-file align_krakenuniq_nt/${i%.fastq}_not_mpa.report  --output align_krakenuniq_nt/${i%.fastq}_not_mpa.txt ${i} ; done 
