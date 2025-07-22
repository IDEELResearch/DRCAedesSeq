###########################################################################################################################################################################
#Programmer: Wenqiao He
#Last update: July 2025
#Purpose: Aedes mosquito viral metagenomic analysis
#This is a bash script for phylogenetic analysis
###########################################################################################################################################################################


#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --time=140:00:00
#SBATCH --mem=100G

# 1. select best-fit model using IQ-TREE   ###########################################################################################################################################

## Set working directory 
cd /path_with_raw_reads/ 

## Load or install model
module load iqtree #version 2.4.0
iqtree -s DENV2_aligned.fasta -bb 1000 -nt 8
iqtree -s DENV4_aligned.fasta -bb 1000 -nt 8
iqtree -s Pan_DENV_aligned.fasta -bb 1000 -nt 8
iqtree -s Bat_faecal_associated_dicistrovirus_4_aligned.fasta -bb 1000 -nt 8
iqtree -s Pegivirus_aligned.fasta -bb 1000 -nt 8
iqtree -s human_blood_associated_dicistrovirus_K4_aligned.fasta -bb 1000 -nt 8
iqtree -s human_blood_associated_dicistrovirus_K5_1_aligned.fasta -bb 1000 -nt 8
iqtree -s human_blood_associated_dicistrovirus_K5_aligned.fasta -bb 1000 -nt 8

## Check the *.iqtree files and select the best-fitting models based on BIC values for use in RAxML-NG.

# 2. Generate phylogenetic trees using RAxML-NG #####################################################################################################################################
## Load or install RAxML-NG by using the pre-compiled binary (version 1.2.2) and generate the phylogenetic trees
/path_with_the_raxml-ng_file/raxml-ng --all --msa DENV2_aligned.fasta --outgroup KR011349.2 --model GTR+F+I+G4 --bs-trees 1000 --threads 8 --prefix DENV2_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa DENV4_aligned.fasta --outgroup EU179857.1 --model GTR+F+I+G4 --bs-trees 1000 --threads 8 --prefix DENV4_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa Pan_DENV_aligned.fasta --outgroup NC012532.1 --model GTR+F+I+G4 --bs-trees 1000 --threads 8 --prefix Pan_DENV_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa Pan_DENV_aligned.fasta --outgroup NC012532.1 --model GTR+F+I+G4 --bs-trees 1000 --threads 8 --prefix Pan_DENV_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa Bat_faecal_associated_dicistrovirus_4_aligned.fasta --model TPM3+F+R2 --bs-trees 1000 --prefix Bat_faecal_associated_dicistrovirus_4_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa Pegivirus_aligned.fasta --model HKY+F+G4 --outgroup AF176573.1,OQ832071.1 --bs-trees 1000 --prefix Pegivirus_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa human_blood_associated_dicistrovirus_K4_aligned.fasta --model TPM3+F+I+G4 --bs-trees 1000 --prefix human_blood_associated_dicistrovirus_K4_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa human_blood_associated_dicistrovirus_K5_1_aligned.fasta --model TPM3+F+R3 --bs-trees 1000 --prefix human_blood_associated_dicistrovirus_K5_1_aligned_raxml_ng
/path_with_the_raxml-ng_file/raxml-ng --all --msa human_blood_associated_dicistrovirus_K5_aligned.fasta --model TPM3+F+R2 --bs-trees 1000 --prefix human_blood_associated_dicistrovirus_K5_aligned_raxml_ng
