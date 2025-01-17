The code in this project accompanies an effort to investigate the virome of wild-caught Aedes mosquitoes from Kimpese, the Democratic Republic of the Congo (DRC). The DRC is the second largest country in Africa and has an estimated population of over 100 million. Its tropical climate is conducive to mosquito survival. Considerable attention is devoted to the high burden of malaria in the DRC, but arboviruses remain neglected, with limited studies of humans and mosquitoes to-date. To improve our understanding of the role of Aedes mosquito vectors in arboviral transmission in the DRC, we collected and sequenced the virome of mosquito pools from three areas of Kimpese, a region near the Angola border which has experienced recent arboviral outbreaks.

Full details of this project can be found within this pre-print: Evidence of dengue virus transmission and a diverse Aedes mosquito virome on the Democratic Republic of Congo-Angola border (DOI: https://doi.org/10.1101/2025.01.16.633031).


Overview of the scripts found here:


**Bash scripts:

a. Metagenomic analysis.sh # Bash script for metagenomic data analysis to explore virome in Aedes mosquitoes in the DRC. This file contains the needed code for paired-end reads merging, adaptor removal, reads filtering, and taxonomic classification.

b. Nanopore targeted sequencing analysis.sh # Bash script for nanopore data analysis to explore mosquito species and blood meal. This file contains the needed code for basecalling, quality-filtering, and blast analysis.

*# Bash scripts rely on functioning installations of BBMap (version 38.96), Trimmomatic (version 0.36), bwa-mem2 (version 2.2.1), SAMtools (version 1.21), SPAdes (version 4.0.0), KrakenUniq (version 1.0.4), Guppy (version 6.5.7), and Blast (version 2.14.1). Installing all tools could take 1–1.5 hours. Use package managers (e.g., conda, brew, or apt) for faster installations.

#These scripts can be run on the sequencing data available at NCBI SRA: BioProject IDs PRJNA1200724 and PRJNA1200731. To run the bash scripts, please download the datasets to your local computer or compute clusters. Running of the "metagenomic data analysis.sh" usually takes days to weeks but may vary depending on access to computer clusters resources and use of the SLURM system. Due to large size of the publicly available dataset, analysis may not be feasible on a conventional desktop computer without substantial memory enhancements. Running of the "Nanopore targeted sequencing analysis.sh" usually takes less than 2 days when appropriate computing cluster resources and the SLURM system are utilized.


**R scripts:

a. create_maps.R # R script to create maps for sample collection sites. Data of the DRC administrative boundaries (.shp files) are available at https://gadm.org/ and https://data.humdata.org/.

b. Metagenomic data visualization.R # Visualization of Aedes mosquito viral metagenomic data and nanopore targeted sequencing data analysis. This script can be executed using the example data files.

*# R scripts rely on functioning installations of R software (version 4.2.0) and RStudio (version 2022.02.2). Installing R and RStudio typically takes about 10–30 minutes, depending on your computer's speed, operating system, and familiarity with the process. The required packages are specified within the R scripts. Running of the "create_maps.R" and "Metagenomic data visualization.R" on the example data files takes up to 30 minutes.
