The code in this project accompanies an effort to investigate the virome of wild-caught Aedes mosquitoes from Kimpese, the Democratic Republic of the Congo (DRC). The DRC is the second largest country in Africa and has an estimated population of over 100 million. Its tropical climate is conducive to mosquito survival. Considerable attention is devoted to the high burden of malaria in the DRC, but arboviruses remain neglected, with limited studies of humans and mosquitoes to-date. To improve our understanding of the role of Aedes mosquito vectors in arboviral transmission in the DRC, we collected and sequenced the virome of mosquito pools from three areas of Kimpese, a region near the Angola border which has experienced recent arboviral outbreaks.

Full details of this project can be found within this pre-print: Evidence of dengue virus transmission and a diverse Aedes mosquito virome on the Democratic Republic of Congo-Angola border [doi]

Overview of the scripts found here:


**Bash scripts: 
a. Metagenomic analysis.sh   # Bash script for metagenomic data analysis to explore virome in Aedes mosquitoes in the DRC. This file contains the needed code for paired-end reads merging, adaptor removal, reads filtering, and taxonomic classification.
b. Nanopore targeted sequencing analysis.sh    # Bash script for nanopore data analysis to explore mosquito species and blood meal. This file contains the needed code for basecalling, quality-filtering, and blast analysis.

#Bash scripts rely on functioning installations of BBMerge (version 38.96), Trimmomatic (version 0.36), bwa-mem2 (version 2.2.1), SAMtools (version 1.21), SPAdes (version 4.0.0), KrakenUniq (version 1.0.4), Guppy (version 6.5.7), and Blast (version 2.14.1). 
#These scripts can be run on the sequencing data at NCBI SRA: BioProject IDs PRJNA1200724 and PRJNA1200731. Running of the "metagenomic data analysis.sh" takes more than 7 days but may vary depending on access to computer clusters resources and use of the SLURM system. Due to large size of the publicly-available dataset, analysis may not be feasible on a conventional desktop computer without substantial memory enhancements. Running of the "Nanopore targeted sequencing analysis.sh" takes less than 2 days.


**R scripts:
a. create_maps.R   # R script to create maps for sample collection sites. Data of DRC administrative boundaries (.shp file) are available at https://gadm.org/ and https://data.humdata.org/.
b. Metagenomic data visualization.R   # Visualization of Aedes mosquito viral metagenomic data and nanopore targeted sequencing data analysis. This scripts can be run on the example data files.


 


