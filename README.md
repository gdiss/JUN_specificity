This repository contains all the code used to reproduce the analyses an results presented in The genetic architecture of protein interaction affinity and specificity.

Requirements:
cutadapt 1.1861
PEAR 0.9.1162 
R version 4.3.2 (2023-10-31)
perl 5.16.3 
DiMsum 1.3
MoCHI 1.0


R packages:
ggplot2_3.5.0               
gplots_3.1.3.1              
flowCore_2.14.2             
readxl_1.4.3                
ShortRead_1.60.0            
BiocParallel_1.36.0         
BiocGenerics_0.48.1         
rgl_1.3.1                   
stringr_1.5.1              
stringdist_0.9.12 


To run the analyses:
- Clone the repository or copy the files, create folder 001-merged_reads, 005-Mochi_output, 007-figures and tmp
- Download the fastq files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245326 and place them in the 000-data folder
- Merge barcode-variant association paired-end reads by entering ./001-merge_reads.sh in a bash terminal
- Determine barcode-variant association by copy-pasting the code from 002-bar_var_association.R inside a R console
- Process deepPCA data and build the count tables by copy-pasting the code from 003-binding_Score.R inside a R console
	This is the longest step and will take a few hours. Numbers of threads can be adjusted by changing the Ncore variable at line 119
- Compute binding scores with dimsum by entering ./004-dimsum.sh in a bash terminal
- Fit the thermodynamic model with Mochi by entering ./005-run_mochi.sh in a bash terminal
- Gather and merge data by copy-pasting the code from 006-merge_data.R inside a R console
- Perform all analyses and create all figures by copy-pasting the code from 007-analyses_and_figures.R inside a R console


