#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="ExtractHaseResults"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load squashfs/4.4

# NB! Following will extract all the genes and SNPs! In case of full mapping, 
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=[path to folder where Nextflow executable is]

input_folder=[folder with input .parquet files]
output_file=[path to tab-separated output file]
nr_covariates=[Number of covariates in the model]
snp_ref=bin/hase/data/1000G-30x.ref.gz

#Optional arguments:
#--PThresh               P-value threshold for filtering the results. Defaults None.
#--PhenFilter            Filter to include phenotypes into the association results. One of two: first, file with the list of phenotypes to include. In that case, no header expected. Second, individual phenotype IDs specified and separated by space. Defaults to None.
#--SnpFilter             Filter to include SNPs/variants into the association results.  One of two: first, file with the list of SNPs/variants to include. In that case, no header expected. Second, individual SNP IDs specified and separated by space. Defaults to None.
#--Chunks                Number of chunks which was used for original meta-analysis. Defaults to 100.

NXF_VER=20.10.0 ${nextflow_path}/nextflow run ExtractHaseResults.nf \
--inputfolder ${input_folder} \
--outputfile ${output_file} \
--NumberOfCovariates ${nr_covariates} \
--SnpRef ${snp_ref} \
--PThresh 5e-8 \
-resume \
-profile slurm,singularity

# Compress the output
gzip ${output_file}