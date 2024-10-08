#!/bin/bash
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=200GB
#SBATCH --partition=basic
#SBATCH --job-name=dio_rnaseq
#SBATCH --time=3-00:00:00

module load subread
module load star
module load snakemake
module load cutadapt

snakemake -s Snakemake.rnaseq run --cores 8
