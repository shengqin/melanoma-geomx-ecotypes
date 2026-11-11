#!/bin/bash
#SBATCH --job-name=Ecotyper
#SBATCH --time=02-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --partition=emoding

module load R/4.2.0

Rscript EcoTyper_discovery_bulk.R -c 2026-01-28_TCGA_SKCM_discovery.yaml