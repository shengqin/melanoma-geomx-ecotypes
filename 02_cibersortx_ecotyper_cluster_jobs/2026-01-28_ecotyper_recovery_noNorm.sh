#!/bin/bash
#
#SBATCH --job-name=Ecotyper
#
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --partition=emoding

module load R/4.2.0

Rscript EcoTyper_recovery_bulk.R -d 2026-01-28_TCGA_SKCM_discovery \
    -m /oak/stanford/groups/emoding/analysis/brian/Ecotyper/datasets/2025/2025-01-24_GeoMx_expression_data_noNorm.txt \
    -a /oak/stanford/groups/emoding/analysis/brian/Ecotyper/datasets/2025/2025-01-27_GeoMxmetadata.txt \
    -c Segment_Tags,Pathophysiology,Path_sub_type \
    -o 2026-01-28_GeoMx_recovery_noNorm
