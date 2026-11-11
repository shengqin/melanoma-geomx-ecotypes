#!/bin/bash
#
#SBATCH --job-name=Ciber
#
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --partition=emoding

mkdir -p /oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800/
singularity exec -c -B /oak/stanford/groups/emoding/analysis/brian/Cibersort/bulk_input/:/src/data \
  -B /oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800/:/src/outdir \
  /oak/stanford/groups/emoding/analysis/brian/Cibersort/fractions.sif /src/CIBERSORTxFractions \
  --username YOUR_CIBERSORTX_EMAIL --token YOUR_CIBERSORTX_TOKEN \
  --refsample /src/data/2026-01-28_melanoma_obj_feature_4000.txt \
  --single_cell TRUE \
  --fraction 0 \
  --filter FALSE \
  --G.min 400 \
  --G.max 800