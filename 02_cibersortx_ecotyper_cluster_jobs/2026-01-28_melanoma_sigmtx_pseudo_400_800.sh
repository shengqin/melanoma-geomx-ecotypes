#!/bin/bash

# Define arrays for your configurations. Each array element corresponds to one job.
# Ensure that each array has the same number of elements and that corresponding elements are related to the same job.

outdirNomode=(
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_B_cell_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_CAF_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_Endothelial_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_Macrophage_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_Melanoma_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_Mural_cell_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_NK_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_T-CD4_pseudo/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_sigmtx_400_800_T-CD8_pseudo/"
)

mixture=(
  "/src/data/2026-01-28_pseudobulk_melanoma_B_cell_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_CAF_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_Endothelial_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_Macrophage_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_Melanoma_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_Mural_cell_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_NK_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_T-CD4_cpm.txt"
  "/src/data/2026-01-28_pseudobulk_melanoma_T-CD8_cpm.txt"
)

sigmatrix="/src/data/2026-01-28_melanoma_sigmtx_400_800.txt"

# Common parameters
username="YOUR_CIBERSORTX_EMAIL"
token="YOUR_CIBERSORTX_TOKEN"


# Loop through the arrays
for ((i = 0; i < ${#outdirNomode[@]}; i++)); do
  # Submit job
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_${i}
#SBATCH --time=0-00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --partition=emoding

mkdir -p ${outdirNomode[$i]}
singularity exec -c -B /oak/stanford/groups/emoding/analysis/brian/Cibersort/bulk_input/:/src/data \
  -B ${outdirNomode[$i]}:/src/outdir \
  /oak/stanford/groups/emoding/analysis/brian/Cibersort/fractions.sif /src/CIBERSORTxFractions \
  --username ${username} --token ${token} \
  --sigmatrix ${sigmatrix} \
  --mixture ${mixture[$i]} \
  --verbose TRUE
EOF
done
