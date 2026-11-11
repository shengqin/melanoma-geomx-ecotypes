#!/bin/bash

# Define arrays for your configurations. Each array element corresponds to one job.
# Ensure that each array has the same number of elements and that corresponding elements are related to the same job.

outdir=(
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_TCGA_hires/"
  "/oak/stanford/groups/emoding/analysis/brian/Cibersort/2026-01-28_melanoma_GeoMx_hires/"
)

mixture=(
  "/src/data/2025-01-29_SKCM_coding_symbol_TPM.txt"
  "/src/data/2025-01-24_GeoMx_expression_data.txt"
)

sigmatrix=(
  "/src/data/2026-01-28_melanoma_sigmtx_400_800.txt"
  "/src/data/2026-01-28_melanoma_sigmtx_400_800.txt"
)

refsample=(
  "/src/data/2026-01-28_melanoma_obj_feature_4000.txt"
  "/src/data/2026-01-28_melanoma_obj_feature_4000.txt"
)
# Common parameters
username="YOUR_CIBERSORTX_EMAIL"
token="YOUR_CIBERSORTX_TOKEN"
genelist="/src/data/2025-03-26_coding_gene_list.txt"

# Loop through the arrays
for ((i = 0; i < ${#outdir[@]}; i++)); do
  # Submit job
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_${i}
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --partition=emoding

mkdir -p ${outdir[$i]}
singularity exec -c -B /oak/stanford/groups/emoding/analysis/brian/Cibersort/bulk_input/:/src/data \
  -B ${outdir[$i]}:/src/outdir \
  /oak/stanford/groups/emoding/analysis/brian/Cibersort/hires.sif /src/CIBERSORTxHiRes \
  --username ${username} --token ${token} \
  --sigmatrix ${sigmatrix[$i]} \
  --mixture ${mixture[$i]} \
  --refsample ${refsample[$i]} \
  --subsetgenes ${genelist} \
  --rmbatchSmode TRUE
EOF
done
