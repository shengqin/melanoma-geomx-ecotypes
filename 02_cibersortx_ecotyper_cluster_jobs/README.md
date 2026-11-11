# 02. CIBERSORTx and EcoTyper cluster jobs

SLURM jobs run on Stanford's Oak cluster. CIBERSORTx runs from its Singularity images
(`fractions.sif`, `hires.sif`); EcoTyper runs with R 4.2.

## Run order
| Step | Script | What it does | Main output |
|---|---|---|---|
| 1 | `2026-01-28_melanoma_sigmtx_400_800.sh` | Builds the signature matrix from the single-cell reference (`--single_cell TRUE`, 400–800 genes per population, `--fraction 0`, `--filter FALSE`). | Signature matrix, 2,260 genes (Fig. 5b) |
| 2 | `2026-01-28_melanoma_sigmtx_pseudo_400_800.sh` | Deconvolves the nine sets of pseudobulk mixtures with that matrix (one job per population). | `CIBERSORTx_Results.txt` per population (Fig. 5c) |
| 3 | `2026-01-28_melanoma_hires.sh` | CIBERSORTx high-resolution mode on TCGA-SKCM (protein-coding TPM) and on the GeoMx matrix, with S-mode batch correction and protein-coding genes only. | Cell-type-specific expression per sample (`CIBERSORTxHiRes_NA_<cell type>_Window36.txt`), used for Fig. 8 |
| 4 | `2026-01-28_TCGA_SKCM_discovery.sh` + `2026-01-28_TCGA_SKCM_discovery.yaml` | EcoTyper discovery on TCGA-SKCM, with the precomputed cell-type fractions. | Cell states and ecotypes (Fig. 6) |
| 5 | `2026-01-28_ecotyper_recovery_noNorm.sh` | Applies the TCGA-SKCM ecotype model to the GeoMx AOIs (normal skin removed), annotated by segment, pathophysiology and subtype. | Per-AOI state abundances and ecotype calls (Fig. 7) |

The cell-type fraction tables that steps 4 and 5 and the figures read (`2026-01-28_TCGA_SKCM_fraction.txt`,
`2026-01-28_cibersort_fraction.txt`) are CIBERSORTx fractions obtained with the same signature matrix.

## Discovery settings used in the paper
The paper's discovery run used a cophenetic-coefficient cutoff of **0.997** and **2–4 states** per cell type
(output folder `2026-01-28_TCGA_SKCM_Discovery_997_4_4`), with 60 NMF restarts and at least 4 states per ecotype.
The YAML here still has the values of an earlier exploratory run (cutoff 0.95, up to 6 states, output folder
`..._95_6_4`). To reproduce the paper, set `Cophenetic coefficient cutoff : 0.997` and
`Maximum number of states per cell type : 4`.

## Before running
- Replace `YOUR_CIBERSORTX_EMAIL` and `YOUR_CIBERSORTX_TOKEN` with your own CIBERSORTx credentials
  (https://cibersortx.stanford.edu), and `YOUR_EMAIL` in the `#SBATCH --mail-user` lines.
- Change the `/oak/...` paths, the partition name (`emoding`) and the Singularity image paths to your cluster.
- EcoTyper (`EcoTyper_discovery_bulk.R`, `EcoTyper_recovery_bulk.R`) comes from https://github.com/digitalcytometry/ecotyper;
  run steps 4 and 5 from the EcoTyper folder.
