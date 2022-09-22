#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ccle_mxs"
#SBATCH --output=drug_response/out/01_ccle_generate_vc_exp_matrices.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/src/processing/01_ccle_generate_vc_exp_matrices.py
