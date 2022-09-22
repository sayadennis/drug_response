#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="corr_dpmp_gdsc"
#SBATCH --output=drug_response/out/pearson_correlation_depmap_ccle_gdsc.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/exploratory/pearson_correlation_depmap_ccle_gdsc.py
