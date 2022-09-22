#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=80G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ol_ccle_gdsc"
#SBATCH --output=drug_response/out/overlap_ccle_gdsc_pdx.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/exploratory/overlap_ccle_gdsc_pdx.py
