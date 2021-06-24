#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="cellmeta_sum"
#SBATCH --output=drug_response/out/cell_lines_meta_summary.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/src/data_summary/cell_lines_meta_summary.py
