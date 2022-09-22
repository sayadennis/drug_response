#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ctrp_target"
#SBATCH --output=drug_response/out/03_ctrp_generate_target.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/scripts/03_ctrp_generate_target.py
