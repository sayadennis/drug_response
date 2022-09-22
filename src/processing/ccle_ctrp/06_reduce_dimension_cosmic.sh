#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="reduce_cosmic"
#SBATCH --output=drug_response/out/06_reduce_dimension_cosmic.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/src/processing/06_reduce_dimension_cosmic.py
