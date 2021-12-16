#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="runrecsys"
#SBATCH --output=drug_response/out/run_recsys.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

python drug_response/src/modeling/expression_knn.py
# python drug_response/src/modeling/expression_ubcf.py
python drug_response/src/modeling/vc_knn.py
# python drug_response/src/modeling/vc_ubcf.py
