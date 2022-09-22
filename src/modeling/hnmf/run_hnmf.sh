#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="hnmf"
#SBATCH --output=drug_response/out/run_hnmf.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate cancerdrugenv

din='/projects/b1131/saya/drug_response/depmap/01_cleaned'

python drug_response/src/modeling/hnmf/run_hnmf.py \
    --X1 $din/cleaned_vc_mx.csv \
    --X2 $din/cleaned_exp_mx.csv \
    --X1loss 'kl' \
    --X2loss 'l2' \
    --max_k 500 \
    > drug_response/model_performance/hnmf/results_hnmf.txt
#
