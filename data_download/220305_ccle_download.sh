#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ccle_dl"
#SBATCH --output=drug_response/out/220305_ccle_download.out

"""
Downloaded DepMap Public 22Q1 data on 2022-03-07 from https://depmap.org/portal/download/
"""

cd /projects/b1131/saya/drug_response/depmap/

# README.txt
wget https://ndownloader.figshare.com/files/34008500 -O README.txt
# Achilles_gene_effect.csv
wget https://ndownloader.figshare.com/files/34008383 -O Achilles_gene_effect.csv
# Achilles_gene_effect_uncorrected.csv
wget https://ndownloader.figshare.com/files/34008380 -O Achilles_gene_effect_uncorrected.csv
# Achilles_gene_dependency.csv
wget https://ndownloader.figshare.com/files/34008377 -O Achilles_gene_dependency.csv
# Achilles_common_essentials.csv
wget https://ndownloader.figshare.com/files/34008353 -O Achilles_common_essentials.csv
# Achilles_guide_efficacy.csv
wget https://ndownloader.figshare.com/files/34008359 -O Achilles_guide_efficacy.csv
# Achilles_cell_line_efficacy.csv
wget https://ndownloader.figshare.com/files/34008347 -O Achilles_cell_line_efficacy.csv
# Achilles_cell_line_growth_rate.csv
wget https://ndownloader.figshare.com/files/34008356 -O Achilles_cell_line_growth_rate.csv
# CRISPR_dataset_sources.csv
wget https://ndownloader.figshare.com/files/34008476 -O CRISPR_dataset_sources.csv
# CRISPR_gene_effect.csv
wget https://ndownloader.figshare.com/files/34008491 -O CRISPR_gene_effect.csv
# CRISPR_gene_dependency.csv
wget https://ndownloader.figshare.com/files/34008485 -O CRISPR_gene_dependency.csv
# CRISPR_common_essentials.csv
wget https://ndownloader.figshare.com/files/34008473 -O CRISPR_common_essentials.csv
# common_essentials.csv
wget https://ndownloader.figshare.com/files/34008470 -O common_essentials.csv
# nonessentials.csv
wget https://ndownloader.figshare.com/files/34008497 -O nonessentials.csv
# Achilles_raw_readcounts.csv
wget https://ndownloader.figshare.com/files/34008395 -O Achilles_raw_readcounts.csv
# Achilles_raw_readcounts_failures.csv
wget https://ndownloader.figshare.com/files/34008392 -O Achilles_raw_readcounts_failures.csv
# Achilles_logfold_change.csv
wget https://ndownloader.figshare.com/files/34008443 -O Achilles_logfold_change.csv
# Achilles_logfold_change_failures.csv
wget https://ndownloader.figshare.com/files/34008389 -O Achilles_logfold_change_failures.csv
# Achilles_guide_map.csv
wget https://ndownloader.figshare.com/files/34008362 -O Achilles_guide_map.csv
# Achilles_replicate_map.csv
wget https://ndownloader.figshare.com/files/34008398 -O Achilles_replicate_map.csv
# Achilles_replicate_QC_report_failing.csv
wget https://ndownloader.figshare.com/files/34008401 -O Achilles_replicate_QC_report_failing.csv
# Achilles_dropped_guides.csv
wget https://ndownloader.figshare.com/files/34008350 -O Achilles_dropped_guides.csv
# Achilles_high_variance_genes.csv
wget https://ndownloader.figshare.com/files/34008365 -O Achilles_high_variance_genes.csv
# CCLE_RNAseq_reads.csv
wget https://ndownloader.figshare.com/files/34008455 -O CCLE_RNAseq_reads.csv
# CCLE_expression_full.csv
wget https://ndownloader.figshare.com/files/34008407 -O CCLE_expression_full.csv
# CCLE_expression.csv
wget https://ndownloader.figshare.com/files/34008404 -O CCLE_expression.csv
# CCLE_expression_transcripts_expected_count.csv
wget https://ndownloader.figshare.com/files/34008419 -O CCLE_expression_transcripts_expected_count.csv
# CCLE_expression_proteincoding_genes_expected_count.csv
wget https://ndownloader.figshare.com/files/34008410 -O CCLE_expression_proteincoding_genes_expected_count.csv
# CCLE_RNAseq_transcripts.csv
wget https://ndownloader.figshare.com/files/34008506 -O CCLE_RNAseq_transcripts.csv
# CCLE_segment_cn.csv
wget https://ndownloader.figshare.com/files/34008464 -O CCLE_segment_cn.csv
# CCLE_wes_segment_cn.csv
wget https://ndownloader.figshare.com/files/34008467 -O CCLE_wes_segment_cn.csv
# CCLE_gene_cn.csv
wget https://ndownloader.figshare.com/files/34008428 -O CCLE_gene_cn.csv
# CCLE_wes_gene_cn.csv
wget https://ndownloader.figshare.com/files/34008482 -O CCLE_wes_gene_cn.csv
# CCLE_fusions.csv
wget https://ndownloader.figshare.com/files/34008422 -O CCLE_fusions.csv
# CCLE_fusions_unfiltered.csv
wget https://ndownloader.figshare.com/files/34008425 -O CCLE_fusions_unfiltered.csv
# CCLE_mutations.csv
wget https://ndownloader.figshare.com/files/34008434 -O CCLE_mutations.csv
# CCLE_mutations_bool_hotspot.csv
wget https://ndownloader.figshare.com/files/34008446 -O CCLE_mutations_bool_hotspot.csv
# CCLE_mutations_bool_damaging.csv
wget https://ndownloader.figshare.com/files/34008440 -O CCLE_mutations_bool_damaging.csv
# CCLE_mutations_bool_nonconserving.csv
wget https://ndownloader.figshare.com/files/34008449 -O CCLE_mutations_bool_nonconserving.csv
# CCLE_mutations_bool_otherconserving.csv
wget https://ndownloader.figshare.com/files/34008452 -O CCLE_mutations_bool_otherconserving.csv
# sample_info.csv
wget https://ndownloader.figshare.com/files/34008503 -O sample_info.csv
# Achilles_metadata.csv
wget https://ndownloader.figshare.com/files/34008386 -O Achilles_metadata.csv
