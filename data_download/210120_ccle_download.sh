#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ccle_dl"
#SBATCH --output=drug_response/210120_ccle_download.out

curl https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_maf_20180718.txt -o /projects/b1131/saya/ccle/CCLE_DepMap_18q3_maf_20180718.txt
curl https://data.broadinstitute.org/ccle/CCLE_translocations_SvABA_20181221.xlsx -o /projects/b1131/saya/ccle/CCLE_translocations_SvABA_20181221.xlsx
curl https://data.broadinstitute.org/ccle/gencode.v19.genes.v7_model.patched_contigs.gtf.gz -o /projects/b1131/saya/ccle/gencode.v19.genes.v7_model.patched_contigs.gtf.gz
curl https://data.broadinstitute.org/ccle/CCLE_ABSOLUTE_combined_20181227.xlsx -o /projects/b1131/saya/ccle/CCLE_ABSOLUTE_combined_20181227.xlsx
curl https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt -o /projects/b1131/saya/ccle/Cell_lines_annotations_20181226.txt


curl https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz -o /projects/b1131/saya/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz

curl https://data.broadinstitute.org/ccle/CCLE_RRBS_TSS1kb_20181022.txt.gz -o /projects/b1131/saya/ccle/CCLE_RRBS_TSS1kb_20181022.txt.gz
curl https://data.broadinstitute.org/ccle/CCLE_RRBS_cgi_CpG_clusters_20181119.txt.gz -o /projects/b1131/saya/ccle/CCLE_RRBS_cgi_CpG_clusters_20181119.txt.gz
curl https://data.broadinstitute.org/ccle/CCLE_ABSOLUTE_combined_20181227.xlsx -o /projects/b1131/saya/ccle/CCLE_ABSOLUTE_combined_20181227.xlsx
