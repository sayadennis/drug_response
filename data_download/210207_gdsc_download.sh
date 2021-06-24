#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gdsc_dl"
#SBATCH --output=drug_response/210207_gdsc_download.out

curl https://cog.sanger.ac.uk/cmp/download/cell_line_wes_20180620.zip -o /projects/b1131/saya/gdsc/cell_lines_genomic/cell_line_wes_20180620.zip
curl https://cog.sanger.ac.uk/cmp/download/mutations_20191101.zip -o /projects/b1131/saya/gdsc/cell_lines_genomic/mutations_20191101.zip
curl https://cog.sanger.ac.uk/cmp/download/rnaseq_20191101.zip -o /projects/b1131/saya/gdsc/cell_lines_genomic/rnaseq_20191101.zip
curl https://cog.sanger.ac.uk/cmp/download/cnv_20191101.zip -o /projects/b1131/saya/gdsc/cell_lines_genomic/cnv_20191101.zip
