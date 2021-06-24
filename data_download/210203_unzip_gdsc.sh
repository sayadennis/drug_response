#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="unzip_gdsc"
#SBATCH --output=drug_response/210203_unzip_gdsc.out

unzip /projects/b1131/saya/gdsc/cell_lines_genomic/cell_line_wes_20180620.zip -d /projects/b1131/saya/gdsc/cell_lines_genomic
unzip /projects/b1131/saya/gdsc/cell_lines_genomic/cnv_20191101.zip -d /projects/b1131/saya/gdsc/cell_lines_genomic
unzip /projects/b1131/saya/gdsc/cell_lines_genomic/mutations_20191101.zip -d /projects/b1131/saya/gdsc/cell_lines_genomic
unzip /projects/b1131/saya/gdsc/cell_lines_genomic/rnaseq_20191101.zip -d /projects/b1131/saya/gdsc/cell_lines_genomic
