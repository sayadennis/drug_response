import sys

sys.path.append("drug_response/src/processing")
import ExpressionMatrix

ExpressionMatrix.write_mx(
    fin = "/projects/b1131/saya/drug_response/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct", 
    fout = "/projects/b1131/saya/drug_response/ccle/expression_mx.csv"
)
