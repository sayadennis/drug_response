import sys

sys.path.append("drug_response/src/processing")
import VCMatrix

VCMatrix.write_mx(
    fin = "/projects/b1131/saya/drug_response/ccle/CCLE_DepMap_18q3_maf_20180718.txt", 
    fout = "/projects/b1131/saya/drug_response/ccle/vc_mx_excludesilent.csv"
)
