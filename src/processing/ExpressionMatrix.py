import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime

def write_mx(fin, fout):
    """
    This function allows you to create a RNA expression matrix based on the CCLE expression file.
    """
    exp = pd.read_csv(fin, sep="\t", index_col=1, skiprows=2) # rows = genes // columns = tissue or sample
    exp = exp.T.drop("Name", axis=0) # dropping "ENSG<num*>" gene names
    exp.to_csv(fout, header=True, index=True)
    return
