"""
Rdkit only works in certain versions of Python. 
Although other codes are written for the "cancerdrugenv" conda environment, this scripts needs to be run under the "rdkitenv" for this reason. 
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

din='/projects/b1131/saya/drug_response/ctrp'
dout='/projects/b1131/saya/drug_response/processed_data'

## Load compound data
cpd = pd.read_csv(f'{din}/v20.meta.per_compound.txt', sep='\t')

## Iteratively convert SMILES to Morgan fingerprints of different lengths 

morgan128 = pd.DataFrame(None, index=cpd.master_cpd_id, columns=np.arange(128))
morgan64 = pd.DataFrame(None, index=cpd.master_cpd_id, columns=np.arange(64))

# for 128-bit vector
for i in cpd.index:
    cpd_id = cpd.loc[i,'master_cpd_id']
    mol = Chem.MolFromSmiles(cpd.loc[i,'cpd_smiles'])
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=2, nBits=128)
    vec = np.array(fp)
    morgan128.loc[cpd_id] = vec


# for 64-bit vector 
for i in cpd.index:
    cpd_id = cpd.loc[i,'master_cpd_id']
    mol = Chem.MolFromSmiles(cpd.loc[i,'cpd_smiles'])
    fp = Chem.AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=2, nBits=64)
    vec = np.array(fp)
    morgan64.loc[cpd_id] = vec

## Save encoded fingerprints 
morgan64.to_csv(f'{dout}/morgan_fingerprint_64.csv', index=True, header=True)
morgan128.to_csv(f'{dout}/morgan_fingerprint_128.csv', index=True, header=True)
