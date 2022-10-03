# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:45:41 2020

@author: ppradeep

----------------------------------------------------------------------------------------
This script runs different machine learning algorithms on the dataset to find the best 
model as validated using a 5-fold cross validation on a 80% training set and 20% test 
set.
Input:
    
----------------------------------------------------------------------------------------
"""
import pandas as pd
import pubchempy as pcp

from rdkit import DataStructs, Chem
from rdkit.Chem import AllChem

path = 'C:/Users/ppradeep/Desktop/PFASsimilarity/'

      
#%%
###########################################################################
## 1. Load the test chemicals
###########################################################################
dawsonAll = pd.read_excel(path+'data/Dawson_SimEndChems_frommodel_DED05082020.xlsx', index_col='CASRN') 

###########################################################################
## 2. Replace 'F' by 'H' in the test chemicals SMILES
###########################################################################
dawsonAll['QSAR_READY_SMILES_H'] = dawsonAll['QSAR_READY_SMILES'].str.replace('F', 'H') 

#%%
###########################################################################
## 3. Generate PubChem and Morgan fingerprints for the test list
###########################################################################
for cas in dawsonAll.index:
    smile = dawsonAll.loc[cas,'QSAR_READY_SMILES_H']
    # PubChem FP
    try:
        cmpd = pcp.get_compounds(smile, 'smiles')[0]
        #pfp = DataStructs.CreateFromBitString(cmpd.cactvs_fingerprint)
        pfp = cmpd.cactvs_fingerprint
        dawsonAll.loc[cas,'PubChemFP'] = pfp
    except:
        print("PubChem FP not calculated for %s" %cas)
    # Morgan FP
    try:
        smile_m=smile.replace("(",'').replace(")",'').replace("=",'').replace("H",'')
        m = Chem.MolFromSmiles(smile_m)
        #m = Chem.MolFromSmiles(smile)
        dawsonAll.loc[cas,'MorganFP'] = m
    except:
        print("Morgan FP not calculated for %s" %cas)
  
#%%
###########################################################################
## 4. Load the EPA ALL PFAS List with PubChem and Morgan fingerprints.
## These were calculated earlier in the same way as calculated for the
## test chemicals
###########################################################################
endoMetaAll = pd.read_csv(path+'data/EPA-PFAS_DSSTox_wFPs.csv', index_col='CASRN') 

#%%
###########################################################################
## 5. Calculate Tanimoto similarity (using PubChem and Morgan fingerprints) 
## between the EPA Endo list and the test list
###########################################################################      
similarity = pd.DataFrame(columns=['CASRN', 'DTXSID', 'ENDO CASRN (EPA-PFAS)',\
                                   'ENDO DTXSID (EPA-PFAS)', 'Tanimoto Score (PubChem)', 'Tanimoto Score (Morgan)']) 
i = 0
for cas1 in dawsonAll.index:
    try:
        pfp1 = DataStructs.CreateFromBitString(dawsonAll.loc[cas1, 'PubChemFP'])
    except:
        pass          
    mfp1 = AllChem.GetMorganFingerprintAsBitVect(dawsonAll.loc[cas1, 'MorganFP'],2,nBits=1024)
    for cas2 in endoMetaAll.index:
        i=i+1            
        similarity.loc[i,'CASRN'] = cas1
        similarity.loc[i,'DTXSID'] = dawsonAll.loc[cas1, 'DTXSID']
        similarity.loc[i,'ENDO CASRN (EPA-PFAS)'] = cas2
        similarity.loc[i, 'ENDO DTXSID (EPA-PFAS)'] = endoMetaAll.loc[cas2, 'DTXSID']
        try:
            pfp2 = DataStructs.CreateFromBitString(endoMetaAll.loc[cas2, 'PubChemFP'])
            pscore = DataStructs.FingerprintSimilarity(pfp1, pfp2)
            similarity.loc[i,'Tanimoto Score (PubChem)'] = pscore
        except:
            pass#print("PubChem FP Similarity not calculated for %s-%s" %(cas1,cas2))
        try:
            mfp2 = AllChem.GetMorganFingerprintAsBitVect(endoMetaAll.loc[cas2, 'MorganFP'],2,nBits=1024)
            mscore = DataStructs.FingerprintSimilarity(mfp1, mfp2)
            similarity.loc[i,'Tanimoto Score (Morgan)'] = mscore
        except:
            pass#print("Morgan FP Similarity not calculated for %s-%s" %(cas1,cas2))

#%%
###########################################################################
## 6. Save the similarity file to disk
###########################################################################  
similarity.to_csv(path+'output\EndoSubset-EPA_PFAS_DSSTox-Similarity.csv', index_label = 'Row Number')
