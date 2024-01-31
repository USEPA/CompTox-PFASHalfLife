# -*- coding: utf-8 -*-
"""

@author: rtornero-velez

----------------------------------------------------------------------------------------
This script generates fingerprints for the chemicals
    
----------------------------------------------------------------------------------------
"""
import pandas as pd
import pubchempy as pcp


from rdkit import Chem
from rdkit.Chem import AllChem


#path = 'C:/Users/ppradeep/Desktop/PFASsimilarity/'
#path = 'C:/Users/rtorne02/OneDrive - Environmental Protection Agency (EPA)/HTTK/PFAS_HL_QSAR/'
path = 'C:/Users/rtorne02/OneDrive - Environmental Protection Agency (EPA)/PFAS Project Share/PFAS_Similarity/'
      
#%%
###########################################################################
## 1. Load chemicals candidates to compare with endogenous chemicals
###########################################################################

#candidates = pd.read_excel(path+'input/final_list_OPPT_QSAR-ready_smi_test.xlsx', index_col='CASRN') 
candidates = pd.read_excel(path+'input/Patlewicz500_112923_QSAR-ready_smi.xlsx', index_col='CASRN')
#candidates = pd.read_excel(path+'input/Patlewicz_112923_QSAR-ready_smi.xlsx', index_col='CASRN')


#%%
###########################################################################
## 2. Generate PubChem and Morgan fingerprints for the test list
###########################################################################
import time
start = time.time()
print("hello")

for i in candidates.index:
    #smile = dawsonAll.loc[i,'QSAR_READY_SMILES_H']
    smile = candidates.loc[i,'QSAR_READY_SMILES']
    # PubChem FP
    try:
        cmpd = pcp.get_compounds(smile, 'smiles')[0]
        #pfp = DataStructs.CreateFromBitString(cmpd.cactvs_fingerprint)
        pfp = cmpd.cactvs_fingerprint
        candidates.loc[i,'PubChemFP'] = pfp
    except:
        print("PubChem FP not calculated for %s" %i)
    # Morgan FP
    try:
        #smile_m=smile.replace("(",'').replace(")",'').replace("=",'').replace("H",'')
        #m = Chem.MolFromSmiles(smile_m)        
        molobj= Chem.MolFromSmiles(smile)
        morFP = AllChem.GetMorganFingerprintAsBitVect(molobj, radius=2, nBits=1024)  # Adjust radius and nBits as needed
        bit_string = morFP.ToBitString()
        candidates.loc[i,'MorganFP'] = bit_string
    except:
        print("Morgan FP not calculated for %s" %i)
  
#%%
        
end = time.time()
print("end of calcs")
print(end - start)

###########################################################################
## 4. Save the finger print file to disk
###########################################################################  

#similarity.to_csv(path+'output\EndoSubset-EPA_PFAS_DSSTox-Similarity.csv', index_label = 'Row Number')
#candidates.to_csv(path+'output/final_list_OPPT_QSAR-ready_smi_500_wFP.csv', index_label = 'CASRN')
#candidates.to_csv(path+'input/Wambaugh_110123_QSAR-ready_smi_wFP.csv', index_label = 'CASRN')
#candidates.to_csv(path+'input/TEST_112823_QSAR-ready_smi_wFP.csv', index_label = 'DTXSID')
candidates.to_csv(path+'input/Patlewicz500_112923_QSAR-ready_smi_wFP_test.csv', index_label = 'CASRN')
#candidates.to_csv(path+'input/Patlewicz_112923_QSAR-ready_smi_wFP.csv', index_label = 'CASRN')