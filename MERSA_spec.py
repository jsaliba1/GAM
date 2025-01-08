import pandas as pd
import numpy as np
import scipy.stats as stats
import math

drugs = pd.read_csv('RSdata13.csv', index_col = ['ID'])
mut = pd.read_csv('snp_calls.csv' , sep = "\t")

drug_names = ['GEN','PEN','MET','FUS','TEI','VAN','Ery','CLI','LIN','CIP','RIF','TET','TMP']
n_drugs = len(drug_names)

start = int(input('start: '))
stop = int(input('stop: '))
drug = str(input('drug: '))
q = list(range(start,stop))

mut2 = mut[mut.POS.isin(q)]
mut2 = mut2.drop(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
mut2 = mut2.set_index('ID').T
mut2 = mut2[mut2.index.isin(list(drugs.index))]

cnt = mut2.T.sum()
bad = []
for i in range(len(cnt)):
    if cnt[i] > 1:
        bad.append(cnt.index[i])

mut2 = mut2.drop(bad, axis=0).sort_index()
drugs2 = drugs.drop(bad, axis=0).sort_index()

mut1 = mut2.replace([0,1],[1,2])
drugs1 = drugs2.replace([0,1],[1,-1])

cs_drug = []
z = mut1.T * drugs1[drug]
bon = 0.05/len(z)
for j in range(0,len(z)): 
    N =  z.iloc[j].value_counts()
   
    try:
        a = N[N.index==-2].iloc[0]
    except:
        a = 0
    try:
        b = N[N.index==2].iloc[0]
    except:
        b = 0
    try:
        c = N[N.index==-1].iloc[0]
    except:
        c = 0
    try:
        d = N[N.index==1].iloc[0]
    except:
        d = 0
    
    odd_ratio, p_value = stats.fisher_exact([[a,b],[c,d]])
    
    try: 
        sensitivity = a/(a+c)
    except:
        sensitivity = np.nan
    try:
        specificity = d/(d+b)
    except: 
        specificity = np.nan
    
    n = a+b
    try:
        PPV = a/n
        se = 1.96 * (math.sqrt( (PPV * (1-PPV)) / n ))
        LCI = PPV - se
        UCI = PPV + se
    except:
        PPV = np.nan
        LCI = np.nan
        UCI = np.nan
    
    cs_drug.append((a,b,c,d,sensitivity,specificity,PPV,LCI,UCI,odd_ratio,p_value,n,N.name))
    
csd = pd.DataFrame(cs_drug)