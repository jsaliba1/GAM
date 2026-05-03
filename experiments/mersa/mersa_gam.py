import pandas as pd
import numpy as np
from collections import Counter
import scipy.stats as stats


mut = pd.read_csv('snp_calls.csv' , sep = "\t")
mut = mut.drop(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
mut = mut.set_index('ID').T
#mut = pd.read_csv('gpa.csv', index_col = ['gene']).T #gene    
    
drugs = pd.read_csv('RSdata13.csv', index_col = ['ID'])

mut = mut[mut.index.isin(list(drugs.index))]

drug_names = ['GEN','PEN','MET','FUS','TEI','VAN','Ery','CLI','LIN','CIP','RIF','TET','TMP']
n_drugs = len(drug_names)

y = []
for i in range(len(drugs)):
    x = drugs.iloc[i,:] #taking resitence profile for each strain
    y.append(str(list(x))) #profile added to list
   
counts = Counter(y) #count the number of time each profile occurs 
r = pd.Series(counts)
r = r[r>=1]
indexs = r.index.sort_values() #list of unique resistence profile and sort


mut1 = mut.replace([0,1],[1,2])
drugs1 = drugs.replace([0,1],[1,-1])

cs = []

for i in range(n_drugs):
    cs_drug = []
    z = mut1.T * drugs1[drug_names[i]]
    bon = 0.05/len(z)
    for j in range(len(z)):
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
        #if p_value < bon:
        cs_drug.append(N.name)
    cs.append(cs_drug)


pval = pd.DataFrame()

for i in range(len(indexs)): #filters for each drug and places them into approprate group 
    cs_g = []
    S = indexs[i]
    for j in range(n_drugs):
        if j==0:
            fix = drugs[drugs[drug_names[j]] == int(S[1])]
            if int(S[1]) == 1:
                cs_g.append(cs[j])
        else:
            fix = fix[fix[drug_names[j]] == int(S[1+j*3])]
            if int(S[1+j*3]) == 1:
                cs_g.append(cs[j])
    
    cs_g = {k for lst in cs_g for k in lst}
    listFix = list(fix.index)
    mutS = mut[mut.index.isin(listFix)]
    mutS = mutS.T[mutS.columns.isin(cs_g)].T
    Ms = mutS.sum()
       
    pval = pd.concat([pval,Ms],axis=1).replace(np.nan, 0)

pval.columns = list(indexs)

MD_pval = pd.DataFrame(np.nan, index = pval.index, columns = drug_names)

for i in range(len(pval)):
    A = [0] * n_drugs
    B = [0] * n_drugs
    C = [0] * n_drugs
    D = [0] * n_drugs
    mut_pval = pval.iloc[i] 

    for j in range(len(mut_pval)):
         drug_set = mut_pval.index[j]
         
         if mut_pval[j] > 0: 
            for k in range(1, len(drug_set), 3):
                if drug_set[k] == '0': 
                    C[int((k-1)/3)] += 1 
                elif drug_set[k] == '1': 
                    A[int((k-1)/3)] += 1 
                  
         else:
            for k in range(1, len(drug_set), 3):
                if drug_set[k] == '0': 
                    D[int((k-1)/3)] += 1 
                elif drug_set[k] == '1': 
                    B[int((k-1)/3)] += 1 
            
    for q in range(n_drugs):
        data = [ [A[q],B[q]], [C[q],D[q]] ]
        odd_ratio, p_value = stats.fisher_exact(data)
        
        if odd_ratio > 1:
            MD_pval.iloc[i,q] = p_value
        else:
            MD_pval.iloc[i,q] = 1
        
   
MDbon = 0.05/len(MD_pval)
MD_pval = MD_pval[MD_pval <= MDbon]
MD_pval.dropna(axis=0, how ='all', inplace=True)
MD_pval = MD_pval.set_axis(drug_names,axis=1)

        
    
  


