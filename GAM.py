import pandas as pd
import numpy as np
from collections import Counter
import scipy.stats as stats


mut = pd.read_csv('MUTATIONS.csv')
mut = mut[mut.IS_SYNONYMOUS == False ]

ENA_drugs = pd.read_csv('ENA_drug.csv', index_col = ['UNIQUEID']) 

drug_names = ['AMIKACIN', 'BEDAQUILINE', 'CLOFAZIMINE', 'DELAMANID','ETHAMBUTOL', 'ETHIONAMIDE', 'ISONIAZID',\
              'KANAMYCIN', 'LEVOFLOXACIN', 'LINEZOLID', 'MOXIFLOXACIN', 'RIFAMPICIN', 'RIFABUTIN']
n_drugs = len(drug_names) #number of drugs 

for P in [0]: #drop rate
     for remove_n in [0]: #how much do you want removed 
            for seed in range(1): #how many replicates 
                np.random.seed(seed)
                bugs = pd.read_csv('bugs.csv')
                
                drop_indices = np.random.choice(bugs.index, remove_n, replace=False)
                bugs_subset = bugs.drop(drop_indices)  
                
                bugNames = list(bugs_subset['ID'])
                drugs = ENA_drugs[ENA_drugs.index.isin(bugNames)] #only strains in file input
                matrix = np.random.choice([0, 1], size=(len(drugs), 13), p=[P, 1-P])
                drugs = drugs*matrix
                drugs = drugs.replace('','a') 
               
                y = []
                for i in range(len(drugs)):
                    x = drugs.iloc[i,:] #taking resitence profile for each strain
                    y.append(str(list(x))) #profile added to list
                   
                counts = Counter(y) #count the number of time each profile occurs 
                r = pd.Series(counts)
                
                r = r.sort_values(ascending = False) #sort largest to smallest
                indexs = r.index #list of unique resistence profile 
                
                pval = {}
                for i in range(len(indexs)): #filters for each drug and places them into approprate group 
                    S = indexs[i]
                    fix = drugs[drugs.AMIKACIN == S[2] ]
                    fix = fix[fix.BEDAQUILINE == S[7] ]
                    fix = fix[fix.CLOFAZIMINE == S[12] ]
                    fix = fix[fix.DELAMANID == S[17] ]
                    fix = fix[fix.ETHAMBUTOL == S[22] ]
                    fix = fix[fix.ETHIONAMIDE == S[27] ]
                    fix = fix[fix.ISONIAZID == S[32] ]
                    fix = fix[fix.KANAMYCIN == S[37] ]
                    fix = fix[fix.LEVOFLOXACIN == S[42] ]
                    fix = fix[fix.LINEZOLID == S[47] ]
                    fix = fix[fix.MOXIFLOXACIN == S[52] ]
                    fix = fix[fix.RIFAMPICIN == S[57] ]
                    fix = fix[fix.RIFABUTIN == S[62] ]
                    
                    
                    if len(fix) > 1:
                        listFix = list(fix.index)
                        
                        if S == "["+"'S', "*(n_drugs-1)+"'S']": #control
                            lenS = len(fix)
                            mutS = mut[mut.UNIQUEID.isin(listFix)]  
                            
                            ENA_snps = {}
                            for k in range(lenS):
                                folderK = listFix[k]
                                snps = mutS[mutS.UNIQUEID.isin([folderK])]
                                ENA_snps[folderK] = snps.groupby(['MUTATION','GENE'])['MUTATION'].count()
                        
                            ENA_snps = pd.DataFrame(ENA_snps)
                            SMs = ENA_snps.sum(axis=1) # numebr of time mutation occurs in control 
                        
                        else: #test groups
                            lenX = len(fix)
                            mutX = mut[mut.UNIQUEID.isin(listFix)]
                            
                            ENA_snps = {}
                            for k in range(lenX):
                                folderK = listFix[k]
                                snps = mutX[mutX.UNIQUEID.isin([folderK])]
                                ENA_snps[folderK] = snps.groupby(['MUTATION','GENE'])['MUTATION'].count()
                            
                            ENA_snps = pd.DataFrame(ENA_snps)
                            RMs = ENA_snps.sum(axis=1) # numebr of time mutation occurs in test group
                            
                            
                            df = pd.concat([RMs,SMs],axis=1).replace(np.nan, 0)
                            df = df[df[0] > 0]
                            for j in range(0,len(df)):
                                data = [[ df.iloc[j,0], df.iloc[j,1]] , [ lenX-df.iloc[j,0], lenS-df.iloc[j,1] ]]
                                odd_ratio, p_value = stats.fisher_exact(data)
                                df.iloc[j,0] = odd_ratio
                                df.iloc[j,1] = p_value
                         
                                
                            if len(df) > 0:
                                bon = 0.05/len(df)
                                df = df[df.iloc[:,1] <= bon]  #filter out all non-significant mutations
                                df = df[df.iloc[:,0] > 1] #filter out mutations that are correlated negativly
                                
                                pval[S] = df.iloc[:,1]
                            
                                
                if len(pval) != 0:
                    pval = pd.DataFrame(pval).replace(np.nan, 1)
                    MD_pval = pval.iloc[:,0:n_drugs] #to be writen over
                    
                    for i in range(len(pval)):
                        A = [0] * n_drugs
                        B = [0] * n_drugs
                        C = [0] * n_drugs
                        D = [0] * n_drugs
                        mut_pval = pval.iloc[i] 
                    
                        for j in range(len(mut_pval)):
                             drug_set = mut_pval.index[j]
                             
                             if mut_pval[j] < 0.05: 
                                for k in range(2, len(drug_set), 5):
                                    if drug_set[k] == 'S': 
                                        C[int((k-2)/5)] += 1
                                    elif drug_set[k] == 'R': 
                                        A[int((k-2)/5)] += 1
                                      
                             else:
                                for k in range(2, len(drug_set), 5):
                                    if drug_set[k] == 'S': 
                                        D[int((k-2)/5)] += 1   
                                    elif drug_set[k] == 'R': 
                                        B[int((k-2)/5)] += 1
                                
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
                    MD_pval.to_csv(str(seed)+str(P*100)+str(len(bugs)-remove_n)+'.csv', index = True)
                
                else:
                    print('fail:'+str(seed)+str(P*100)+str(len(bugs)-remove_n))
       
             
