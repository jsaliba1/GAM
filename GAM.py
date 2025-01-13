import pandas as pd
import numpy as np
from collections import Counter
import scipy.stats as stats

class GAM:
    def __init__(self, mutation_file, drug_file, bug_file):
        self.mut = pd.read_csv(mutation_file,low_memory=False)
        self.mut = self.mut[self.mut.IS_SYNONYMOUS == False]
        self.ENA_drugs = pd.read_csv(drug_file, index_col=['UNIQUEID'])
        self.bugs = pd.read_csv(bug_file)
        self.drug_names = self.ENA_drugs.columns.tolist()
        self.n_drugs = len(self.drug_names)
        self.pval = {}

    def group_bugs(self, P=0, remove_n=0, seed=0):
        np.random.seed(seed)
        drop_indices = np.random.choice(self.bugs.index, remove_n, replace=False)
        bugs_subset = self.bugs.drop(drop_indices)  
        
        bug_names = list(bugs_subset['ID'])
        drugs = self.ENA_drugs[self.ENA_drugs.index.isin(bug_names)]
        
        matrix = np.random.choice([0, 1], size=(len(drugs), 13), p=[P, 1-P])
        drugs = drugs * matrix
        self.drugs = drugs.replace('', 'a') 

        y = [str(list(self.drugs.iloc[i,:])) for i in range(len(self.drugs))]
        counts = Counter(y)
        r = pd.Series(counts)
        r = r.sort_values(ascending=False)
        self.indexs = r.index
        return 
        
    def analyze_resistance(self):
        for i in range(len(self.indexs)):
            S = self.indexs[i]
            fix = self.drugs
            for j in range(self.n_drugs):
                fix = fix[fix[self.drug_names[j]] == S[j * 5 + 2]]

            if len(fix) > 1:
                listFix = list(fix.index) #list of bug names in a group
                if S == "[" + "'S', " * (self.n_drugs - 1) + "'S']":  # control
                    lenS = len(fix)
                    mutS = self.mut[self.mut.UNIQUEID.isin(listFix)]
                    ENA_snps = {}
                    for k in range(lenS):
                        folderK = listFix[k]
                        snps = mutS[mutS.UNIQUEID.isin([folderK])]
                        ENA_snps[folderK] = snps.groupby(['MUTATION', 'GENE'])['MUTATION'].count()
                    ENA_snps = pd.DataFrame(ENA_snps)
                    SMs = ENA_snps.sum(axis=1) #how often a mut is in a group
                else:  # test groups
                    lenX = len(fix)
                    mutX = self.mut[self.mut.UNIQUEID.isin(listFix)]
                    ENA_snps = {}
                    for k in range(lenX):
                        folderK = listFix[k]
                        snps = mutX[mutX.UNIQUEID.isin([folderK])]
                        ENA_snps[folderK] = snps.groupby(['MUTATION', 'GENE'])['MUTATION'].count()
                    ENA_snps = pd.DataFrame(ENA_snps)
                    RMs = ENA_snps.sum(axis=1)

                    df = pd.concat([RMs, SMs], axis=1).replace(np.nan, 0)
                    df = df[df[0] > 0]
                    for j in range(0, len(df)):
                        data = [[df.iloc[j, 0], df.iloc[j, 1]], [lenX - df.iloc[j, 0], lenS - df.iloc[j, 1]]]
                        odd_ratio, p_value = stats.fisher_exact(data)
                        df.iloc[j, 0] = odd_ratio
                        df.iloc[j, 1] = p_value

                    if len(df) > 0:
                        bon = 0.05 / len(df)
                        df = df[df.iloc[:, 1] <= bon]  # filter out all non-significant mutations
                        df = df[df.iloc[:, 0] > 1]  # filter out mutations that are correlated negatively
                        self.pval[S] = df.iloc[:, 1]
        return

    def score_resistance(self):
        if len(self.pval) != 0:
            pval = pd.DataFrame(self.pval).replace(np.nan, 1)
            MD_pval = pval.iloc[:, 0:self.n_drugs]
            for i in range(len(pval)):
                A = [0] * self.n_drugs
                B = [0] * self.n_drugs
                C = [0] * self.n_drugs
                D = [0] * self.n_drugs
                mut_pval = pval.iloc[i]
                for j in range(len(mut_pval)):
                    drug_set = mut_pval.index[j]
                    if mut_pval[j] < 0.05:
                        for k in range(2, len(drug_set), 5):
                            if drug_set[k] == 'S':
                                C[int((k - 2) / 5)] += 1
                            elif drug_set[k] == 'R':
                                A[int((k - 2) / 5)] += 1
                    else:
                        for k in range(2, len(drug_set), 5):
                            if drug_set[k] == 'S':
                                D[int((k - 2) / 5)] += 1
                            elif drug_set[k] == 'R':
                                B[int((k - 2) / 5)] += 1

                for q in range(self.n_drugs):
                    data = [[A[q], B[q]], [C[q], D[q]]]
                    odd_ratio, p_value = stats.fisher_exact(data)
                    if odd_ratio > 1:
                        MD_pval.iloc[i, q] = p_value
                    else:
                        MD_pval.iloc[i, q] = 1

            MDbon = 0.05 / len(MD_pval)
            MD_pval = MD_pval[MD_pval <= MDbon]
            MD_pval.dropna(axis=0, how='all', inplace=True)
            MD_pval = MD_pval.set_axis(self.drug_names, axis=1)
        return MD_pval

if __name__ == "__main__":
    # Example usage:
    analysis = GAM('MUTATIONS.csv', 'ENA_drug.csv', 'bugs.csv')
    analysis.group_bugs(0, 0, 0)  # Example parameters for analysis
    analysis.analyze_resistance()
    MD_pval = analysis.score_resistance()