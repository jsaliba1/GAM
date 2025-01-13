import pandas as pd
import numpy as np
import collections
import scipy.stats as stats
import math
import os

class FineAnalysis:
    def __init__(self, mutation_file, drug_file, input_file, reg_folder_path, run_num):
        self.mut = pd.read_csv(mutation_file)
        self.ENA_drugs = pd.read_csv(drug_file, index_col=['UNIQUEID']).replace(['S', 'R', 'I'], [1, -1, 0])
        self.folder_names = list(pd.read_csv(input_file)['ID'])
        self.reg_folder = reg_folder_path
        self.un2 = ['rrs', 'rpoB', 'katG', 'embB', 'fabG1', 'gyrA']
        self.run = run_num
        self.ssp = {} 

    def analyze_snps(self):
        ENA_snps = {}
        for k in range(len(self.folder_names)):
            folderK = self.folder_names[k]
            snps = self.mut[self.mut.UNIQUEID.isin([folderK])]
            snps = snps[snps.GENE.isin(self.un2)] #filter out non gene targets 
            ENA_snps[folderK] = snps.groupby(['MUTATION', 'GENE'])['MUTATION'].count()
            if self.run != 0:
                un1 = [el for el, count in collections.Counter(list(snps.GENE)).items() if count > 1] 
                if un1 == self.un2[5]:
                    un1 = [el for el, count in collections.Counter(list(snps.GENE)).items() if count > 2] 
                snps = snps[snps.GENE.isin(un1)]
                x = snps.groupby(['MUTATION', 'GENE'])['MUTATION'].count()
                for i in x.index:
                    ENA_snps[folderK][i] = 0
        
        ENA_snps = pd.DataFrame(ENA_snps).replace([1, np.nan], [2, 1])
        self.ENA_snps = ENA_snps.reindex(sorted(ENA_snps.columns), axis=1)
        return
    
    def score_snps(self):
        names = ['AMIKACIN', 'ETHAMBUTOL', 'ETHIONAMIDE', 'ISONIAZID', 'KANAMYCIN', 'LEVOFLOXACIN', 'MOXIFLOXACIN', 'RIFAMPICIN']
        for i in range(len(names)):
            meta = {}
            z = self.ENA_snps * self.ENA_drugs[names[i]]
            reg = pd.read_csv(os.path.join(self.reg_folder, names[i]+'.csv'), index_col=['MUTATION', 'GENE'])
            z = z[z.index.isin(reg.index)]

            for j in range(len(z)):
                N = z.iloc[j].value_counts()
                try:
                    a = N[N.index == -2].iloc[0]
                except:
                    a = 0
                try:
                    b = N[N.index == 2].iloc[0]
                except:
                    b = 0
                try:
                    c = N[N.index == -1].iloc[0]
                except:
                    c = 0
                try:
                    d = N[N.index == 1].iloc[0]
                except:
                    d = 0

                odd_ratio, p_value = stats.fisher_exact([[a, b], [c, d]])

                try:
                    sensitivity = a / (a + c)
                except:
                    sensitivity = np.nan
                try:
                    specificity = d / (d + b)
                except:
                    specificity = np.nan
                n = a + b
                try:
                    PPV = a / n
                    se = 1.96 * (math.sqrt((PPV * (1 - PPV)) / n))
                    LCI = PPV - se
                    UCI = PPV + se
                except:
                    PPV = np.nan
                    LCI = np.nan
                    UCI = np.nan

                meta[z.index[j]] = {'a': a, 'b': b, 'c': c, 'd': d, 'sensitivity': sensitivity, 'specificity': specificity, 'PPV': PPV,\
                                    'LCI': LCI, 'UCI': UCI,'OD': odd_ratio, 'pval': p_value, 'n1': n}

            meta = pd.DataFrame(meta).T
            meta['Q'] = meta.n1 / sum(meta.n1) * 100  # freqency
            meta = meta[meta.n1 >= 5]
            meta = meta[meta.LCI >= .25]
            meta = meta[meta.OD >= 1]
            if len(meta) > 0:
                bon = .05 / len(meta)
                meta = meta[meta.pval <= bon]
            else:
                print(names[i])

            v = reg[reg.index.isin(meta.index)]
            a = 0
            b = 0
            c = 0
            d = 0
            e = 0
            for j in v:
                if sum(v[j]) == len(v[j]):
                    d += 1
                elif sum(v[j]) == -len(v[j]):
                    c += 1
                elif sum(v[j]) > len(v[j]):
                    b += 1
                elif sum(v[j]) < -len(v[j]):
                    a += 1
                else:
                    e += 1
                    
            n = a + c
            try:
                sensitivity = a / n
                se = 1.96 * (math.sqrt((sensitivity * (1 - sensitivity)) / n))
                SLCI = sensitivity - se
                SUCI = sensitivity + se
            except:
                sensitivity = np.nan
                SLCI = np.nan
                SUCI = np.nan
            n = b + d
            try:
                specificity = d / n
                se = 1.96 * (math.sqrt((specificity * (1 - specificity)) / n))
                fLCI = specificity - se
                fUCI = specificity + se
            except:
                specificity = np.nan
                fLCI = np.nan
                fUCI = np.nan
            n = a + b
            try:
                PPV = a / n
                se = 1.96 * (math.sqrt((PPV * (1 - PPV)) / n))
                LCI = PPV - se
                UCI = PPV + se
            except:
                PPV = np.nan
                LCI = np.nan
                UCI = np.nan

            self.ssp[names[i]] = {'a': a, 'b': b, 'c': c, 'd': d, 'sensitivity': sensitivity, 'sen_LCI': SLCI, 'sen_UCI': SUCI,\
                                  'specificity': specificity,'spec_LCI': fLCI, 'spec_UCI': fUCI, 'PPV': PPV, 'LCI': LCI, 'UCI': UCI}
        return self.ssp

if __name__ == "__main__":
    # Example usage:
    analysis = FineAnalysis('MUTATIONS.csv', 'ENA_drug.csv', 'input_max.csv', 'Reg/')
    analysis.analyze_snps()  # Example parameters for analysis
    ssp = analysis.score_snps()