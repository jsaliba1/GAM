import pandas as pd
import numpy as np
from sklearn.utils import resample
from pysnptools.snpreader import SnpData
from fastlmm.association import single_snp

mut = pd.read_csv('MUTATIONS.csv',low_memory=False).iloc[:,:3]

bugs = pd.read_csv('bugs.csv')
bugNames1 = list(bugs['ID'])

bugNames = list(bugs['ID'])[:int(len(bugs)/2)]
mutation_data = mut[mut.UNIQUEID.isin(bugNames)]

mutation_data['combo'] = mutation_data['GENE'] +"_"+ mutation_data['MUTATION']
mutation_data['value'] = 1
matrix1 = mutation_data.pivot_table(index='combo', columns='UNIQUEID', values='value')

bugNames = list(bugs['ID'])[int(len(bugs)/2):]
mutation_data = mut[mut.UNIQUEID.isin(bugNames)]

mutation_data['combo'] = mutation_data['GENE'] +"_"+ mutation_data['MUTATION']
mutation_data['value'] = 1
matrix2 = mutation_data.pivot_table(index='combo', columns='UNIQUEID', values='value')

mutation_matrix = matrix1.join(matrix2)
mutation_matrix[mutation_matrix > 0] = 1
mutation_matrix = mutation_matrix.replace(np.nan,0).astype('int8')  # Convert to boolean for Jaccard

mutation_array = mutation_matrix.values

# Compute pairwise intersections using matrix multiplication
intersection_matrix = np.dot(mutation_array.T, mutation_array)

# Compute pairwise unions using the inclusion-exclusion principle
row_sums = mutation_array.sum(axis=0)  # Sum of each column
union_matrix = np.add.outer(row_sums, row_sums) - intersection_matrix


union_matrix[union_matrix == 0] = 1

# Compute the Jaccard similarity
kinship_matrix = intersection_matrix / union_matrix
np.fill_diagonal(kinship_matrix, 1)

samples = mutation_matrix.columns #IDs
kinship_matrix1 = pd.DataFrame(kinship_matrix, index=samples, columns=samples).astype('float16')


#Convert data to SnpData formate
def MakeSnpData(data):
    iid = [["fam0", col] for col in data.index]
    sid = data.columns.tolist()
    val = data.values
    return SnpData(iid=iid, sid=sid, val=val)

# Load the data and match index
phenotypes1 = pd.read_csv('ENA_drug.csv', index_col=0).replace(['S','R'], [0,1])
phenotypes1 = phenotypes1[phenotypes1.index.isin(bugNames1)].astype('int8')
lineage_data1 = pd.read_csv('DST_SAMPLES.csv',index_col='UNIQUEID').replace(np.nan, 'Unknown')  # Lineage, country, and site data
lineage_data1 = lineage_data1[lineage_data1.index.isin(bugNames1)]


for size in [179,529,879,1579,2279,2979,4329,5779,7179]:
    for seed in range(10):
        print(size, seed)
        bugNames = resample(bugNames1, n_samples=size, replace=False, random_state=seed) #size of datasets
        
        # Fix data size and formate
        phenotypes = phenotypes1[phenotypes1.index.isin(bugNames)]  # Convert resistance/sensitivity
        phenotypes = MakeSnpData(phenotypes)
        
        snp_matrix = mutation_matrix[bugNames].T  # Filter columns based on UNIQUEID (bugNames)
        snp_matrix = MakeSnpData(snp_matrix)
        
        kinship_matrix = kinship_matrix1[bugNames].loc[bugNames] # Filter kinship matrix to match the bugNames
        kinship_matrix = MakeSnpData(kinship_matrix)
        
        lineage_data = lineage_data1[lineage_data1.index.isin(bugNames)]  # Filter
        
        # Make Dummies
        lineage_dummies = pd.get_dummies(lineage_data['LINEAGE'], drop_first=True)
        site_dummies = pd.get_dummies(lineage_data['SITE'], drop_first=True)
        country_dummies = pd.get_dummies(lineage_data['COUNTRY'], drop_first=True)
        lineage_data = pd.concat([lineage_dummies, site_dummies, country_dummies], axis=1)
        lineage_data = MakeSnpData(lineage_data)
        
        results_df = single_snp(test_snps=snp_matrix , pheno=phenotypes ,\
                                K0=kinship_matrix, covar=lineage_data,\
                                count_A1=False, leave_out_one_chrom=False)
        
        results_df.to_csv('LMM/'+str(size)+'_lmm_'+str(seed)+'.csv')
