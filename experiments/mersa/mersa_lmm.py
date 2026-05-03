import pandas as pd
from pysnptools.snpreader import SnpData
from fastlmm.association import single_snp
from scipy.sparse import csr_matrix, lil_matrix

def sparse_block_dot_product(sparse_array, block_size):
    n = sparse_array.shape[0]
    result = lil_matrix((n, n))  # Use LIL for efficient incremental updates
    for i in range(0, n, block_size):
        for j in range(0, n, block_size):
            block_result = sparse_array[i:i+block_size].dot(sparse_array[j:j+block_size].T)
            result[i:i+block_size, j:j+block_size] = block_result
    return result.tocsr()  # Convert to CSR for fast access

def block_union_matrix(sparse_array, intersection_matrix, block_size):
    n_rows = sparse_array.shape[0]
    row_sums = sparse_array.sum(axis=1).A1  # Convert row sums to 1D dense array
    union_matrix = lil_matrix((n_rows, n_rows))  # Initialize as LIL for efficient updates

    for i in range(0, n_rows, block_size):
        for j in range(0, n_rows, block_size):
            row_block_end = min(i + block_size, n_rows)
            col_block_end = min(j + block_size, n_rows)

            row_block = row_sums[i:row_block_end]
            col_block = row_sums[j:col_block_end]
            intersection_block = intersection_matrix[i:row_block_end, j:col_block_end]

            # Compute union for the block
            union_block = row_block[:, None] + col_block - intersection_block
            union_matrix[i:row_block_end, j:col_block_end] = union_block

    return union_matrix.tocsr()

phenotypes1 = pd.read_csv('RSdata13.csv', index_col=0).astype('int8') #gene
mut = pd.read_csv('snp_calls.csv', sep="\t") #snp 
mut = mut.drop(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
mut = mut.set_index('ID').T
#mut = pd.read_csv('gpa.csv', index_col=0).T
mutation_matrix = mut[mut.index.isin(list(phenotypes1.index))].astype('int8')

mutation_array = csr_matrix(mutation_matrix.values)

# Compute pairwise intersections
intersection_matrix = sparse_block_dot_product(mutation_array, block_size=1000)

union_matrix = block_union_matrix(mutation_array, intersection_matrix, block_size=1000)

# Avoid division by zero
union_matrix.data[union_matrix.data == 0] = 1

# Compute Jaccard similarity
kinship_matrix = intersection_matrix.multiply(union_matrix.power(-1)) # Element-wise division
kinship_matrix.setdiag(1)

samples = mutation_matrix.index  # IDs
kinship_matrix_df = pd.DataFrame.sparse.from_spmatrix(kinship_matrix, index=samples, columns=samples)



#Convert data to SnpData formate
def MakeSnpData(data):
    iid = [["fam0", col] for col in data.index]
    sid = data.columns.tolist()
    val = data.values
    return SnpData(iid=iid, sid=sid, val=val)


phenotypes = MakeSnpData(pd.DataFrame(phenotypes1))
snp_matrix = MakeSnpData(pd.DataFrame(mutation_matrix))
kinship_matrix = MakeSnpData(kinship_matrix_df)

results_df = single_snp(test_snps=snp_matrix , pheno=phenotypes ,\
                        K0=kinship_matrix,count_A1=False, leave_out_one_chrom=False)
