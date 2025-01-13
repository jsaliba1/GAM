import pandas as pd

truth = pd.read_csv("SNP_correct.csv")  # Columns: Drug, Mutation
lmm_results = pd.read_csv("LMM/7179_lmm_0.csv")  # Columns: PHENO (Drug), SNP (Mutation)

# Define the significance threshold
P_VALUE_THRESHOLD = 0.05/len(lmm_results)
lmm_results_filtered = lmm_results[lmm_results['PValue'] < P_VALUE_THRESHOLD]

truth.rename(columns={'DRUG': 'Drug', 'MUT': 'Mutation'}, inplace=True)
lmm_results_filtered.rename(columns={'Pheno': 'Drug', 'SNP': 'Mutation'}, inplace=True)

comparison_results = {}

unique_drugs = truth['Drug'].unique()

for drug in unique_drugs:
    true_mutations = set(truth[truth['Drug'] == drug]['Mutation'])
    
    # LMM identified mutations for the current drug
    predicted_mutations = set(lmm_results_filtered[lmm_results_filtered['Drug'] == drug.upper()]['Mutation'])
    
    true_positives = true_mutations & predicted_mutations  # Intersection
    false_positives = predicted_mutations - true_mutations  # In LMM but not in truth
    false_negatives = true_mutations - predicted_mutations  # In truth but not in LMM
    
    
    comparison_results[drug] = {
        "True Positives": len(true_positives),
        "False Positives": len(false_positives),
        "False Negatives": len(false_negatives), 
    }

results_df = pd.DataFrame.from_dict(comparison_results, orient='index')

print(results_df)
