import os
import pandas as pd

truth = pd.read_csv("SNP_correct.csv")  # Columns: Drug, Mutation
lmm_folder = "LMM/"  # Update with your folder path

truth.rename(columns={'DRUG': 'Drug', 'MUT': 'Mutation'}, inplace=True)

summary_results = []

# Iterate over all LMM result files in the folder
for lmm_file in os.listdir(lmm_folder):
    if lmm_file.endswith(".csv"):  # Process only CSV files
        
        lmm_results = pd.read_csv(os.path.join(lmm_folder, lmm_file)) # Load the current LMM result file
        
        P_VALUE_THRESHOLD = 0.05 / len(lmm_results)
        lmm_results_filtered = lmm_results[lmm_results['PValue'] < P_VALUE_THRESHOLD]
        
        # Standardize column names for consistency
        lmm_results_filtered.rename(columns={'Pheno': 'Drug', 'SNP': 'Mutation'}, inplace=True)
        
        comparison_results = []
        
        # Get the unique drugs in the truth
        unique_drugs = truth['Drug'].unique()
        
        for drug in unique_drugs:
            true_mutations = set(truth[truth['Drug'] == drug]['Mutation']) # mutations for the current drug
            
            # LMM identified mutations for the current drug
            predicted_mutations = set(lmm_results_filtered[lmm_results_filtered['Drug'] == drug.upper()]['Mutation'])
            
            true_positives = true_mutations & predicted_mutations  # Intersection
            false_positives = predicted_mutations - true_mutations  # In LMM but not in truth
            false_negatives = true_mutations - predicted_mutations  # In truth but not in LMM
            
            comparison_results.append({
                "File": lmm_file,
                "Drug": drug,
                "True Positives": len(true_positives),
                "False Positives": len(false_positives),
            })
        
        summary_results.extend(comparison_results)


summary_df = pd.DataFrame(summary_results)
print(summary_df)
