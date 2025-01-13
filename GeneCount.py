import os
import pandas as pd

lmm_folder = "LMM/"  # Replace with the folder path
truth_file = "SNP_correct.csv"  # Replace with the truth file path

truth = pd.read_csv(truth_file)  # Columns: DRUG, MUT, GENE

summary_results = []
aggregated_results = []

for lmm_file in os.listdir(lmm_folder):
    if lmm_file.endswith(".csv"):  # Process only CSV files
        lmm_results = pd.read_csv(os.path.join(lmm_folder, lmm_file))
        
        P_VALUE_THRESHOLD = 0.05 / len(lmm_results)# Define the significance threshold
        
        # Filter LMM results based on the threshold
        lmm_results_filtered = lmm_results[lmm_results['PValue'] < P_VALUE_THRESHOLD]
        lmm_results_filtered['GENE'] = lmm_results_filtered['SNP'].str.split('_').str[0]
        
        valid_lmm_results = lmm_results_filtered[
            lmm_results_filtered.apply(
                lambda row: (
                    (row['SNP'] in truth['MUT'].values) and
                    (row['GENE'] in truth['GENE'].values)
                ), axis=1
            )
        ]
        
        comparison_results = {}
        unique_drugs = ['Amikacin','Ethambutol','Ethionamide','Isoniazid','Kanamycin','Levofloxacin','Moxifloxacin','Rifampicin']
        
        total_true_positives = 0
        total_false_positives = 0
        
        
        for drug in unique_drugs:
            # mutations and genes for the current drug
            true_genes = set(truth[truth['DRUG'] == drug]['GENE'])
            
            # LMM identified genes for the current drug
            predicted_genes = set(valid_lmm_results[valid_lmm_results['Pheno'] == drug.upper()]['GENE'])
            
            true_positive_genes = true_genes & predicted_genes  # Intersection
            false_positive_genes = predicted_genes - true_genes  # In LMM but not in truth
            
            total_true_positives += len(true_positive_genes)
            total_false_positives += len(false_positive_genes)
            
            comparison_results[drug] = {
                "True Positives": len(true_positive_genes),
                "False Positives": len(false_positive_genes)
            }
        
        # Add per-drug results to summary
        for drug, metrics in comparison_results.items():
            summary_results.append({
                "File": lmm_file,
                "Drug": drug,
                **metrics
            })
        
summary_df = pd.DataFrame(summary_results)
print(summary_df)
