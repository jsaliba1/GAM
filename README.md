# GAM: Group Association Model

GAM is a statistical pipeline designed to detect genotype-phenotype associations involved in microbial drug resistance. This tool is particularly effective for datasets with complex relationships, such as genomic data, and helps researchers determine significant features for predictive modeling.

**This repository holds the input data and code to analyse the results produced by GAM.**

Note: this pipeline was developed for a study on  *M. tuberculosis*  and  *S. aureus*, but in principle can be used on any microbe. 
  
  
## Code and Packages 

**Code was written and tested using Python (version 3.8.1) and Jupyter Notebook (version 7.3.2)**

Additional Packages: 
Pandas (version 1.3.5),
Numpy (version 1.18.2),
Scipy (version 1.6.2),
Pysnptools (version 0.0.2),
Fastlmm (version 0.0.1), and 
Sklearn (version 1.6.1)
  
  
## Recommendations

GAM is best suited for datasets with the following characteristics:

1) **High-Quality Phenotype Data:** Datasets should include accurate and reliable phenotype data to ensure robust analysis.

2) **Diverse Resistance Profiles:** A broad range of resistance profiles is essential, with multiple entries for each profile to capture variations across all drugs. This helps avoid incorporating resistance-associated mutations into the control group.

3) **Distributed Resistance:** Resistance to a specific drug should be represented across groups with varying resistance profiles to     improve statistical power.

4) **Multi-Drug Resistance:** Most isolates in the dataset should exhibit resistance to multiple drugs, especially outside a single drug family, to prevent categorization bias.

5) **Adequate Sample Size:** Sufficient numbers of isolates per resistance type are needed, depending on the complexity and prevalence of resistance profiles in the population.

By following these recommendations, GAM can effectively analyze datasets and provide meaningful insights into drug resistance patterns. However, GAM is not limited to datasets with these characteristics. 
