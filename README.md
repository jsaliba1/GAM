# GAM: Group Association Model

GAM is a statistical pipeline designed to detect genotype-phenotype associations involved in microbial drug resistance. This tool is particularly effective for datasets with complex relationships, such as genomic data, and helps researchers determine significant features for predictive modeling.

**This repository holds the input data and code to analyse the results produced by GAM.**

Note: this pipeline was developed for a study on  *M. tuberculosis*  and  *S. aureus*, but in principle can be used on any microbe. 


## Data availability 

All data files to reproduce GAM results for *M. tuberculosis*  and  *S. aureus* are available on [Figshare](https://doi.org/10.6084/m9.figshare.28191398) in Data.zip. 

**Must have at least 30GB free memory**

Note: Data.zip should be unzipped in same path as code being executed. Some files must be unzipped twice. 


## Code and Packages 

**Code was written and tested using Python (version 3.8.1) and Jupyter Notebook (version 7.3.2)**

Additional Packages: 
Pandas (version 1.3.5),
Numpy (version 1.18.2),
Scipy (version 1.6.2),
Pysnptools (version 0.0.2),
Fastlmm (version 0.0.1), and 
Sklearn (version 1.6.1)
  
Install these libraries using pip.

  
## Recommendations

GAM is best suited for datasets with the following characteristics:

1) **High-Quality Phenotype Data:** Datasets should include accurate and reliable phenotype data to ensure robust analysis.

2) **Diverse Resistance Profiles:** A broad range of resistance profiles is essential, with multiple entries for each profile to capture variations across all drugs. This helps avoid incorporating resistance-associated mutations into the control group.

3) **Distributed Resistance:** Resistance to a specific drug should be represented across groups with varying resistance profiles to     improve statistical power.

4) **Multi-Drug Resistance:** Most isolates in the dataset should exhibit resistance to multiple drugs, especially outside a single drug family, to prevent categorization bias.

5) **Adequate Sample Size:** Sufficient numbers of isolates per resistance type are needed, depending on the complexity and prevalence of resistance profiles in the population.

By following these recommendations, GAM can effectively analyze datasets and provide meaningful insights into drug resistance patterns. However, GAM is not limited to datasets with these characteristics. 


## Features

### GAM
- Groups bacterial strains based on resistance patterns.
- Performs resistance analysis and identifies significant mutations.
- Scores mutations for their association with drug resistance using statistical metrics.

### FineAnalysis
- Focuses on key genes associated with resistance.
- Calculates sensitivity, specificity, and positive predictive value (PPV) for mutations.

## Example usage

### GAM.py
**Initialization:** Load the necessary input files by specifying their paths during the creation of a GAM object.
`analysis = GAM('MUTATIONS.csv', 'ENA_drug.csv', 'bugs.csv')`

**Grouping:** group and create a modified resistance matrix.
`analysis.group_bugs(P=0, remove_n=0, seed=0)`
- P: Probability of assigning a drug resistance state (default is 0 for full susceptibility).
- remove_n: Number of isolates to exclude from the analysis (default is 0).
- seed: Random seed for reproducibility (default is 0).

**Resistance Analysis:** Perform resistance analysis to calculate odds ratios and p-values for mutations.
`analysis.analyze_resistance()`

**Scoring Resistance:** Calculate mutation-drug associations and return a filtered matrix of p-values.
`df = analysis.score_resistance()`

Note: if control group is small consider using MERSA_gam.py instead. 

### FineAnalysis.py
**Initialize with input files**
`analysis = FineAnalysis('MUTATIONS.csv', 'ENA_drug.csv', 'input_max.csv', 'Reg/', run_num=1)`

**Analyze SNPs and score them**
`analysis.analyze_snps()` 

**score SNPs**
`results = analysis.score_snps()`

### All other files either use similar logic or can be executed without user input


## Output
The `score_resistance` method outputs a DataFrame where:

Rows represent significant mutations.
Columns represent drugs.
Values are p-values after statistical correction.

## Data Summary 

**MUTATIONS.csv** - Contains mutation data for all isolates collected by the Cryptic consortium. This file is used as an input for GAM. UNIQUEID (string): Identifier for the isolate. `MUTATION` (string): Mutation label. `GENE` (string): Gene where the mutation occurs. `IS_SYNONYMOUS` (Boolean): Indicates if the mutation is synonymous.
###
**ENA_drug.csv** - Stores resistance phenotypes for different drugs. It is used to define resistance groupings in GAM. `UNIQUEID` (string): Identifier for the isolate Drug columns (binary: Resistant/Susceptible): Resistance status for each drug.
### 
**bugs.csv** - Contains a list of isolates to be included in the study, ensuring they belong to a group relevant to the GAM analysis. `ID` (string): Identifier for the isolate. Group (string): Group assignment.
### 
**DST_SAMPLES.csv** - Metadata file providing lineage, site, and country of origin for each sample. Used to categorize isolates in the analysis. `UNIQUEID` (string): Identifier for the isolate. `LINEAGE` (string): Lineage classification. `SITE` (string): Sampling site. `COUNTRY` (string): Country of origin.
### 
**GAM_train.csv** - The main input file for the machine learning (ML) model, indicating the presence or absence of significant mutations identified by GAM. `UNIQUEID` (string): Identifier for the isolate. Mutation columns (binary: 1/0): Presence or absence of significant mutations.
### 
**Y9_train.csv** - Label file for the training dataset used in the ML model. Corresponds to GAM_train.csv. `UNIQUEID` (string): Identifier for the isolate.
### 
**LMM/** (Folder) - Stores results from the LMM analysis across different sample sizes. These results are used as input for calculating summary statistics and filtering significant mutations.

## Notes
Ensure all input files are correctly formatted and placed in the working directory.
Results depend on the quality and completeness of input datasets.
For questions or issues, please open an issue in the repository.

## License
This project is licensed under the MIT License. See `LICENSE` for details.
