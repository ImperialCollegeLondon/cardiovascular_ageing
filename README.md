# Environmental and genetic predictors of human cardiovascular ageing

## Abstract
Cardiovascular ageing is a complex interaction of physiological processes leading to accumulated irreversible damage at cellular, tissue and organ level. Here, we use machine learning cardiac imaging analysis to predict cardiovascular age in UK Biobank participants to identify environmental modifiers and associated common genetic variants to provide insights into this process and provide targets to attenuate ageing. 

## Content

The steps involved in this process are described in this repository.

1. **Cardiac Image analysis**
- Automated pipeline for image segmentation and motion analysis [code](https://github.com/baiwenjia/ukbb_cardiac/tree/2b6d6371be9a666a41627926324030c31897f877)   

2. **Predicting cardiovascular age** 
- CatBoost model to predict age from a dataset comprising cardiac MRI derived phenotypes [code](https://github.com/ImperialCollegeLondon/cardiovascular_ageing/tree/main/predicting%20cardiac%20age)
- Subsequent bias correction on these predictions and cardiovascular age gap ("cardiovascular age delta") is generated for each subject. 

3. **Phenotype analysis**
- Assessing associations between cardiovascular age delta and categorical (e.g. presence of disease) / continuous risk factors (e.g. amount of alcohol consumed per day) [code](https://github.com/ImperialCollegeLondon/cardiovascular_ageing/tree/main/phenotype%20analysis)
- Assessing associations between cardiovascular age delta and mortality [NEEDS TO BE UPLOADED- johannamielke]
- Medication analysis [code](https://github.com/ImperialCollegeLondon/cardiovascular_ageing/tree/main/self-rep-med-analysis)

4. **Genetic analysis**
- Assessing associations between cardiovascular age delta and genetic common variants. [COMMON VARIANT ANALYSIS NEEDS TO BE UPLOADED- chbender]
- Assessing associations between cardiovascular age delta and genetic rare variants [code](https://github.com/ImperialCollegeLondon/cardiovascular_ageing/tree/main/genetic%20analysis) (update once neater structure)
- Computation of Polygenic Risk Score (PRS) [code](https://github.com/ImperialCollegeLondon/cardiovascular_ageing/tree/main/genetic%20analysis/prs).

5. **PheWAS**
- Assessing associations between cardiovascular age delta PRS and an unbiased selection of phenotypes from the UK Biobank dataset. [NEEDS TO BE UPLOADED- seanzheng33]

## License
Distributed under an MIT LICENSE license. See LICENSE for more information.

## Citation
TBC
