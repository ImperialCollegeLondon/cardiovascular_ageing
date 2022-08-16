---
title: "Cardiovascular age-delta GWAS and Burden tests"
author: changlubio@github
---

# High quality common SNPs set
1. Liftover array data to GRCh38 (https://github.com/dnanexus-rnd/liftover_plink_beds)
2. Quality check with Plink

## 1. Liftover
```
$java -jar dxCompiler-2.10.2.jar compile liftover_plink_beds/liftover_plink_beds.wdl 
=> workflow-GFZYPGQJ4fX339BGG4KPxB31 
$dx download /Derived/Genotypes/grch38/liftover_input.json
$dx run workflow-GFZYPGQJ4fX339BGG4KPxB31 -f liftover_input.json --brief --destination /Derived/Genotypes/grch38/
=> analysis-GFb4zy0J4fX50p4Q9zPxjjb6
```
Analysis took 22h 18m on RAP. Merged:

Total genotyping rate is 0.969609. 803700 variants and 488377 people pass filters and QC. 

ukb_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed  
ukb_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim  
ukb_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam  
ukb_GRCh38_full_analysis_set_plus_decoy_hla_merged.log  
ukb_GRCh38_full_analysis_set_plus_decoy_hla_merged.nosex 

chr1-22:

ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed  
ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim  
ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam  
ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.log 

## 2. Quality check
```
RAP-> Swiss Army Knife-> command line->  
plink2 
--bfile /mnt/project/Derived/Genotypes/grch38/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged 
--out 500K_array_snps_GRCh38_qc_pass 
--mac 100 
--maf 0.01 
--hwe 1e-15 
--mind 0.1 
--geno 0.1 
--write-snplist 
--write-samples 
--no-id-header 
--threads 8 
```

# REGENIE Step 1
Using `eid_andcardiacAGE_someadj_May2022.csv` the unadjusted column of aged-delta
See also the DNANexus gitbook [here](https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/gwas-ex#regenie-step-1).
```
RAP-> Swiss Army Knife-> command line->  

regenie 
--step 1 
--bed /mnt/project/Derived/Genotypes/grch38/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged 
--phenoFile /mnt/project/Derived/Phenotypes/cardiac_age_cmri/regenie_Jul2022/eid_andcardiacAGE_someadj_May2022.txt 
--covarFile /mnt/project/Derived/Phenotypes/cardiac_age_cmri/regenie_Jul2022/covariates.cov 
--phenoCol catb_delta_with_t1_bc_cole 
--covarCol Sex,Age,Age2,Array,PC{1:10} 
--extract /mnt/project/Derived/Phenotypes/cardiac_age_cmri/regenie_Aug2022/500K_array_snps_GRCh38_qc_pass.snplist 
--keep /mnt/project/Derived/Phenotypes/cardiac_age_cmri/regenie_Aug2022/500K_array_snps_GRCh38_qc_pass.id 
--bsize 1000 
--lowmem 
--lowmem-prefix tmp_regnie_preds 
--threads 16 
--gz 
--out regenie_step1_agedelta  
```

# REGENIE Step 2 - GWAS
```
sh 02_automate_regenie_gwas.sh
```

# REGENIE Step 2 - Burden test
Memory error and gene processing errors were randomly encountered. In order to resolve this, break each chromosomes into chunks, this got me through to processing whole-genome wide gene sets.

1. Generate chunks
```
python3 geneset_chunks.py [size_of_chunk] [, separated list of chroms] [out_dir] 
```
2. Run burden test
```
sh 03_automate_regenie_burden.sh
```
