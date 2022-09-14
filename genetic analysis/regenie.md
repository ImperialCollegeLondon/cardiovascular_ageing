# "Cardiovascular age-delta GWAS and Burden tests on the UKB Exome data"
@changlubio

This is a step-by-step tutorial to perform genome-wide association and gene burden test on the UKB exomes using a pre-generated phenotype 'age-delta', using [Regenie](https://rgcgithub.github.io/regenie/) and the [UKB-RAP](https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/gwas-ex#regenie-step-1). Codes were run using the dx command line tools and the Regenie in the Swiss Army Knife app. Major steps include:

1. Generation of high quality common (HQ-common) SNPs set.
2. Run Regenie step 1 on the HQ-common SNPs to fit a whole genome regression model that captures a good fraction of the phenotype variance attributable to genetic effects.
3. Run Regenie step 2 to test for association with the phenotype conditional upon the prediction from the step 1. 

## HQ common SNPs set
To acquire a set of HQ common SNPs, we first converted the array data to GRCh38. Plink2 was then used to do quality check.

### (1) Liftover
Liftover was performed on RAP following https://github.com/dnanexus-rnd/liftover_plink_beds. It took 22.5h. 803700 variants and 488377 people were retained in the final set; total genotyping rate is 0.969609.

```
$java -jar dxCompiler-2.10.2.jar compile liftover_plink_beds/liftover_plink_beds.wdl 
=> workflow-GFZYPGQJ4fX339BGG4KPxB31 
$dx download /Derived/Genotypes/grch38/liftover_input.json
$dx run workflow-GFZYPGQJ4fX339BGG4KPxB31 -f liftover_input.json --brief --destination /Derived/Genotypes/grch38/
=> analysis-GFb4zy0J4fX50p4Q9zPxjjb6
```
Outputs are in Derived/Genotypes/grch38/:

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

### (2) Quality check

For SNPs in Chromosome 1-22, require: MAC > 100, MAF > 1%, HWE p-value > 1e-15, genotyping missingness < 10%; also require sample missingness < 10%.

```
RAP-> Swiss Army Knife-> command line->  
plink2 
--bfile /mnt/project/Derived/Genotypes/grch38/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged 
--out 500K_array_snps_GRCh38_qc_pass 
--mac 100 --maf 0.01 --hwe 1e-15 
--mind 0.1 --geno 0.1 
--write-snplist --write-samples 
--no-id-header 
--threads 8 
```

## REGENIE Step 1
Using `eid_andcardiacAGE_someadj_May2022.csv` the unadjusted column of aged-delta

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

## REGENIE Step 2 - GWAS
```
for chrom in c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 cX cY
do 
    chrom_b0_v1=$chrom\_b0_v1
    cmd="
    regenie \
    --step 2 \
    --pred $destination_dir/regenie_step1_agedelta_pred.list\
    --bgen $genotype_prefix/ukb23150_$chrom_b0_v1.bgen --ref-first \
    --sample $genotype_prefix/ukb23150_$chrom_b0_v1.sample \
    --keep $destination_dir/500K_array_snps_GRCh38_qc_pass.id \
    --phenoFile $phenotype_prefix/eid_andcardiacAGE_someadj_May2022.txt \
    --covarFile $phenotype_prefix/covariates.cov \
    --phenoCol catb_delta_with_t1_bc_cole --covarColList Sex,Age,Age2,Array,PC{1:10} \
    --firth --approx \
    --pThresh 0.01 \
    --bsize 1000 \
    --out agedelta_firth_wstep1_$chrom_b0_v1
    "
    # echo $cmd
    dx run swiss-army-knife -icmd="$cmd" \
        --destination "$destination_dir" \
        -y --brief \
        --tag agedelta \
        --instance-type mem1_hdd1_v2_x8
done
```

## REGENIE Step 2 - Burden test
Memory error and gene processing errors were randomly encountered. In order to resolve this, break each chromosomes into smaller chunks. Finally it got me through to processing all gene sets across the whole genome.

1. Generate chunks
```
python3 geneset_chunks.py [size_of_chunk] [, separated list of chroms] [out_dir] 
```
2. Run burden test
```
for chrom in c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 cX
do
    for n in 0 1 2
    do
        chrom_b0_v1=$chrom\_b0_v1
        helperid=$chrom\_n$n
        if exists_in_list "$list_of_set" " " $helperid;   then
            cmd="
            regenie \
            --step 2 \
            --pred $destination_dir/regenie_step1_agedelta_pred.list \
            --bgen $genotype_prefix/ukb23150_$chrom_b0_v1.bgen --ref-first \
            --sample $genotype_prefix/ukb23150_$chrom_b0_v1.sample \
            --keep $destination_dir/500K_array_snps_GRCh38_qc_pass.id \
            --phenoFile $phenotype_prefix/eid_andcardiacAGE_someadj_May2022.txt \
            --covarFile $phenotype_prefix/covariates.cov \
            --phenoCol catb_delta_with_t1_bc_cole --covarColList Sex,Age,Age2,Array,PC{1:10} \
            --set-list $path_to_custom_helper_files/ukb23149_450k_OQFE_$helperid.sets.tsv  \
            --anno-file $path_to_450kwes_helper_files/ukb23149_450k_OQFE.annotations.txt.gz \
            --mask-def $destination_dir/custome_mask.txt \
            --aaf-bins 0.01,0.001 --bsize 200 \
            --out agedelta_regenie_burden_02_$chrom_b0_v1.chunk_$n
            "
            # echo $cmd
            # --set-list $path_to_450kwes_helper_files/ukb23149_450k_OQFE.sets.txt.gz  \
            # --extract-setlist \"GRHL3(ENSG00000158055)\" \
            # --exclude $path_to_450kwes_helper_files/ukb23149_450k_OQFE.90pct10dp_qc_variants.txt \
            dx run swiss-army-knife -icmd="$cmd" \
                --destination "$destination_dir" \
                -y --brief \
                --tag agedelta \
                --instance-type mem1_ssd1_v2_x8
        fi
    done
done
```
