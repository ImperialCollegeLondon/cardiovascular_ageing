#!/usr/bin/env bash
# https://documentation.dnanexus.com/getting-started/cli-quickstart
# Usage: <script_name.sh> 

# input
genotype_prefix="/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ BGEN\ format\ -\ interim\ 450k\ release/"
phenotype_prefix="/mnt/project/Derived/Phenotypes/cardiac_age_cmri/regenie_Jul2022/"
path_to_450kwes_helper_files="/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ interim\ 450k\ release/helper_files"

# output
destination_dir="/mnt/project/Derived/Phenotypes/cardiac_age_cmri/regenie_Aug2022/"
# c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 cX cY
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
