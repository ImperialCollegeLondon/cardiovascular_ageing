#!/bin/bash
# plink-clump 0.0.1
# DNA nexus applet source code

main() {
    
    set -exo pipefail
    
    echo "bgen_input: '$bgen_input_prefix'"

    dx download --lightweight "$summary_stats" -o summary_stats.txt
    awk 'NR!=1{print $3}' summary_stats.txt > rsids.include

    dx download --lightweight "$bgen_sample" -o bgen_sample.sample
    dx download --lightweight "$bgen_input" -o bgen_input.bgen
    
    docker load --input /docker/plink_bgen.tar
    docker run --detach --interactive --tty --volume "${PWD}":/tmp --name plink_bgen plink_bgen
    PLINK="docker exec -w /tmp plink_bgen plink"
    PLINK2="docker exec -w /tmp plink_bgen plink2"
    BGENIX="docker exec -w /tmp plink_bgen bgenix"
    
    $BGENIX -g bgen_input.bgen -index
    $BGENIX -g bgen_input.bgen -incl-rsids rsids.include > tmp.bgen
    $PLINK2 --bgen tmp.bgen ref-first \
        --sample bgen_sample.sample \
        --real-ref-alleles --rm-dup exclude-all \
        --make-bed --out tmp
    rm bgen_input.bgen bgen_input.bgen.bgi
        
    $PLINK --bfile tmp \
        --clump summary_stats.txt \
        --clump-p1 1 \
        --clump-r2 $clump_r2 \
        --clump-kb 250 \
        --clump-snp-field ID \
        --clump-field P \
        --out $bgen_input_prefix
        
    awk 'NR!=1{print $3}' $bgen_input_prefix.clumped > $bgen_input_prefix.snps
    
    $PLINK2 --bfile tmp --extract $bgen_input_prefix.snps --out $bgen_input_prefix.clumped --make-bed

    geno_bed=$(dx upload $bgen_input_prefix.clumped.bed --brief)
    geno_bim=$(dx upload $bgen_input_prefix.clumped.bim --brief)
    geno_fam=$(dx upload $bgen_input_prefix.clumped.fam --brief)
    plink_log=$(dx upload $bgen_input_prefix.clumped.log --brief)
    snp_table=$(dx upload $bgen_input_prefix.clumped --brief)

    dx-jobutil-add-output geno_bed "$geno_bed" --class=file
    dx-jobutil-add-output geno_bim "$geno_bim" --class=file
    dx-jobutil-add-output geno_fam "$geno_fam" --class=file
    dx-jobutil-add-output plink_log "$plink_log" --class=file
    dx-jobutil-add-output snp_table "$snp_table" --class=file
}
