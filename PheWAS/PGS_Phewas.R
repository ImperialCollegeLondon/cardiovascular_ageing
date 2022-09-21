# Code to run PheWAS for polygenic score on UKB traits
# Author: Sean L Zheng, Date: 21-Sep 2022

library(data.table)
library(PheWAS)

# Run if HPC capable of multicore processing
ncores <- detectCores()

# Set output file names and plot title
output_file <- 'lms-ware-analysis/live/sean/phewas/cardiac_aging/results/phewas.result.agedeltas.adjusted.150922.txt'
png_file <- 'lms-ware-analysis/live/sean/phewas/cardiac_aging/results/phewas.result.agedeltas.adjusted.150922.pdf'
plot_title <- 'Cardiac Age Delta PheWAS'

# Load in PRS scores, formating file to PheWAS format (ID, PRS score)
geno = read.csv("lms-ware-analysis/live/sean/phewas/cardiac_aging/age_deltas.csv")
link = fread('lms-ware-analysis/live/sean/Sample_Phenotype_Files/ukbb_eid_link_file.txt')
geno1 = merge(geno, link, by='eid_40616')
c = geno1[,c('eid_47602','catb_delta_with_t1_bc_cole')]
colnames(c) <- c('id','variant')
genotypes=as.data.frame((c),stringsAsFactors=FALSE)
colnames(genotypes)[1]<-c("id")
genotypes$id<-as.integer(genotypes$id)
genotypes <- genotypes[,c('id','variant')]

# Load in clinical data and covariates and join to form final clinical data frame
phenotypes<-fread("lms-ware-analysis/live/sean/phewas/ICD_wide.txt",header=TRUE,sep="\t") #read data summary file
covariates <- fread('lms-ware-analysis/live/sean/Sample_Phenotype_Files/ukbb_genetic_covariates_488k.txt')
colnames(covariates)[1] = 'id'
covariates$age2 = (covariates$age)^2
data=inner_join(inner_join(phenotypes,genotypes) ,covariates)

# Run PheWAS
results=phewas(data=data,phenotypes=names(phenotypes),genotypes="variant",covariates=c('age','age2','sex','pc1','pc2','pc3','pc4','pc5','pc6','pc7','pc8','pc9','pc10'), additive.genotypes=FALSE, cores=ncores,min.records=10)
print('PheWAS analysis complete')
write.table(results,output_file,row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t") #output

# Plotting
pdf(png_file, width = 12, height = 8, useDingbats = F)
phewasManhattan(results, title=plot_title, point.size=1, OR.direction=T, annotate.phenotype.description = T, annotate.size = 2)
dev.off()
