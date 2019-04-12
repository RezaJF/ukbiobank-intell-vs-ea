#Neil Davies 15/03/17
#This runs the MR-Base scripts for the education MR paper
#I load the Hill and Lee GWAS and restrict then to common SNPs available in the June 2017 release of UK Biobank (i.e. HRC SNPs only)

#Load MR base packages
library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(foreign)
#Load TwoSampleMR package
library(TwoSampleMR)

#Install MR base instruments
devtools::install_github("MRCIEU/MRInstruments")

library(MRInstruments)

ao<-available_outcomes()

#Load Hill et al cognition results 7,710,315 SNPs
file=paste(path2, "/gwas_results/Intelligence - David Hill 171023.txt",sep="")
hill_intell_exp_dat<-read.table(file, header = TRUE)

#Load Lee et al education 10,101,242 SNPs
file=paste(path2, "/gwas_results/EduYears_Main.txt",sep="")
okbay_eduyears_exp_dat<-read.table(file, header = TRUE)

#Load UK Biobank June 2017 11,554,957 SNPs (HRC)
file=paste(path2, "/rawdata/qc_rsid.txt",sep="")
ukbb_snps<-read.table(file, header = FALSE)

#Restrict both GWAS to a common set of 7,238,401 SNPs
hill_coefficients_common<-hill_intell_exp_dat[hill_intell_exp_dat$SNP_ID %in% lee_eduyears_exp_dat$MarkerName,]
lee_coefficients_common<-lee_eduyears_exp_dat[lee_eduyears_exp_dat$MarkerName %in% hill_intell_exp_dat$SNP_ID,]

#Restrict both GWAS to UK Biobank 7,233,234 SNPs
hill_coefficients_common<-hill_coefficients_common[hill_coefficients_common$SNP_ID %in% ukbb_snps$V1,]
lee_coefficients_common<-lee_coefficients_common[lee_coefficients_common$MarkerName %in% ukbb_snps$V1,]

#rm(ukbb_snps,file,hill_intell_exp_dat,lee_eduyears_exp_dat)

#Format Hill and Lee data
hill_coefficients_fm<-format_data(hill_coefficients_common,snp_col = "SNP_ID", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",type = "cog")
lee_coefficients_fm<-format_data(lee_coefficients_common,snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval",type= "EA")
rm(hill_coefficients_common,lee_coefficients_common)

#Clump the Hill data - this results in 197 genome-wide significant SNPs 
hill_coefficients<-clump_data(subset(hill_coefficients_fm, pval.cog<5e-08),clump_r2 = 0.01)

#Clump the Lee data to give 505 SNPs
lee_coefficients<-clump_data(subset(lee_coefficients_fm, pval.EA<5e-08),clump_r2 = 0.01)

#Create list of SNPs to extract from Biobank
#Select RS IDs
file=paste(path2, "workingdata/hill_intelligence_snps_coef.txt",sep="")
write.table(hill_coefficients, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/hill_intelligence_snps.txt",sep="")
write.table(subset(hill_coefficients,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

file=paste(path2, "workingdata/lee_eduyears_snps_coef.txt",sep="")
write.table(lee_coefficients, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/lee_eduyears_snps.txt",sep="")
write.table(subset(lee_coefficients,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Create combined list of SNPs detected in either GWAS to extract from UK Biobank

#The following code clumps the EA + cognition results:

#Combining the EA and cognition hits
#Clumping at r2<0.01 results in 79 EA hits and 197 IQ SNPs (6 common)
combined_snps<-unique(rbind(subset(lee_coefficients,select=c(SNP)),subset(hill_coefficients,select=c(SNP))))
head(combined_snps)
#Output for extraction of 270 from UKBB
file=paste(path2, "workingdata/lee_hill_all_681_rsid.txt",sep="")
write.table(combined_snps, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Extract the combined list of 238 SNPs (38 were due to overlapping signals for EA and cognition)
snps<-combined_snps$SNP
hill_coefficients_all<-hill_intell_exp_dat[hill_intell_exp_dat$SNP_ID %in% snps,]
lee_coefficients_all<-lee_eduyears_exp_dat[lee_eduyears_exp_dat$MarkerName %in% snps,]

#Clump the resulting lists - results in 219 EA SNPs and 222 cognition SNPs
hill_coefficients_all<-clump_data(format_data(hill_coefficients_all,snp_col = "SNP_ID", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"),clump_r2 = 0.01)
lee_coefficients_all<-clump_data(format_data(lee_coefficients_all,snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval"),clump_r2 = 0.01)

#Clean up work space
rm(combined_snps,hill_coefficients,hill_intell_exp_dat,lee_coefficients,lee_eduyears_exp_dat,file,snps)

#Output the summary data
file=paste(path2, "workingdata/lee_hill_intelligence_snps_coef_all.txt",sep="")
write.table(hill_coefficients_all, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/lee_hill_intelligence_snps_all.txt",sep="")
write.table(subset(hill_coefficients_all,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)
file=paste(path2, "workingdata/lee_eduyears_snps_coef_all.txt",sep="")
write.table(lee_coefficients_all, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/lee_eduyears_snps_all.txt",sep="")
write.table(subset(lee_coefficients_all,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#List of all 598 SNPs to extract from UKBB. This is the list of SNPs that are selected in either the cognition or EA clumping
combined_snps<-unique(rbind(subset(lee_coefficients_all,select=c(SNP)),subset(hill_coefficients_all,select=c(SNP))))
head(combined_snps)

#Output for extraction from UKBB
file=paste(path2, "workingdata/lee_eduyears_snps_all_rsid.txt",sep="")
write.table(combined_snps, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Output rsid effect and other allele
file<-paste(path2, "rawdata/qc_rsid2.txt",sep="")
snpstats<-read.table(file)

#Select 598 SNPs from the two GWAS
snpstats<-snpstats[snpstats$V1 %in% combined_snps$SNP,]
file<-paste(path2, "rawdata/598_snpstats.txt",sep="")
write.table(snpstats,file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)