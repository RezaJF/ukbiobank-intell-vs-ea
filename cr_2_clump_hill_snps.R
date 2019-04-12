#Neil Davies 15/03/17
#This runs the MR-Base scripts for the education MR paper
#I load the Hill and Okbay GWAS and restrict then to common SNPs available in the June 2017 release of UK Biobank (i.e. HRC SNPs only)

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

#Load Okbay et al education 8,146,840 SNPs
file=paste(path2, "/gwas_results/EduYears_Main.txt",sep="")
okbay_eduyears_exp_dat<-read.table(file, header = TRUE)

#Load UK Biobank June 2017 11,554,957 SNPs (HRC)
file=paste(path2, "/rawdata/qc_rsid.txt",sep="")
ukbb_snps<-read.table(file, header = FALSE)

#Restrict both GWAS to a common set of 7,710,315 SNPs
hill_coefficients_common<-hill_intell_exp_dat[hill_intell_exp_dat$SNP_ID %in% okbay_eduyears_exp_dat$MarkerName,]
okbay_coefficients_common<-okbay_eduyears_exp_dat[okbay_eduyears_exp_dat$MarkerName %in% hill_intell_exp_dat$SNP_ID,]

#Restrict both GWAS to UK Biobank 7,303,122 SNPs
hill_coefficients_common<-hill_coefficients_common[hill_coefficients_common$SNP_ID %in% ukbb_snps$V1,]
okbay_coefficients_common<-okbay_coefficients_common[okbay_coefficients_common$MarkerName %in% ukbb_snps$V1,]

rm(ukbb_snps,file,hill_intell_exp_dat,okbay_eduyears_exp_dat)

#Format Hill and Okbay data
hill_coefficients_fm<-format_data(hill_coefficients_common,snp_col = "SNP_ID", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",type = "cog")
okbay_coefficients_fm<-format_data(okbay_coefficients_common,snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval",type= "EA")
rm(hill_coefficients_common,okbay_coefficients_common)
clump_dat
#Clump the Hill data - this results in 197 genome-wide significant SNPs 
hill_coefficients<-clump_data(subset(hill_coefficients_fm, pval.cog<5e-08),clump_r2 = 0.01)

#Clump the Okbay data to give 79 SNPs
okbay_coefficients<-clump_data(subset(okbay_coefficients_fm, pval.EA<5e-08),clump_r2 = 0.01)

#Create list of SNPs to extract from Biobank
#Select RS IDs
file=paste(path2, "workingdata/hill_intelligence_snps_coef.txt",sep="")
write.table(hill_coefficients, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/hill_intelligence_snps.txt",sep="")
write.table(subset(hill_coefficients,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

file=paste(path2, "workingdata/okbay_eduyears_snps_coef.txt",sep="")
write.table(okbay_coefficients, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/okbay_eduyears_snps.txt",sep="")
write.table(subset(okbay_coefficients,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Create combined list of SNPs detected in either GWAS to extract from UK Biobank

#The following code clumps the EA + cognition results:

#Combining the EA and cognition hits
#Clumping at r2<0.01 results in 79 EA hits and 197 IQ SNPs (6 common)
combined_snps<-unique(rbind(subset(okbay_coefficients,select=c(SNP)),subset(hill_coefficients,select=c(SNP))))
head(combined_snps)
#Output for extraction of 270 from UKBB
file=paste(path2, "workingdata/okbay_hill_all_270_rsid.txt",sep="")
write.table(combined_snps, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Extract the combined list of 238 SNPs (38 were due to overlapping signals for EA and cognition)
snps<-combined_snps$SNP
hill_coefficients_all<-hill_intell_exp_dat[hill_intell_exp_dat$SNP_ID %in% snps,]
okbay_coefficients_all<-okbay_eduyears_exp_dat[okbay_eduyears_exp_dat$MarkerName %in% snps,]

#Clump the resulting lists - results in 219 EA SNPs and 222 cognition SNPs
hill_coefficients_all<-clump_data(format_data(hill_coefficients_all,snp_col = "SNP_ID", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"),clump_r2 = 0.01)
okbay_coefficients_all<-clump_data(format_data(okbay_coefficients_all,snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval"),clump_r2 = 0.01)

#Clean up work space
rm(combined_snps,hill_coefficients,hill_intell_exp_dat,okbay_coefficients,okbay_eduyears_exp_dat,file,snps)

#Output the summary data
file=paste(path2, "workingdata/hill_intelligence_snps_coef_all.txt",sep="")
write.table(hill_coefficients_all, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/hill_intelligence_snps_all.txt",sep="")
write.table(subset(hill_coefficients_all,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

file=paste(path2, "workingdata/okbay_eduyears_snps_coef_all.txt",sep="")
write.table(okbay_coefficients_all, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
file=paste(path2, "workingdata/okbay_eduyears_snps_all.txt",sep="")
write.table(subset(okbay_coefficients_all,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#List of all 238 SNPs to extract from UKBB. This is the list of SNPs that are selected in either the cognition or EA clumping
combined_snps<-unique(rbind(subset(okbay_coefficients_all,select=c(SNP)),subset(hill_coefficients_all,select=c(SNP))))
head(combined_snps)

#Output for extraction from UKBB
file=paste(path2, "workingdata/okbay_eduyears_snps_all_rsid.txt",sep="")
write.table(combined_snps, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Check the allele frequency of selected SNPs in UKBB
file<-paste(path2, "rawdata/qc_rsid.txt",sep="")
ukbb_snps<-read.table(file, header = FALSE)
ukbb_snps_select<-ukbb_snps[ukbb_snps$V1 %in% combined_snps$SNP,]

max(ukbb_snps_select$V2)
min(ukbb_snps_select$V2)

#Max and min MAF are 0.499835 and 0.0235191

#Output rsid effect and other allele
file<-paste(path2, "rawdata/qc_rsid2.txt",sep="")
snpstats<-read.table(file)

#Select 238 SNPs from the two GWAS
snpstats<-snpstats[snpstats$V1 %in% combined_snps$SNP,]
