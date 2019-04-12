#Neil Davies 12/07/18
#This extracts the coefficients for the 17 SNPs reported in the Snieker paper


#Load MR base packages
library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(foreign)
#Load TwoSampleMR package
library(TwoSampleMR)



infile=paste(path2, "gwas_results/iq - sniekers - full gwas.txt",sep="")
snieker_disc_exp_dat<-read.table(infile, header = TRUE)
outfile<-paste(path2, "gwas_results/iq - sniekers - full gwas.rda",sep="")
saveRDS(snieker_disc_exp_dat,outfile)

#12,104,294 SNPs in Snieker GWAS
snieker_format<-format_data(snieker_disc_exp_dat, type="exposure", snp_col = "rsid", beta_col = "Beta", se_col = "SE", eaf_col = "MAF", effect_allele_col = "ref", other_allele_col = "alt", pval_col = "p_value")

#Restrict to 9,199,317 SNPs in UK Biobank
file<-paste(path2, "rawdata/qc_rsid2.txt",sep="")
ukbb_snps<-read.table(file, header = FALSE)
snieker_coefficients_common<-snieker_format[snieker_format$SNP %in% ukbb_snps$V1,]

#Clump using r2=0.01 & p<5E-08 results in 16 SNPs.
snieker_clump<-clump_data(subset(snieker_coefficients_common, pval.exposure<5e-08),clump_r2 = 0.01)

#Output list to extract from UKBB
file<-paste(path2, "rawdata/snieker_intelligence_snps.txt",sep="")
write.table(subset(snieker_clump,select=c(SNP)), file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)

#Output coefficient data to construct scores:
file<-paste(path2, "rawdata/snieker_intelligence_snps_coef_all.txt",sep="")
write.table(snieker_clump, file, sep="\t" , quote=FALSE, row.names = FALSE, col.names = TRUE)
