#!/bin/bash

#module add languages/gcc-7.1.0
#module add apps/qctool-2.0 

source("scripts/filepaths_r.r")

bgenix="${path5}/bin/gavinband-bgen-798eca81c0fa/build/apps/bgenix"

#genfile="{path6}/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen"
genfile="{path6}/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/dosage_bgen"

# Extract SNPs of interest using bgenix
for chr in {1..22}
	do
		echo "CHR ${chr}"
		${bgenix} -g  ${genfile}/ukb_imp_chr${chr}_v2.bgen -i ${genfile}/ukb_bgi_chr${chr}_v2.bgi -incl-rsids ~/uk_biobank/snplists/hill_okbay_262_snps.txt > ~/uk_biobank/extracted_snps/hill_okbay_262_snps_${chr}.bgen
	done 




# Convert SNPs into additive format using plink2
for chr in {1..22}
	do
		echo "CHR ${chr}"
		 plink2 --bgen ~/uk_biobank/extracted_snps/hill_okbay_262_snps_${chr}.bgen \
		 --sample ~/uk_biobank/raw_data/ukb878_imp_chr10_v2_s487406.sample \
		 --export A --out ~/uk_biobank/extracted_snps/okbay_hill_${chr}
	done 

