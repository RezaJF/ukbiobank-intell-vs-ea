//Neil Davies 24/07/17
//This imports the individual SNP data into Stata and cleans it into additive format

cap prog drop import_snps
prog def import_snps
args file

import delimited "rawdata/`file'", delimiter(tab, collapse) encoding(ISO-8859-1) clear
rename fid n_eid
drop iid pat mat phenotype
compress

save "rawdata/`file'.dta",replace 
end

fs "rawdata/okbay_eduyears_snps_all_270*.raw"
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}


use "rawdata/okbay_eduyears_snps_all_27001.raw.dta",clear
forvalues k=2(1)22{
	if `k'<10{
		local k="0`k'"
		}
	cap:joinby n_eid using "rawdata/okbay_eduyears_snps_all_270`k'.raw.dta",unmatched(both)
	cap:drop _m
	}
gen n=_n
compress
save "workingdata/okbay_eduyears_snps_all_270_clean.dta",replace

//270 SNPs
import delimited "workingdata/okbay_hill_all_270_rsid.txt",clear
rename v1 rsid
save "workingdata/270_snps",replace

//Add in effect and other allele
import delimited rawdata/qc_rsid.txt, delimiter(space) clear 
rename v1 rsid
rename v2 A_allele
rename v3 B_allele
rename v4 A_allele_F
rename v5 MAF

joinby rsid using "workingdata/270_snps",
gen n=_n
save "workingdata/snp_data",replace



//Match in the UKBB effect and other allele from snp-stats file
//Import SNP-stats file into stat
import delimited "rawdata/270_snpstats.txt",clear


foreach i in alternateids	rsid	chromosome	position	alleleA	alleleB	comment	HWexactpvalue	HWlrtpvalue	alleleAcount	alleleBcount	alleleAfrequency	alleleBfrequency	minorallelefrequency	minorallele	majorallele	info	imputeinfo	missingproportion	A	B	AA	AB	BB	NULL	total{
	local k =`k'+1
	rename v`k' `i'
	}
keep rsid alleleA alleleB
gen n=_n
save "workingdata/alleles",replace

use "workingdata/okbay_hill_snps_clean.dta",clear

joinby n using "workingdata/alleles",unmatched(both) _merge(X)

//Rename each SNP
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local allelea=alleleA[`i']
	local alleleb=alleleB[`i']
	
	cap: rename `rsid'_ `rsid'_`allelea'_`alleleb'
	}
compress
save "workingdata/okbay_hill_snps_clean.dta",replace

//****************************************
//Repeat for Snieker intelligence SNPs
//****************************************

fs "rawdata/snieker_intelligence_snps_*.raw"
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}
	
use "rawdata/snieker_intelligence_snps_01.raw.dta",clear
forvalues k=2(1)22{
	if `k'<10{
		local k="0`k'"
		}
	cap:joinby n_eid using "rawdata/snieker_intelligence_snps_`k'.raw.dta",unmatched(both)
	cap:drop _m
	}
gen n=_n
compress
save "workingdata/snieker_intelligence_snps_clean.dta",replace

//16 SNPs
//Add in effect and other allele

//Match in the UKBB effect and other allele from snp-stats file
//Import SNP-stats file into stat
import delimited "rawdata/snieker_snpstats.txt",clear

foreach i in alternateids	rsid	chromosome	position	alleleA	alleleB	comment	HWexactpvalue	HWlrtpvalue	alleleAcount	alleleBcount	alleleAfrequency	alleleBfrequency	minorallelefrequency	minorallele	majorallele	info	imputeinfo	missingproportion	A	B	AA	AB	BB	NULL	total{
	local k =`k'+1
	rename v`k' `i'
	}
keep rsid alleleA alleleB
gen n=_n
save "workingdata/snieker_alleles",replace

use "workingdata/snieker_intelligence_snps_clean.dta",clear

joinby n using "workingdata/snieker_alleles",unmatched(both) _merge(X)

//Rename each SNP
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local allelea=alleleA[`i']
	local alleleb=alleleB[`i']
	
	cap: rename `rsid'_ `rsid'_`allelea'_`alleleb'
	}
compress
save "workingdata/snieker_snps_clean.dta",replace

