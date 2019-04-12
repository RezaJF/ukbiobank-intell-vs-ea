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

fs "rawdata/ea3_598_*.raw"
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}


	
use "rawdata/ea3_598_SNPS_01.raw.dta",clear
forvalues k=2(1)22{
	if `k'<10{
		local k="0`k'"
		}
	cap:joinby n_eid using "rawdata/ea3_598_SNPS_`k'.raw.dta",unmatched(both)
	cap:drop _m
	}
gen n=_n
compress
save "workingdata/ea3_598_clean.dta",replace

//598 SNPs
import delimited "workingdata/lee_eduyears_snps_all_rsid.txt",clear
rename v1 rsid
save "workingdata/598_snps",replace

//Add in effect and other allele
import delimited rawdata/qc_rsid.txt, delimiter(space) clear 
rename v1 rsid
rename v2 A_allele
rename v3 B_allele
rename v4 A_allele_F
rename v5 MAF

joinby rsid using "workingdata/598_snps",
gen n=_n
save "workingdata/snp_data_598",replace


//Match in the UKBB effect and other allele from snp-stats file
//Import SNP-stats file into stata
import delimited "rawdata/598_snpstats.txt",clear
rename v1 rsid
rename v2 A_allele
rename v3 B_allele
rename v4 A_allele_F
rename v5 MAF
gen n=_n
save "workingdata/alleles_598",replace

use "workingdata/ea3_598_clean.dta",clear

joinby n using "workingdata/alleles_598",unmatched(both) _merge(X)

//Rename each SNP
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local allelea=A_allele[`i']
	local alleleb=B_allele[`i']
	
	cap: rename `rsid'_ `rsid'_`allelea'_`alleleb'
	}
compress
drop n-MAF
save "workingdata/ea3_598_clean_final.dta",replace
