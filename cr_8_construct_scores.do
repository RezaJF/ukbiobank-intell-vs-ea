//Neil Davies 12/07/18
//This creates a IQ and education score in UKBB data to conduct bivariate MR on IPD

import delimited "rawdata/hill_intelligence_snps_coef.txt", encoding(ISO-8859-1)clear
rename snp hill_rsid
rename effect hill_effect
rename other hill_other
rename pvalc hill_pval
rename beta hill_beta
keep hill*
gen n=_n
save "workingdata/hill_coef",replace

import delimited "rawdata/okbay_eduyears_snps_coef_all.txt", encoding(ISO-8859-1)clear
rename snp okbay_rsid
rename effect okbay_effect
rename other okbay_other
rename pvalex okbay_pval
rename beta okbay_beta
keep okbay*
gen n=_n
rename okbay_rsid rsid
joinby rsid using "workingdata/rsid_info",unmatched(master)
keep if okbay_pval <5E-08
rename rsid okbay_rsid
keep okbay_rsid okbay_effect okbay_other okbay_beta n
save "workingdata/okbay_coef",replace

import delimited "rawdata/snieker_intelligence_snps_coef_all.txt", encoding(ISO-8859-1)clear
rename snp snieker_rsid
rename effect snieker_effect
rename other snieker_other
rename pvalexposure snieker_pval
rename beta snieker_beta
keep snieker*
gen n=_n
save "workingdata/snieker_coef",replace

use "workingdata/analysis_dataset_interim",clear

drop _m
//Merge in the Sniekers SNPs
joinby n_eid using "workingdata/snieker_snps_clean.dta",unmatched(master)

drop n
gen n=_n
//Merge in coefficients from Okbay and Hill
drop _m
joinby n using "workingdata/hill_coef",unmatched(master)
drop _m
joinby n using "workingdata/okbay_coef",unmatched(master)
drop _m
joinby n using "workingdata/snieker_coef",unmatched(master)
drop _m

//Construct the scores:
gen allele_score_cog=0
gen allele_score_ea=0
gen allele_score_sniekers=0

	
forvalues i=1(1)197{
	local rsid=hill_rsid[`i']
	local effect=hill_effect[`i']
	local other=hill_other[`i']
	local beta=hill_beta[`i']
	cap:ds `rsid'_`effect'_`other'
	if _rc==111{
		cap:ds `rsid'_`other'_`effect'
		if _rc!=111{
			di "Reversed `rsid' effect=`effect' other=`other' beta=`beta'"
			replace allele_score_cog=`rsid'_`other'_`effect'*`beta'*-1+allele_score_cog
			}
		else{
			di "Error `rsid' effect=`effect' other=`other' beta=`beta'"
			}
		}
	else{
		di "Not reversed `rsid' effect=`effect' other=`other' beta=`beta'"
		replace allele_score_cog=`rsid'_`effect'_`other'*`beta'+allele_score_cog
		}
	}
	
forvalues i=1(1)75{
	local rsid=okbay_rsid[`i']
	local effect=okbay_effect[`i']
	local other=okbay_other[`i']
	local beta=okbay_beta[`i']
	cap:ds `rsid'_`effect'_`other'
	if _rc==111{
		cap:ds `rsid'_`other'_`effect'
		if _rc!=111{
			di "Not reversed `rsid' effect=`effect' other=`other' beta=`beta'"
			replace allele_score_ea=`rsid'_`other'_`effect'*`beta'+allele_score_ea
			}
		else{
			di "Error `rsid' effect=`effect' other=`other' beta=`beta'"
			}
		}
	else{
		di "Reversed `rsid' effect=`effect' other=`other' beta=`beta'"
		replace allele_score_ea=`rsid'_`effect'_`other'*`beta'*-1+allele_score_ea
		}
	}

forvalues i=1(1)16{
	local rsid=snieker_rsid[`i']
	local effect=snieker_effect[`i']
	local other=snieker_other[`i']
	local beta=snieker_beta[`i']
	cap:ds `rsid'_`effect'_`other'
	if _rc==111{
		cap:ds `rsid'_`other'_`effect'
		if _rc!=111{
			di "Not reversed `rsid' effect=`effect' other=`other' beta=`beta'"
			replace allele_score_sniekers=`rsid'_`other'_`effect'*`beta'+allele_score_sniekers
			}
		else{
			di "Error `rsid' effect=`effect' other=`other' beta=`beta'"
			}
		}
	else{
		di "Reversed `rsid' effect=`effect' other=`other' beta=`beta'"
		replace allele_score_sniekers=`rsid'_`effect'_`other'*`beta'*-1+allele_score_sniekers
		}
	}
//Unweighted Hill score:
gen allele_score_unweight_cog=0
forvalues i=1(1)197{
	local rsid=hill_rsid[`i']
	local effect=hill_effect[`i']
	local other=hill_other[`i']
	local beta=hill_beta[`i']
	cap:ds `rsid'_`effect'_`other'
	if _rc==111{
		cap:ds `rsid'_`other'_`effect'
		if _rc!=111{
			if `beta'>0{
				replace allele_score_unweight_cog=`rsid'_`other'_`effect'+allele_score_unweight_cog
				}
			if `beta'<0{
				replace allele_score_unweight_cog=`rsid'_`other'_`effect'*-1+allele_score_unweight_cog
				}
			}
		else{
			di "Error `rsid' effect=`effect' other=`other' beta=`beta'"
			}
		}
	else{
		di "Not reversed `rsid' effect=`effect' other=`other' beta=`beta'"
			if `beta'>0{
				replace allele_score_unweight_cog=`rsid'_`effect'_`other'*-1+allele_score_unweight_cog
				}
			if `beta'<0{
				replace allele_score_unweight_cog=`rsid'_`effect'_`other'+allele_score_unweight_cog
			}
		}
	}	
compress
save "workingdata/analysis_dataset_okbay_hill_score",replace

keep n_eid allele_score_unweight_cog
save "workingdata/unweighted_hill_score",replace
	
