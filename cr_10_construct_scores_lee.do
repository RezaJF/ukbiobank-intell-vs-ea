//Neil Davies 12/07/18
//This creates a IQ and education score in UKBB data to conduct bivariate MR on IPD

import delimited "workingdata/lee_eduyears_snps_coef_all.txt", encoding(ISO-8859-1)clear
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
save "workingdata/lee_coef",replace


//Construct EA3 score

use "workingdata/analysis_dataset_interim",clear

drop _m
//Merge in the lee SNPs
joinby n_eid using "workingdata/ea3_598_clean_final.dta",unmatched(master)

drop n
gen n=_n
//Merge in coefficients from Okbay and Hill
drop _m

joinby n using "workingdata/lee_coef",unmatched(master)

//Construct the scores:
gen allele_score_ea3=0

forvalues i=1(1)486{
	local rsid=okbay_rsid[`i']
	local effect=okbay_effect[`i']
	local other=okbay_other[`i']
	local beta=okbay_beta[`i']
	cap:ds `rsid'_`effect'_`other'
	if _rc==111{
		cap:ds `rsid'_`other'_`effect'
		if _rc!=111{
			di "Not reversed `rsid' effect=`effect' other=`other' beta=`beta'"
			if `beta'>0{
				replace allele_score_ea=`rsid'_`other'_`effect'+allele_score_ea
				}
			if `beta'<0{
				replace allele_score_ea=`rsid'_`other'_`effect'*-1+allele_score_ea
				}
			}
		else{
			di "Error `rsid' effect=`effect' other=`other' beta=`beta'"
			}
		}
	else{
		di "Reversed `rsid' effect=`effect' other=`other' beta=`beta'"
		if `beta'>0{
			replace allele_score_ea=`rsid'_`effect'_`other'*-1+allele_score_ea
			}
		if `beta'<0{
			replace allele_score_ea=`rsid'_`effect'_`other'+allele_score_ea
		}		
		}
	}
egen x=std(allele_score_ea3)
replace allele_score=x
keep n_eid allele_score
	
compress
save "workingdata/lee_score",replace
	
