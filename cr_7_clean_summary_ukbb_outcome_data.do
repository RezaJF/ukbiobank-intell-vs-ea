//Neil Davies 12/07/18
//This cleans the summary results file from the UK Biobank data

//Save education + intelligence SNPs 
//Create dta with rsid and indicators for Hill and Okbay studies

//197 intelligence SNPs
import delimited "workingdata/hill_intelligence_snps.txt", delimiter(space) clear 
rename v1 rsid
gen n=_n
save "workingdata/hill_snps",replace

//79 Education SNPs
import delimited "workingdata/okbay_eduyears_snps.txt", delimiter(space) clear 
rename v1 rsid
gen n=_n
save "workingdata/okbay_snps",replace

//Clumped list of 219 EA SNPs and 222 cognition SNPs
import delimited "workingdata/hill_intelligence_snps_all.txt", delimiter(space) clear 
rename v1 rsid
gen n=_n
save "workingdata/clumped_hill_snps",replace

import delimited "workingdata/okbay_eduyears_snps_all.txt", delimiter(space) clear 
rename v1 rsid
gen n=_n
save "workingdata/clumped_okbay_snps",replace

//Create merged file of all 270 SNPs
use "workingdata/okbay_snps",clear
gen mstudy_okbay=1
joinby rsid using "workingdata/hill_snps",unmatched(both)
gen mstudy_hill=(_m==3|_m==2)
keep rsid mstudy*

joinby rsid using "workingdata/clumped_hill_snps",unmatched(both)
gen mstudy_hill_clump=(_m==3)
drop _m

joinby rsid using "workingdata/clumped_okbay_snps",unmatched(both)
gen mstudy_okbay_clump=(_m==3)
drop _m
replace mstudy_okbay=0 if mstudy_okbay==.
compress
save "workingdata/rsid_info",replace

//Clean the results from the UKBB analysis
//Clean results and run summary data MR
use results/outcome_snp2.dta,clear
duplicates drop
drop if depvar=="cov_male"

bys var :gen n=_N
drop if n==11

gen rsid=word(subinstr(var,"_"," ",4),1)
gen ukbb_effect=word(subinstr(var,"_"," ",4),2)
gen ukbb_other=word(subinstr(var,"_"," ",4),3)
duplicates drop


//We need to create a single dataset with one line per SNP, and a column for the coefficient and SE of each of the 22 outcomes:
levels depvar
foreach i in `r(levels)'{
	preserve
	keep if depvar=="`i'"

	rename coef ukbb_beta_`i'
	rename stderr ukbb_se_`i'
	rename pval ukbb_pval_`i'
	keep rsid ukbb_*

	save "workingdata/outcome_`i'",replace
	restore
	}
	
//Clean Hill and Okbay results
use results/outcome_snp2.dta,clear
drop if depvar=="out_happiness"
drop if depvar=="cov_male"

bys var :gen n=_N
drop if n==11

//We need to create a single dataset with one line per SNP, and a column for the coefficient and SE of each of the 23 outcomes:
levels depvar
foreach i in `r(levels)'{
	preserve
	keep if depvar=="`i'"
	gen rsid=word(subinstr(var,"_"," ",4),1)
	gen ukbb_effect=word(subinstr(var,"_"," ",4),2)
	gen ukbb_other=word(subinstr(var,"_"," ",4),3)
	rename coef ukbb_beta_`i'
	rename stderr ukbb_se_`i'
	rename pval ukbb_pval_`i'
	keep rsid ukbb_*

	save "workingdata/outcome_`i'",replace
	restore
	}
//Save 270 education + intelligence SNPs 
gen rsid=word(subinstr(var,"_"," ",4),1)
keep rsid 
duplicates drop
save "workingdata/snplist",replace

//Create file with Okbay education weights
//Can't use the full Okbay results, because the outcome results estimated on the non-interim release
import delimited rawdata/okbay_educ_discovery_replication.txt, clear 
rename markername rsid
joinby rsid using "workingdata/snplist",unmatched(using)
drop chr pos _m
rename a1 okbay_a1
rename a2 okbay_a2
rename eaf okbay_eaf
rename beta okbay_beta
rename se okbay_se
rename pval okbay_p

compress
save "workingdata/okbay",replace

//Merge with Hill SNPs

import delimited "gwas_results/Intelligence - David Hill 171023.txt", encoding(ISO-8859-1)clear
rename snp_id rsid
joinby rsid using "workingdata/okbay",unmatched(using)
drop _m

rename a1 hill_a1
rename a2 hill_a2

rename beta hill_beta
rename se hill_se
rename p hill_p

save "workingdata/hill_okbay_results",replace
	
	
use "workingdata/hill_okbay_results",clear

//Hill and Okbay results are harmonized. Correlation of 0.84 between the coefficients.
//Merge with the outcome results
#delimit ;
foreach i in out_alcohol
out_bmi
out_cancer
out_dead
out_depression2
out_dia_bp
out_diabetes
out_exsmoker
out_gripstrength
out_heartattack
out_height
out_highbloodpressure
out_income_over_100k
out_income_over_31k
out_income_over_52k
out_income_under_18k
out_phys_m_act
out_phys_v_act
out_sedentary
out_smoker
out_stroke
out_sys_bp{;
	joinby rsid using "workingdata/outcome_`i'",unmatched(both);
	tab _m;
	drop _m;
	gen ukbb_aw_`i'=1/ukbb_se_`i'^2;
	};
#delimit cr

//Harmonize alleles in Hill GWAS to be cognition increasing
gen flip=(hill_beta<0)

//Flip Hill results
gen Xhill_a1=hill_a1 if flip==0
replace Xhill_a1=hill_a2 if flip==1

gen Xhill_a2=hill_a2 if flip==0
replace Xhill_a2=hill_a1 if flip==1

gen Xhill_beta=hill_beta if flip==0
replace Xhill_beta=hill_beta*-1 if flip==1 

//Flip Okbay results
gen Xokbay_a1=okbay_a1 if flip==0
replace Xokbay_a1=okbay_a2 if flip==1

gen Xokbay_a2=okbay_a2 if flip==0
replace Xokbay_a2=okbay_a1 if flip==1

gen Xokbay_beta=okbay_beta if flip==0
replace Xokbay_beta=okbay_beta*-1 if flip==1 

replace okbay_eaf =1-okbay_eaf  if flip==1 

drop hill_a1 hill_a2 hill_beta okbay_a1 okbay_a2 okbay_beta 

foreach i in a1 a2 beta{
	rename Xhill_`i' hill_`i'
	rename Xokbay_`i' okbay_`i'
	}

//Check for 24 palindromic SNPs
gen palindromic=((ukbb_effect=="A"&ukbb_other=="T")|(ukbb_effect=="G"&ukbb_other=="C")|(ukbb_effect=="T"&ukbb_other=="A")|(ukbb_effect=="C"&ukbb_other=="G"))

//Get the UKBB allele frequencies

preserve

use "$path2/workingdata/full_merged_dataset.dta",clear

joinby n_eid using "workingdata/okbay_hill_snps_clean.dta",unmatched(master)
drop _m
joinby n_eid using "workingdata/interim.dta",unmatched(master)

drop if _m!=3
//rename rsid Rsid
rename rsid Rsid
collapse (mean) rs*

gen rsid =""
gen ukbb_eaf=.
set obs 270
rename rsid Rsid
ds rs*
foreach i in `r(varlist)'{
	local j=`j'+1
	di "`i'"
	di "`j'"
	replace ukbb_eaf=`i'[1] in `j'
	replace Rsid ="`i'" in `j'
	}
gen Xrsid=word(subinstr(Rsid,"_"," ",4),1)
gen ukbb_effect=word(subinstr(Rsid,"_"," ",4),2)
gen ukbb_other=word(subinstr(Rsid,"_"," ",4),3)

replace ukbb_eaf=ukbb_eaf/2

drop Rsid
rename Xrsid rsid
	
save "workingdata/ukbb_eaf",replace
restore

//rename Rsid rsid
joinby rsid using "workingdata/ukbb_eaf",unmatched(master)
drop _m

gen palindromic_prob=((ukbb_eaf >0.3&ukbb_eaf <0.7) & palindromic ==1)

//Harmonize UKBB to Hill filled SNPs
gen flip_ukbb=(hill_a1==ukbb_other & hill_a2==ukbb_effect)

gen ukbb_flip_effect=ukbb_effect if flip_ukbb==0
replace ukbb_flip_effect=ukbb_other if flip_ukbb==1
gen ukbb_flip_other=ukbb_other if flip_ukbb==0
replace ukbb_flip_other=ukbb_effect if flip_ukbb==1

drop ukbb_effect ukbb_other
rename ukbb_flip_effect ukbb_effect
rename ukbb_flip_other ukbb_other
replace ukbb_eaf=1-ukbb_eaf if flip_ukbb==1

//drop X
#delimit ;
foreach i in out_alcohol
out_bmi
out_cancer
out_dead
out_depression
out_dia_bp
out_diabetes
out_exsmoker
out_gripstrength
out_heartattack
out_height
out_highbloodpressure
out_income_over_100k
out_income_over_31k
out_income_over_52k
out_income_under_18k
out_phys_m_act
out_phys_v_act
out_sedentary
out_smoker
out_stroke
out_sys_bp{;
	gen X=ukbb_beta_`i' if flip_ukbb==0;
	replace X=ukbb_beta_`i'*-1 if flip_ukbb==1;
	drop ukbb_beta_`i';
	rename X ukbb_beta_`i';
	};
#delimit cr
//Match in the SNP selection info
joinby rsid using "workingdata/rsid_info",unmatched(both)
keep if _m==3

compress
save "workingdata/summary_datafile",replace
