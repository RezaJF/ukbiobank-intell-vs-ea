//Neil Davies 10/01/18
//This estimates the associations of the educaiton and cognition SNPs on the outcomes in the sample 
//of individuals who did not take the intelligence battery


use "workingdata/analysis_dataset_interim",clear

tabstat cov_male if exclude==0 & (within_fam_id ==1|within_fam_id ==.),c(s) stats(N sum mean) 
tabstat yob if exclude==0,c(s) stats(N mean sd) 

tab mob, gen(imob_)

tabstat out_highbloodpressure out_diab out_stroke out_heart out_depression out_cancer out_dead out_exsmoker out_smoker out_income_under_18k out_income_over_31k ///
	out_income_over_52k out_income_over_100k if exclude==0 & (within_fam_id ==1|within_fam_id ==.),c(s) stats(N sum mean) 

tabstat out_gripstrength out_height out_bmi out_dia_bp out_sys_bp out_alcohol out_sedentary   out_phys_v_act  out_phys_m_act ///
	 if exclude==0 & (within_fam_id ==1|within_fam_id ==.),c(s) stats(N mean sd ) 

//Estimate the association of each outcome with education and cognition SNPs

//We can't include arterial stiffness, happiness or intelligence, because these were all measured at the same times as intelligence.
order out_happiness out_intell out_arterial_stiffness 


reg cov_male imob_* yob_* pc_* cov_male  if exclude==0 & (within_fam_id ==1|within_fam_id ==.),ro
regsave imob_1 using "results/outcome_snp2", detail(all) pval ci replace 

//ds out_gripstrength -out_highbloodpressure
rename rsid Rsid
ds out_alcohol out_depression2  out_exsmoker out_income_over_100k out_income_over_31k out_income_over_52k out_income_under_18k out_phys_m_act out_phys_v_act out_sedentary out_smoker
foreach j in  `r(varlist)'{
	di "`j'"
	}
foreach j in  `r(varlist)'{
	ds rs*
	foreach i in `r(varlist)'{
		di "`j' `i'"
		reg `j' `i'  imob_* yob_* pc_* cov_male  if exclude==0 & (within_fam_id ==1|within_fam_id ==.),ro
		regsave `i' using "results/outcome_snp2", detail(all) pval ci append 
		}

	}
ds out_bmi out_cancer  out_dead out_dia_bp out_diabetes out_gripstrength out_heartattack out_height out_highbloodpressure out_stroke out_sys_bp
foreach j in  `r(varlist)'{
	di "`j'"
	}
foreach j in  `r(varlist)'{
	ds rs*
	foreach i in `r(varlist)'{
		reg `j' `i'  imob_* yob_* pc_* cov_male  if exclude==0 & (within_fam_id ==1|within_fam_id ==.),ro
		regsave `i' using "results/outcome_snp2", detail(all) pval ci append 
		}

	}


/*
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
//Can use the full Okbay results, because the outcome results estimated on the non-interim release
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


/*
use "workingdata/hill_okbay_results",clear

//Hill and Okbay results are harmonized. Correlation of 0.84 between the coefficients.
//Merge with the outcome results
#delimit ;
foreach i in out_alcohol
out_arterial_stiffness
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
	joinby rsid using "workingdata/outcome_`i'",unmatched(both);
	tab _m;
	drop _m;
	gen ukbb_aw_`i'=1/ukbb_se_`i'^2;
	};
#delimit cr
drop if chr==.
//Identify the Hill SNPs
joinby rsid using  "workingdata/combinedresults_final",unmatched(master) _merge(XX)
drop XX-ukbb_pval


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

drop hill_a1 hill_a2 hill_beta okbay_a1 okbay_a2 okbay_beta 

foreach i in a1 a2 beta{
	rename Xhill_`i' hill_`i'
	rename Xokbay_`i' okbay_`i'
	}


//Harmonize UKBB to Hill filled SNPs
gen flip_ukbb=(hill_a1==ukbb_other)

gen ukbb_flip_effect=ukbb_effect if flip_ukbb==0
replace ukbb_flip_effect=ukbb_other if flip_ukbb==1
gen ukbb_flip_other=ukbb_other if flip_ukbb==0
replace ukbb_flip_other=ukbb_effect if flip_ukbb==1

drop ukbb_effect ukbb_other
rename ukbb_flip_effect ukbb_effect
rename ukbb_flip_other ukbb_other

#delimit ;
foreach i in out_alcohol
out_arterial_stiffness
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

//Create new results file
mregger ukbb_beta_out_alcohol hill_beta [aweight=ukbb_aw_out_alcohol] if hill_snp==1 & hill_p<5E-08 ,gxse(hill_se) ivw  heter
regsave hill_beta using "results/hill_okbay_univariate_results",detail(all) pval replace addlabel(depvar, "hill",outcome,"out_alcohol")

cap prog drop univariate_analysis
prog def univariate_analysis
args outcome exposure hill_snp

//Run univariate summary data analysis
//Using all 186 SNPs

mregger ukbb_beta_`outcome' `exposure'_beta [aweight=ukbb_aw_`outcome'] if hill_snp==`hill_snp' ,gxse(`exposure'_se) ivw  heter
regsave `exposure'_beta using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'")

//MR-Egger
mregger ukbb_beta_`outcome' `exposure'_beta [aweight=ukbb_aw_`outcome'] if hill_snp==`hill_snp' ,gxse(`exposure'_se) 
regsave using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'") 	
	
//Weighted median
mrmedian ukbb_beta_`outcome' ukbb_se_`outcome' `exposure'_beta `exposure'_se if hill_snp==`hill_snp',  w
regsave using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'") 	
	
//Modal
mrmodal ukbb_beta_`outcome' ukbb_se_`outcome' `exposure'_beta `exposure'_se if hill_snp==`hill_snp',  weight
regsave using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'") 

/*
//Plot results
mreggerplot ukbb_beta_`outcome' ukbb_se hill_beta hill_se if hill_snp==1 & low_MAF==0, gpci
mrmedian ukbb_beta_`outcome' ukbb_se hill_beta hill_se if hill_snp==1 & low_MAF==0, weighted
addplot : function _b[beta]*x if hill_snp==1 & low_MAF==0, range(0.01 0.05) lc(gs0) lp(shortdash) lw(vthin)
mrmodal ukbb_beta_`outcome' ukbb_se hill_beta hill_se if hill_snp==1 & low_MAF==0,  weight
addplot : function _b[beta]*x if hill_snp==1 & low_MAF==0, range(0.01 0.05) lc(gs0) lp(longdash) ///
	legend(order(3 "MR-Egger" 2 "MR-Egger 95% CI" 7 "Weighted median" 8 "Modal") rows(1) si(vsmall) symx(*.5)) ///
	graphregion(color(white))  plotregion(lc(white)) title("Effect of cognition on education") ///
	xtitle("SNP-intelligence association") ytitle("SNP-education association")
graph export "results/supp_figure_1_effect_cognition_education.eps", as(pdf) replace fontface("Calibri")	
*/

end

rename ukbb_beta_out_arterial_stiffness ukbb_beta_arterial_stiffness
rename ukbb_aw_out_arterial_stiffness ukbb_aw_arterial_stiffness
rename ukbb_se_out_arterial_stiffness ukbb_se_arterial_stiffness

rename ukbb_beta_out_highbloodpressure ukbb_beta_highbloodpressure
rename ukbb_aw_out_highbloodpressure ukbb_aw_highbloodpressure
rename ukbb_se_out_highbloodpressure ukbb_se_highbloodpressure


#delimit ;
foreach i in out_alcohol
arterial_stiffness
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
highbloodpressure
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
	univariate_analysis `i' hill 1;
	univariate_analysis `i' okbay 0;
	};

reg ukbb_beta_arterial_stiffness hill_beta okbay_beta [aweight=ukbb_aw_arterial_stiffness],ro nocons
regsave hill_beta okbay_beta using "results/hill_okbay_bivariate_results",detail(all) pval replace addlabel(outcome,"arterial_stiffness")

#delimit ;
foreach i in out_alcohol
arterial_stiffness
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
highbloodpressure
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
	di "`i'";
	reg ukbb_beta_`i' hill_beta okbay_beta [aweight=ukbb_aw_`i'],ro;
	regsave hill_beta okbay_beta using "results/hill_okbay_bivariate_results",detail(all) pval append addlabel(outcome,"`i'");
	};

use results/hill_okbay_bivariate_results.dta,clear

*/
