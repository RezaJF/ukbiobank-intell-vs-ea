//Neil Davies 21/12/17
//This merges the cognition SNPs with the 75 education SNPs

//Estimate association with educational attainment in sample that:
//A) Did not take the fluid intelligence battery
//B) Was not in the interim release

use "workingdata/analysis_dataset_interim",clear
tab mob, gen(imob_)
//************************************************
//Estimate the effect of intelligence on education
//************************************************

//Need to restrict the sampe to the sample without cognition and that was not in the interim release.
//Identify the Hill SNPs
drop n _m rsid alleleA alleleB
gen n=_n
joinby n using "workingdata/clumped_hill_snps",unmatched(master)

forvalues i=1(1)222{
	local rsid=rsid[`i']
	order `rsid'
	}
egen z_eduyears=std(eduyears2)
//They replicate on the cognition outcome
//Estimate the association of these variants and educational attainment in the remaining sample
reg z_eduyears rs55659265_G_A  imob_* yob_* pc_* cov_male if exclude==0,ro
regsave rs55659265_G_A using "results/hill_snps_educ", detail(all) pval ci replace 

ds rs9307617_C_G-rs62036613_A_G /*rs2992037_A_G-rs8138438_A_G*/
foreach i in `r(varlist)'{
	reg z_eduyears `i'  imob_* yob_* pc_* cov_male  if exclude==0,ro
	regsave `i' using "results/hill_snps_educ", detail(all) pval ci append 
	}

//Next we need to clean the Hill coefficients on cognition for the 270 SNPs (197 cognition SNPs and 73 education SNPs + 1 shared.)
use	"results/hill_snps_educ",clear
gen rsid=word(subinstr(var,"_"," ",4),1)
gen ukbb_other=word(subinstr(var,"_"," ",4),2)
gen ukbb_effect=word(subinstr(var,"_"," ",4),3)
rename coef ukbb_beta
rename stderr ukbb_se
rename pval ukbb_pval
keep rsid ukbb_*
duplicates drop
save "workingdata/hill_educ_snps",replace

//Open Hill results and extract the 270 SNPs
import delimited "gwas_results/Intelligence - David Hill 171023.txt", encoding(ISO-8859-1)clear

rename snp_id rsid

joinby rsid using "workingdata/hill_educ_snps",unmatched(using)

//Harmonize alleles in Hill GWAS to be cognition increasing
gen flip=(beta<0)
gen hill_effect=a1 if flip==0
replace hill_effect=a2 if flip==1
gen hill_other=a2 if flip==0
replace hill_other=a1 if flip==1
gen hill_beta=beta if flip==0
replace hill_beta=beta*-1 if flip==1 
rename se hill_se
drop flip
rename p hill_pval

//Harmonize UKBB to Hill filled SNPs
gen flip_ukbb=(hill_effect==ukbb_other)
gen ukbb_flip_effect=ukbb_effect if flip_ukbb==0
replace ukbb_flip_effect=ukbb_other if flip_ukbb==1
gen ukbb_flip_other=ukbb_other if flip_ukbb==0
replace ukbb_flip_other=ukbb_effect if flip_ukbb==1
gen ukbb_flip_beta=ukbb_beta if flip_ukbb==0
replace ukbb_flip_beta=ukbb_beta*-1 if flip_ukbb==1

drop a1 a2 beta _merge ukbb_beta ukbb_effect ukbb_other

rename ukbb_flip_effect ukbb_effect
rename ukbb_flip_other ukbb_other
rename ukbb_flip_beta ukbb_beta 

order rsid chr base_pair_position hill_effect hill_other hill_beta hill_se hill_pval ukbb_effect ukbb_other ukbb_beta ukbb_se  
drop flip_ukbb
compress
save "workingdata/combinedresults",replace

gen hill_snp=1

save "workingdata/combinedresults_final",replace
use "workingdata/combinedresults_final",clear

//drop _m
//Add indicators for study origin
joinby rsid using "workingdata/clumped_hill_snps",unmatched(both)
gen mstudy_hill_clump=(_m==3)
drop _m

joinby rsid using "workingdata/clumped_okbay_snps",unmatched(both)
gen mstudy_okbay_clump=(_m==3)

//Run summary data regressions:

gen aw=1/ukbb_se^2

//Using all 194 SNPs

mregger ukbb_beta hill_beta [aweight=aw] if  mstudy_hill_clump==1 & hill_pval<5e-08,gxse(hill_se) ivw  heter
regsave hill_beta using "results/mregger_ivw_outcome",detail(all) pval replace addlabel(depvar, "education")

//MR-Egger
mregger ukbb_beta hill_beta [aweight=aw] if  mstudy_hill_clump==1 & hill_pval<5e-08,gxse(hill_se) 
regsave using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "education") 	
	
//Weighted median
mrmedian ukbb_beta ukbb_se hill_beta hill_se if  mstudy_hill_clump==1 & hill_pval<5e-08,  w
regsave using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "education") 	
	
//Modal
mrmodal ukbb_beta ukbb_se hill_beta hill_se if  mstudy_hill_clump==1 & hill_pval<5e-08,  weight
regsave using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "education") 


//Plot results
mreggerplot ukbb_beta ukbb_se hill_beta hill_se if  mstudy_hill_clump==1 & hill_pval<5e-08 , gpci
mrmedian ukbb_beta ukbb_se hill_beta hill_se if  mstudy_hill_clump==1 & hill_pval<5e-08 , weighted
addplot : function _b[beta]*x if  mstudy_hill_clump==1 & hill_pval<5e-08, range(0.01 0.05) lc(gs0) lp(shortdash) lw(vthin)
mrmodal ukbb_beta ukbb_se hill_beta hill_se if  mstudy_hill_clump==1 & hill_pval<5e-08,  weight
addplot : function _b[beta]*x if  mstudy_hill_clump==1 & hill_pval<5e-08 , range(0.01 0.05) lc(gs0) lp(longdash) ///
	legend(order(3 "MR-Egger" 2 "MR-Egger 95% CI" 7 "Weighted median" 8 "Modal") rows(1) si(vsmall) symx(*.5)) ///
	graphregion(color(white))  plotregion(lc(white))  ///
	xtitle("SNP-intelligence association (SD)") ytitle("SNP-education association (SD)")
graph export "results/supp_figure_1_effect_cognition_education.eps", as(pdf) replace fontface("Calibri")	

//Data for export
use "results/mregger_ivw_outcome",clear
gen lower_ci=coef-1.96*stderr
gen upper_ci=coef+1.96*stderr

gen double pval2=pval
replace pval2=2*normal(-abs(coef/stderr)) if pval==0
drop pval
rename pval pval

order cmdline cmd var coef stderr lower_ci upper_ci pval k 

gen order=1 in 5
replace order=2 in 3
replace order=3 in 4
replace order=4 in 2
replace order=5 in 1

sort order



//************************************************
//Estimate the effect of education on intelligence
//************************************************

//No restrictions


use "workingdata/analysis_dataset_interim",clear
tab mob, gen(imob_)
//Normalise intelligence 
egen Xout_intell=std(out_intell)

replace out_intell=Xout_intell
drop Xout_intell

drop n _m rsid alleleA alleleB
gen n=_n
joinby n using "workingdata/clumped_okbay_snps",unmatched(master)

forvalues i=1(1)219{
	local rsid=rsid[`i']
	order `rsid'
	}

egen z_intell=std(out_intell)
reg z_intell rs6568547_C_T  imob_* yob_* pc_* cov_male,ro
regsave rs6568547_C_T using "results/hill_snps_intell", detail(all) pval ci replace 

ds rs9648380_A_G-rs7319102_G_A
foreach i in `r(varlist)'{
	reg  z_intell `i'  imob_* yob_* pc_* cov_male ,ro
	regsave `i' using "results/hill_snps_intell", detail(all) pval ci append 
	}

//Next we need to clean the Okbay coefficients on EA for the (194 cognition SNPs and 75 education SNPs + 1 shared.)
use	"results/hill_snps_intell",clear
gen rsid=word(subinstr(var,"_"," ",4),1)
gen ukbb_other=word(subinstr(var,"_"," ",4),2)
gen ukbb_effect=word(subinstr(var,"_"," ",4),3)
rename coef ukbb_beta
rename stderr ukbb_se
rename pval ukbb_pval
keep rsid ukbb_*
duplicates drop
save "workingdata/hill_intell_snps",replace	

//Open Okbay results and extract the 270 SNPs
import delimited "gwas_results/EduYears_Main.txt", encoding(ISO-8859-1)clear

rename markername rsid

joinby rsid using "workingdata/hill_intell_snps",unmatched(using)

//Harmonize alleles in Okbay GWAS to be EA increasing
gen flip=(beta<0)
gen okbay_effect=a1 if flip==0
replace okbay_effect=a2 if flip==1
gen okbay_other=a2 if flip==0
replace okbay_other=a1 if flip==1
gen okbay_beta=beta if flip==0
replace okbay_beta=beta*-1 if flip==1 
rename se okbay_se
drop flip
rename pval okbay_pval

//Harmonize UKBB to Hill filled SNPs
gen flip_ukbb=(okbay_effect==ukbb_other)
gen ukbb_flip_effect=ukbb_effect if flip_ukbb==0
replace ukbb_flip_effect=ukbb_other if flip_ukbb==1
gen ukbb_flip_other=ukbb_other if flip_ukbb==0
replace ukbb_flip_other=ukbb_effect if flip_ukbb==1
gen ukbb_flip_beta=ukbb_beta if flip_ukbb==0
replace ukbb_flip_beta=ukbb_beta*-1 if flip_ukbb==1

drop a1 a2 beta _merge ukbb_beta ukbb_effect ukbb_other

rename ukbb_flip_effect ukbb_effect
rename ukbb_flip_other ukbb_other
rename ukbb_flip_beta ukbb_beta 

order rsid chr pos okbay_effect okbay_other okbay_beta okbay_se okbay_pval ukbb_effect ukbb_other ukbb_beta ukbb_se  
drop flip_ukbb
compress
save "workingdata/combinedresults_intell",replace	
use "workingdata/combinedresults_intell",clear

//Add indicators for study origin
joinby rsid using "workingdata/clumped_hill_snps",unmatched(both)
gen mstudy_hill_clump=(_m==3)
drop _m

joinby rsid using "workingdata/clumped_okbay_snps",unmatched(both)
gen mstudy_okbay_clump=(_m==3)


//Run summary data regressions:

gen aw=1/ukbb_se^2

//Using all 75 SNPs

mregger ukbb_beta okbay_beta [aweight=aw] if mstudy_okbay_clump==1 & okbay_pval<5e-08,gxse(okbay_se) ivw  heter
regsave okbay_beta using "results/mregger_ivw_outcome_EA_on_cog",detail(all) pval replace addlabel(depvar, "education")

//MR-Egger
mregger ukbb_beta okbay_beta [aweight=aw] if mstudy_okbay_clump==1 & okbay_pval<5e-08,gxse(okbay_se) 
regsave using "results/mregger_ivw_outcome_EA_on_cog",detail(all) pval append addlabel(depvar, "education") 	
	
//Weighted median
mrmedian ukbb_beta ukbb_se okbay_beta okbay_se if mstudy_okbay_clump==1 & okbay_pval<5e-08,  w
regsave using "results/mregger_ivw_outcome_EA_on_cog",detail(all) pval append addlabel(depvar, "education") 	
	
//Modal
mrmodal ukbb_beta ukbb_se okbay_beta okbay_se if mstudy_okbay_clump==1 & okbay_pval<5e-08,  weight
regsave using "results/mregger_ivw_outcome_EA_on_cog",detail(all) pval append addlabel(depvar, "education") 


//Plot results
mreggerplot ukbb_beta ukbb_se okbay_beta okbay_se if mstudy_okbay_clump==1 & okbay_pval<5e-08, gpci
mrmedian ukbb_beta ukbb_se okbay_beta okbay_se if mstudy_okbay_clump==1 & okbay_pval<5e-08, weighted
addplot : function _b[beta]*x , range(0.01 0.05) lc(gs0) lp(shortdash) lw(vthin)
mrmodal ukbb_beta ukbb_se okbay_beta okbay_se if mstudy_okbay_clump==1 & okbay_pval<5e-08,  weight
addplot : function _b[beta]*x, range(0.01 0.05) lc(gs0) lp(longdash) ///
	legend(order(3 "MR-Egger" 2 "MR-Egger 95% CI" 7 "Weighted median" 8 "Modal") rows(1) si(vsmall) symx(*.5)) ///
	graphregion(color(white))  plotregion(lc(white)) ///
	xtitle("SNP-education association (SD)") ytitle("SNP-intelligence association (SD)") 
graph export "results/supp_figure_1b_effect_education_on_cognition.eps", as(pdf) replace fontface("Calibri")	
	
//Data for export
use "results/mregger_ivw_outcome_EA_on_cog",clear
gen lower_ci=coef-1.96*stderr
gen upper_ci=coef+1.96*stderr

gen double pval2=pval
replace pval2=2*normal(-abs(coef/stderr)) if pval==0
drop pval
rename pval pval


order cmdline cmd var coef stderr lower_ci upper_ci pval k 

gen order=1 in 5
replace order=2 in 3
replace order=3 in 4
replace order=4 in 2
replace order=5 in 1

sort order

//Check how many SNPs overlap:
//Intelligence
use "workingdata/combinedresults_final",clear


//Add indicators for study origin
joinby rsid using "workingdata/clumped_hill_snps",unmatched(both)
gen mstudy_hill_clump=(_m==3)
drop _m

joinby rsid using "workingdata/clumped_okbay_snps",unmatched(both)
gen mstudy_okbay_clump=(_m==3)

count if mstudy_hill_clump==1 & hill_pval<5e-08

joinby rsid using "workingdata/combinedresults_intell", unmatched(master) _merge(XX)

count if mstudy_hill_clump==1 & hill_pval<5e-08 & XX==3 & okbay_pval<5e-08

//Education
use "workingdata/combinedresults_intell",clear
//Add indicators for study origin
joinby rsid using "workingdata/clumped_hill_snps",unmatched(both)
gen mstudy_hill_clump=(_m==3)
drop _m

joinby rsid using "workingdata/clumped_okbay_snps",unmatched(both)
gen mstudy_okbay_clump=(_m==3)

count if mstudy_okbay_clump==1 & okbay_pval<5e-08

joinby rsid using "workingdata/combinedresults_final", unmatched(master) _merge(XX)

count if mstudy_okbay_clump==1 & okbay_pval<5e-08 & XX==3 & hill_pval<5e-08
