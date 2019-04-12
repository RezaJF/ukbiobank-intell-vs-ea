//Neil Davies 14/03/19
//This creates bias component plots for the Sniekers and Okbay scores

//Bias test as per Davies et al. (2017)

cap prog drop hetero_test
prog def hetero_test

args outcome exposure iv 

macro shift 3
local covar="`*'"

cap drop _const
cap gen _const  = 1

di "outcome=`outcome'"
di "exposure=`exposure'"
di "instrumen=`iv'"
di "covariates=`covar'"

gmm (`outcome' - {xb1:`exposure' `covar' _const})  ///
	(`outcome' - {xb2:`exposure' `covar' _const}) , ///
	instruments(1:`exposure' `covar') ///
	instruments(2:`iv' `covar') ///
	winit(unadjusted,independent) onestep  ///
	vce(cluster mob) ///
	deriv(1/xb1 = -1) ///
	deriv(2/xb2 = -1)
drop _const

local outcome2=substr("`outcome'",1,16)
est sto results_`outcome2'

lincom _b[xb1:`exposure']-_b[xb2:`exposure']
local het_p=2*(1-normal(abs(r(estimate)/r(se))))

regsave `exposure' using "results/bias_plots_basic_adjusted_`outcome'_`iv'", detail(all) pval ci replace addvar(het_p,`het_p')

end


use "workingdata/analysis_dataset_okbay_hill_score",clear

egen z_eduyears2=std(eduyears2)
egen z_out_intell=std(out_intell)
tab mob,gen(imob_)

replace z_eduyears2=. if interim==1

//Keep the unrelated individuals
keep if (within_fam_id ==1|within_fam_id ==.)
drop if interim==1

//Get SD data
tabstat cov_comp_bodysize10 cov_comp_height10 cov_birthweight cov_dist_lon cov_birth_location_imd cov_birth_location_easting cov_birth_location_northing ///
	cov_breastfed cov_father_alive cov_mother_alive cov_num_sisters cov_num_brothers cov_matsmoking ,stats(sd) c(s) varw(30)


ds cov_comp_bodysize10 cov_comp_height10 cov_birthweight cov_dist_lon cov_birth_location_imd cov_birth_location_easting cov_birth_location_northing  ///
		cov_breastfed cov_father_alive cov_mother_alive cov_num_sisters cov_num_brothers cov_matsmoking
foreach i in `r(varlist)'{
	hetero_test `i' z_eduyears2 allele_score_ea imob_* cov_male pc_*
	hetero_test `i'  z_out_intell allele_score_sniekers imob_* cov_male pc_* 
	}



//Clean the results and construct coefficient plots

use "results/bias_plots_basic_adjusted_cov_mother_alive_allele_score_sniekers",clear

#delimit ;
foreach i in 
cov_comp_bodysize10 cov_comp_height10 cov_birthweight cov_dist_lon cov_birth_location_imd cov_birth_location_easting cov_birth_location_northing  ///
		cov_breastfed cov_father_alive cov_mother_alive cov_num_sisters cov_num_brothers cov_matsmoking{;
	append using "results/bias_plots_basic_adjusted_`i'_allele_score_sniekers";
	};
#delimit cr

gen outcome=substr(word(cmdline ,2),2,.)
gen instrument="allele_score_sniekers" if substr(var,1,3)=="xb2"
replace instrument="z_out_intell" if substr(var,1,3)=="xb1"

gen het_test=coef[_n+1] if substr(var,1,3)=="xb2"
order outcome instrument coef stderr het_test
drop if var=="het_p"

keep outcome instrument coef stderr het_test cmdline
compress


cap  prog drop sd_adjust
prog def sd_adjust
args variable sd
replace coef=coef/`sd' if outcome=="`variable'"
replace stderr=stderr/`sd' if outcome=="`variable'"
end

sd_adjust cov_comp_bodysize10	.6760717
sd_adjust cov_comp_height10	.6756012
sd_adjust cov_birthweight	.6602215
sd_adjust cov_dist_lon	146.8861
sd_adjust cov_birth_location_imd	9711.428
sd_adjust cov_birth_location_easting	92650.99
sd_adjust cov_birth_location_northing	166613.1
sd_adjust cov_breastfed	.4530773
sd_adjust cov_father_alive	.4224312
sd_adjust cov_mother_alive	.4889007
sd_adjust cov_num_sisters	1.091811
sd_adjust cov_num_brothers	1.142767
sd_adjust cov_matsmoking	.458508


save "results/merged_allele_score_sniekers",replace


//Non-genetic covariates
//Merge in labels
joinby outcome using "workingdata/labels2.dta",unmatched(master) update
sort order
drop if order==.

drop if outcome=="cov_current_location_easting"
drop if outcome=="cov_current_location_northing"

mkmat coef stderr if instrument=="allele_score_sniekers" & substr(outcome,1,2)!="al" & outcome!="cov_male",matrix(results) rownames(order)
matrix results_as_ea=results'

mkmat coef stderr if instrument=="z_out_intell"& substr(outcome,1,2)!="al"& outcome!="cov_male",matrix(results) rownames(order)
matrix results_out_ea=results'

#delimit ;
coefplot (matrix(results_out_ea) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0))
		 (matrix(results_as_ea) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(-0.1))
	     , graphregion(color(white))  xtitle("Bias in SD units",size(small)) plotregion(color(white)) grid(none) xline(0) ysize(5.345) xsize(7.27) ylabel(,labsize(tiny)) xlabel(,labsize(tiny)) 
	     xlabel(-1(0.5)1,noticks)
		 legend(order( 2 "Measured intelligence" 4 "Sniekers et al allele score" ) row(1) size(tiny))
		 headings(1 = "{bf:Birth location}"
			 6 = "{bf:Current location}"
			 9 = "{bf:Early life}"
			 14 = "{bf:Family characteristics}"
			 , labsize(tiny)) 
	coeflabels(1="Easting"
2="Northing"
3="Index of Multiple Deprivation"
5="Distance from London"
9="Birthweight"
10="Breastfed"
11="Mother smoked in pregnancy"
12="Comparative body size age 10"
13="Comparative height age 10"
14="Father alive"
15="Mother alive"
16="Number of brothers"
17="Number of sisters");


#delimit cr
graph export "results/figure_supp_X2_balance_phenotypes_intell.eps", as(pdf) replace 


use "results/bias_plots_basic_adjusted_cov_mother_alive_allele_score_ea",clear
#delimit ;
foreach i in 
cov_comp_bodysize10 cov_comp_height10 cov_birthweight cov_dist_lon cov_birth_location_imd cov_birth_location_easting cov_birth_location_northing  ///
		cov_breastfed cov_father_alive cov_mother_alive cov_num_sisters cov_num_brothers cov_matsmoking{;
	append using "results/bias_plots_basic_adjusted_`i'_allele_score_ea";
	};
#delimit cr

gen outcome=substr(word(cmdline ,2),2,.)
gen instrument="allele_score_ea" if substr(var,1,3)=="xb2"
replace instrument="z_eduyears2" if substr(var,1,3)=="xb1"

gen het_test=coef[_n+1] if substr(var,1,3)=="xb2"
order outcome instrument coef stderr het_test
drop if var=="het_p"

keep outcome instrument coef stderr het_test cmdline
compress

cap  prog drop sd_adjust
prog def sd_adjust
args variable sd
replace coef=coef/`sd' if outcome=="`variable'"
replace stderr=stderr/`sd' if outcome=="`variable'"
end

sd_adjust cov_comp_bodysize10	.6760717
sd_adjust cov_comp_height10	.6756012
sd_adjust cov_birthweight	.6602215
sd_adjust cov_dist_lon	146.8861
sd_adjust cov_birth_location_imd	9711.428
sd_adjust cov_birth_location_easting	92650.99
sd_adjust cov_birth_location_northing	166613.1
sd_adjust cov_breastfed	.4530773
sd_adjust cov_father_alive	.4224312
sd_adjust cov_mother_alive	.4889007
sd_adjust cov_num_sisters	1.091811
sd_adjust cov_num_brothers	1.142767
sd_adjust cov_matsmoking	.458508

save "results/merged_allele_score_ea",replace

//Non-genetic covariates
//Merge in labels
joinby outcome using "workingdata/labels2.dta",unmatched(master) update
sort order
drop if order==.

drop if outcome=="cov_current_location_easting"
drop if outcome=="cov_current_location_northing"

mkmat coef stderr if instrument=="allele_score_ea" & substr(outcome,1,2)!="al" & outcome!="cov_male",matrix(results) rownames(order)
matrix results_as_ea=results'

mkmat coef stderr if instrument=="z_eduyears2"& substr(outcome,1,2)!="al"& outcome!="cov_male",matrix(results) rownames(order)
matrix results_out_ea=results'

#delimit ;
coefplot (matrix(results_out_ea) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0))
		 (matrix(results_as_ea) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(-0.1))
		 , graphregion(color(white))  xtitle("Bias in SD units",size(small)) plotregion(color(white)) grid(none) xline(0) ysize(5.345) xsize(7.27) ylabel(,labsize(tiny)) xlabel(,labsize(tiny)) 
	     xlabel(-1(0.5)1,noticks)
		 legend(order( 2 "Measured educational attainment" 4 "Okbay et al allele score"  ) row(1) size(tiny))
		 headings(1 = "{bf:Birth location}"
			 6 = "{bf:Current location}"
			 9 = "{bf:Early life}"
			 14 = "{bf:Family characteristics}"
			 , labsize(tiny)) 
	coeflabels(1="Easting"
2="Northing"
3="Index of Multiple Deprivation"
5="Distance from London"
9="Birthweight"
10="Breastfed"
11="Mother smoked in pregnancy"
12="Comparative body size age 10"
13="Comparative height age 10"
14="Father alive"
15="Mother alive"
16="Number of brothers"
17="Number of sisters");

graph export "results/figure_supp_X2_balance_phenotypes_ea.eps", as(pdf) replace ;
#delimit cr
//Create Excel spreadsheets of results
save "workingdata/temp",replace
use "workingdata/temp",clear
drop if outcome=="cov_male"
keep if instrument=="bw12"
gen genetic=(substr(outcome,1,2)!="al")
rename coef bw12_coef
gen bw12_lci=bw12_coef-1.96*stderr
gen bw12_uci=bw12_coef+1.96*stderr
rename het_test bw12_het_pval
format %9.3f bw12_* 
keep outcome bw12_* order genetic
sort genetic order
order  outcome bw12_coef bw12_lci bw12_uci bw12_het_pval
save "results/bias_test_bw12",replace

use "workingdata/temp",clear
drop if outcome=="cov_male"
keep if instrument=="allele_score_21"
gen genetic=(substr(outcome,1,2)=="al")

rename coef prs_coef
gen prs_lci=prs_coef-1.96*stderr
gen prs_uci=prs_coef+1.96*stderr
rename het_test prs_het_pval
format %9.3f prs_* 
keep outcome prs_* order genetic
sort genetic order
order  outcome prs_coef prs_lci prs_uci prs_het_pval
save "results/bias_test_prs",replace

use "workingdata/temp",clear
drop if outcome=="cov_male"
keep if instrument=="more_educ_15" & analysis=="allele_score_21"

gen genetic=(substr(outcome,1,2)!="al")

rename coef educ_coef
gen educ_lci=educ_coef-1.96*stderr
gen educ_uci=educ_coef+1.96*stderr
format %9.3f educ_* 
keep outcome educ_* order genetic
sort genetic order
order  outcome educ_coef educ_lci educ_uci 
save "results/bias_test_educ",replace	
	
joinby outcome using "results/bias_test_bw12"
joinby outcome using "results/bias_test_prs"
drop order genetic

format %9.2e *het_pval

//Finally calculating the association between each of the covariates and the outcomes:

use "workingdata/full_merged_dataset",clear

//Gen updated weights from Hughes et al.
gen weight2=34.29 if weight!=1
replace weight2=12.37 if weight==1

drop allele_score_6
reg out_alcohol allele_score_1 imob_* cov_male pc_*[pweight=weight2],ro
regsave allele_score_1 using "results/covariate_outcome",replace detail(all)
ds out_*
foreach j in `r(varlist)'{
	ds allele_score_* cov_*
	foreach i in `r(varlist)'{
		di "`i' `j'"
		reg `j' `i' imob_* cov_male pc_*[pweight=weight2],ro
		regsave `i' using "results/covariate_outcome",append detail(all)
		}
	}
use "results/covariate_outcome",clear

egen pval=2*(1-norm(abs(coef/stderr)))

keep var depvar coef pval
levels depvar 
foreach i in `r(levels)'{
	preserve
	keep if depvar=="`i'"
	gen n=_n
	local j=`j'+1
	save "workingdata/temp`j'"
	restore
	}
use "workingdata/temp1",clear
rename coef coef_1
rename pval pval_1
forvalues i=2(1)25{
	joinby n using "workingdata/temp`i'",
	rename coef coef_`i'
	rename pval pval_`i'
	}
rename var outcome	
joinby outcome using "workingdata/labels.dta",unmatched(master)	
drop _m
joinby outcome using "workingdata/labels2.dta",unmatched(master) update

forvalues i=1(1)25{
	replace coef_`i'=. if pval_`i'>0.05
	}

drop n pval*
replace order=order+48 if _m==3
sort order
