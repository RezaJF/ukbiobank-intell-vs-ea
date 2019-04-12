//Neil Davies 15/01/18
//Clean run summary data MR

use "workingdata/summary_datafile",clear

drop mstudy_hill mstudy_okbay

//Create new results file
mregger ukbb_beta_out_alcohol hill_beta [aweight=ukbb_aw_out_alcohol] if mstudy_okbay_clump==1 & okbay_p <5e-08 ,gxse(hill_se) ivw  heter
regsave hill_beta using "results/hill_okbay_univariate_results",detail(all) pval replace addlabel(exposure, "hill",outcome,"out_alcohol")

cap prog drop univariate_analysis
prog def univariate_analysis
args outcome exposure study

//Run univariate summary data analysis
//Using all 262 SNPs

mregger ukbb_beta_`outcome' `exposure'_beta [aweight=ukbb_aw_`outcome'] if mstudy_`study'_clump==1 & `study'_p <5e-08 ,gxse(`exposure'_se) ivw  heter
regsave `exposure'_beta using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'")

//MR-Egger
mregger ukbb_beta_`outcome' `exposure'_beta [aweight=ukbb_aw_`outcome'] if mstudy_`study'_clump==1 & `study'_p <5e-08,gxse(`exposure'_se) 
regsave using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'") 	
	
//Weighted median
mrmedian ukbb_beta_`outcome' ukbb_se_`outcome' `exposure'_beta `exposure'_se  if mstudy_`study'_clump==1 & `study'_p <5e-08,  w
regsave using "results/hill_okbay_univariate_results",detail(all) pval append addlabel(exposure, "`exposure'",outcome,"`outcome'") 	
	
//Modal
mrmodal ukbb_beta_`outcome' ukbb_se_`outcome' `exposure'_beta `exposure'_se if mstudy_`study'_clump==1 & `study'_p <5e-08,  weight
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

rename ukbb_beta_out_highbloodpressure ukbb_beta_highbloodpressure
rename ukbb_aw_out_highbloodpressure ukbb_aw_highbloodpressure
rename ukbb_se_out_highbloodpressure ukbb_se_highbloodpressure

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
	univariate_analysis `i' hill hill;
	univariate_analysis `i' okbay okbay;
	};
	
	
//Repeat for the bivariate results
#delimit ;
reg ukbb_beta_out_bmi hill_beta okbay_beta [aweight=ukbb_aw_out_bmi] if mstudy_hill_clump==1,ro nocons;
regsave hill_beta okbay_beta using "results/bivariate_results_hill",detail(all) pval replace addlabel(outcome,"BMI");

reg ukbb_beta_out_bmi hill_beta okbay_beta [aweight=ukbb_aw_out_bmi]  if mstudy_okbay_clump==1,ro nocons;
regsave hill_beta okbay_beta using "results/bivariate_results_okbay",detail(all) pval replace addlabel(outcome,"BMI");


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
	reg ukbb_beta_`i' hill_beta okbay_beta [aweight=ukbb_aw_`i']  if mstudy_hill_clump==1,ro nocons;
	regsave hill_beta okbay_beta using "results/bivariate_results_hill",detail(all) pval append addlabel(outcome,"`i'");
	
	reg ukbb_beta_`i' hill_beta okbay_beta [aweight=ukbb_aw_`i'] if mstudy_okbay_clump==1,ro nocons;
	regsave hill_beta okbay_beta using "results/bivariate_results_okbay",detail(all) pval append addlabel(outcome,"`i'");
	};
#delimit cr

//Clean univariate results
use "results/hill_okbay_univariate_results",clear	
duplicates drop

replace coef=coef*-1

gen lower_ci=coef-stderr*1.96
gen upper_ci=coef+stderr*1.96
order cmd outcome coef stderr lower_ci upper_ci pval
gen order=1 if strpos(cmdline,"ivw")!=0
replace order=2 if strpos(var,"slope")!=0
replace order=3 if strpos(var,"_cons")!=0
replace order=4 if cmd=="mrmedian"
replace order=5 if cmd=="mrmodal"

bys outcome order:gen X=_n
drop if X==3
sort outcome exposure order 

gen estimator="IVW" if order==1
replace estimator="MR-Egger:slope" if order==2
replace estimator="MR-Egger:constant" if order==3
replace estimator="Weighted median" if order==4
replace estimator="Weighted mode" if order==5

joinby outcome using workingdata/results_order,unmatched(master)
drop _m

sort exposure nx order estimator
 

gen cont=0
#delimit ;
foreach i in 
out_gripstrength
out_arterial_stiffness
out_height
out_bmi
out_dia_bp
out_sys_bp
out_intell
out_happiness
out_alcohol
out_sedentary
out_phys_m_act
out_phys_v_act{;
	replace cont =1 if outcome=="`i'";
	};
#delimit cr

//Plot results using coefplot
//Genetic scores
mkmat coef stderr if strpos(cmdline,"ivw")!=0 & exposure=="hill" & cont==0,matrix(results_hill) rownames(nx)
matrix results_hill_iv=results_hill'

mkmat coef stderr if strpos(cmdline,"ivw")!=0 & exposure!="hill" & cont==0,matrix(results_okbay) rownames(nx)
matrix results_okbay_iv=results_okbay'

#delimit ;
coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.15))
		 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(-0.15)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100",size(small))
	/*xlabel(-20(5)20,noticks)*/ rescale(100)  legend(off)
	headings(1 = "{bf:Morbidity}"
			 7 = "{bf:Mortality}"
			 8 = "{bf:Health behaviours}"
			 10 = "{bf:Income}"
			 , labsize(small)) 
	coeflabels(1="Hypertension"
2="Diabetes"
3="Stroke"
4="Heart attack"
5="Depression"
6="Cancer"
7="Died"
8="Ever smoked"
9="Currently smoke"
10="Income over £18k"
11="Income over £31k"
12="Income over £52k"
13=".                            Income over £100k", wrap(51));
graph export "results/figure_2_univariate_ivreg_bin.pdf", as(pdf) replace fontface("Calibri (Body)");
graph save "results/figure_2_univariate_ivreg_bin", replace;
#delimit cr

//Continious variables
mkmat coef stderr if strpos(cmdline,"ivw")!=0 & exposure=="hill" & cont==1,matrix(results_hill) rownames(nx)
matrix results_hill_iv=results_hill'

mkmat coef stderr if strpos(cmdline,"ivw")!=0 & exposure!="hill" & cont==1,matrix(results_okbay) rownames(nx)
matrix results_okbay_iv=results_okbay'


#delimit ;
coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.15))
		 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(-0.15)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
	legend(order(2 "Intelligence" 4 "Education" ) row(1) size(vsmall))
	xlabel(-3(1)3,noticks)  xtick(none) ytick(none) 
	headings(14 = "{bf:Indicator of aging}"
			 15 = "{bf:Anthropometry}"
			 17 = "{bf:Blood pressure}"
			 19 = "{bf:Health behaviours}"
			 , labsize(small)) 
	coeflabels(14="Gripstrength (kg)"
15="Height (cm)"
16="BMI (kg/m2)"
17="Diastolic (mmHg)"
18="Systolic (mmHg)"
19="Alcohol consumption (0 low to 5 high)"
20="Hours watching television per day"
21="Vigorous physical activity (days/week)"
22="Moderate physical activity (days/week)"
, wrap(40));	
#delimit cr

graph save "results/figure_2_univariate_ivreg_cont", replace 

graph combine "results/figure_2_univariate_ivreg_bin" "results/figure_2_univariate_ivreg_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/figure_2_univariate_ivreg.eps", as(pdf) replace fontface("Calibri")



//Clean bivariate results
//Hill 
foreach x in hill /* okbay*/{
	use results/bivariate_results_`x'.dta,clear
	replace coef=coef*-1
	joinby outcome using workingdata/results_order,unmatched(both)

	sort var nx

	gen lower_ci=coef-1.96*stderr
	gen upper_ci=coef+1.96*stderr
	order outcome coef stderr lower upper  pva

	gen cont=0
	#delimit ;
	foreach i in 
	out_gripstrength
	out_arterial_stiffness
	out_height
	out_bmi
	out_dia_bp
	out_sys_bp
	out_intell
	out_happiness
	out_alcohol
	out_sedentary
	out_phys_m_act
	out_phys_v_act{;
		replace cont =1 if outcome=="`i'";
		};
	#delimit cr

	drop if nx==.
	//Plot results using coefplot
	//Genetic scores
	mkmat coef stderr if var=="hill_beta" & cont==0,matrix(results_hill) rownames(nx)
	matrix results_hill_iv=results_hill'

	mkmat coef stderr if  var=="okbay_beta" & cont==0,matrix(results_okbay) rownames(nx)
	matrix results_okbay_iv=results_okbay'

	#delimit ;
	coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.15))
			 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(-0.15)) 
		, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100",size(small))
		/*xlabel(-20(5)20,noticks)*/ rescale(100)  legend(off)
		headings(1 = "{bf:Morbidity}"
				 7 = "{bf:Mortality}"
				 8 = "{bf:Health behaviours}"
				 10 = "{bf:Income}"
				 , labsize(small)) 
		coeflabels(1="Hypertension"
	2="Diabetes"
	3="Stroke"
	4="Heart attack"
	5="Depression"
	6="Cancer"
	7="Died"
	8="Ever smoked"
	9="Currently smoke"
	10="Income over £18k"
	11="Income over £31k"
	12="Income over £52k"
	13=".                            Income over £100k", wrap(51));
	graph export "results/figure_2_bivariate_ivreg_bin.pdf", as(pdf) replace fontface("Calibri (Body)");
	graph save "results/figure_2_bivariate_ivreg_bin", replace;
	#delimit cr

	//Continious variables
	mkmat coef stderr if var=="hill_beta" & cont==1,matrix(results_hill) rownames(nx)
	matrix results_hill_iv=results_hill'

	mkmat coef stderr if  var=="okbay_beta" & cont==1,matrix(results_okbay) rownames(nx)
	matrix results_okbay_iv=results_okbay'


	#delimit ;
	coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.15))
			 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(-0.15)) 
		, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
		legend(order(2 "Intelligence" 4 "Education" ) row(1) size(vsmall))
		xlabel(-3(1)3,noticks)  xtick(none) ytick(none) 
		headings(14 = "{bf:Indicator of aging}"
				 15 = "{bf:Anthropometry}"
				 17 = "{bf:Blood pressure}"
				 19 = "{bf:Health behaviours}"
				 , labsize(small)) 
		coeflabels(14="Gripstrength (kg)"
	15="Height (cm)"
	16="BMI (kg/m2)"
	17="Diastolic (mmHg)"
	18="Systolic (mmHg)"
	19="Alcohol consumption (0 low to 5 high)"
	20="Hours watching television per day"
	21="Vigorous physical activity (days/week)"
	22="Moderate physical activity (days/week)"
	, wrap(40));	
	#delimit cr

	graph save "results/figure_2_bivariate_ivreg_cont", replace 

	graph combine "results/figure_2_bivariate_ivreg_bin" "results/figure_2_bivariate_ivreg_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
	graph export "results/figure_2_bivariate_`x'_ivreg.eps", as(pdf) replace fontface("Calibri")
	}
 
