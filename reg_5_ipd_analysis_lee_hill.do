//Neil Davies 12/07/18
//This runs the bivariate MR of education and multiple outcomes on IPD
//It uses the Okbay discovery sample and the cognition score at p<5E-08

use "workingdata/analysis_dataset_okbay_hill_score",clear

joinby n_eid using "workingdata/lee_score",unmatched(master)
drop _m
joinby n_eid using "workingdata/unweighted_hill_score",unmatched(master)


egen z_eduyears2=std(eduyears2)
egen z_out_intell=std(out_intell)
tab mob,gen(imob_)

ivreg2 cov_male (z_eduyears2 z_out_intell=allele_score_ea3 allele_score_unweight_cog)  imob_* yob_* pc_* cov_male if _n<_N*.01,ro first cluster(mobi) endog(z_eduyears2 z_out_intell) partial(cov_male imob_* yob_* pc_*)
matrix X=e(first)
local ea_f=X[8,1]
local intell_f=X[8,2]

di "`ea_f'"
regsave using "results/ivreg2_outcome_lee_hill", detail(all) pval ci replace  addlabel(outcome,"`i'",intell_f,`intell_f',ea_f,`ea_f') 


ds out_alcohol out_depression out_exsmoker out_income_over_100k out_income_over_31k out_income_over_52k out_income_under_18k out_phys_m_act out_phys_v_act out_sedentary out_smoker
foreach j in  `r(varlist)'{
	di "`j'"

	ivreg2 `j' (z_eduyears2 z_out_intell=allele_score_ea3 allele_score_unweight_cog)  imob_* yob_* pc_* cov_male,ro first cluster(mobi) endog(z_eduyears2 z_out_intell)  partial(cov_male imob_* yob_* pc_*)
	matrix X=e(first)
	local ea_f=X[8,1]
	local intell_f=X[8,2]
	regsave z_eduyears2 z_out_intell using "results/ivreg2_outcome_lee_hill", detail(all) pval ci append  addlabel(outcome,"`j'",intell_f,`intell_f',ea_f,`ea_f') 
	}
	
	
	
ds out_bmi out_cancer  out_dead out_dia_bp out_diabetes out_gripstrength out_heartattack out_height out_highbloodpressure out_stroke out_sys_bp  /*out_depress*/
foreach j in  `r(varlist)'{
	di "`j'"
	ivreg2 `j' (z_eduyears2 z_out_intell=allele_score_ea3 allele_score_unweight_cog)  imob_* yob_* pc_* cov_male,ro first endog(z_eduyears2 z_out_intell) partial(cov_male imob_* yob_* pc_*) cluster(mobi) 
	matrix X=e(first)
	local ea_f=X[8,1]
	local intell_f=X[8,2]
	regsave z_eduyears2 z_out_intell  using "results/ivreg2_outcome_lee_hill", detail(all) pval ci append addlabel(outcome,"`j'",intell_f,`intell_f',ea_f,`ea_f') 	
	}
	

use "results/ivreg2_outcome_lee_hill",clear


//Clean bivariate results
	replace outcome="out_depression" if  outcome=="out_depression2"
	replace outcome="highbloodpressure" if outcome=="out_highbloodpressure"
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
	drop if Fdf2 ==.
	drop if nx==.
	//Plot results using coefplot
	//Genetic scores
	mkmat coef stderr if var=="z_out_intell" & cont==0,matrix(results_hill) rownames(nx)
	matrix results_hill_iv=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==0,matrix(results_okbay) rownames(nx)
	matrix results_okbay_iv=results_okbay'

	#delimit ;
	coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
			 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
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
	graph export "results/IPD_figure_3_bivariate_ivreg_bin.pdf", as(pdf) replace fontface("Calibri (Body)");
	graph save "results/IPD_figure_3_bivariate_ivreg_bin", replace;
	#delimit cr

	//Continious variables
	mkmat coef stderr if var=="z_out_intell" & cont==1,matrix(results_hill) rownames(nx)
	matrix results_hill_iv=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==1,matrix(results_okbay) rownames(nx)
	matrix results_okbay_iv=results_okbay'


	#delimit ;
	coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
			 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
		, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
		legend(order(2 "Cognition" 4 "Education" ) row(1) size(vsmall))
		xlabel(,noticks)  xtick(none) ytick(none) 
		headings(14 = "{bf:Indicators of aging}"
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
	21="Vigorous exercise (days/week)"
	22="Moderate exercise (days/week)"
	, wrap(40));	
	#delimit cr

	graph save "results/IPD_figure_3_bivariate_ivreg_cont", replace 

	graph combine "results/IPD_figure_3_bivariate_ivreg_bin" "results/IPD_figure_3_bivariate_ivreg_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
	graph export "results/IPD_lee_hill_figure_3_bivariate_`x'_ivreg.eps", as(pdf) replace fontface("Calibri")

 //Clean the results for export to Excel
 
 
