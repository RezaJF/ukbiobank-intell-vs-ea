//Neil Davies 12/07/18
//This runs the bivariate MR of education and multiple outcomes on IPD
//It uses the Okbay discovery sample and the cognition score at p<5E-08

use "workingdata/analysis_dataset_okbay_hill_score",clear

egen z_eduyears2=std(eduyears2)
egen z_out_intell=std(out_intell)
tab mob,gen(imob_)

replace z_eduyears2=. if interim==1

xtivreg cov_male (z_eduyears2 =allele_score_ea )  imob_* yob_* cov_male ,i(famid) fe  vce(cluster famid )
regsave using "results/xtivreg2_outcome_education", detail(all) pval ci replace  

xtivreg cov_male (z_out_intell =allele_score_sniekers )  imob_* yob_* cov_male ,i(famid) fe  vce(cluster famid )
regsave using "results/xtivreg2_outcome_intell", detail(all) pval ci replace  


ds out_alcohol out_depression out_exsmoker out_income_over_100k out_income_over_31k out_income_over_52k out_income_under_18k out_phys_m_act out_phys_v_act out_sedentary out_smoker
foreach j in  `r(varlist)'{
	di "`j'"

	xtivreg `j' (z_eduyears2 =allele_score_ea )  imob_* yob_* cov_male ,i(famid) fe  vce(cluster famid )
	regsave using "results/xtivreg2_outcome_education", detail(all) pval ci append  

	xtivreg `j' (z_out_intell =allele_score_sniekers )  imob_* yob_* cov_male ,i(famid) fe  vce(cluster famid )
	regsave using "results/xtivreg2_outcome_intell", detail(all) pval ci append  
	}
	
ds out_bmi out_cancer  out_dead out_dia_bp out_diabetes out_gripstrength out_heartattack out_height out_highbloodpressure out_stroke out_sys_bp out_depress
foreach j in  `r(varlist)'{
	di "`j'"
	xtivreg `j' (z_eduyears2 =allele_score_ea )  imob_* yob_* cov_male ,i(famid) fe  vce(cluster famid )
	regsave using "results/xtivreg2_outcome_education", detail(all) pval ci append  

	xtivreg `j' (z_out_intell =allele_score_sniekers )  imob_* yob_* cov_male ,i(famid) fe  vce(cluster famid )
	regsave using "results/xtivreg2_outcome_intell", detail(all) pval ci append  
	}
	

use "results/xtivreg2_outcome_education",clear
append using "results/xtivreg2_outcome_intell",
gen outcome=depvar
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

	drop if !(var=="z_out_intell"| var=="z_eduyears2")|outcome=="cov_male"
	
	//Plot results using coefplot
	//Genetic scores
	mkmat coef stderr if var=="z_out_intell" & cont==0,matrix(results_hill) rownames(nx)
	matrix results_hill_iv=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==0,matrix(results_okbay) rownames(nx)
	matrix results_okbay_iv=results_okbay'

	#delimit ;
	coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.15))
			 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(-0.15)) 
		, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100 (95%CI)",size(vsmall))
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
	graph export "results/IPD_figure_S2_within_family_ivreg_bin.pdf", as(pdf) replace fontface("Calibri (Body)");
	graph save "results/IPD_figure_S2_within_family_ivreg_bin", replace;
	#delimit cr

	//Continious variables
	mkmat coef stderr if var=="z_out_intell" & cont==1,matrix(results_hill) rownames(nx)
	matrix results_hill_iv=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==1,matrix(results_okbay) rownames(nx)
	matrix results_okbay_iv=results_okbay'


	#delimit ;
	coefplot (matrix(results_hill_iv) , se(2) ms(C) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
			 (matrix(results_okbay_iv) , se(2) ms(T) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(-0.15)) 
		, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference (95%CI)",size(vsmall))
		legend(order(2 "Intelligence" 4 "Education" ) row(1) size(vsmall))
		xlabel(,noticks)  xtick(none) ytick(none) 
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

	graph save "results/IPD_figure_S2_within_family_ivreg_cont", replace 

	graph combine "results/IPD_figure_S2_within_family_ivreg_bin" "results/IPD_figure_S2_within_family_ivreg_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
	graph export "results/IPD_figure_S2_within_family_`x'_ivreg.eps", as(pdf) replace fontface("Calibri")

 //Add in the IVW results to create a plot:
 save "results/within_family",replace
 use "results/within_family",clear
 gen analysis="within_family"
  duplicates drop 
 append using "results/hill_okbay_univariate_results"
 duplicates drop 
 
 drop if lower_ci==. & strpos(cmdline,"ivw")==0
 //Drop regression used to create the summary data MR results file
 drop if _n==89

replace coef=coef*-1  if lower_ci==.
replace lower_ci=coef-stderr*1.96 if lower_ci==.
replace upper_ci=coef+stderr*1.96  if upper_ci==.
drop _m
joinby outcome using workingdata/results_order,unmatched(master) update

replace analysis="univariate IVW" if analysis==""
replace var="z_eduyears2" if substr(var,-10,10)=="okbay_beta"
replace var="z_out_intell" if  substr(var,-9,10)=="hill_beta"

sort nx outcome var analysis
order nx outcome var analysis coef stderr lower_ci upper_ci pval

//Create plots:
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

	replace cont=0 if cont==.
//Plot results using coefplot
	//Genetic scores
	mkmat coef stderr if var=="z_out_intell" & cont==0 & analysis=="within_family",matrix(results_hill) rownames(nx)
	matrix results_hill_within=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==0 & analysis=="within_family",matrix(results_okbay) rownames(nx)
	matrix results_okbay_within=results_okbay'
	
	mkmat coef stderr if var=="z_out_intell" & cont==0 & analysis=="univariate IVW",matrix(results_hill) rownames(nx)
	matrix results_hill_ivw=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==0 & analysis=="univariate IVW",matrix(results_okbay) rownames(nx)
	matrix results_okbay_ivw=results_okbay'
	

	#delimit ;
	coefplot (matrix(results_hill_within) , se(2) ms(C) msize(vsmall) mc(dkgreen) ciopts(lc(dkgreen) lwidth(vthin)) offset(+0.3))
			 (matrix(results_hill_ivw) , se(2) ms(C) msize(vsmall) mc(eltgreen) ciopts(lc(eltgreen) lwidth(vthin)) offset(+0.15))
			  (matrix(results_okbay_within) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(-0.15))
			   (matrix(results_okbay_ivw) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.3))
		, graphregion(color(white))  transform(* = min(max(@,-10),30))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100 (95%CI)",size(vsmall))
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
	graph export "results/IPD_figure_S2_within_family_ivreg_bin.pdf", as(pdf) replace fontface("Calibri (Body)");
	graph save "results/IPD_figure_S2_within_family_ivreg_bin", replace;
	#delimit cr

	//Continious variables
	mkmat coef stderr if var=="z_out_intell" & cont==1 & analysis=="within_family",matrix(results_hill) rownames(nx)
	matrix results_hill_within=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==1 & analysis=="within_family",matrix(results_okbay) rownames(nx)
	matrix results_okbay_within=results_okbay'
	
	mkmat coef stderr if var=="z_out_intell" & cont==1 & analysis=="univariate IVW",matrix(results_hill) rownames(nx)
	matrix results_hill_ivw=results_hill'

	mkmat coef stderr if  var=="z_eduyears2" & cont==1 & analysis=="univariate IVW",matrix(results_okbay) rownames(nx)
	matrix results_okbay_ivw=results_okbay'
	

	#delimit ;
	coefplot (matrix(results_hill_within) , se(2) ms(C) msize(vsmall) mc(dkgreen) ciopts(lc(dkgreen) lwidth(vthin)) offset(+0.3))
			 (matrix(results_hill_ivw) , se(2) ms(C) msize(vsmall) mc(eltgreen) ciopts(lc(eltgreen) lwidth(vthin)) offset(+0.15))
			  (matrix(results_okbay_within) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(-0.15))
			   (matrix(results_okbay_ivw) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.3))
		, graphregion(color(white))  plotregion(lc(white)) transform(* = min(max(@,-4),3))  grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference (95%CI)",size(vsmall))
		legend(order(2 "Intelligence" 4 "Education" ) row(1) size(vsmall))
		xlabel(,noticks)  xtick(none) ytick(none) 
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

	graph save "results/IPD_figure_S2_within_family_ivreg_cont", replace 

	graph combine "results/IPD_figure_S2_within_family_ivreg_bin" "results/IPD_figure_S2_within_family_ivreg_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
	graph export "results/IPD_figure_S2_within_family_`x'_ivreg.eps", as(pdf) replace fontface("Calibri")

