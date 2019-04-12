//Neil Davies 03/07/15
//This creates an indicator for birth month in the UK Biobank data


use "$path1/raw_data/biobank_phenotypes_nmd_150417.dta", clear

//Gen indicator for months and year of birth
gen year_month_birth=100*n_34_0_0+n_52_0_0
tab year_month_birth
replace year_month_birth=100*n_34_0_0+n_52_0_0

//Program to combine the fields 0 and 1 for each variable
cap prog drop merge_var
prog def merge_var
replace `1'_0_0=`1'_1_0 if (`1'_0_0==.|`1'_0_0<0 )& (`1'_1_0>0 & `1'_1_0!=.)
end

cap prog drop merge_var2
prog def merge_var2
replace `1'_0_0=`1'_1_`2' if (`1'_0_`2'==.|`1'_0_`2'<0 )& (`1'_1_`2'>0 & `1'_1_`2'!=.)
end

cap prog drop merge_svar
prog def merge_svar
replace `1'_0_0=`1'_1_0 if (`1'_0_0=="")& (`1'_1_0!="")
end

//Clean years of full time education
//Use Okbay method for defining years of education
forval i = 0/1 {
	forval j = 0/5 {
		g EA_`i'_`j' = 20 if n_6138_0_0 == 1
		replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 2
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 3
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 4
		replace EA_`i'_`j' = 19 if n_6138_`i'_`j' == 5
		replace EA_`i'_`j' = 15 if n_6138_`i'_`j' == 6
		replace EA_`i'_`j' = 7 if n_6138_`i'_`j' == -7
		replace EA_`i'_`j' = . if n_6138_`i'_`j' == -3
		}
	}
// take max 
egen eduyears2 = rmax(EA_*_*)
drop EA_*


//Identify individuals who were not born in England:
merge_var n_1647
gen born_english=(n_1647_0_0==1 & n_20115_0_0==.)


//Generate dob variable
gen yob=n_34_0_0
gen mob=n_52_0_0
gen dob=ym(yob,mob)


//Clean the outcome data
//Physical Exercise
/*
Maximum	7
Decile 9	7
Decile 8	6
Decile 7	5
Decile 6	4
Median	3
Decile 4	3
Decile 3	2
Decile 2	1
Decile 1	0
Minimum	0
2304 items have value -3 (Prefer not to answer)
24680 items have value -1 (Do not know)
*/

merge_var n_904
merge_var n_884
gen out_phys_v_act=n_904_0_0 if n_904_0_0 >=0 &n_904_0_0!=.
gen out_phys_m_act=n_884_0_0 if n_884_0_0 >=0 & n_884_0_0!=.

//Sedentary activity
/*
Maximum	24
Decile 9	5
Decile 8	4
Decile 7	4
Decile 6	3
Median	3
Decile 4	2
Decile 3	2
Decile 2	2
Decile 1	1
Minimum	0

24259 items have value -10 (Less than an hour a day)
767 items have value -3 (Prefer not to answer)
3838 items have value -1 (Do not know)
*/

merge_var n_1070
gen out_sedentary=n_1070_0_0 if n_1070_0_0>0 & n_1070_0_0
replace out_sedentary=0 if n_1070_0_0==-10 & n_1070_0_0

//Income
//See http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=738

merge_var n_738
gen out_income_under_18k=(n_738_0_0>1) if n_738_0_0>0 &n_738_0_0!=.
gen out_income_over_31k=(n_738_0_0>2) if n_738_0_0>0 &n_738_0_0!=.
gen out_income_over_52k=(n_738_0_0>3) if n_738_0_0>0 &n_738_0_0!=.
gen out_income_over_100k=(n_738_0_0>4) if n_738_0_0>0 &n_738_0_0!=.

//Smoking

merge_var n_20116
gen out_smoker=(n_20116_0_0==2) if n_20116_0_0>=0 & n_20116_0_0!=.
gen out_exsmoker=(n_20116_0_0==2|n_20116_0_0==1) if n_20116_0_0>=0 & n_20116_0_0!=.

//Alcohol consumption

merge_var n_1558
gen out_alcohol=6-n_1558_0_0 if n_1558_0_0>0 & n_1558_0_0!=.

//Depression
gen out_depression2=(n_2090_0_0 ==1|n_2090_1_0==1) if (n_2090_0_0 >=0&n_2090_0_0!=.) | (n_2090_1_0 >=0  &n_2090_1_0!=.)

gen out_hes_depress=0
ds s_41202* s_41204_0_*
foreach i in `r(varlist)'{
	di "`i'"
	replace out_hes_depress =1 if substr(`i',1,3)=="F32"|substr(`i',1,3)=="F33"|substr(`i',1,3)=="F34"|substr(`i',1,3)=="F38"|substr(`i',1,3)=="F39"
	}
replace out_depression2=1 if out_hes_depress==1	

//Happiness
merge_var n_4526
gen out_happiness=(6-n_4526_0_0) if n_4526_0_0>0 & n_4526_0_0!=.

//Cognition
merge_var n_20016
gen out_intell=n_20016_0_0

//Blood pressure
merge_var n_4080
merge_var n_4079
egen out_sys_bp=rowmean(n_4080_0_1 n_4080_0_0)
egen out_dia_bp=rowmean(n_4079_0_1 n_4079_0_0)

//Anthropometry
merge_var n_21001
merge_var n_50
gen out_bmi=n_21001_0_0
gen out_height=n_50_0_0

//Arterial Stiffness
merge_var n_21021
merge_svar s_4206

xi:reg n_21021_0_0 i.s_4206_0_0 
predict out_arterial_stiffness,res

//Grip strength
merge_var n_46
merge_var n_47
merge_svar s_38
egen X=rowmean(n_46_0_0 n_47_0_0)
xi:reg X i.s_38_0_0
predict out_gripstrength if X!=.,res
drop X

//Mortality
gen out_dead=(n_40018_0_0!=.)

//Diagnosed with cancer
gen out_cancer=(n_40008_0_0>=30 & n_40008_0_0!=.) if n_40008_0_0>=30 | n_40008_0_0==.

//Had heart attack or stroke
merge_var n_6150

merge_var2 n_6150 1
merge_var2 n_6150 2
merge_var2 n_6150 3

gen out_heartattack=(n_6150_0_0==1|n_6150_0_1==1|n_6150_0_2==1|n_6150_0_3==1) if n_6150_0_0!=-3 & n_6150_0_0!=.
gen out_stroke=(n_6150_0_0==3|n_6150_0_1==3|n_6150_0_2==3|n_6150_0_3==3) if n_6150_0_0!=-3 & n_6150_0_0!=.

//Diagnosed with diabetes
merge_var n_2443
merge_var n_2976
gen out_diabetes=(n_2443_0_0==1) if n_2443_0_0>=0 
replace out_diabete=. if n_2976_0_0<=21 & n_2976_0_0!=.

//Diagnosed with hypertension
merge_var n_2966
gen out_highbloodpressure=(n_2966_0_0>0 &  n_2966_0_0!=.) if n_2966_0_0>0

//Gender
gen male=(n_31_0_0==1)

keep n_845_0_0 n_6138_* out_* male year_month_birth  eduyears2      born_english  yob mob dob  n_eid  
compress

//Limit the sample just to English born individuals
//drop if born_english!=1
joinby n_eid using "$path2/workingdata/covariates",unmatched(master) _merge(XXX)
drop if XXX!=3

//Create month of birth variable

gen mobi=100*yob+mob

compress



save "$path2/workingdata/cleaned_biobank_outcomes_ENGLISH",replace

