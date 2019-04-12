//Neil Davies 17/07/18
//This creates table 1 a description of the SNPs used in the analysis

//Extract EAF from UKBB:

use "workingdata/analysis_dataset_interim",clear

replace rsid=""
gen EAF=.
rename rsid Rsid
ds rs*
foreach i in `r(varlist)'{
	local l=`l'+1
	replace Rsid ="`i'" in `l'
	tabstat `i',save
	local X=el(r(StatTotal),1,1)/2
	di "`X'"
	replace EAF=`X' in `l'
	}
keep Rsid EAF
drop if Rsid==""
rename Rsid rsid
save "workingdata/EAF",replace	
use "workingdata/EAF",clear
replace rsid=word(subinstr(rsid,"_"," ",4),1)
save "workingdata/EAF",replace	

//194 Hill SNPs
use "workingdata/combinedresults_final",clear


//Add indicators for study origin
joinby rsid using "workingdata/clumped_hill_snps",unmatched(both)
gen mstudy_hill_clump=(_m==3)
drop _m

drop if hill_pval>5e-08
drop if mstudy_hill_clump !=1

joinby rsid using "workingdata/EAF",unmatched(both)
order rsid chr base_pa hill_e hill_o EAF
sort chr base_pair_position 

//75 Education SNPs
use "workingdata/combinedresults_intell",clear

//Add indicators for study origin
joinby rsid using "workingdata/clumped_okbay_snps",unmatched(both)
gen mstudy_okbay_clump=(_m==3)


drop if okbay_pval>5e-08
drop if mstudy_okbay_clump !=1
sort chr pos 

drop _m
joinby rsid using "workingdata/EAF",unmatched(master)
order rsid chr pos okbay_e okbay_o EAF
sort chr pos

