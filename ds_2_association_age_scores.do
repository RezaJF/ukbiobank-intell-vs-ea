//Neil Davies 14/03/19
//This checks whether the scores were associated with age at assessment centre visit

use "$path1/raw_data/biobank_phenotypes_nmd_150417.dta", clear
tab mob,gen(imob_)
//Gen age at assessment centre visit
gen dob=mdy(n_52_0_0,15,n_34_0_0)
gen cov_age= ts_53_0_0-dob 

keep n_eid cov_age
save "workingdata/cov_age",replace


use "workingdata/analysis_dataset_okbay_hill_score",clear
joinby n_eid using "workingdata/cov_age",unmatched(master)
tab mob,gen(imob_)
egen z_as_ea=std(allele_score_ea)
egen z_as_intell=std(allele_score_sniekers)


replace cov_age =cov_age /365.25
reg  cov_age z_as_ea if (within_fam_id ==1|within_fam_id ==.) & eduyears2!=.& out_intell!=.,ro cluster(mobi) 
reg  cov_age z_as_intell if (within_fam_id ==1|within_fam_id ==.) & eduyears2!=.& out_intell!=.,ro cluster(mobi) 

reg  cov_age z_as_ea  pc_* cov_male if (within_fam_id ==1|within_fam_id ==.) & eduyears2!=.& out_intell!=.,ro cluster(mobi) 
reg  cov_age z_as_intell pc_* cov_male if (within_fam_id ==1|within_fam_id ==.) & eduyears2!=.& out_intell!=.,ro cluster(mobi) 
