//Neil Davies 14/03/19
//This cleans the association estimates

use "workingdata/summary_datafile",clear

#delimit ;
foreach i in out_phys_m_act
out_phys_v_act
out_sedentary
out_alcohol
out_sys_bp
out_dia_bp
out_bmi
out_height
out_gripstrength
out_income_over_100k
out_income_over_52k
out_income_over_31k
out_income_under_18k
out_smoker
out_exsmoker
out_dead
out_cancer
out_depression
out_heartattack
out_stroke
out_diabetes
out_highbloodpressure{;
	order ukbb_beta_`i' ukbb_se_`i' ukbb_pval_`i';
	};

order Rsid chr base_pair_position okbay_beta okbay_se okbay_p okbay_eaf hill_beta hill_se hill_p
sort chr base_pair_position

drop ukbb_aw_out_alcohol-mstudy_okbay_clump
compress
save "results/snp-phenotype_associations",replace
