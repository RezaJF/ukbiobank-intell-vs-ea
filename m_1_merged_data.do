//Neil Davies 08/05/17
//This merges the allele scores to the ROSLA dataset and defines the exclusions for cryptic relatedness and ethnicity



//Clean exclusions list non-europeans
import delimited "rawdata/data.non_europeans_exclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/non-europeans",replace

//Europeans
import delimited "rawdata/data.europeans_inclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/europeans",replace

//Sex mismatches, sex chromosome aneuploidy, and excess heterogeneity
import delimited "rawdata/meta.recommended_exclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/exclusions",replace

//Relateds
import delimited "rawdata/data.relateds_exclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/relateds",replace

//Generate indicator for interim release
use "$path3/data.4688.dta",clear
gen interim=(n_22000_0_0!=2000 & n_22000_0_0!=.)
keep n_eid interim
compress
save "workingdata/interim",replace

//PCs
import delimited "rawdata/data.pcas1-40.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
forvalues i=3(1)42{
	local j=`i'-2
	rename v`i' pc_`j'
	}
save "workingdata/PCs",replace

//Standard covariates
import delimited "rawdata/data.covariates.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2 m
save "workingdata/sample",replace

//Load up the allele scores
use "workingdata/full_sample_allele_scores",clear

//Recomended exclusions
joinby n_eid using  "workingdata/exclusions",unmatched(master)
drop if _m==3
drop _m

//Non-europeans
joinby n_eid using  "workingdata/non-europeans",unmatched(master)
drop if _m==3
drop _m

//Europeans
joinby n_eid using  "workingdata/europeans",unmatched(master)
drop if _m!=3
drop _m

//Relateds
joinby n_eid using   "workingdata/relateds",unmatched(master)
//Generate indicator for relateds
gen relateds=(_m==3)

//Match in family ID (generated in cr_0_within_family_ids.do)
drop _m
joinby n_eid using "workingdata/family_ids",unmatched(master)

//drop if _m==3
drop _m

//PCs
joinby n_eid using "workingdata/PCs",unmatched(master)
drop _m

//Join phenotype data
joinby n_eid using "workingdata/cleaned_biobank_outcomes_ENGLISH.dta",unmatched(master)
drop _m

compress

//Merge in additional phenotypic baseline variables

joinby n_eid using  "workingdata/cov_dist_long",unmatched(master)
drop _m
joinby n_eid using   "workingdata/birth_location_imd_rural_urban",unmatched(master)
rename imd cov_birth_location_imd
drop dist
rename cov_urban cov_birth_location_urban
drop _merge
joinby n_eid using   "workingdata/birth_location",unmatched(master)
drop _m
rename n_129_0_0 cov_birth_location_northing
rename n_130_0_0 cov_birth_location_easting

joinby n_eid using  "workingdata/current_location",unmatched(master)
drop _m
rename n_20074_0_0 cov_current_location_easting
rename n_20075_0_0 cov_current_location_northing

drop data_error
compress


//Final exclusion list
joinby n_eid using "rawdata/exclusions_170726.dta", unmatched(master)
drop if _m==3
drop _m
compress

tab yob, gen(yob_)
gen yob_sex=cov_male*yob
tab yob_sex, gen(yob_sex_)
drop allele_score_1-allele_score_52 XXX
compress
save "workingdata/full_merged_dataset",replace


use "workingdata/full_merged_dataset",clear
joinby n_eid using "workingdata/okbay_hill_snps_clean.dta",unmatched(master) _merge(XX)
tab XX
drop if XX!=3
drop XX
//Restrict to the interim dataset
joinby n_eid using "workingdata/interim.dta",unmatched(master)
drop if _m!=3

//Generate indicator for took fluid intelligence test
gen took_intell=(out_intell!=.)

gen exclude=(took_intell==1|interim==1|born_english==0) 
compress
save "workingdata/analysis_dataset_interim",replace

