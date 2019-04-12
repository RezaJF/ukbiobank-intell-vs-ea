//Neil Davies 12/03/2019
//This simulates a multivariable MR study with measurement error on the exposures

clear
set obs 250000

//Generate errors
gen v=rnormal()
gen u=rnormal()
gen e=rnormal()
gen c=rnormal()

//Generate correlated instruments
gen rg=rnormal()
gen z_x1=0.3*rnormal()+0.7*rg
gen z_x2=0.3*rnormal()+0.7*rg

//Generate the exposures
gen x1=z_x1+c+v
gen x2=z_x2+c+u

//Generate the outcome
gen y=x1+x2+c+e

//Run IV regression
ivreg2 y (x1 x2 =z*),ro

//Repeat with measurement error
//Generate measurement error
gen me1=rnormal()
gen me2=rnormal()

gen x1_measured=x1+me1
gen x2_measured=x2+me1

//Run IV regression
ivreg2 y (x1_measured x2_measured =z*),ro
//Remains unbiased

//Run IV regression with one mismeasured and one measured accurately variable
ivreg2 y (x1 x2_measured =z*),ro


