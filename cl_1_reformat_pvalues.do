//Neil Davies 12/03/19
//This cleans the p-values 


cap prog drop pvalue_conv_3
prog def pvalue_conv_3
args pvalue

gen pvalue_conv=""

replace pvalue=substr(`pvalue',1,3)+"x10"+substr(`pvalue',-3,3) if real(substr(`pvalue',-2,2))>3
replace pvalue="0.00"+string(round(real(substr(`pvalue',1,3)))) if real(substr(`pvalue',-2,2))==3
replace pvalue="0.0"+string(round(real(substr(`pvalue',1,3)))) if real(substr(`pvalue',-2,2))==2
replace pvalue="0."+string(round(10*real(substr(`pvalue',1,4)))) if real(substr(`pvalue',-2,2))==1
replace pvalue="1.00" if real(substr(`pvalue',-2,2))==0

end


cap prog drop pvalue_conv_2
prog def pvalue_conv_2
args pvalue

gen pvalue_conv=""
replace pvalue=substr(`pvalue',1,3)+"x10"+substr(`pvalue',-3,3) if real(substr(`pvalue',-2,2))>3
replace pvalue="0.00"+string(round(real(substr(`pvalue',1,3)))) if real(substr(`pvalue',-2,2))==3
replace pvalue="0.0"+string(round(real(substr(`pvalue',1,3)))) if real(substr(`pvalue',-2,2))==2
replace pvalue="0."+string(round(10*real(substr(`pvalue',1,4)))) if real(substr(`pvalue',-2,2))==1
replace pvalue="1.00" if real(substr(`pvalue',-2,2))==0

end
