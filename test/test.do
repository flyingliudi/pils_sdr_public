cscript

do compile

local seed 12345671
set seed `seed'
set sortseed `seed'

cap log close

log using tmp, text replace


dgp_two_sample, 		///
	nobs_main(2000) 	///
	nobs_aux(1000)		///
	d3_fv(4)		///
	x3vars(x3_c1 x3_fv1)	///
	d2(2)			///
	dta_main(sample_main)	///
	dta_aux(sample_aux)

use sample_main

pils y x1_* (x2_c* i.(x2_fv*)= x3_c1 i.x3_fv1) 		///
 	using sample_aux 				///
 	, extram(x3_c* i.(x3_fv*)) 			///
 	sdrc(sir) 					///
 	sdrfv(logit)


log close
