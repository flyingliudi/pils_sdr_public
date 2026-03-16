cscript

do compile

local seed 12345671
set seed `seed'
set sortseed `seed'

dgp_two_sample, 		///
	nobs_main(1000) 	///
	nobs_aux(1000)		///
	d3_fv(4)		///
	x3vars(x3_c1 x3_fv1)	///
	d2(2)			///
	dta_main(sample_main)	///
	dta_aux(sample_aux)

use sample_main

local spec y x1_* (x2_c* i.(x2_fv*)= x3_c1 i.x3_fv1) using sample_aux 			
local opts extram(x3_c* i.(x3_fv*)) 			

foreach sdrc in sir pls pir {
	foreach sdrfv in sir logit probit{
		pils `spec', 		///
			`opts'		///
			sdrc(`sdrc') 	///
			sdrfv(`sdrfv')

//mkassert eclass, saving(true/true_`sdrc'_`sdrfv'.do, replace)
		
		qui do true/true_`sdrc'_`sdrfv'
		di "ok: `sdrc' + `sdrfv'
	}
}


