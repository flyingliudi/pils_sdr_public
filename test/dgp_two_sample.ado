*! version 1.0.0  01feb2026
/*
DGP for two sample data
=======================

* Sample A

	- We observe {Y, x1, x3}
	- Y = b0 + x1'b1 + x2'b2 + x_{3I}'b3 + u
	- x2 are missing
	- x3 = (x_{3I}', x_{3E}')' are use to match x2

* Sample B

	- We observe {x2, x3}
*/

program dgp_two_sample, sclass
	syntax [, b0(real 1)		///
		b1(real 1)		///
		b2(real 1)		///
		b3(real 1)		///
		d1(int 1)		///
		d2(int 2)		///
		x3vars(string)		///
		d3_cont(int 6)		///
		d3_fv(int 4)		///
		rho_main(real 0.4)	///
		rho_aux(real 0.4)	///
		rho_fv_main(real 0.4)	///
		rho_fv_aux(real 0.4)	///
		nobs_main(int 1000)	///
		nobs_aux(int 2000)	///
		dta_main(string)	///
		dta_aux(string)]
	
	if (`"`x3vars'"' == "") {
		local x3vars x3_c1 x3_fv1
	}

	if (`"`dta_main'"' == "") {
		local dta_main sample_main
	}

	if (`"`dta_aux'"' == "") {
		local dta_aux sample_aux
	}

	/* ---------------------------------------------------------- */
	// create frame for two dataset
	tempname main aux
	frame create `main'
	frame create `aux'

	/* ---------------------------------------------------------- */
	// DGP in two sample
	foreach fr in main aux {
		frame ``fr'' {
			GenOneSample, 			///
				nobs(`nobs_`fr'')	///
				rho(`rho_`fr'')		///
				d3_cont(`d3_cont')	///
				rho_fv(`rho_fv_`fr'')	///
				d3_fv(`d3_fv')		///
				d2(`d2')		///
				x3vars(`x3vars')	///
				d1(`d1')		///
				b0(`b0')		///
				b1(`b1')		///
				b2(`b2')		///
				b3(`b3')
		}
	}

	/* ---------------------------------------------------------- */
	// save main sample
	frame `main' {
		save `dta_main', replace
	}

	/* ---------------------------------------------------------- */
	// save auxiliary sample
	frame `aux' {
		drop y x1*
		save `dta_aux', replace
	}

	sret local dta_main `dta_main'
	sret local dta_aux `dta_aux'
end

					//----------------------------//
					// generate one sample
					//----------------------------//
program GenOneSample

	syntax, nobs(string)	///
		rho(string)	///
		d3_cont(string)	///
		rho_fv(string)	///
		d3_fv(string)	///
		d2(string)	///
		x3vars(string)	///
		d1(string)	///
		b0(string)	///
		b1(string)	///
		b2(string)	///
		b3(string)
		
				// set obs
	clear
	set obs `nobs'
				// X3 continous
	GenX3Continous, 	///
		nobs(`nobs') 	///
		rho(`rho') 	///
		d3_cont(`d3_cont')
	local x3list_cont `s(x3list_cont)'
	
				// X3 binary
	GenX3Binary,			///
		nobs(`nobs')		///
		rho_fv(`rho_fv')	///
		d3_fv(`d3_fv')
	local x3list_fv `s(x3list_fv)'
	
				// true index
	GenTrueIndex,				///
		d2(`d2')			///	
		x3list_cont(`x3list_cont')	///
		x3list_fv(`x3list_fv')	
	local zlist `s(zlist)'
		
				// x2
	GenX2, zlist(`zlist')	///
		d2(`d2')
	local x2list `s(x2list)'
	
				// x1
	GenX1, x3vars(`x3vars')	///
		d1(`d1')
	local x1list `s(x1list)'
	
				// y
	GenY, x1list(`x1list')		///
		x2list(`x2list')	///
		x3list(`x3vars')	///
		b0(`b0')		///
		b1(`b1')		///
		b2(`b2')		///
		b3(`b3')	
end
					//----------------------------//
					// generate continous x3
					//----------------------------//
program GenX3Continous, sclass
	syntax, 		///
		nobs(string)	///
		rho(string)	///
		d3_cont(string)
	
	qui forvalues i=1/`d3_cont' {
		gen double x3_c`i' = .
		local vlist `vlist' x3_c`i'
	}

	mata: draw_mvnormal_trunc(	///
		`"`vlist'"',		///
		`nobs',			///
		`d3_cont', 		///
		`rho')
	
	sret local x3list_cont `vlist'
end

					//----------------------------//
					// generate binary x3
					//----------------------------//
program GenX3Binary, sclass
	syntax, nobs(string)	///
		rho_fv(string)	///
		d3_fv(string)
	
	if (`d3_fv' == 0) {
		exit
		// NotReached
	}

	qui forvalues i=1/`d3_fv' {
		gen byte x3_fv`i' = .
		local vlist `vlist' x3_fv`i'
	}

	mata: draw_rademacher(	///
		"`vlist'",	///
		`nobs',		///
		`rho_fv',	///
		`d3_fv')
	
	sret local x3list_fv `vlist'
end

					//----------------------------//
					// generate true index
					//----------------------------//
program GenTrueIndex, sclas
	syntax, d2(string)		///
		x3list_cont(string)	///
		[x3list_fv(string)]
	
	qui forvalues i=1/`d2' {
		gen double z`i' = .
		local vlist `vlist' z`i'
	}

	mata: get_z(			///
		`"`vlist'"', 		///
		`d2', 			///
		`"`x3list_cont'"', 	///
		`"`x3list_fv'"')
	
	sret local zlist `vlist'
end
					//----------------------------//
					// generate X1
					//----------------------------//
program GenX1, sclass
	syntax, x3vars(string)	///
		d1(string)

	qui forvalues i=1/`d1' {
		tempvar eta
		gen double x1_`i' = rnormal()

		foreach var of local x3vars {
			replace x1_`i' = x1_`i' + `var'
		}
		local vlist `vlist' x1_`i'
	}

	sret local x1list `vlist'
end
					//----------------------------//
					// generate X2
					//----------------------------//
program GenX2, sclass
	syntax, zlist(string)	///
		d2(string)

	local k_c = 1
	local k_fv = 1
	qui forvalues i=1/`d2' {
		gen double x2_`i' = .

		local zvar : word `i' of `zlist'

		if (`i' <= ceil(`d2'/2)) {
			replace x2_`i' = `zvar' + 9*normalden(`zvar')	///
			 	+ rnormal()
			rename x2_`i'  x2_c`k_c'
			local vlist `vlist' x2_c`k_c'
			local k_c = `k_c' + 1
		}
		else {
			replace x2_`i' = `zvar' + rnormal() > 0
			rename x2_`i' x2_fv`k_fv'
			local vlist `vlist' x2_fv`k_fv'
			local k_fv = `k_fv' + 1
		}
	}

	sret local x2list `vlist'
end
					//----------------------------//
					// gen Y
					//----------------------------//
program GenY
	syntax, x1list(string)	///
		x2list(string)	///
		x3list(string)	///
		b0(string)	///
		b1(string)	///
		b2(string)	///
		b3(string)	
	
	gen double y = rnormal() + `b0'

	qui forvalues i=1/3 {
		foreach var in `x`i'list' {
			replace y = y + `var'*`b`i''
		}
	}
end


/*-----------------------------------------------------------------------------
	mata utility	
-----------------------------------------------------------------------------*/
mata:
mata set matastrict on

					//----------------------------//
					// get V matrix
					//----------------------------//
real matrix get_V(			///
	real scalar	d3_cont,	///
	real scalar	rho)
{
	real matrix	V
	real scalar	i, j

	V = J(d3_cont, d3_cont, .)

	for (i=1; i<=d3_cont; i++) {
		for (j=1; j<=d3_cont; j++) {
			V[i, j] = rho^(abs(i-j))
		}
	}

	return(V)
}

					//----------------------------//
					// draw truncated mvnormal 
					//----------------------------//
void draw_mvnormal_trunc(		///
	string scalar	vlist,		///
	real scalar	nobs,		///
	real scalar	d3_cont,	///
	real scalar	rho)
{
	real matrix 	V, L, X3
	real scalar	i
	real rowvector	z, x

	V = get_V(d3_cont, rho)
	L = cholesky(V)
	Vinv = invsym(V)
	X3 = J(nobs, d3_cont, .)

	i = 1
	while (i<=nobs) {
		z = rnormal(1, d3_cont, 0, 1)
		x = z*L

		tmp = x*Vinv*x'
		if (tmp <= 9) {
			X3[i, .] = x
			i++
		}
	}


	st_store(., tokens(vlist), X3)
}
					//----------------------------//
					// draw rademacher
					//----------------------------//
void draw_rademacher(		///
	string scalar	vlist,	///	
	real scalar	nobs,	///
	real scalar	rho,	///
	real scalar	d3_fv)
{
	real matrix	V, L, X3
	real scalar	i
	real rowvector	z, x

	V = get_V(d3_fv, rho)
	L = cholesky(V)
	X3 = J(nobs, d3_fv, .)

	for (i=1; i<=nobs; i++) {
		z = rnormal(1, d3_fv, 0, 1)
		x = z*L:>=0
		X3[i, .] = x
	}

	st_store(., tokens(vlist), X3)
}
					//----------------------------//
					// get z
					//----------------------------//
void get_z(				///
	string scalar	vlist,		///
	real scalar	d2,		///
	string scalar	st_x3list_cont,	///
	string scalar	st_x3list_fv)
{
	real matrix	X3_cont, X3_fv, Z
	real colvector	theta_c, theta_f, z
	real scalar	i

	st_view(X3_cont, ., st_x3list_cont)

	if (st_x3list_fv != "") {
		st_view(X3_fv, ., st_x3list_fv)
	}
	else {
		X3_fv = 0
	}

	theta_c = J(cols(X3_cont), 1, .)
	for (i=1; i<=length(theta_c); i++) {
		if (mod(i, 3) == 1) {
			theta_c[i] = 1
		}
		else if (mod(i, 3) == 2) {
			theta_c[i] = 0
		}
		else {
			theta_c[i] = -1
		}
	}

	theta_f = J(cols(X3_fv), 1, .)
	for (i = 1; i<= length(theta_f); i++) {
		if (mod(i, 3) == 1) {
			theta_f[i] = 0.5
		}
		else if (mod(i, 3) == 2) {
			theta_f[i] = 0
		}
		else {
			theta_f[i] = -0.5
		}
	}

	Z = J(st_nobs(), d2, .)
	for (i=1; i<=d2; i++) {
		if (mod(i, 2) == 0) {
			z = X3_cont*theta_c :- X3_fv*theta_f
		}
		else {
			z = X3_cont*theta_c :+ X3_fv*theta_f
		}

		Z[., i] = z
	}

	st_store(., tokens(vlist), Z)
}


end

