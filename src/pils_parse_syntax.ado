*! version 1.0.0  06feb2026
program pils_parse_syntax
	version 16.0

	_on_colon_parse `0'
	local before `s(before)'
	local after `s(after)'

	/* ---------------------------------------------------------- */
	// get touse1 and touse2
	local 0 `before'
	syntax , touse1(passthru) 	///
		touse2(passthru)	

	/* ---------------------------------------------------------- */
	// syntax signature
	local 0 `after'
	syntax anything(equalok name=eq)	///
		using				///
		[if]				///
		[in]				///
		[, extram(passthru)		///
		sdrc(string) 			///
		sdrfv(string) 			///
		q_pir(passthru)			///
		q_pls(passthru)]
	
	/* ---------------------------------------------------------- */
	// create object
	pils_obj create

	/* ---------------------------------------------------------- */
	// parse the main syntax
					// append two samples
	AppendSample `if' `in' `using', `touse1' `touse2' 
					// parse eq
	ParseEq `eq', `extram' `touse1' `touse2'  
					// parse sdrc
	ParseSdrc, `sdrc'
					// parse sdrfv
	ParseSdrFv, `sdrfv'
					// parse q_pir and q_pls
	ParseQ, `q_pir' `q_pls'	
end

					//----------------------------//
					// append two samples
					//----------------------------//
program AppendSample 
	syntax [if] [in] using, ///
		touse1(string)	///
		touse2(string)	

	tempvar src
	qui append `using', generate(`src')

	marksample mytouse
	qui gen byte `touse1' = 0
	qui replace `touse1' = `mytouse' if `src' == 0

	qui gen byte `touse2' = 0
	qui replace `touse2' = `mytouse' if `src' == 1

	pils_obj getname
	mata: `OBJ'._touse1 = `"`touse1'"'
	mata: `OBJ'._touse2 = `"`touse2'"'
end
					//----------------------------//
					// parse equation
					//----------------------------//
program ParseEq
	syntax anything(equalok name=eq)	///
		, touse1(string)		///
		touse2(string)			///
		[extram(varlist numeric fv)]

	_iv_parse `eq'
	local st_y `s(lhs)'
	local st_x1 `s(exog)'
	local st_x2 `s(endog)'
	local st_x3i `s(inst)'

	CheckEqSyntax, depvar(`st_y') x2(`st_x2') x3i(`st_x3i')

	_rmcoll `st_x1' if `touse1', expand
	local st_x1 `r(varlist)'

	_rmcoll `st_x2' if `touse2', expand
	local st_x2 `r(varlist)'

	_rmcoll `st_x3i', expand
	local st_x3i `r(varlist)'

	_rmcoll `extram', expand
	local extram `r(varlist)'
						// get x3e
	local x3 : list st_x3i | extram
	local st_x3e : list x3 - st_x3i
						// depvar
	_fv_check_depvar `st_y'
						// x2
	CheckVarSource `st_x2', touse(`touse2') i(2) 
						// x3i and x3e
	CheckVarSource `st_x3i' `st_x3e', touse(`touse1') i(1) 
	CheckVarSource `st_x3i' `st_x3e', touse(`touse2') i(2) 
						// markout variables
	markout `touse1' `st_y' `st_x1' `st_x3i' `st_x3e' 
	markout `touse2' `st_x3i' `st_x3e' `st_x2'

	sum `touse1' if `touse1', meanonly
	if (r(N) == 0)  error 2000

	sum `touse2' if `touse2', meanonly
	if (r(N) == 0)  error 2000

	SplitFvContinous `st_x2'
	local st_x2c `s(vars_c)'
	local st_x2fv `s(vars_fv)'
	local st_x2fv_raw `s(vars_fv_raw)'

	SplitFvContinous `st_x3i' `st_x3e'
	local st_x3c `s(vars_c)'
	local st_x3fv `s(vars_fv)'

	SplitFvContinous `st_x3i'
	local st_x3i `s(vars_c)' `s(vars_fv)'

	pils_obj getname
	mata: `OBJ'._st_y = `"`st_y'"'
	mata: `OBJ'._st_x1 = `"`st_x1'"'
	mata: `OBJ'._st_x2c = `"`st_x2c'"'
	mata: `OBJ'._st_x2fv = `"`st_x2fv'"'
	mata: `OBJ'._st_x3i = `"`st_x3i'"'
	mata: `OBJ'._st_x3c = `"`st_x3c'"'
	mata: `OBJ'._st_x3fv = `"`st_x3fv'"'
	mata: `OBJ'._st_x2fv_raw = `"`st_x2fv_raw'"'
end
					//----------------------------//
					// check equation syntax
					//----------------------------//
program CheckEqSyntax
	syntax [, depvar(string)	///
		x2(string)		///
		x3i(string)]
	
	if (`"`depvar'"' == "") {
		di as err "must specify the dependent variable"
		exit 198
	}
	
	if (`"`x2'"' == "") {
		di as err "must specify the matched variables {it:varlist_x2}"
		exit 198
	}

	if (`"`x3i'"' == "") {
		di as err "must specify the common variables {it:varlist_x3i}"
		exit 198
	}
end
					//----------------------------//
					// Check vars source
					//----------------------------//
program CheckVarSource
	syntax varlist(numeric fv),	///
		touse(string)		///
		i(string)		

	sum `varlist' if `touse', meanonly
	if (r(N) == 0) {
		di as err "variables {it:`varlist'} are missing "	///
			"in sample `i'"
		exit 2000
	}
end

					//----------------------------//
					// parse sdr
					//----------------------------//
program ParseSdrc
	cap syntax [, sir		///
		pls			///
		pir			///
		slice(integer 10)]
	
	local errmsg di as err "option {bf:sdrc()} allows only one of "	///
			"{bf:sir}, {bf:pls}, or {bf:pir}"
	if _rc {
		`errmsg'
		exit 198
		// NotReached
	}
	
	local sdrc `sir' `pls' `pir'

	local nsdr : list sizeof sdrc
	if (`nsdr' == 0 ) {
		local sdrc sir
	}
	else if (`nsdr' > 1) {
		`errmsg'
		exit 198
		// NotReached
	}

	pils_obj getname
	mata: `OBJ'._sdrc = `"`sdrc'"'
	mata: `OBJ'._sdrc_slice = `slice'
end
					//----------------------------//
					// parse sdr fv
					//----------------------------//
program ParseSdrFv
	cap syntax [, sir		///
		logit			///
		probit			///
		slice(integer 2)]
	
	local errmsg di as err "option {bf:sdrfv()} allows only one of " ///
			"{bf:sir}, {bf:logit}, or {bf:probit}"
	if _rc {
		`errmsg'
		exit 198
		// NotReached
	}
	
	local sdrfv `sir' `logit' `probit'

	local nsdr : list sizeof sdrfv
	if (`nsdr' == 0 ) {
		local sdrfv sir
	}
	else if (`nsdr' > 1) {
		`errmsg'
		exit 198
		// NotReached
	}

	pils_obj getname
	mata: `OBJ'._sdrfv = `"`sdrfv'"'
	mata: `OBJ'._sdrfv_slice = `slice'

end

					//----------------------------//
					// split fv and continous
					//----------------------------//
program SplitFvContinous, sclass
	syntax anything(name=myvars)

	foreach var of local myvars {
		_ms_parse_parts `var'

		if ("`r(type)'"' == "factor") {
			if (r(omit) == 0) {
				local vars_fv `vars_fv' `var'
				local vars_fv_raw `vars_fv_raw' `r(name)'
			}
		}
		else {
			local vars_c `vars_c' `var'
		}
	}

	sret local vars_fv `vars_fv'	
	sret local vars_fv_raw `vars_fv_raw'	
	sret local vars_c `vars_c'
end
					//----------------------------//
					// parse Q
					//----------------------------//
program ParseQ
	syntax [, q_pls(integer 0)	///
		q_pir(integer 0) ]
	

	pils_obj getname
	mata: `OBJ'._q_pls = `q_pls'
	mata: `OBJ'._q_pir = `q_pir'
end
