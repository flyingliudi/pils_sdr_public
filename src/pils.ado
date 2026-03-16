*! version 1.0.0  06feb2026
/*
	plugin least squares using sufficient dimension reduction
*/
program pils
	version 16.0

	cap noi Compute `0'
	local rc = _rc

	SafeExit, rc(`rc')
end

					//----------------------------//
					// Compute
					//----------------------------//
program Compute

	/* ---------------------------------------------------------- */	
	// preserve start
	preserve

	/* ---------------------------------------------------------- */
	// parse syntax
	tempvar touse1 touse2
 	pils_parse_syntax, touse1(`touse1') touse2(`touse2') : `0'

	/* ---------------------------------------------------------- */
	// compute
	pils_obj getname
	mata: `OBJ'.compute()

	/* ---------------------------------------------------------- */
	// display
	Display

	/* ---------------------------------------------------------- */
	// preserve end
	restore
end

					//----------------------------//
					// safe exit
					//----------------------------//
program SafeExit
	syntax, rc(string)

	if `rc' {
		nobreak pils_obj rm
		exit `rc'
		// NotReached
	}

	pils_obj rm
end

					//----------------------------//
					// Display
					//----------------------------//
program Display

					// header
	Header
					// coefficient table
	_coef_table
end

					//----------------------------//
					// Header
					//----------------------------//
program Header
					// title
	local title as txt `"`e(title)'"'

	local col1 = 39
	local col2 = 69
					// Number of obs
	local n_obs1 as txt _col(`col1') `"Number of obs in sample 1"'  ///
		_col(67) "=" _col(`col2') as res %10.0fc e(main_nobs)

	local n_obs2 as txt _col(`col1') 		///
		`"Number of obs in sample 2"'	///
		_col(67) "=" _col(`col2') as res %10.0fc e(aux_nobs)
	
					// k matched 
	local k_match as txt _col(`col1') "Number of matched vars." ///
		_col(67) "=" _col(`col2') as res %10.0fc e(k_matched)

					// k_overlap
	local k_overlap as txt _col(`col1') 	///
		"Number of overlapping vars." ///
		_col(67) "=" _col(`col2') as res %10.0fc e(k_overlap)

					// sdrc
	local sdrc as txt "Continuous SDR:" _col(17) "{bf:`e(sdrc)'}"

					// sdrfv
	local sdrfv as txt "Binary SDR:" _col(17) "{bf:`e(sdrfv)'}"


	di
	di `title' `n_obs1'
	di `n_obs2'
	di `sdrfv' `k_match'
	di `sdrc' `k_overlap'
	di
end
