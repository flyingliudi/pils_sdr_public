*! version 1.0.0  11feb2026
*! DoNotDistribute

findfile "pils_macdefs.matah", path(".")
include `"`r(fn)'"'

version 16.0

mata:
mata set matastrict on

/*-----------------------------------------------------------------------------
	class for pils-sdr estimator				
-----------------------------------------------------------------------------*/
class `PILS' {

public:
	string scalar	_touse1
	string scalar	_touse2

	string scalar	_st_y

	string scalar	_st_x1

	string scalar	_st_x2c
	string scalar	_st_x2fv
	string scalar	_st_x2fv_raw
	string vector	_st_x2fv_vec

	string scalar	_st_x3i

	string scalar	_st_x3c
	string scalar	_st_x3fv
	string scalar	_st_x3

	string scalar	_sdrc
	string scalar	_sdrfv

	real scalar	_sdrc_slice
	real scalar	_sdrfv_slice

	real scalar	_q_pir
	real scalar	_q_pls

	void		compute()

protected:

	void		set_data()
	void		sdr_aux()
	real colvector	sdr_c()
	real colvector	sdr_fv()
	real colvector	sir()
	real colvector	sir_fv()
	real colvector	pls()
	real colvector	pir()
	real matrix	get_Xi()
	real colvector	run_logit_probit()
	real matrix	normalize_theta()
	void		sdr_predict_main()
	void		kernel_link_aux()
	real colvector	kernel_func()
	real colvector	compute_glink()
	real colvector	compute_glink_loo()
	real scalar	get_kernel_bandwidth()
	void		run_regress()
	real colvector	run_regress_wrk()
	void		post_result()
	void		compute_variance()
	real matrix	compute_Omega2()
	real matrix	compute_phi_21()
	real matrix	compute_phi_21_l()
	real matrix	compute_phi_22()
	real matrix	compute_phi_22_l()
	real matrix	compute_psi_index()
	real matrix	compute_psi_logit()
	real matrix	compute_psi_probit()
	void		load_st_sdr()

	real colvector	_main_y
	real matrix	_main_x1
	real matrix	_main_x3
	real matrix	_main_x3i
	real scalar	_main_nobs
	real scalar	_main_ones
	real matrix	_main_x2c_z
	real matrix	_main_x2fv_z
	real matrix	_main_x2_z
	real matrix	_main_x2c_glink
	real matrix	_main_x2fv_glink
	real matrix	_main_x2_glink
	real rowvector	_main_x3_bar
	real matrix	_main_Xhat
	real rowvector	_main_x3i_bar
	string scalar	_st_x3_demean

	real matrix	_aux_x2c
	real matrix	_aux_x2fv
	real matrix	_aux_x2
	real matrix	_aux_x3
	real scalar	_aux_nobs
	real scalar	_aux_ones
	real matrix	_aux_x2c_z
	real matrix	_aux_x2fv_z
	real matrix	_aux_x2_z
	real rowvector	_aux_x3_bar

	real matrix	_theta_sdr_x2c
	real matrix	_theta_sdr_x2fv

	real colvector	_beta
	real colvector	_beta_x2
	string matrix	_cs
	real matrix	_V

	real scalar	_d2c
	real scalar	_d2fv
	real scalar	_d2
	real scalar	_d1
	real scalar	_d3i
	real scalar	_d3
}
					//----------------------------//
					// compute
					//----------------------------//
void `PILS'::compute()
{
	set_data()

	sdr_aux()

	sdr_predict_main()

	kernel_link_aux()

	run_regress()

	compute_variance()

	_beta[length(_beta)] = _beta[length(_beta)]	///
		- _main_x3i_bar*_beta[_d1 + _d2 + 1.._d1+_d2+_d3i]

	post_result()
}

					//----------------------------//
					// set data
					//----------------------------//
void `PILS'::set_data()
{
	_st_x3 = _st_x3c + " " + _st_x3fv
	_st_x2fv_vec = tokens(_st_x2fv_raw)

	/* ---------------------------------------------------------- */	
	// main sample data
	_main_y = st_data(., _st_y, _touse1)
	_main_x1 = st_data(., _st_x1, _touse1)
	_main_nobs = rows(_main_y)
	_main_x3 = st_data(., _st_x3, _touse1)
	_main_ones = J(_main_nobs, 1, 1)

	// demean x3 in the main sample
	_main_x3_bar = mean(_main_x3)
	_main_x3 = _main_x3 :- _main_x3_bar

	// extract x3i before demean
	_main_x3i = st_data(., _st_x3i, _touse1)
	_main_x3i_bar = mean(_main_x3i)
	_main_x3i = _main_x3i :- _main_x3i_bar

	real rowvector	idx_x3_demean
	_st_x3_demean = st_tempname(cols(_main_x3))
	idx_x3_demean = st_addvar("double", _st_x3_demean)
	st_store(., idx_x3_demean, _touse1, _main_x3)
	_st_x3_demean = invtokens(_st_x3_demean)

	/* ---------------------------------------------------------- */
	// auxiliary sample
	_aux_x3 = st_data(., _st_x3, _touse2)
	_aux_x2c = st_data(., _st_x2c, _touse2)
	if (_st_x2fv != "") {
		_aux_x2fv = st_data(., _st_x2fv, _touse2)
	}
	else {
		_aux_x2fv = J(rows(_aux_x3),0, .)
	}
	_aux_x2 = (_aux_x2c, _aux_x2fv)
	_aux_nobs = rows(_aux_x3)
	_aux_ones = J(_aux_nobs, 1, 1)

	// demean x3 in auxiliary sample
	_aux_x3_bar = mean(_aux_x3)
	_aux_x3 = _aux_x3 :- _aux_x3_bar
	st_store(., idx_x3_demean, _touse2, _aux_x3)

	/* ---------------------------------------------------------- */
	// sdr index coefficient matrix
	_theta_sdr_x2c = J(cols(_aux_x3), cols(_aux_x2c), .)
	if (_sdrfv == "logit" | _sdrfv == "probit") {
		_theta_sdr_x2fv = J(cols(_aux_x3)+1, cols(_aux_x2fv), .)
	}
	else {
		_theta_sdr_x2fv = J(cols(_aux_x3), cols(_aux_x2fv), .)
	}

	/* ---------------------------------------------------------- */
	// dimension of parameters
	_d2c = cols(_aux_x2c)
	_d2fv = cols(_aux_x2fv)
	_d2 = _d2c + _d2fv
	_d1 = cols(_main_x1)
	_d3i = cols(_main_x3i)
	_d3 = cols(_main_x3)

	/* ---------------------------------------------------------- */
	// predicted index
	_main_x2c_z = J(_main_nobs, _d2c, .)	
	_main_x2fv_z = J(_main_nobs, _d2fv, .)	

	_aux_x2c_z = J(_aux_nobs, _d2c, .)	
	_aux_x2fv_z = J(_aux_nobs, _d2fv, .)	

	/* ---------------------------------------------------------- */
	// predicted glink
	_main_x2c_glink = J(_main_nobs, _d2c, .)	
	_main_x2fv_glink = J(_main_nobs, _d2fv, .)	
}

					//----------------------------//
					// compute 
					//----------------------------//
void `PILS'::sdr_aux()
{
	real scalar	i	

	for (i=1; i<=_d2c;  i++) {
		_theta_sdr_x2c[., i] = sdr_c(_aux_x2c[., i],  _aux_x3)
	}

	for (i=1; i<=_d2fv;  i++) {
		_theta_sdr_x2fv[., i] = sdr_fv(	///
			_aux_x2fv[., i],  	///
			_aux_x3,		///
			_st_x2fv_vec[i])
	}
}

					//----------------------------//
					// sdr for continous variables
					//----------------------------//
real colvector `PILS'::sdr_c(	///
	real colvector	x2,	///
	real matrix	x3)
{
	real matrix	myx3
	real colvector	theta

	myx3 = x3:-mean(x3)

	if (_sdrc == "sir") {
		theta = sir(x2, myx3, _sdrc_slice)	
	}
	else if (_sdrc == "pls") {
		theta = pls(x2, myx3, ., _q_pls)
	}
	else if (_sdrc == "pir") {
		theta = pir(x2, myx3, _sdrc_slice, _q_pir)
	}

	theta = normalize_theta(theta)

	return(theta)
}

					//----------------------------//
					// sliced inverse regression
					//----------------------------//
real colvector `PILS'::sir(	///
	real colvector	x2,	///
	real matrix	x3,	///
	real scalar	H,	///
	|real colvector	Xi_xy)
{
	real matrix	X, x3bar_H, Sigma, Sigma_H, tmp, V
	real scalar	k, s0, s1, nobs_k, nobs, p, mk, d3
	real rowvector	x3bar_h, x3bar
	real colvector	mk_H, theta, w, b

	X = (x2, x3)
	nobs = rows(X)
	p = cols(X)
	d3 = p - 1

	_sort(X, 1)
	
	nobs_k = floor(nobs/H)
	x3bar_H = J(H, d3, .)
	mk_H = J(H, 1, .)

	/* ---------------------------------------------------------- */
	// get mean of x3 in each slice
	for (k=1; k<=H; k++) {
		s0 = (k-1)*nobs_k + 1
		s1 = min((k*nobs_k, nobs))

		nobs_k = s1 - s0 + 1

		x3bar_h = colsum(X[s0..s1, 2..p])/nobs_k

		x3bar_H[k, .] = x3bar_h
		mk_H[k] = nobs_k
	}

	/* ---------------------------------------------------------- */
	// compute Sigma_H
	Sigma_H = J(d3, d3, 0)
	x3bar = mean(X[., 2..p])

	for (k=1; k<= H; k++) {
		tmp = x3bar_H[k, .] :- x3bar
		tmp = mk_H[k]*(tmp'*tmp)

		Sigma_H = Sigma_H + tmp
	}
	Sigma_H = Sigma_H/nobs

	tmp = x3:-x3bar
	Sigma = cross(tmp, tmp)/nobs

	eigensystem(invsym(Sigma)*Sigma_H, V=., .)

	theta = Re(V[., 1])

	Xi_xy = Sigma*theta

	return(theta)
}

					//----------------------------//
					// partial least squares
					//----------------------------//
real colvector `PILS'::pls(		///
	real colvector	x2,		///
	real matrix	x3,		///
	|real colvector	mysigma_l,	///
	real scalar	myq)
{
	real matrix	Sigma_l, Sigma, Xi
	real rowvector	x3bar, x2bar
	real scalar	nobs
	real colvector	theta

	x3bar = mean(x3)
	x2bar = mean(x2)
	nobs = rows(x2)

	if ((args() == 2) | (mysigma_l == .)) {
		Sigma_l = cross(x3:-x3bar, x2 :- x2bar)/nobs
	}
	else {
		Sigma_l = mysigma_l
	}
	Sigma = cross(x3:-x3bar, x3:-x3bar)/nobs

	Xi = get_Xi(Sigma, Sigma_l, myq)

	theta = Xi*invsym(Xi'*Sigma*Xi)*Xi'*Sigma_l
	return(theta)
}

					//----------------------------//
					// get Xi matrix
					//----------------------------//
real matrix `PILS'::get_Xi(		///
	real matrix	Sigma,		///
	real colvector	Sigma_l,	///
	real scalar	myq)
{
	real scalar	d3, i, q, a, tol
	real matrix	Xi, tmp
	real rowvector	L

	d3 = rows(Sigma)
	Xi = Sigma_l

	tmp = 1
	for (i=2; i<= d3; i++) {
		tmp = Sigma*tmp
		Xi = (Xi, tmp*Sigma_l)
	}

	eigensystem(Xi*(Xi'), ., L=.)

	L = Re(L)
	q = 0
	a = 1.5
	tol = 1e-8
	for (i=1; i<= d3-1; i++) {
		if (L[i]/L[i+1] > a & L[i+1] > tol) {
			q = q + 1
		}
	}

	if (myq > 0) {
		Xi = Xi[., 1..myq]
	}
	else {
		Xi = Xi[., 1..q]
	}


	return(Xi)
}

					//----------------------------//
					// partial inverse regression
					//----------------------------//
real colvector	`PILS'::pir(	///
	real colvector	x2,	///
	real matrix	x3,	///
	real scalar	H,	///
	real scalar	myq)
{
	real colvector	Xi_xy, theta 
	real matrix	Sigma

	(void) sir(x2, x3, H, Xi_xy)
	theta = pls(x2, x3, Xi_xy, myq)

	return(theta)
}
					//----------------------------//
					// sdr for factor variable
					//----------------------------//
real colvector `PILS'::sdr_fv(	///
	real colvector	x2,	///
	real matrix	x3,	///
	|string scalar	st_x2)
{
	real matrix	myx3
	real colvector	theta

	myx3 = x3:-mean(x3)

	if (_sdrfv == "sir") {
		theta = sir_fv(x2, myx3, _sdrfv_slice)
		theta = normalize_theta(theta)
	}
	else if (_sdrfv == "logit" | _sdrfv == "probit" ) {
		theta = run_logit_probit(_sdrfv, st_x2, _st_x3_demean, _touse2)
	}	

	return(theta)
}

					//----------------------------//
					// ran logit
					//----------------------------//
real colvector `PILS'::run_logit_probit(	///
	string scalar	cmd,			///
	string scalar	st_x2,			///
	string scalar	st_x3,			///
	string scalar	touse)
{
	string scalar	ss
	real colvector	theta

	ss = sprintf("qui %s  %s %s if %s", cmd, st_x2, st_x3, touse)
	stata(ss)

	theta = st_matrix("e(b)")'
	return(theta)
}

					//----------------------------//
					// normalized theta make theta[1] = 1
					//----------------------------//
real matrix `PILS'::normalize_theta(
	real matrix	theta)
{
	return(theta:/theta[1, .])
}

					//----------------------------//
					// sdr predict in the main sample
					//----------------------------//
void `PILS'::sdr_predict_main()
{
	real scalar	i
	real matrix	tmp_main, tmp_aux


	for (i=1; i<= _d2c; i++) {
		_aux_x2c_z[., i] = _aux_x3*_theta_sdr_x2c[., i]
		_main_x2c_z[., i] = _main_x3*_theta_sdr_x2c[., i]
	}

	if (_sdrfv == "logit" | _sdrfv == "probit") {
		tmp_main = (_main_x3, _main_ones)
		tmp_aux = (_aux_x3, _aux_ones)
	}
	else {
		tmp_main = _main_x3
		tmp_aux = _aux_x3
	}

	for (i=1; i<= _d2fv; i++) {
		_main_x2fv_z[., i] = tmp_main*_theta_sdr_x2fv[., i]
		_aux_x2fv_z[., i] = tmp_aux*_theta_sdr_x2fv[., i]
	}

	_aux_x2_z = (_aux_x2c_z, _aux_x2fv_z)
	_main_x2_z = (_main_x2c_z, _main_x2fv_z)
}

					//----------------------------//
					// kernel link aux
					//----------------------------//
void `PILS'::kernel_link_aux()
{
	real scalar	i

	for (i=1; i<= _d2c; i++) {
		_main_x2c_glink[., i] = compute_glink(	///
			_aux_x2c[., i],			///
			_aux_x2c_z[.,i], 		///
			_main_x2c_z[., i])
	}

	for (i=1; i<= _d2fv; i++) {
		_main_x2fv_glink[., i] = compute_glink(	///
			_aux_x2fv[., i],		///
			_aux_x2fv_z[.,i], 		///
			_main_x2fv_z[., i])
	}

	_main_x2_glink = (_main_x2c_glink, _main_x2fv_glink)
}

					//----------------------------//
					// kernel function
					//----------------------------//
real colvector	`PILS'::kernel_func(
	real colvector	u)
{
	real colvector	ku

	ku = 2*exp(-u:^2/2)/sqrt(2*pi()) - exp(-u:^2/4)/sqrt(4*pi())

	return(ku)
}

					//----------------------------//
					// compute glink
					//----------------------------//
real colvector	`PILS'::compute_glink(	///
	real colvector	aux_x2,		///
	real colvector	aux_z,		///
	real colvector	main_z)
{
	real scalar	h, i
	real colvector	glink, u, ku

	h = get_kernel_bandwidth(aux_z, rows(aux_z))	

	glink = J(rows(main_z), 1, .)

	for (i=1; i<= rows(main_z); i++) {
		u = (aux_z:-main_z[i])/h
		ku = kernel_func(u)
		glink[i] = mean(aux_x2:*ku)/mean(ku)
	}

	return(glink)
}

					//----------------------------//
					// get kernel bandwidth
					//----------------------------//
real scalar `PILS'::get_kernel_bandwidth(	///
	real colvector	z,			///
	real scalar	m)
{
	real scalar	h

	h = 0.5*sqrt(variance(z))*(log(m)/m)^(3/20)

	return(h)
}

					//----------------------------//
					// run regress
					//----------------------------//
void `PILS'::run_regress()
{
	real colvector	y

	_main_Xhat = ( _main_x1, 	///
		_main_x2c_glink, 	///
		_main_x2fv_glink, 	///
		_main_x3i,		///
		_main_ones)
	
	_cs = tokens(_st_x1+ " " 	///
		+ _st_x2c + " " 	///
		+ _st_x2fv + " "	///
		+ _st_x3i + " "		///
		+ "_cons")'
	
	_cs = J(rows(_cs), 1, " "), _cs
	
	_beta = run_regress_wrk(_main_y, _main_Xhat)
	_beta_x2 = _beta[cols(_main_x1)+1..cols(_main_x1)+_d2]
}

					//----------------------------//
					// run regress wrk
					//----------------------------//
real colvector	`PILS'::run_regress_wrk(	///
	real colvector	y,			///
	real matrix	X)
{
	real colvector	b
	real matrix XpX, Xpy

	XpX = quadcross(X, X)
	Xpy = quadcross(X, y)

	b = invsym(XpX)*Xpy

	return(b)
}

					//----------------------------//
					// post result
					//----------------------------//
void `PILS'::post_result()
{
	string scalar	bs, Vs, ss

	/* ---------------------------------------------------------- */
	// b and V
	bs = st_tempname()

	st_matrix(bs, _beta')
	st_matrixcolstripe(bs, _cs)

	Vs = st_tempname()
	st_matrix(Vs, _V)
	st_matrixcolstripe(Vs, _cs)
	st_matrixrowstripe(Vs, _cs)

	ss = sprintf("ereturn post %s %s", bs, Vs)
	stata(ss)

	/* ---------------------------------------------------------- */
	// scalars
					// main nobs
	st_numscalar("e(main_nobs)", _main_nobs)
					// aux nobs
	st_numscalar("e(aux_nobs)", _aux_nobs)
					// k_matches
	st_numscalar("e(k_matched)", _d2)
					// k_matches
	st_numscalar("e(k_overlap)", _d3)
	/* ---------------------------------------------------------- */
	// macros
					// sdrc
	st_global("e(sdrc)", strupper(_sdrc))
					// sdrfv
	if (_sdrfv == "sir") {
		_sdrfv = "SIR"
	}
	else {
		_sdrfv = strproper(_sdrfv)
	}
	st_global("e(sdrfv)", _sdrfv)
					// title
	st_global("e(title)", "Two sample regression with SDR")
					// cmd
	st_global("e(cmd)", "pils")	


}

					//----------------------------//
					// compute variance
					//----------------------------//
void `PILS'::compute_variance()
{
	real matrix	Qx, Qx_inv, Omega_1, Omega_2
	real colvector	res


	Qx = cross(_main_Xhat, _main_Xhat)/_main_nobs
	Qx_inv = invsym(Qx)

	res = _main_y - _main_Xhat*_beta
	Omega_1 = cross(_main_Xhat:*res, _main_Xhat:*res)/_main_nobs

	Omega_2 = compute_Omega2()


	_V = Qx_inv*(Omega_1 + _main_nobs/_aux_nobs*Omega_2)*Qx_inv/_main_nobs

}

					//----------------------------//
					// compute Omega2
					//----------------------------//
real matrix `PILS'::compute_Omega2()
{
	real matrix	phi_21, phi_22, phi
	real matrix	Omega_2

	phi_21 = compute_phi_21()
	phi_22 = compute_phi_22()
	phi = phi_21 + phi_22
	Omega_2 = cross(phi, phi)/_aux_nobs
	return(Omega_2)
}

					//----------------------------//
					// compute phi_21
					//----------------------------//
real matrix `PILS'::compute_phi_21()
{
	real scalar	l
	real matrix	phi_21

	/*-------------------------------------*/
	// first fix l, the index in x2

	phi_21 = 0

	for (l = 1; l<= _d2; l++) {
		phi_21 =  phi_21 :+ compute_phi_21_l(l)
	}

	return(phi_21)
}
					//----------------------------//
					// compute phi_21 for one l
					//----------------------------//
real matrix `PILS'::compute_phi_21_l(real scalar l)
{
	real colvector	eta, aux_x2_hat
	real scalar	b2
	real matrix	phi_21_l
					// eta = x2_l - glink(z_l)
	// notice we use leave-one-out glink prediction
	aux_x2_hat = compute_glink_loo(	///
		_aux_x2[., l], 		///
		_aux_x2_z[., l], 	///
		_aux_x2_z[., l])
	eta = _aux_x2[., l]  - aux_x2_hat
					// b2 = b_x2[l]
	b2 = _beta_x2[l]

					// E(main_x1|main_z_l = aux_z_l)	
	real scalar	k
	real matrix	aux_x1_hat	

	aux_x1_hat = J(_aux_nobs, _d1, .)
	for (k=1; k<=_d1; k++) {
		aux_x1_hat[., k] = compute_glink(	///
			_main_x1[., k],			///
			_main_x2_z[., l],		///
			_aux_x2_z[., l])
	}

					// E(main_x2_glink|main_z_l = aux_z_l)
	real matrix	aux_glink_hat
	aux_glink_hat = J(_aux_nobs, _d2, .)
	for (k=1; k<=_d2; k++) {
	    	aux_glink_hat[., k] = compute_glink(	///
	    		_main_x2_glink[., k],		///
	    		_main_x2_z[., l],		///
	    		_aux_x2_z[., l])
	}

					// E(x3^A|z_l^A = Z_{l, j}^B)
	real matrix	aux_x3i_hat
	real colvector	theta
	aux_x3i_hat = J(_aux_nobs, _d3i, .)
	for (k=1; k<= _d3i; k++) {
		theta = run_regress_wrk(	///
			_main_x3i[., k],	///
			(_main_x2_z[., l] ,_main_ones) )
		aux_x3i_hat[., k] = (_aux_x2_z[., l], _aux_ones)*theta

//!!		aux_x3i_hat[., k] = compute_glink(	///
//!!			_main_x3i[., k],		///
//!!			_main_x2_z[., l],		///
//!!			_aux_x2_z[., l])
	}

	phi_21_l = (aux_x1_hat, aux_glink_hat, aux_x3i_hat, _aux_ones):*eta*b2

	return(phi_21_l)
}
					//----------------------------//
					// compute phi_22
					//----------------------------//
real matrix `PILS'::compute_phi_22()
{
	real scalar	l	
	real matrix	phi_22

	phi_22 = 0

	for (l=1; l<=_d2; l++) {
		phi_22 = phi_22 :+ compute_phi_22_l(l)
	}

	return(phi_22)
}

					//----------------------------//
					// compute phi_22_l
					//----------------------------//
real matrix `PILS'::compute_phi_22_l(	///
	real scalar	l)
{
	real scalar	b2, k
	real matrix	g3hat, g3res, phi_22_l
	real colvector	theta, g_deriv

	b2 = _beta_x2[l]

	// g3^B(Z^A)
	g3hat = J(_main_nobs, _d3, .)
	for (k=1; k<=_d3; k++) {
		theta  = run_regress_wrk(	///
			_aux_x3[., k],		///
		 	(_aux_x2_z[., l], _aux_ones) )
		g3hat[., k] = (_main_x2_z[., l], _main_ones)*theta

//!!		g3hat[., k] = compute_glink(	///
//!!			_aux_x3[., k], 		///
//!!			_aux_x2_z[., l],	///
//!!			_main_x2_z[., l])
	}

	// g3 residuals
	g3res = _main_x3 - g3hat
	
	// derivative of g2l
	real scalar h	
	real colvector	gph, gmh

	h = get_kernel_bandwidth(_aux_x2_z[., l], _aux_nobs)
	gph = compute_glink(_aux_x2[., l], _aux_x2_z[., l], _main_x2_z[., l]:+h)
	gmh = compute_glink(_aux_x2[., l], _aux_x2_z[., l], _main_x2_z[., l]:-h)
	g_deriv = (gph - gmh)/(2*h)

	real matrix	tmp
	tmp = cross(_main_Xhat, (g3res:*g_deriv))/_main_nobs

	// influence function for index coefficient
	real colvector	psi_l

	if ( (l <= _d2c) | (l > _d2c & _sdrfv == "sir") ) {
		psi_l = compute_psi_index(l)
	}
	else if (_sdrfv == "logit") {
		psi_l = compute_psi_logit(l)
	}
	else if (_sdrfv == "probit") {
		psi_l = compute_psi_probit(l)
	}
	phi_22_l = psi_l*tmp'*b2

	return(phi_22_l)
}
					//----------------------------//
					// compute influence funtion for index
					//----------------------------//
real matrix `PILS'::compute_psi_index(	///
	real scalar	l)	
{
	real scalar	j, k
	real matrix	theta_jk, id, psi_index
	real colvector	idx
	real rowvector	theta_bar

	theta_jk = J(_aux_nobs, _d3, .)
	id = range(1, _aux_nobs, 1)

	for (j=1; j<=_aux_nobs; j++) {
		idx = selectindex(id:!=j)
		
		if (l <= _d2c) {
			theta_jk[j, .] = sdr_c(		///
				_aux_x2[idx, l], 	///
				_aux_x3[idx, .])'
		}
		else {
			k = l - _d2c
			theta_jk[j, .] = sdr_fv(	///
				_aux_x2[idx, l],	///
				_aux_x3[idx, .],	///
				_st_x2fv_vec[k])'
		}
	}

	theta_bar = mean(theta_jk)

	psi_index = (theta_bar:-theta_jk)*(_aux_nobs - 1)

	return(psi_index)
}


					//----------------------------//
					// sir fv
					//----------------------------//
real colvector	`PILS'::sir_fv(	///
	real colvector	x2,	///
	real matrix	x3,	///
	real scalar	H)
{
	real matrix	X, x3bar_H, Sigma, Sigma_H, tmp, V
	real scalar	k, nobs_k, nobs, p, mk, d3
	real rowvector	x3bar_h, x3bar
	real colvector	mk_H, theta, w, b, idx, uniq_x2

	X = (x2, x3)
	nobs = rows(X)
	p = cols(X)
	d3 = p - 1

	_sort(X, 1)
	x3bar_H = J(H, d3, .)
	mk_H = J(H, 1, .)

	uniq_x2 = uniqrows(X[., 1])

	/* ---------------------------------------------------------- */
	// get mean of x3 in each slice
	for (k=1; k<=H; k++) {
		idx = selectindex(X[., 1]:== uniq_x2[k])
		nobs_k = length(idx)
		x3bar_h = colsum(X[idx, 2..p])/nobs_k
		x3bar_H[k, .] = x3bar_h
		mk_H[k] = nobs_k
	}

	/* ---------------------------------------------------------- */
	// compute Sigma_H
	Sigma_H = J(d3, d3, 0)
	x3bar = mean(X[., 2..p])

	for (k=1; k<= H; k++) {
		tmp = x3bar_H[k, .] :- x3bar
		tmp = mk_H[k]*(tmp'*tmp)

		Sigma_H = Sigma_H + tmp
	}
	Sigma_H = Sigma_H/nobs

	tmp = x3:-x3bar
	Sigma = cross(tmp, tmp)/nobs

	eigensystem(invsym(Sigma)*Sigma_H, V=., .)

	theta = Re(V[., 1])

	return(theta)
	
}
					//----------------------------//
					// compute glink loo
					//----------------------------//
real colvector	`PILS'::compute_glink_loo(	///
	real colvector	aux_x2,			///
	real colvector	aux_z,			///
	real colvector	main_z)
{
	real scalar	h, i, nobs
	real colvector	glink, u, ku, idx, oi

	h = get_kernel_bandwidth(aux_z, rows(aux_z))	

	glink = J(rows(main_z), 1, .)
	nobs = rows(aux_z)
	oi = rangen(1, nobs, 1) 

	for (i=1; i<= rows(main_z); i++) {
		idx = selectindex( oi:!= i)
		u = (aux_z[idx]:-main_z[i])/h
		ku = kernel_func(u)
		glink[i] = mean(aux_x2[idx]:*ku)/mean(ku)
	}

	return(glink)
}

					//----------------------------//
					// compute psi for logit or probit
					//----------------------------//
real matrix `PILS'::compute_psi_logit(	///
	real scalar l)
{
	real matrix	psi, x3, sc, H, Hinv
	real colvector	z, y
	real scalar	N, p

	x3 = (_aux_x3, _aux_ones)
	z = _aux_x2_z[., l]

	y = _aux_x2[., l]
	sc = y:/(1:+exp(z)) :- (1:-y):/(1 :+exp(-z))
	sc = sc:*x3
	N = _aux_nobs

	H = -exp(-z):/(1 :+ exp(-z)):^2
	H =cross(x3, H, x3)
	H = H/N
	Hinv = qrinv(H)

	psi = -sc*Hinv

	p = cols(psi)
	psi = psi[., 1..p-1]

	return(psi)
}
					//----------------------------//
					// compute psi for probit
					//----------------------------//
real matrix `PILS'::compute_psi_probit(	///
	real scalar l)
{
	real matrix	psi, x3, sc, H, Hinv
	real colvector	z, y, ncdf, one_m_ncdf, npdf
	real scalar	N, p

	x3 = (_aux_x3, _aux_ones)
	z = _aux_x2_z[., l]

	y = _aux_x2[., l]
	N = _aux_nobs

	ncdf = normal(z)
	one_m_ncdf = 1:-ncdf
	npdf = normalden(z)

	sc = npdf:*(y:/ncdf - (1:-y):/one_m_ncdf)
	sc = sc:*x3

	H = -npdf:*( y:*(npdf:+z:*ncdf):/(ncdf:^2)	///
		:+(1:-y):*(npdf:-z:*one_m_ncdf):/(one_m_ncdf:^2))

	H = cross(x3, H, x3)
	H = H/N
	Hinv = qrinv(H)

	psi = -sc*Hinv

	p = cols(psi)
	psi = psi[., 1..p-1]

	return(psi)
}

end
