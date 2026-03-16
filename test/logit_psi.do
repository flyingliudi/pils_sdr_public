
cscript

do compile

local seed 12345671
set seed `seed'
set sortseed `seed'

cap log close

log using tmp, text replace


dgp_two_sample, 		///
	nobs_main(300) 		///
	nobs_aux(100)		///
	d3_fv(4)		///
	x3vars(x3_c1 x3_fv1)	///
	d2(2)			///
	dta_main(sample_main)	///
	dta_aux(sample_aux)

use sample_main

//!! pils `spec' y x1_* (x2_c* i.(x2_fv*)= x3_c1 i.x3_fv1) 	///
//!! 	using sample_aux 				///
//!! 	, extram(x3_c* i.(x3_fv*)) 			///
//!! 	sdrc(sir) 					///
//!! 	sdrfv(logit)

logit x2_fv1 x3_c* i.(x3_fv*), vce(robust)

fvexpand x3_c* i.(x3_fv*)
local x3 `r(varlist)'

predict double mysc , score

mata:
x3 = st_data(., `"`x3'"')
N = rows(x3)
ones = J(N, 1, 1)
x3 = (x3, ones)
b = st_matrix("e(b)")
z = x3*(b')

y = st_data(., "x2_fv1")

sc = y:/(1:+exp(z)) :- (1:-y):/(1 :+exp(-z))

mysc = st_data(., "mysc")
assert(mreldif(sc, mysc)<=1e-8)

sc = sc:*x3

H = -exp(-z):/(1 :+ exp(-z)):^2
H =cross(x3, H, x3)
H = H/N
Hinv = qrinv(H)

psi = -sc*Hinv
qc = N/(N-1)
myV = cross(psi, psi)/N^2*qc
V = st_matrix("e(V)") 

mreldif(myV, V)


end
log close
