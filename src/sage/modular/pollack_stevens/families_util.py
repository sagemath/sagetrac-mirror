def ps_normalize(f, p, p_prec):
	"""reduces all of the coefficients of the power series modulo p^N"""
	v = Sequence(f)
	v = [v[a] % (p ** p_prec) for a in range(len(v))]
	S = f.parent()
	return S(v)

def logp_fcn(p, p_prec, z):
	"""this is the *function* on Z_p^* which sends z to log_p(z) using a power series truncated at p_prec terms"""
	R = Qp(p, 2 * p_prec)
	z = z / R.teichmuller(z)
    return sum([((-1) ** (m - 1)) * ((z - 1) ** m) / m for m in range(1, p_prec)])

def logpp(p, p_prec):
	"""returns the (integral) power series for log_p(1+p*z) -- extra p here!"""
	SS = PolynomialRing(QQ, 'y')
    y = SS.gen()
    return sum([((-1) ** (m - 1)) * ((p * y) ** m) / m for m in range(1, p_prec)])

def logpp_gam(p, p_prec):
	"""returns the (integral) power series log_p(1+p*z)*(1/log_p(1+p)) where the denominator is computed with some accuracy"""
	L = logpp(p, p_prec)
	loggam = ZZ(logp_fcn(p, p_prec * (p ** 2), 1 + p))
	return ps_normalize(L / loggam, p, p_prec)

@cached_function
def logpp_binom(n, p, p_prec):
	"""returns the (integral) power series p^n*(log_p(1+p*z)/log_p(1+p) choose n)"""
	#prod=1+0*z
	L = logpp_gam(p, p_prec)
    ans = prod([(L - j) for j in range(n)])
	#for j in range(0,n):
	#	prod=prod*(L-j)
	ans *= (p ** n) / factorial(n)
	
	return ps_normalize(ans.truncate(p_prec), p, p_prec)

#@cached_function
def automorphy_factor_matrix(p, a, c, k, chi, p_prec, var_prec, R):
	S = PolynomialRing(R, 'z')
    z = S.gens()[0]
    w = R.gen()
    aut = S(1)
    for n in range(1, var_prec):
        LB = logpp_binom(n, p, p_prec)
        ta = ZZ(Qp(p, 2 * max(p_prec, var_prec)).teichmuller(a))
        arg = (a / ta - 1) / p + c / (p * ta) * z
        aut += LB(arg) * (w ** n)
    aut *= (ta ** k)
    if not (chi is None):
        aut *= chi(a)    
    aut = Sequence(aut)
    if len(aut) < p_prec:
        aut += [0] * (p_prec - len(aut))
    return aut