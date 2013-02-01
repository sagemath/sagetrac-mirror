from sage.structure.sage_object import SageObject

def ps_normalize(f,p,N):
	"""reduces all of the coefficients of the power series modulo p^N"""
	v=Sequence(f)
	v=[v[a]%(p^N) for a in range(len(v))]
	S=f.parent()
	f=S(v)

	return f
		
def logp_fcn(p,N,z):
	"""this is the *function* on Z_p^* which sends z to log_p(z) using a power series truncated at N terms"""
	R=pAdicField(p,2*N)
	z=z/R.teichmuller(z)
	ans=0
	for m in range(1,N):
		ans=ans+(-1)**(m-1)*(z-1)**m/m

	return ans

def logpp(p,N,z):
	"""returns the (integral) power series for log_p(1+p*z) -- extra p here!"""
	ans=0
	for m in range(1,N):
		ans=ans+(-1)^(m-1)*(p*z)^m/m

	return ans

def explog_s_ell(p, ell, N, deg, w, a = 1):
	"""returns the power series in `w` for `\langle ell\rangle^{as}`. here, must substitute `w = ((1+p)^s - 1)/p)`. the power series is cut off at degree deg and the coefficients are to precision N."""
	R = w.parent()
	if ell.is_zero():
		return R.zero_element()
	if ell.is_one():
		return R(1)
	QQp = Qp(p, prec = max(20, N))
	teich = QQp.teichmuller(ell)
	principal_part = ell / teich
	r = a * principal_part.log(0) / QQp(1+p).log(0)
	b = QQp(1)
	R = w.parent()
	ret = R(1)
	pN = p ** N
	for i in range(1, deg + 1):
		b = b * p * (r - (i - 1)) / i
		ret += (b.lift() % pN) * (w ** i)
	return ret

def ell_power_fam(p, ell, k, N, deg, w):
	r"""
	Returns the power series in `w` for `\ell^{k-1}\langle\ell\rangle^{s+2-k}`. Here, must substitute `w=((1+p)^s-1)/p). The power series is cut off at degree deg and the coefficients are to precision N.
	"""
	R = w.parent()
	if ell.is_zero():
		return R.zero_element()
	if ell.is_one():
		return R(1)
	QQp = Qp(p, prec = max(20, N))
	teich = (QQp.teichmuller(ell) ** (k - 2)).lift()
	ret = (teich * ell * explog_s_ell(p, ell, N, deg, w)).list()
	pN = p ** N
	ret = [c % pN for c in ret]
	return R(ret)

def logpp_gam(p,N,z):
	"""returns the (integral) power series log_p(1+p*z)*(1/log_p(1+p)) where the denominator is computed with some accuracy"""
	L=logpp(p,N,z)
	loggam=ZZ(logp_fcn(p,N*p^2,1+p))
	return ps_normalize(L/loggam,p,N)

@cached_function
def logpp_binom(n,p,N,z):
	"""returns the (integral) power series p^n*(log_p(1+p*z)/log_p(1+p) choose n)"""
	prod=1+0*z
	L=logpp_gam(p,N,z)
	for j in range(0,n):
		prod=prod*(L-j)
	prod=p^n*prod/factorial(n)

	return ps_normalize(prod.truncate(N),p,N)
	
## This function returns the power series K_{a,c,r}(z,w) in the notes.  The answer is returned as a vector whose j-th coordinate is the coefficient of z^j (which itself is a power series in w). 
@cached_function
def aut(p,Mw,Mm,a,c,r,chi,w):
	R=w.parent()
	S=PolynomialRing(R,'zz')
	SS=PolynomialRing(QQ,'yy')
	yy=SS.gen()
	zz=S.gens()[0]
	
	ans=1+0*zz
	for n in range(1,Mw): 
		LB=logpp_binom(n,p,Mm,yy)
		ta=teich(a,p,2*max(Mw,Mm))
		v=(a/ta-1)/p+c/(p*ta)*zz
		ans=ans+w^n*LB(v)

	ans = ans * ta^r * chi(a)

	v=Sequence(ans)
	while len(v)<Mm:
		v=v+[0*w]
	return v

def derivative_in_weight(F, k, p, M):
    r"""Given a power series, return its derivative with respect to the weight. Recall the given `F(w)`, the function of the weight is `f(k)=F(((1+p)^k-1)/p)`. Here, the weight is normalized such that elliptic curves are weight 0."""
    return F.derivative().substitute((((1 + p) ** k) - 1)/p) * ((1 + p) ** k) * logp_fcn(p, M, 1 + p) / p

def derivative_in_weight_fam(F, p, M):
    r"""
    Given a power series F, return the power series such that if it were specialized at weight k, i.e. at `w=((1+p)^k-1)/p`, one would get the derivative_in_weight(F, k, p, M).
    """
    w = F.parent().gen()
    return F.derivative() * (1 + p * w) * logp_fcn(p, M, 1 + p) / p