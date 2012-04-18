from sage.rings.padics.all import pAdicField
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.big_oh import O
from sage.rings.arith import binomial, gcd, kronecker
from sage.rings.infinity import Infinity
from sage.matrix.all import Matrix

# @cached_function
def teich(a,p,M):
    r"""
    Returns the Teichmuller representative of `a` in `Qp` as an element of `Z` 
    (`Qp` is created with precision `O(p^M)`)

    INPUT:
        - ``a`` -- element of `Qp`
        - ``p`` -- prime
        - ``M`` -- precision, `O(p^M)`

    OUTPUT:
        
    The Teichmuller representative of a as an element of `Z`

    EXAMPLES:

        sage: teich(2,5,10)
        6139557
        sage: teich(3,7,10)
        146507973
    """
    R = pAdicField(p,M)
    return ZZ(R.teichmuller(a))

# @cached_function
def logp(p,z,M):
    r"""
    Returns the truncation `\sum_{j=1}^{M-1} (-1)^j/j z^j` of `log_p(1+z)`
        
    INPUT:
        - ``p`` -- prime 
        - ``z`` -- variable
        - ``M`` -- precision

    OUTPUT:
    
    The truncated power series expansion `\sum_{j=1}^{M-1} (-1)^j/j z^j` 
    of `log_p(1+z)`

    EXAMPLES:

        sage: R.<z> = QQ['z']
        sage: logp(5,z,10)
        -1/9*z^9 + 1/8*z^8 - 1/7*z^7 + 1/6*z^6 - 1/5*z^5 + 1/4*z^4 - 1/3*z^3 + 1/2*z^2 - z

        sage: R.<z> = QQ['z']
        sage: logp(7,z,6)
        -1/5*z^5 + 1/4*z^4 - 1/3*z^3 + 1/2*z^2 - z          
    """
    ans = 0
    one = ZZ(1)
    for j in range(1,M):
        ans=ans+((-one)**j)/j*(z**j)
    return ans

# @cached_function
def loggam_binom(p,gam,z,n,M):
    r"""
    Returns the list of coefficients in the power series
    expansion (up to precision `M`) of `{\log_p(z)/\log_p(\gamma) \choose n}`

    INPUT:

        - ``p`` --  prime
        - ``gam`` -- topological generator e.g., `1+p`
        - ``z`` -- variable
        - ``n`` -- nonnegative integer
        - ``M`` -- precision

    OUTPUT:

    The list of coefficients in the power series expansion of 
    `{\log_p(z)/\log_p(\gamma) \choose n}`
    
    EXAMPLES:

        sage: R.<z> = QQ['z']  
        sage: loggam_binom(5,1+5,z,2,4)
        [0, -3/205, 651/84050, -223/42025]
        sage: loggam_binom(5,1+5,z,3,4)
        [0, 2/205, -223/42025, 95228/25845375]
    """
    L = logp(p,z,M)
    logpgam = L.substitute(z = (gam-1)) #log base p of gamma
    loggam = L/logpgam                  #log base gamma 
    return z.parent()(binomial(loggam,n)).truncate(M).list()
        
# @cached_function
def phi_on_Da(Phi,a,D):
    """
    Returns `\Phi_{\chi}` where `\chi` is a character of conductor `D` 

    INPUT:
        - ``Phi`` -- overconvergent `U_p`-eigensymbol
        - ``a`` -- integer in [0..p-1]
        - ``D`` -- conductor of the quadratic twist `\chi`

    OUTPUT:

    `\Phi_{\chi}`

    EXAMPLES:
    
    """
    p = Phi.p()
    ans = Phi.zero_elt()
    for b in range(1,abs(D)+1):
        if gcd(b,D)==1:        
            M1 = Matrix(2,2,[1,b/abs(D),0,1])
            ans=ans+Phi.eval(M1*Matrix(2,2,[a,1,p,0])).act_right(M1).scale(kronecker(D,b)).normalize()
    return ans.normalize()

# @cached_function
def basic_integral(Phi,a,j,ap,D):
    """
    Returns `\int_{a+pZ_p} (z-{a})^j d\Phi(0-infty)` 
         -- see formula [Pollack-Stevens, sec 9.2] 

    INPUT:

        - ``Phi`` -- overconvergnt `U_p`-eigensymbol
        - ``a`` -- integer in [0..p-1]
        - ``j`` -- positive integer
        - ``ap`` -- Hecke eigenvalue?
        - ``D`` -- conductor of the quadratic twist `\chi`

    OUTPUT:

    `\int_{a+pZ_p} (z-{a})^j d\Phi(0-\infty)` 

    EXAMPLES:

    """
    M = Phi.num_moments()
    p = Phi.p()
    ap = ap*kronecker(D,p)
    ans = 0
    for r in range(j+1):
        ans = ans+binomial(j,r)*((a-teich(a,p,M))**(j-r))*(p**r)*phi_on_Da(Phi,a,D).moment(r)
    return ans/ap

def pLfunction_coef(Phi,ap,n,D,gam,error=None):
    """
    Returns the nth coefficient of the `p`-adic `L`-function in the 
    `T`-variable of a quadratic twist of `\Phi`.  
    
    If error is not specified, then the correct error bound is computed 
    and the answer is return modulo that accuracy.

    INPUT:
        - ``Phi`` -- overconvergent Hecke-eigensymbol
        - ``ap`` -- eigenvalue of `U_p`
        - ``n`` -- index of desired coefficient
        - ``D`` -- discriminant of quadratic twist
        - ``gam`` -- topological generator of `1 + pZ_p`
    
    OUTPUT:

    The `n`th coefficient of the `p`-adic `L`-function in the `T`-variable
    of a quadratic twist of `\Phi`

    EXAMPLES:
    
    """
    S = QQ[['z']]
    z = S.gen()
    p = Phi.p()
    M = Phi.num_moments()
    R = pAdicField(p,M)
    lb = loggam_binom(p,gam,z,n,2*M)
    dn = 0
    if n == 0:
        err = M
    else:
        if (error == None):
            err = min([j+lb[j].valuation(p) for j in range(M,len(lb))])
        else:
            err = error
        lb = [lb[a] for a in range(M)]
        
    for j in range(len(lb)):
        cjn = lb[j]
        temp = 0
        for a in range(1,p):
            temp = temp + (teich(a,p,M)**(-j))*basic_integral(Phi,a,j,ap,D)
        dn = dn + cjn*temp
    return dn + O(p**err)

def pLfunction(Phi,ap,quad_twist=None):
    """
    Returns the `p`-adic `L`-function in the `T`-variable of a quadratic 
    twist of `\Phi`

    INPUT:
        - ``Phi`` -- overconvergent Hecke-eigensymbol
        - ``ap`` -- eigenvalue at `p`
        - ``quad_twist`` -- conductor of quadratic character

    OUTPUT:

    Returns the `p`-adic `L`-function of a quadratic twist of `\Phi`

    EXAMPLES:
    
    """
    if quad_twist == None:
        D = 1
    else:
        D = quad_twist
    M = Phi.num_moments()
    p = Phi.p()
    gam = 1+p
    for a in range(1,p):
        for j in range(M):
            basic_integral(Phi,a,j,ap,D)
    SS = QQ[['T']]
    T = SS.gen()
    ans = pLfunction_coef(Phi,ap,0,D,gam)+0*T
    S = QQ[['z']]
    z = S.gen()
    err = Infinity
    n = 1
    while (err > 0) and (n < M):
        lb = loggam_binom(p,gam,z,n,2*M)
        err = min([j+lb[j].valuation(p) for j in range(M,len(lb))])
        if err > 0:
            dn = pLfunction_coef(Phi,ap,n,D,gam,error=err)
            ans = ans+dn*(T**n)
        n = n+1
    return ans

def lambda_inv(L):
    """
    Returns the `\lambda` invariant of `L`
    """
    v = L.list()
    vals = [v[a].valuation() for a in range(len(v))]
    return vals.index(min(vals))
