"""
EXAMPLES::

    sage: var('a b c d e f k m n p s t y')
    (a, b, c, d, e, f, k, m, n, p, s, t, y)
    sage: w0=SR.wild(0); w1=SR.wild(1); w2=SR.wild(2); w3=SR.wild(3)
    sage: w4=SR.wild(4); w5=SR.wild(5); w6=SR.wild(6); w7=SR.wild(7)

    sage: (sqrt(c+x)*c).match(w0 * sqrt(w0+w1))
    {$0: c, $1: x}
    sage: (sqrt(c+x)*c).match(w1 * sqrt(w0+w1))
    {$0: x, $1: c}
    sage: (sqrt(c+x)*x).match(w0 * sqrt(w0+w1))
    {$0: x, $1: c}
    sage: (sqrt(c+x)*x).match(w1 * sqrt(w0+w1))
    {$0: c, $1: x}

    sage: def check(s, p):
    ....:     m = s.match(p)
    ....:     if (m and (p.subs(m, flags=1) - s).is_trivial_zero()):
    ....:         return True
    ....:     return False
    
    sage: assert check((a+b)*(a+c), (w0+w1)*(w0+w2))
    sage: assert check((a+b)*(a+c), (w0+w1)*(w1+w2))
    sage: assert check((c+b)*(a+c), (w0+w1)*(w0+w2))
    sage: assert check((c+b)*(a+c), (w0+w1)*(w1+w2))
    
    sage: assert check((a+b)^(a+c), (w0+w1)^(w0+w2)) # known bug
    sage: assert check((a+b)^(a+c), (w0+w1)^(w1+w2))
    sage: assert check((c+b)^(a+c), (w0+w1)^(w0+w2))
    sage: assert check((c+b)^(a+c), (w0+w1)^(w1+w2))
    
    sage: assert check((a^b)^(a+c), (w0^w1)^(w0+w2))
    sage: assert not check((a^b)^(a+c), (w0^w1)^(w1+w2))
    sage: assert check((c^b)^(a+c), (w0^w1)^(w0+w2))
    sage: assert not check((c^b)^(a+c), (w0^w1)^(w1+w2))

    sage: assert check (hypergeometric((a+c,b+c), (c,d), x),
    ....: hypergeometric((w0+w2,w1+w2), (w2,w3), x))
    sage: assert not check (hypergeometric((a+c,b+c), (c,d), x),
    ....: hypergeometric((w0+w1,w1+w2), (w2,w3), x))

    sage: assert check(((a+b+1)*(a+b+2))*b, w1*((w0+w1+w2)*(w0+w1+w3)))
    sage: assert check(((a+b+1)*(a+b+2))*b, w0*((w0+w1+w2)*(w0+w1+w3)))
    sage: assert check(((a+b+1)*(a+b+2))*a, w1*((w0+w1+w2)*(w0+w1+w3)))
    sage: assert check(((a+b+1)*(a+b+2))*a, w0*((w0+w1+w2)*(w0+w1+w3)))
    sage: assert check(((a+b+1)*(a+b+2))*b, w4*((w4+w1+w2)*(w4+w1+w3)))
    sage: assert check(((a+b+1)*(a+b+2))*a, w4*((w4+w1+w2)*(w4+w1+w3)))

    
    sage: assert check((a*b*c+b*c*d+c*d*e)*(b*c*e+c*d*f+d*e*f),
    ....: (w0*w2*w4+w4*w3*w1+w2*w4*w3)*(w4*w1*w2+w3*w1*w5+w3*w4*w5))
    sage: assert check((a*b*d+b*c*d+c*d*e)*(b*c*e+c*d*f+d*e*f),
    ....: (w0*w2*w4+w4*w3*w1+w2*w4*w3)*(w2*w3*w1+w3*w4*w5+w4*w1*w5))
    sage: assert check((a*b*d+b*c*d+c*d*e)*(b*c*a+b*d*f+d*a*f),
    ....: (w0*w2*w4+w4*w3*w1+w2*w4*w3)*(w2*w3*w1+w3*w4*w5+w4*w1*w5))
"""
