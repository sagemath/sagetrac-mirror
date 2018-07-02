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
    ....:     if not m:
    ....:         return 0
    ....:     if not (p.subs(m, flags=1) - s).is_trivial_zero():
    ....:         raise AithmeticError
    ....:     m = s.match(p, all=True)
    ....:     if not m:
    ....:         raise AithmeticError
    ....:     n = 0
    ....:     for mm in m:
    ....:         if (p.subs(mm, flags=1) - s).is_trivial_zero():
    ....:             n = n+1
    ....:         else:
    ....:             raise AithmeticError
    ....:     return n
    
    sage: check(sqrt(c+x)*c, w0 * sqrt(w0+w1))
    1
    sage: check(sqrt(c+x)*c, w1 * sqrt(w0+w1))
    1
    sage: check(sqrt(c+x)*x, w0 * sqrt(w0+w1))
    1
    sage: check(sqrt(c+x)*x, w1 * sqrt(w0+w1))
    1

    sage: check(sqrt(c+x)*x, w0 * sqrt(w1+w2))
    2
    sage: check(sqrt(a+b)*c, w0 * sqrt(w1+w2))
    2

    sage: check((a+b)*(a+c), (w0+w1)*(w0+w2))
    2
    sage: check((a+b)*(a+c), (w0+w1)*(w1+w2))
    2
    sage: check((c+b)*(a+c), (w0+w1)*(w0+w2))
    2
    sage: check((c+b)*(a+c), (w0+w1)*(w1+w2))
    2
    
    sage: check((a+b)*(a+c), (w0+w1)*(w3+w2))
    8
    sage: check((a+b)*(d+c), (w0+w1)*(w3+w2))
    8

    sage: check((a+b)^(a+c), (w0+w1)^(w0+w2))
    1
    sage: check((a+b)^(a+c), (w0+w1)^(w1+w2))
    1
    sage: check((c+b)^(a+c), (w0+w1)^(w0+w2))
    1
    sage: check((c+b)^(a+c), (w0+w1)^(w1+w2))
    1
    
    sage: check((a^b)^(a+c), (w0^w1)^(w0+w2))
    1
    sage: check((a^b)^(a+c), (w0^w1)^(w1+w2))
    0
    sage: check((c^b)^(a+c), (w0^w1)^(w0+w2))
    1
    sage: check((c^b)^(a+c), (w0^w1)^(w1+w2))
    0

    sage: check(a+b+c, w0+w1+w2)
    6
    sage: check(a+b+c, w0+w1)
    6
    sage: check(a+b+c+d, w0+w1+w2+w3)
    24
    sage: check(a+b+c+d, w0+w1+w2)
    36

    sage: check (hypergeometric((a+c,b+c), (c,d), x),
    ....: hypergeometric((w0+w2,w1+w2), (w2,w3), x))
    1
    sage: check (hypergeometric((a+c,b+c), (c,d), x),
    ....: hypergeometric((w0+w1,w1+w2), (w2,w3), x))
    0

    sage: check(((a+b+1)*(a+b+2))*b, w1*((w0+w1+w2)*(w0+w1+w3)))
    2
    sage: check(((a+b+1)*(a+b+2))*b, w0*((w0+w1+w2)*(w0+w1+w3)))
    2
    sage: check(((a+b+1)*(a+b+2))*a, w1*((w0+w1+w2)*(w0+w1+w3)))
    2
    sage: check(((a+b+1)*(a+b+2))*a, w0*((w0+w1+w2)*(w0+w1+w3)))
    2
    sage: check(((a+b+1)*(a+b+2))*b, w4*((w4+w1+w2)*(w4+w1+w3)))
    2
    sage: check(((a+b+1)*(a+b+2))*a, w4*((w4+w1+w2)*(w4+w1+w3)))
    2

    sage: check(((a+b+2)*(d+e+2))*2, w0*((w0+w1+w2)*(w0+w4+w3)))
    8
    
    sage: check((a*b*c+b*c*d+c*d*e)*(b*c*e+c*d*f+d*e*f),
    ....: (w0*w2*w4+w4*w3*w1+w2*w4*w3)*(w4*w1*w2+w3*w1*w5+w3*w4*w5))
    1
    sage: check((a*b*d+b*c*d+c*d*e)*(b*c*e+c*d*f+d*e*f),
    ....: (w0*w2*w4+w4*w3*w1+w2*w4*w3)*(w2*w3*w1+w3*w4*w5+w4*w1*w5))
    1
    sage: check((a*b*d+b*c*d+c*d*e)*(b*c*a+b*d*f+d*a*f),
    ....: (w0*w2*w4+w4*w3*w1+w2*w4*w3)*(w2*w3*w1+w3*w4*w5+w4*w1*w5))
    1

    sage: check(-3*(b*x + a)^(5/2)*a/b^3, w0*w1^w2*w3^w4)
    2
    sage: check(3*a*b^2*x^(m + 2),        w0*w1^w2*w3^w4)
    2
    sage: check(2*(d*x)^(m + 6)*a*b,    w0*w1*(w1*w3)^w4)
    0
    sage: check(-sqrt(x*y)*sqrt(a*b), w0*(w1*w5)^w2*(w3*w4)^w2)
    8
    sage: check(-sqrt(-(a*x + 1)*(a*x - 1))*(a*x + 1)*
    ....: ((x+1)*(a*x - 1))^(1/2),    w0*(w1*w5)^w2*(w3*w5)^w4)
    2
"""
