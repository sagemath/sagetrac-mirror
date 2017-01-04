"""
        References
        ==========

        .. [1] Marko Petkovsek, Herbert S. Wilf, Doron Zeilberger, A = B,
               AK Peters, Ltd., Wellesley, MA, USA, 1997, pp. 73--100

TESTS::

    sage: _ = var('a b k m n r')
    sage: SR(1).gosper_sum(k)
    k
    sage: SR(1).gosper_sum(k,0,n)
    n + 1
    sage: SR(1).gosper_sum(k, a, b)
    -a + b + 1
    sage: a.gosper_sum(k)
    a*k
    sage: a.gosper_sum(k,0,n)
    a*(n + 1)
    sage: a.gosper_sum(k,a,b)
    -(a - b - 1)*a
    sage: k.gosper_sum(k)
    1/2*(k - 1)*k
    sage: k.gosper_sum(k,0,n)
    1/2*(n + 1)*n
    sage: k.gosper_sum(k, a, b)
    -1/2*(a + b)*(a - b - 1)
    sage: k.gosper_sum(k, a+1, b)
    -1/2*(a + b + 1)*(a - b)
    sage: k.gosper_sum(k, a, b+1)
    -1/2*(a + b + 1)*(a - b - 2)
    sage: (k^3).gosper_sum(k)
    1/4*(k - 1)^2*k^2
    sage: (k^3).gosper_sum(k, a, b)
    -1/4*(a^2 + b^2 - a + b)*(a + b)*(a - b - 1)

    sage: (1/k).gosper_sum(k)
    Traceback (most recent call last):
    ...
    ValueError: expression not Gosper-summable
    sage: x = (1/k/(k+1)/(k+2)/(k+3)/(k+5)/(k+7)).gosper_sum(k,a,b)
    sage: y = sum(1/k/(k+1)/(k+2)/(k+3)/(k+5)/(k+7) for k in range(5,500))
    sage: assert x.subs(a==5,b==499) == y

The following are from A==B, p.78 ff. Since the resulting expressions
get more and more complicated with many correct representations that
may differ depending on future capabilities, we check correctness by
doing random summations::

    sage: def check(ex, var, val1, val2from, val2to):
    ....:     import random
    ....:     symb = SR.var('symb')
    ....:     s1 = ex.gosper_sum(var, val1, symb)
    ....:     R = random.SystemRandom
    ....:     val2 = R().randint(val2from, val2to)
    ....:     s2 = sum(ex.subs(var==i) for i in range(val1, val2+1))
    ....:     assert s1.subs(symb==val2) == s2

    sage: def check_unsolvable(ex, *args):
    ....:     try:
    ....:         SR(ex).gosper_sum(*args)
    ....:         raise AssertionError
    ....:     except ValueError:
    ....:         pass

    sage: check((4*n+1) * factorial(n)/factorial(2*n+1), n, 0, 10, 20)
    sage: check_unsolvable(factorial(k), k,0,n)
    sage: (k * factorial(k)).gosper_sum(k,1,n)
    n*factorial(n) + factorial(n) - 1
    sage: (k * factorial(k)).gosper_sum(k)
    factorial(k)
    sage: check_unsolvable(binomial(n,k), k,0,n)
    sage: ((-1)^k*binomial(n,k)).gosper_sum(k,0,n)
    0
    sage: ((-1)^k*binomial(n,k)).gosper_sum(k,0,a)
    -(-1)^a*(a - n)*binomial(n, a)/n
    sage: (binomial(1/2,a-k+1)*binomial(1/2,a+k)).gosper_sum(k,1,b)
    (2*a + 2*b - 1)*(a - b + 1)*b*binomial(1/2, a + b)*binomial(1/2, a - b + 1)/((2*a + 1)*a)
    sage: t = (binomial(2*k,k)/4^k).gosper_sum(k,0,n); t
    (2*n + 1)*binomial(2*n, n)/4^n
    sage: t = t.gosper_sum(n,0,n); t
    1/3*(2*n + 3)*(2*n + 1)*binomial(2*n, n)/4^n
    sage: t = t.gosper_sum(n,0,n); t
    1/15*(2*n + 5)*(2*n + 3)*(2*n + 1)*binomial(2*n, n)/4^n
    sage: t = t.gosper_sum(n,0,n); t
    1/105*(2*n + 7)*(2*n + 5)*(2*n + 3)*(2*n + 1)*binomial(2*n, n)/4^n
    sage: (binomial(2*k+2*a,2*a)*binomial(2*k,k)/binomial(k+a,a)/4^k).gosper_sum(k,0,n)
    (2*a + 2*n + 1)*binomial(2*a + 2*n, 2*a)*binomial(2*n, n)/(4^n*(2*a + 1)*binomial(a + n, a))
    sage: (4^k/binomial(2*k,k)).gosper_sum(k,0,n)
    1/3*(2*4^n*n + 2*4^n + binomial(2*n, n))/binomial(2*n, n)

    # The following are from A==B, 5.7 Exercises
    sage: for k in range(1,5): (n^k).gosper_sum(n,0,m)
    1/2*(m + 1)*m
    1/6*(2*m + 1)*(m + 1)*m
    1/4*(m + 1)^2*m^2
    1/30*(3*m^2 + 3*m - 1)*(2*m + 1)*(m + 1)*m
    sage: for k in range(1,4): (n^k*2^n).gosper_sum(n,0,m)
    2*2^m*m - 2*2^m + 2
    2*2^m*m^2 - 4*2^m*m + 6*2^m - 6
    2*2^m*m^3 - 6*2^m*m^2 + 18*2^m*m - 26*2^m + 26
    sage: (1 / (n^2 + sqrt(5)*n - 1)).gosper_sum(n,0,m) # known bug
    ...                                           # TODO: algebraic solutions
    sage: check((n^4 * 4^n / binomial(2*n, n)), n, 0, 10, 20)
    sage: check((factorial(3*n) / (factorial(n) * factorial(n+1) * factorial(n+2) * 27^n)), n,0,10,20)
    sage: (binomial(2*n, n)^2 / (n+1) / 4^(2*n)).gosper_sum(n,0,m)
    (2*m + 1)^2*binomial(2*m, m)^2/(4^(2*m)*(m + 1))
    sage: (((4*n-1) * binomial(2*n, n)^2) / (2*n-1)^2 / 4^(2*n)).gosper_sum(n,0,m)
    -binomial(2*m, m)^2/4^(2*m)
    sage: check(n * factorial(n-1/2)^2 / factorial(n+1)^2, n,0,10,20)

    sage: (n^2 * a^n).gosper_sum(n,0,m)
    (a^2*a^m*m^2 - 2*a*a^m*m^2 - 2*a*a^m*m + a^m*m^2 + a*a^m + 2*a^m*m + a^m - a - 1)*a/(a - 1)^3
    sage: ((n - r/2)*binomial(r, n)).gosper_sum(n,0,m)
    1/2*(m - r)*binomial(r, m)
    sage: x = var('x')
    sage: (factorial(n-1)^2 / factorial(n-x) / factorial(n+x)).gosper_sum(n,1,m)
    (m^2*factorial(m - 1)^2*factorial(x + 1)*factorial(-x + 1) + x^2*factorial(m + x)*factorial(m - x) - factorial(m + x)*factorial(m - x))/(x^2*factorial(m + x)*factorial(m - x)*factorial(x + 1)*factorial(-x + 1))
    sage: ((n*(n+a+b)*a^n*b^n)/factorial(n+a)/factorial(n+b)).gosper_sum(n,1,m) # known bug
    ...   # TODO: take limit if result is NaN

    sage: check_unsolvable(1/n, n,1,m)
    sage: check_unsolvable(1/n^2, n,1,m)
    sage: check_unsolvable(1/n^3, n,1,m)
    sage: ((6*n + 3) / (4*n^4 + 8*n^3 + 8*n^2 + 4*n + 3)).gosper_sum(n,1,m)
    (m + 2)*m/(2*m^2 + 4*m + 3)
    sage: (2^n * (n^2 - 2*n - 1)/(n^2 * (n+1)^2)).gosper_sum(n,1,m)
    -2*(m^2 - 2^m + 2*m + 1)/(m + 1)^2
    sage: ((4^n * n^2)/((n+2) * (n+1))).gosper_sum(n,1,m)
    2/3*(2*4^m*m - 2*4^m + m + 2)/(m + 2)
    sage: check_unsolvable(2^n / (n+1), n,0,m-1)
    sage: check((4*(1-n) * (n^2-2*n-1) / n^2 / (n+1)^2 / (n-2)^2 / (n-3)^2), n, 4, 10, 20)
    sage: check(((n^4-14*n^2-24*n-9) * 2^n / n^2 / (n+1)^2 / (n+2)^2 / (n+3)^2), n, 1, 10, 20)

Exercises 3 (h), (i), (j) require symbolic product support so we leave
them out for now.

::

    sage: _ = var('a b k m n r')
    sage: check_unsolvable(binomial(2*n, n) * a^n, n)
    sage: (binomial(2*n,n)*(1/4)^n).gosper_sum(n)
    2*(1/4)^n*n*binomial(2*n, n)
    sage: ((k-1) / factorial(k)).gosper_sum(k)
    -k/factorial(k)
    sage: F(n, k) = binomial(n, k) / 2^n
    sage: check_unsolvable(F(n, k), k)
    sage: _ = (F(n+1, k) - F(n, k)).gosper_sum(k)
    sage: F(n, k) = binomial(n, k)^2 / binomial(2*n, n)
    sage: check_unsolvable(F(n, k), k)
    sage: _ = (F(n+1, k) - F(n, k)).gosper_sum(k)
    sage: F(n, k) = binomial(n,k) * factorial(n) / factorial(k) / factorial(a-k) / factorial(a+n)
    sage: check_unsolvable(F(n, k), k)
    sage: _ = (F(n+1, k) - F(n, k)).gosper_sum(k)   # known bug

    sage: (1/n/(n+1)/(n+2)/(n+5)).gosper_sum(n)
    1/720*(55*n^5 + 550*n^4 + 1925*n^3 + 2510*n^2 - 1728)/((n + 4)*(n + 3)*(n + 2)*(n + 1)*n)
    sage: (1/n/(n+1)/(n+2)/(n+7)).gosper_sum(n)
    1/1050*(91*n^7 + 1911*n^6 + 15925*n^5 + 66535*n^4 + 142534*n^3 + 132104*n^2 - 54000)/((n + 6)*(n + 5)*(n + 4)*(n + 3)*(n + 2)*(n + 1)*n)
    sage: (1/n/(n+1)/(n+2)/(n+5)/(n+7)).gosper_sum(n)
    1/10080*(133*n^7 + 2793*n^6 + 23275*n^5 + 97755*n^4 + 213472*n^3 + 206892*n^2 - 103680)/((n + 6)*(n + 5)*(n + 4)*(n + 3)*(n + 2)*(n + 1)*n)
"""

