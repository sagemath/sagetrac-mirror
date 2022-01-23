## -*- encoding: utf-8 -*-
"""
This file (./programmation_doctest.sage) was *autogenerated* from ./programmation.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./programmation_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./programmation.tex, line 118::

  sage: 2*3; 3*4; 4*5          # one comment, 3 results
  6
  12
  20

Sage example in ./programmation.tex, line 137::

  sage: 123 + \
  ....: 345
  468

Sage example in ./programmation.tex, line 201::

  sage: import keyword; sorted(keyword.kwlist)
  [...'and', 'as', 'assert', ..., 'class', 'continue', 'def', 'del'...]

Sage example in ./programmation.tex, line 338::

  sage: y = 3; y = 3 * y + 1; y = 3 * y + 1; y
  31

Sage example in ./programmation.tex, line 368::

  sage: a, b = 10, 20 # (a, b) = (10, 20) and [10, 20] are also possible
  sage: a, b = b, a
  sage: a, b
  (20, 10)

Sage example in ./programmation.tex, line 382::

  sage: temp = a; a = b; b = temp    # equivalent to: a, b = b, a

Sage example in ./programmation.tex, line 391::

  sage: x, y = var('x, y'); a = x ; b = y
  sage: a, b
  (x, y)
  sage: a = a + b ; b = a - b ; a = a - b
  sage: a, b
  (y, x)

Sage example in ./programmation.tex, line 416::

  sage: 2 + 2 == 2^2, 3 * 3 == 3^3
  (True, False)

Sage example in ./programmation.tex, line 483::

  sage: for k in [1..5]:
  ....:    print(7*k)  # block containing a single instruction
  7
  14
  21
  28
  35

Sage example in ./programmation.tex, line 687::

  sage: S = 0 ; k = 0         #        The sum S starts to 0
  sage: while e^k <= 10^6:    #            e^13 <= 10^6 < e^14
  ....:     S = S + k^2       #           accumulates the squares k^2
  ....:     k = k + 1
  sage: S
  819

Sage example in ./programmation.tex, line 734::

  sage: x = 10^4; u = 1; n = 0          # invariant: u = 2^n
  sage: while u <= x: n = n+1; u = 2*u  # or n += 1; u *= 2
  sage: n
  14

Sage example in ./programmation.tex, line 880::

  sage: U = 1.0              # or U = 1. or U = 1.000
  sage: for n in [1..20]:
  ....:   U = 1 / (1 + U^2)
  sage: U
  0.682360434761105

Sage example in ./programmation.tex, line 942::

  sage: S = 0 ; n = 10
  sage: for k in [1..n]:
  ....:     S = S + (2*k) * (2*k+1)
  sage: S
  1650

Sage example in ./programmation.tex, line 961::

  sage: n, k = var('n, k') ; res = sum(2*k*(2*k+1), k, 1, n)
  sage: res, factor(res)    # result expanded, factorised
  (4/3*n^3 + 3*n^2 + 5/3*n, 1/3*(4*n + 5)*(n + 1)*n)

Sage example in ./programmation.tex, line 1074::

  sage: U = 2.0; V = 50.0
  sage: while V-U >= 1.0e-6:      # 1.0e-6 stands for 1.0*10^-6
  ....:   temp = U
  ....:   U = 2 * U * V / (U + V)
  ....:   V = (temp + V) / 2
  sage: U, V
  (9.99999999989256, 10.0000000001074)

Sage example in ./programmation.tex, line 1166::

  sage: U = 0.0        # the sum S0 is empty, of value zero
  sage: V = -1.0       # S1 = -1/1^3
  sage: n = 0          # U and V contain S(2n) and S(2n+1)
  sage: while U-V >= 1.0e-6:
  ....:   n = n+1             # n += 1 is equivalent
  ....:   U = V + 1/(2*n)^3   # going from S(2n-1) to S(2n)
  ....:   V = U - 1/(2*n+1)^3 # going from S(2n) to S(2n+1)
  sage: V, U
  (-0.901543155458595, -0.901542184868447)

Sage example in ./programmation.tex, line 1404::

  sage: u = 6 ; n = 0
  sage: while u != 1:
  ....:   if u % 2 == 0: # the operator % yields the remainder
  ....:     u = u//2     # //: Euclidean division quotient
  ....:   else:
  ....:     u = 3*u+1
  ....:   n = n+1
  sage: n
  8

Sage example in ./programmation.tex, line 1508::

  sage: def fct2 (x, y):
  ....:   return x^2 + y^2
  sage: a = var('a')
  sage: fct2 (a, 2*a)
  5*a^2

Sage example in ./programmation.tex, line 1557::

  sage: def foo (u):
  ....:   t = u^2
  ....:   return t*(t+1)
  sage: t = 1 ; u = 2
  sage: foo(3), t, u
  (90, 1, 2)

Sage example in ./programmation.tex, line 1570::

  sage: a = b = 1
  sage: def f(): global a; a = b = 2
  sage: f(); a, b
  (2, 1)

Sage example in ./programmation.tex, line 1591::

  sage: def AHmean (u, v):
  ....:    u, v = min(u, v), max(u, v)
  ....:    while v-u > 2.0e-8:
  ....:       u, v = 2*u*v/(u+v), (u+v)/2
  ....:    return (u+v) / 2

Sage example in ./programmation.tex, line 1604::

  sage: AHmean (1., 2.)
  1.41421356237309
  sage: AHmean                      # corresponds to a function
    <function AHmean at ...>

Sage example in ./programmation.tex, line 1687::

  sage: def fact1 (n):
  ....:    res = 1
  ....:    for k in [1..n]: res = res*k
  ....:    return res

Sage example in ./programmation.tex, line 1693::

  sage: def fact2 (n):
  ....:    if n == 0: return 1
  ....:    else: return n*fact2(n-1)

Sage example in ./programmation.tex, line 1728::

  sage: def fib1 (n):
  ....:   if n == 0 or n == 1: return n
  ....:   else:
  ....:     U = 0 ; V = 1 # the initial terms u0 and u1
  ....:     for k in [2..n]: W = U+V ; U = V ; V = W
  ....:     return V
  sage: fib1(8)
  21

Sage example in ./programmation.tex, line 1769::

  sage: def fib2 (n):
  ....:   if 0 <= n <= 1: return n      # for n = 0 or n = 1
  ....:   else: return fib2(n-1) + fib2(n-2)

Sage example in ./programmation.tex, line 1857::

  sage: a = 2; n = 6; res = 1         # 1 is the product neutral element
  sage: for k in [1..n]: res = res*a
  sage: res                           # the value of res is 2^6
  64

Sage example in ./programmation.tex, line 1958::

  sage: def pow1 (a, n):
  ....:   if n == 0: return 1
  ....:   elif n % 2 == 0: b = pow1 (a, n//2); return b*b
  ....:   else: return a * pow1(a, n-1)

Sage example in ./programmation.tex, line 1968::

  sage: pow1 (2, 11)                 # result is 2^11
  2048

Sage example in ./programmation.tex, line 2010::

  sage: def pow2 (u, k):
  ....:   v = 1
  ....:   while k != 0:
  ....:     if k % 2 == 0: u = u*u ; k = k//2
  ....:     else: v = v*u ; k = k-1
  ....:   return v

Sage example in ./programmation.tex, line 2022::

  sage: pow2 (2, 10)                 # result is 2^10
  1024

Sage example in ./programmation.tex, line 2109::

  sage: def fib3 (n):
  ....:   A = matrix ([[0, 1], [1, 1]]) ; X0 = vector ([0, 1])
  ....:   return (A^n*X0)[0]

Sage example in ./programmation.tex, line 2114::

  sage: def fib4 (n):
  ....:   return (matrix([[0,1], [1,1]])^n * vector([0,1]))[0]

Sage example in ./programmation.tex, line 2195::

  sage: for k in [1..6]: print('%2d^4 = %4d' % (k, k^4))
   1^4 =    1
   2^4 =   16
   3^4 =   81
   4^4 =  256
   5^4 =  625
   6^4 = 1296

Sage example in ./programmation.tex, line 2274::

  sage: L = [10, 20, 30]
  sage: L
  [10, 20, 30]
  sage: []                    # [] is the empty list
  []

Sage example in ./programmation.tex, line 2299::

  sage: L[1], len(L), len([])
  (20, 3, 0)

Sage example in ./programmation.tex, line 2309::

  sage: L[2] = 33
  sage: L
  [10, 20, 33]

Sage example in ./programmation.tex, line 2318::

  sage: L = [11, 22, 33]
  sage: L[-1], L[-2], L[-3]
  (33, 22, 11)

Sage example in ./programmation.tex, line 2336::

  sage: L = [0, 11, 22, 33, 44, 55]
  sage: L[2:4]
  [22, 33]
  sage: L[-4:4]
  [22, 33]
  sage: L[2:-2]
  [22, 33]
  sage: L[:4]
  [0, 11, 22, 33]
  sage: L[4:]
  [44, 55]

Sage example in ./programmation.tex, line 2359::

  sage: L = [0, 11, 22, 33, 44, 55, 66, 77]
  sage: L[2:6] = [12, 13, 14]        # substitutes [22, 33, 44, 55]

Sage example in ./programmation.tex, line 2393::

  sage: L = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
  sage: L[3:len(L)-5] == L[3-len(L):-5]
  True
  sage: [5 in L, 6 in L]
  [True, False]

Sage example in ./programmation.tex, line 2417::

  sage: L = [1, 2, 3] ; L + [10, 20, 30]
  [1, 2, 3, 10, 20, 30]
  sage: 4 * [1, 2, 3]
  [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]

Sage example in ./programmation.tex, line 2441::

  sage: L = 5*[10, 20, 30] ; L[:3]+L[3:] == L
  True

Sage example in ./programmation.tex, line 2459::

  sage: [1..3, 7, 10..13]
  [1, 2, 3, 7, 10, 11, 12, 13]

Sage example in ./programmation.tex, line 2485::

  sage: list(map (cos, [0, pi/6, pi/4, pi/3, pi/2]))
  [1, 1/2*sqrt(3), 1/2*sqrt(2), 1/2, 0]

Sage example in ./programmation.tex, line 2496::

  sage: list(map (lambda t: cos(t), [0, pi/6, pi/4, pi/3, pi/2]))
  [1, 1/2*sqrt(3), 1/2*sqrt(2), 1/2, 0]

Sage example in ./programmation.tex, line 2534::

  sage: list(map (lambda t: N(cos(t)), [0, pi/6, pi/4, pi/3, pi/2]))
  [1.00000000000000, 0.866025403784439, 0.707106781186548,
  0.500000000000000, 0.000000000000000]

Sage example in ./programmation.tex, line 2546::

  sage: list(map (N, map (cos, [0, pi/6, pi/4, pi/3, pi/2])))
  [1.00000000000000, 0.866025403784439, 0.707106781186548,
  0.500000000000000, 0.000000000000000]

Sage example in ./programmation.tex, line 2551::

  sage: list(map (compose(N, cos), [0, pi/6, pi/4, pi/3, pi/2]))
  [1.00000000000000, 0.866025403784439, 0.707106781186548,
  0.500000000000000, 0.000000000000000]

Sage example in ./programmation.tex, line 2564::

  sage: list(filter (is_prime, [1..55]))
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]

Sage example in ./programmation.tex, line 2578::

  sage: p = 37 ; list(filter (lambda n: n^4 % p == 7, [0..p-1]))
  [3, 18, 19, 34]

Sage example in ./programmation.tex, line 2595::

  sage: list(map(lambda n:2*n+1, [0..15]))
  [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
  sage: [2*n+1 for n in [0..15]]
  [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]

Sage example in ./programmation.tex, line 2607::

  sage: list(filter (is_prime, [1..55]))
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
  sage: [p for p in [1..55] if is_prime(p)]
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]

Sage example in ./programmation.tex, line 2621::

  sage: list(filter (is_prime, [4*n+1 for n in [0..20]]))
  [5, 13, 17, 29, 37, 41, 53, 61, 73]
  sage: [n^2 for n in [1..20] if is_prime(n)]
  [4, 9, 25, 49, 121, 169, 289, 361]

Sage example in ./programmation.tex, line 2651::

  sage: reduce (lambda x, y: 10*x+y, [1, 2, 3, 4])
  1234

Sage example in ./programmation.tex, line 2657::

  sage: reduce (lambda x, y: 10*x+y, [9, 8, 7, 6], 1)
  19876

Sage example in ./programmation.tex, line 2667::

  sage: L = [2*n+1 for n in [0..9]]
  sage: reduce (lambda x, y: x*y, L, 1)
  654729075

Sage example in ./programmation.tex, line 2743::

  sage: prod ([2*n+1 for n in [0..9]], 1) # a list with for
  654729075
  sage: prod ( 2*n+1 for n in [0..9])     # without a list
  654729075
  sage: prod (n for n in [0..19] if n%2 == 1)
  654729075

Sage example in ./programmation.tex, line 2764::

  sage: def fct (x): return 4/x == 2
  sage: all (fct(x) for x in [2, 1, 0])
  False
  sage: any (fct(x) for x in [2, 1, 0])
  True

Sage example in ./programmation.tex, line 2793::

  sage: [[x, y] for x in [1..2] for y in [6..8]]
  [[1, 6], [1, 7], [1, 8], [2, 6], [2, 7], [2, 8]]

Sage example in ./programmation.tex, line 2800::

  sage: [[[x, y] for x in [1..2]] for y in [6..8]]
  [[[1, 6], [2, 6]], [[1, 7], [2, 7]], [[1, 8], [2, 8]]]

Sage example in ./programmation.tex, line 2810::

  sage: list(map (lambda x, y: [x, y], [1..3], [6..8]))
  [[1, 6], [2, 7], [3, 8]]

Sage example in ./programmation.tex, line 2829::

  sage: L = [[1, 2, [3]], [4, [5, 6]], [7, [8, [9]]]]
  sage: flatten (L, max_level = 1)
  [1, 2, [3], 4, [5, 6], 7, [8, [9]]]
  sage: flatten (L, max_level = 2)
  [1, 2, 3, 4, 5, 6, 7, 8, [9]]
  sage: flatten (L)             # equivalent to flatten (L, max_level = 3)
  [1, 2, 3, 4, 5, 6, 7, 8, 9]

Sage example in ./programmation.tex, line 2852::

  sage: x = var('x')
  sage: factor(diff(x*exp(x), [x, x]))
  (x + 2)*e^x
  sage: list(map(lambda n: factor(diff(x*exp(x), n*[x])), [0..6]))
  [x*e^x, (x + 1)*e^x, (x + 2)*e^x, (x + 3)*e^x, (x + 4)*e^x,
  (x + 5)*e^x, (x + 6)*e^x]
  sage: [factor (diff (x*exp(x), n*[x])) for n in [0..6]]
  [x*e^x, (x + 1)*e^x, (x + 2)*e^x, (x + 3)*e^x, (x + 4)*e^x,
  (x + 5)*e^x, (x + 6)*e^x]

Sage example in ./programmation.tex, line 2907::

  sage: L = [1, 8, 5, 2, 9] ; L.reverse() ; L
  [9, 2, 5, 8, 1]
  sage: L.sort() ; L
  [1, 2, 5, 8, 9]
  sage: L.sort(reverse = True) ; L
  [9, 8, 5, 2, 1]

Sage example in ./programmation.tex, line 3031::

  sage: def homogLex (P, Q):
  ....:  sp = sum (P) ; sq = sum (Q)
  ....:  if sp < sq: return int(-1)
  ....:  elif sp > sq: return int(1)
  ....:  else: return alpha (P, Q)

Sage example in ./programmation.tex, line 3038::

  sage: homogLex ([2, 3, 4, 6, 4], [2, 3, 4, 5, 6])
  -1

Sage example in ./programmation.tex, line 3121::

  sage: def fct1(L):
  ....:   return [list(filter(lambda n: n % 2 == 0, L)),
  ....:           list(filter(lambda n: n % 2 == 1, L))]

Sage example in ./programmation.tex, line 3126::

  sage: fct1([1..10])
  [[2, 4, 6, 8, 10], [1, 3, 5, 7, 9]]

Sage example in ./programmation.tex, line 3145::

  sage: def fct2 (L):
  ....:   res0 = [] ; res1 = []
  ....:   for k in L:
  ....:     if k%2 == 0: res0.append(k) # or res0[len(res0):] = [k]
  ....:     else: res1.append(k)        # or res1[len(res1):] = [k]
  ....:   return [res0, res1]

Sage example in ./programmation.tex, line 3157::

  sage: def fct3a (L, res0, res1):
  ....:   if L == []: return [res0, res1]
  ....:   elif L[0]%2 == 0: return fct3a(L[1:], res0+[L[0]], res1)
  ....:   else: return fct3a (L[1:], res0, res1+[L[0]])

Sage example in ./programmation.tex, line 3163::

  sage: def fct3 (L): return fct3a (L, [], [])

Sage example in ./programmation.tex, line 3195::

  sage: def subSequences (L):
  ....:   if L == []: return []
  ....:   res = [] ; start = 0 ; k = 1
  ....:   while k < len(L):    # 2 consecutive terms are defined
  ....:     if L[k-1] > L[k]:
  ....:       res.append (L[start:k]) ; start = k
  ....:     k = k+1
  ....:   res.append (L[start:k])
  ....:   return res

Sage example in ./programmation.tex, line 3212::

  sage: subSequences([1, 4, 1, 5])
  [[1, 4], [1, 5]]
  sage: subSequences([4, 1, 5, 1])
  [[4], [1, 5], [1]]

Sage example in ./programmation.tex, line 3255::

  sage: S = 'This is a character string.'

Sage example in ./programmation.tex, line 3277::

  sage: S = 'This is a déjà-vu example.'
  sage: print(S)
  This is a déjà-vu example.

Sage example in ./programmation.tex, line 3322::

  sage: S='one two three four five six seven'; L=S.split(); L
  ['one', 'two', 'three', 'four', 'five', 'six', 'seven']

Sage example in ./programmation.tex, line 3368::

  sage: L1 = [11, 22, 33] ; L2 = L1
  sage: L1[1] = 222 ; L2.sort() ; L1, L2
  ([11, 33, 222], [11, 33, 222])
  sage: L1[2:3] = []; L2[0:0] = [6, 7, 8]
  sage: L1, L2
  ([6, 7, 8, 11, 33], [6, 7, 8, 11, 33])

Sage example in ./programmation.tex, line 3422::

  sage: L1 = [11, 22, 33] ; L2 = L1 ; L3 = L1[:]
  sage: [L1 is L2, L2 is L1, L1 is L3, L1 == L3]
  [True, True, False, True]

Sage example in ./programmation.tex, line 3439::

  sage: La = [1, 2, 3] ; L1 = [1, La] ; L2 = copy(L1)
  sage: L1[1][0] = 5         # [1, [5, 2, 3]] for L1 and L2
  sage: [L1 == L2, L1 is L2, L1[1] is L2[1]]
  [True, False, True]

Sage example in ./programmation.tex, line 3498::

  sage: def lexInverse (P, Q):
  ....:   P1 = copy(P) ; P1.reverse()
  ....:   Q1 = copy(Q) ; Q1.reverse()
  ....:   return - alpha (P1, Q1)

Sage example in ./programmation.tex, line 3593::

  sage: S0 = (); S1 = (1, ); S2 = (1, 2)
  sage: [1 in S1, 1 == (1)]
  [True, True]

Sage example in ./programmation.tex, line 3610::

  sage: S1 = (1, 4, 9, 16, 25); [k for k in S1]
  [1, 4, 9, 16, 25]

Sage example in ./programmation.tex, line 3627::

  sage: L1 = [0..4]; L2 = [5..9]
  sage: list(zip(L1, L2))
  [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]
  sage: list(map(lambda x, y:(x, y), L1, L2))
  [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]

Sage example in ./programmation.tex, line 3656::

  sage: E = Set([1, 2, 4, 8, 2, 2, 2]); F = Set([7, 5, 3, 1]); E, F
  ({8, 1, 2, 4}, {1, 3, 5, 7})

Sage example in ./programmation.tex, line 3678::

  sage: E = Set([1, 2, 4, 8, 2, 2, 2]); F = Set([7, 5, 3, 1])
  sage: 5 in E, 5 in F, E + F == F | E
  (False, True, True)
  sage: E & F, E - F, E ^^ F
  ({1}, {8, 2, 4}, {2, 3, 4, 5, 7, 8})

Sage example in ./programmation.tex, line 3700::

  sage: E = Set([1, 2, 4, 8, 2, 2, 2])
  sage: [E[k] for k in [0..len(E)-1]], [t for t in E]
  ([8, 1, 2, 4], [8, 1, 2, 4])

Sage example in ./programmation.tex, line 3713::

  sage: def included (E, F): return E+F == F

Sage example in ./programmation.tex, line 3728::

  sage: sorted(Set([Set([]), Set([1]), Set([2]), Set([1, 2])]), key=str)
  [{1, 2}, {1}, {2}, {}]
  sage: sorted(Set([    (),     (1, ),    (2, ),    (1, 2) ]))
  [(), (1,), (1, 2), (2,)]

Sage example in ./programmation.tex, line 3744::

  sage: def Parts (EE):
  ....:   if EE == Set([]): return Set([EE])
  ....:   else:
  ....:     return withOrWithout (EE[0], Parts(Set(EE[1:])))

Sage example in ./programmation.tex, line 3754::

  sage: def withOrWithout (a, E):
  ....:   return Set (map (lambda F: Set([a])+F, E)) + E

Sage example in ./programmation.tex, line 3762::

  sage: sorted(Parts(Set([1, 2, 3])), key=str)
  [{1, 2, 3}, {1, 2}, {1, 3}, {1}, {2, 3}, {2}, {3}, {}]

Sage example in ./programmation.tex, line 3804::

  sage: D={}; D['one']=1; D['two']=2; D['three']=3; D['ten']=10
  sage: D['two'] + D['three']
  5

Sage example in ./programmation.tex, line 3857::

  sage: D = {'a0':'b0', 'a1':'b1', 'a2':'b2', 'a3':'b0',\
  ....: 'a4':'b3', 'a5':'b3'}
  sage: E  = Set(D.keys()) ; Imf = Set(D.values())
  sage: Imf == Set(map (lambda t:D[t], E))     # is equivalent
  True

Sage example in ./programmation.tex, line 3894::

  sage: def injective(D):
  ....:   return len(D) == len (Set(D.values()))

"""

