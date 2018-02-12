"""
Miscellaneous generic functions

A collection of functions implementing generic algorithms in arbitrary
groups, including additive and multiplicative groups.

In all cases the group operation is specified by a parameter
'operation', which is a string either one of the set of
multiplication_names or addition_names specified below, or 'other'.
In the latter case, the caller must provide an identity, inverse() and
op() functions.

::

    multiplication_names = ( 'multiplication', 'times', 'product', '*')
    addition_names       = ( 'addition', 'plus', 'sum', '+')


Also included are a generic function for computing multiples (or
powers), and an iterator for general multiples and powers.

EXAMPLES:

Some examples in the multiplicative group of a finite field:

- Discrete logs::

    sage: K = GF(3^6,'b')
    sage: b = K.gen()
    sage: a = b^210
    sage: discrete_log(a, b, K.order()-1)
    210

- Linear relation finder::

    sage: F.<a>=GF(3^6,'a')
    sage: a.multiplicative_order().factor()
    2^3 * 7 * 13
    sage: b=a^7
    sage: c=a^13
    sage: linear_relation(b,c,'*')
    (13, 7)
    sage: b^13==c^7
    True

- Orders of elements::

    sage: k.<a> = GF(5^5)
    sage: b = a^4
    sage: order_from_multiple(b,5^5-1,operation='*')
    781
    sage: order_from_bounds(b,(5^4,5^5),operation='*')
    781

Some examples in the group of points of an elliptic curve over a finite field:

- Discrete logs::

    sage: F=GF(37^2,'a')
    sage: E=EllipticCurve(F,[1,1])
    sage: F.<a>=GF(37^2,'a')
    sage: E=EllipticCurve(F,[1,1])
    sage: P=E(25*a + 16 , 15*a + 7 )
    sage: P.order()
    672
    sage: Q=39*P; Q
    (36*a + 32 : 5*a + 12 : 1)
    sage: discrete_log(Q,P,P.order(),operation='+')
    39

- Linear relation finder::

    sage: F.<a>=GF(3^6,'a')
    sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
    sage: P=E(a^5 + a^4 + a^3 + a^2 + a + 2 , 0)
    sage: Q=E(2*a^3 + 2*a^2 + 2*a , a^3 + 2*a^2 + 1)
    sage: linear_relation(P,Q,'+')
    (1, 2)
    sage: P == 2*Q
    True

- Orders of elements::

    sage: k.<a> = GF(5^5)
    sage: E = EllipticCurve(k,[2,4])
    sage: P = E(3*a^4 + 3*a , 2*a + 1 )
    sage: M = E.cardinality(); M
    3227
    sage: plist = M.prime_factors()
    sage: order_from_multiple(P, M, plist, operation='+')
    3227
    sage: Q = E(0,2)
    sage: order_from_multiple(Q, M, plist, operation='+')
    7
    sage: order_from_bounds(Q, Hasse_bounds(5^5), operation='+')
    7
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                          John Cremona  <john.cremona@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy

from sage.misc.all import prod
import sage.rings.integer_ring as integer_ring
import sage.rings.integer
from sage.arith.srange import xsrange

#
# Lists of names (as strings) which the user may use to identify one
# of the standard operations:
#
multiplication_names = ( 'multiplication', 'times', 'product', '*')
addition_names       = ( 'addition', 'plus', 'sum', '+')


# deprecation(24256)
from sage.structure.element import generic_power as power


def multiple(a, n, operation='*', identity=None, inverse=None, op=None):
    r"""
    Returns either `na` or `a^n`, where `n` is any integer and `a` is
    a Python object on which a group operation such as addition or
    multiplication is defined.  Uses the standard binary algorithm.

    INPUT:  See the documentation for ``discrete_logarithm()``.

    EXAMPLES::

        sage: multiple(2,5)
        32
        sage: multiple(RealField()('2.5'),4)
        39.0625000000000
        sage: multiple(2,-3)
        1/8
        sage: multiple(2,100,'+') == 100*2
        True
        sage: multiple(2,100) == 2**100
        True
        sage: multiple(2,-100,) == 2**-100
        True
        sage: R.<x>=ZZ[]
        sage: multiple(x,100)
        x^100
        sage: multiple(x,100,'+')
        100*x
        sage: multiple(x,-10)
        1/x^10

    Idempotence is detected, making the following fast::

        sage: multiple(1,10^1000)
        1

        sage: E=EllipticCurve('389a1')
        sage: P=E(-1,1)
        sage: multiple(P,10,'+')
        (645656132358737542773209599489/22817025904944891235367494656 : 525532176124281192881231818644174845702936831/3446581505217248068297884384990762467229696 : 1)
        sage: multiple(P,-10,'+')
        (645656132358737542773209599489/22817025904944891235367494656 : -528978757629498440949529703029165608170166527/3446581505217248068297884384990762467229696 : 1)


    """
    from operator import inv, mul, neg, add

    if operation in multiplication_names:
        identity = a.parent()(1)
        inverse  = inv
        op = mul
    elif operation in addition_names:
        identity = a.parent()(0)
        inverse  = neg
        op = add
    else:
        if identity is None or inverse is None or op is None:
            raise ValueError("identity, inverse and operation must all be specified")

    if n == 0:
        return identity

    if n < 0:
        n = -n
        a = inverse(a)

    if n == 1:
        return a

    # check for idempotence, and store the result otherwise
    aa = op(a,a)
    if aa == a:
        return a

    if n == 2:
        return aa

    if n == 3:
        return op(aa,a)

    if n == 4:
        return op(aa,aa)

    # since we've computed a^2, let's start squaring there
    # so, let's keep the least-significant bit around, just
    # in case.
    m = n & 1
    n = n >> 1

    # One multiplication can be saved by starting with
    # the second-smallest power needed rather than with 1
    # we've already squared a, so let's start there.
    apow = aa
    while n&1 == 0:
        apow = op(apow,apow)
        n = n >> 1
    power = apow
    n = n >> 1

    # now multiply that least-significant bit in...
    if m:
        power = op(power,a)

    # and this is straight from the book.
    while n != 0:
        apow = op(apow,apow)
        if n&1 != 0:
            power = op(power,apow)
        n = n >> 1

    return power


#
# Generic iterator for looping through multiples or powers
#

class multiples:
    r"""
    Return an iterator which runs through ``P0+i*P`` for ``i`` in ``range(n)``.

    ``P`` and ``P0`` must be Sage objects in some group; if the operation is
    multiplication then the returned values are instead ``P0*P**i``.

    EXAMPLES::

        sage: list(multiples(1,10))
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: list(multiples(1,10,100))
        [100, 101, 102, 103, 104, 105, 106, 107, 108, 109]

        sage: E=EllipticCurve('389a1')
        sage: P=E(-1,1)
        sage: for Q in multiples(P,5): print((Q, Q.height()/P.height()))
        ((0 : 1 : 0), 0.000000000000000)
        ((-1 : 1 : 1), 1.00000000000000)
        ((10/9 : -35/27 : 1), 4.00000000000000)
        ((26/361 : -5720/6859 : 1), 9.00000000000000)
        ((47503/16641 : 9862190/2146689 : 1), 16.0000000000000)

        sage: R.<x>=ZZ[]
        sage: list(multiples(x,5))
        [0, x, 2*x, 3*x, 4*x]
        sage: list(multiples(x,5,operation='*'))
        [1, x, x^2, x^3, x^4]
        sage: list(multiples(x,5,indexed=True))
        [(0, 0), (1, x), (2, 2*x), (3, 3*x), (4, 4*x)]
        sage: list(multiples(x,5,indexed=True,operation='*'))
        [(0, 1), (1, x), (2, x^2), (3, x^3), (4, x^4)]
        sage: for i,y in multiples(x,5,indexed=True): print("%s  times %s = %s"%(i,x,y))
        0  times x = 0
        1  times x = x
        2  times x = 2*x
        3  times x = 3*x
        4  times x = 4*x

        sage: for i,n in multiples(3,5,indexed=True,operation='*'):  print("3 to the power %s = %s" % (i,n))
        3 to the power 0 = 1
        3 to the power 1 = 3
        3 to the power 2 = 9
        3 to the power 3 = 27
        3 to the power 4 = 81
    """
    def __init__(self,P,n,P0=None,indexed=False, operation='+', op=None):
        """
        Create a multiples iterator

        INPUT:

        - ``P`` - step value: any Sage object on which a binary
                             operation is defined
        - ``n`` - number of multiples: non-negative integer
        - ``P0`` - offset (default 0): Sage object which can be 'added' to P
        - ``indexed`` - boolean (default False)

          If ``indexed==False`` then the iterator delivers ``P0+i*P``
          (if ``operation=='+'``) or ``P0*P**i`` (if
          ``operation=='*'``), for ``i`` in ``range(n)``.

          If ``indexed==True`` then the iterator delivers tuples
          ``(i,P0+i*P)`` or ``(i,P0*P**i)``.

        - ``operation`` - string: '+' (default ) or '*' or `other`.

          If `other`, a function ``op()`` must be supplied (a function
          of 2 arguments) defining the group binary operation; also
          ``P0`` must be supplied.
        """
        if n<0:
            raise ValueError('n cannot be negative in multiples')

        from operator import mul, add

        if operation in multiplication_names:
            if P0 is None: P0 = P.parent()(1)
            self.op = mul
        elif operation in addition_names:
            if P0 is None: P0 = P.parent()(0)
            self.op = add
        else:
            self.op = op
            if P0 is None:
                raise ValueError("P0 must be supplied when operation is neither addition nor multiplication")
            if op is None:
                raise ValueError("op() must both be supplied when operation is neither addition nor multiplication")

        self.P=copy(P)
        self.Q=copy(P0)
        assert self.P is not None and self.Q is not None
        self.i = 0
        self.bound = n
        self.indexed = indexed


    def __next__(self):
        """
        Returns the next item in this multiples iterator.
        """
        if self.i >= self.bound:
            raise StopIteration
        i = self.i
        val = self.Q
        self.i +=1
        self.Q=self.op(self.Q,self.P)
        if self.indexed:
            return (i,val)
        else:
            return val

    next = __next__

    def __iter__(self):
        """
        Standard member function making this class an iterator.
        """
        return self


def bsgs(a, b, bounds, operation='*', identity=None, inverse=None, op=None):
    r"""
    Totally generic discrete baby-step giant-step function.

    Solves `na=b` (or `a^n=b`) with `lb\le n\le ub` where ``bounds==(lb,ub)``,
    raising an error if no such `n` exists.

    `a` and `b` must be elements of some group with given identity,
    inverse of ``x`` given by ``inverse(x)``, and group operation on
    ``x``, ``y`` by ``op(x,y)``.

    If operation is '*' or '+' then the other
    arguments are provided automatically; otherwise they must be
    provided by the caller.

    INPUT:

    - ``a``    - group element
    - ``b``    - group element
    - ``bounds`` - a 2-tuple of integers ``(lower,upper)`` with ``0<=lower<=upper``
    - ``operation`` - string: '*', '+', 'other'
    - ``identity`` - the identity element of the group
    - ``inverse()``  - function of 1 argument ``x`` returning inverse of ``x``
    - ``op()`` - function of 2 arguments ``x``, ``y`` returning ``x*y`` in group

    OUTPUT:

    An integer `n` such that `a^n = b` (or `na = b`).  If no
    such `n` exists, this function raises a ValueError exception.

    NOTE: This is a generalization of discrete logarithm.  One
    situation where this version is useful is to find the order of
    an element in a group where we only have bounds on the group
    order (see the elliptic curve example below).

    ALGORITHM: Baby step giant step.  Time and space are soft
    `O(\sqrt{n})` where `n` is the difference between upper and lower
    bounds.

    EXAMPLES::

        sage: b = Mod(2,37);  a = b^20
        sage: bsgs(b, a, (0,36))
        20

        sage: p=next_prime(10^20)
        sage: a=Mod(2,p); b=a^(10^25)
        sage: bsgs(a, b, (10^25-10^6,10^25+10^6)) == 10^25
        True

        sage: K = GF(3^6,'b')
        sage: a = K.gen()
        sage: b = a^210
        sage: bsgs(a, b, (0,K.order()-1))
        210

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^500
        sage: bsgs(z,w,(0,229))
        40

    An additive example in an elliptic curve group::

        sage: F.<a> = GF(37^5)
        sage: E = EllipticCurve(F, [1,1])
        sage: P = E.lift_x(a); P
        (a : 28*a^4 + 15*a^3 + 14*a^2 + 7 : 1)

    This will return a multiple of the order of P::

        sage: bsgs(P,P.parent()(0),Hasse_bounds(F.order()),operation='+')
        69327408

    AUTHOR:

        - John Cremona (2008-03-15)
    """
    Z = integer_ring.ZZ

    from operator import inv, mul, neg, add

    if operation in multiplication_names:
        identity = a.parent()(1)
        inverse  = inv
        op = mul
    elif operation in addition_names:
        identity = a.parent()(0)
        inverse  = neg
        op = add
    else:
        if identity is None or inverse is None or op is None:
            raise ValueError("identity, inverse and operation must be given")

    lb, ub = bounds
    if lb<0 or ub<lb:
        raise ValueError("bsgs() requires 0<=lb<=ub")

    if a.is_zero() and not b.is_zero():
        raise ValueError("No solution in bsgs()")

    ran = 1 + ub - lb   # the length of the interval

    c = op(inverse(b),multiple(a,lb,operation=operation))

    if ran < 30:    # use simple search for small ranges
        i = lb
        d = c
#        for i,d in multiples(a,ran,c,indexed=True,operation=operation):
        for i0 in range(ran):
            i = lb + i0
            if identity == d:        # identity == b^(-1)*a^i, so return i
                return Z(i)
            d = op(a,d)
        raise ValueError("No solution in bsgs()")

    m = ran.isqrt()+1  # we need sqrt(ran) rounded up
    table = dict()     # will hold pairs (a^(lb+i),lb+i) for i in range(m)

    d=c
    for i0 in xsrange(m):
        i = lb + i0
        if identity==d:        # identity == b^(-1)*a^i, so return i
            return Z(i)
        table[d] = i
        d=op(d,a)

    c = op(c,inverse(d))     # this is now a**(-m)
    d=identity
    for i in xsrange(m):
        j = table.get(d)
        if j is not None:  # then d == b*a**(-i*m) == a**j
            return Z(i*m + j)
        d=op(c,d)

    raise ValueError("Log of %s to the base %s does not exist in %s."%(b,a,bounds))

def discrete_log_rho(a, base, ord=None, operation='*', hash_function=hash):
    """
    Pollard Rho algorithm for computing discrete logarithm in cyclic
    group of prime order.
    If the group order is very small it falls back to the baby step giant step
    algorithm.

    INPUT:

    - ``a`` -- a group element
    - ``base`` -- a group element
    - ``ord`` -- the order of ``base`` or ``None``, in this case we try
      to compute it
    - ``operation`` -- a string (default: ``'*'``) denoting whether we
      are in an additive group or a multiplicative one
    - ``hash_function`` -- having an efficient hash function is critical
      for this algorithm (see examples)

    OUTPUT: an integer `n` such that `a = base^n` (or `a = n*base`)

    ALGORITHM: Pollard rho for discrete logarithm, adapted from the
    article of Edlyn Teske, 'A space efficient algorithm for group
    structure computation'.

    EXAMPLES::

        sage: F.<a> = GF(2^13)
        sage: g = F.gen()
        sage: discrete_log_rho(g^1234, g)
        1234

        sage: F.<a> = GF(37^5)
        sage: E = EllipticCurve(F, [1,1])
        sage: G = (3*31*2^4)*E.lift_x(a)
        sage: discrete_log_rho(12345*G, G, ord=46591, operation='+')
        12345

    It also works with matrices::

        sage: A = matrix(GF(50021),[[10577,23999,28893],[14601,41019,30188],[3081,736,27092]])
        sage: discrete_log_rho(A^1234567, A)
        1234567

    Beware, the order must be prime::

        sage: I = IntegerModRing(171980)
        sage: discrete_log_rho(I(2), I(3))
        Traceback (most recent call last):
        ...
        ValueError: for Pollard rho algorithm the order of the group must be prime

    If it fails to find a suitable logarithm, it raises a ``ValueError``::

        sage: I = IntegerModRing(171980)
        sage: discrete_log_rho(I(31002),I(15501))
        Traceback (most recent call last):
        ...
        ValueError: Pollard rho algorithm failed to find a logarithm

    The main limitation on the hash function is that we don't want to have
    `hash(x*y) = hash(x) + hash(y)`::

        sage: I = IntegerModRing(next_prime(2^23))
        sage: def test():
        ....:     try:
        ....:          discrete_log_rho(I(123456),I(1),operation='+')
        ....:     except Exception:
        ....:          print("FAILURE")
        sage: test()  # random failure
        FAILURE

    If this happens, we can provide a better hash function::

        sage: discrete_log_rho(I(123456),I(1),operation='+', hash_function=lambda x: hash(x*x))
        123456

    AUTHOR:

    - Yann Laigle-Chapuy (2009-09-05)

    """
    from six.moves import range
    from sage.rings.integer import Integer
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from operator import mul, add, pow

    # should be reasonable choices
    partition_size = 20
    memory_size = 4

    if operation in addition_names:
        mult = add
        power = mul
        if ord is None:
            ord = base.additive_order()
    elif operation in multiplication_names:
        mult = mul
        power = pow
        if ord is None:
            ord = base.multiplicative_order()
    else:
        raise ValueError

    ord = Integer(ord)

    if not ord.is_prime():
        raise ValueError("for Pollard rho algorithm the order of the group must be prime")

    # check if we need to set immutable before hashing
    mut = hasattr(base,'set_immutable')

    isqrtord=ord.isqrt()

    if isqrtord < partition_size: #setup to costly, use bsgs
        return bsgs(base,a, bounds=(0,ord), operation=operation)

    reset_bound = 8*isqrtord # we take some margin

    I=IntegerModRing(ord)

    for s in range(10): # to avoid infinite loops
        # random walk function setup
        m=[I.random_element() for i in range(partition_size)]
        n=[I.random_element() for i in range(partition_size)]
        M=[mult(power(base,Integer(m[i])),power(a,Integer(n[i]))) for i in range(partition_size)]

        ax = I.random_element()
        x = power(base,Integer(ax))
        if mut:
            x.set_immutable()

        bx = I(0)

        sigma=[(0,None)]*memory_size
        H={} # memory
        i0=0
        nextsigma = 0
        for i in range(reset_bound):
                    #random walk, we need an efficient hash
            s=hash_function(x) % partition_size
            (x,ax,bx) = (mult(M[s],x), ax+m[s], bx+n[s])
            if mut:
                x.set_immutable()
            # look for collisions
            if x in H:
                ay,by=H[x]
                if bx == by:
                    break
                else:
                    res = sage.rings.integer.Integer((ay-ax)/(bx-by))
                    if power(base,res) == a:
                        return res
                    else:
                        break
            # should we remember this value?
            elif i >= nextsigma:
                if sigma[i0][1] is not None:
                    H.pop(sigma[i0][1])
                sigma[i0]=(i,x)
                i0 = (i0+1) % memory_size
                nextsigma = 3*sigma[i0][0] #3 seems a good choice
                H[x]=(ax,bx)

    raise ValueError("Pollard rho algorithm failed to find a logarithm")

def discrete_log(a, base, ord=None, bounds=None, operation='*', identity=None, inverse=None, op=None):
    r"""
    Totally generic discrete log function.

    INPUT:

    - ``a``    - group element
    - ``base`` - group element (the base)
    - ``ord``  - integer (multiple of order of base, or ``None``)
    - ``bounds`` - a priori bounds on the log
    - ``operation`` - string: '*', '+', 'other'
    - ``identity`` - the group's identity
    - ``inverse()`` - function of 1 argument ``x`` returning inverse of ``x``
    - ``op()`` - function of 2 arguments ``x``, ``y`` returning ``x*y`` in group

    ``a`` and ``base`` must be elements of some group with identity
    given by identity, inverse of ``x`` by ``inverse(x)``, and group
    operation on ``x``, ``y`` by ``op(x,y)``.

    If operation is '*' or '+' then the other
    arguments are provided automatically; otherwise they must be
    provided by the caller.

    OUTPUT: Returns an integer `n` such that `b^n = a` (or `nb = a`),
    assuming that ``ord`` is a multiple of the order of the base `b`.
    If ``ord`` is not specified, an attempt is made to compute it.

    If no such `n` exists, this function raises a ValueError exception.

    .. warning::

       If ``x`` has a log method, it is likely to be vastly faster
       than using this function.  E.g., if ``x`` is an integer modulo
       `n`, use its log method instead!

    ALGORITHM: Pohlig-Hellman and Baby step giant step.

    EXAMPLES::

        sage: b = Mod(2,37);  a = b^20
        sage: discrete_log(a, b)
        20
        sage: b = Mod(2,997);  a = b^20
        sage: discrete_log(a, b)
        20

        sage: K = GF(3^6,'b')
        sage: b = K.gen()
        sage: a = b^210
        sage: discrete_log(a, b, K.order()-1)
        210

        sage: b = Mod(1,37);  x = Mod(2,37)
        sage: discrete_log(x, b)
        Traceback (most recent call last):
        ...
        ValueError: No discrete log of 2 found to base 1
        sage: b = Mod(1,997);  x = Mod(2,997)
        sage: discrete_log(x, b)
        Traceback (most recent call last):
        ...
        ValueError: No discrete log of 2 found to base 1

    See :trac:`2356`::

        sage: F.<w> = GF(121)
        sage: v = w^120
        sage: v.log(w)
        0

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^50
        sage: discrete_log(w,z)
        50

    An example where the order is infinite: note that we must give
    an upper bound here::

        sage: K.<a> = QuadraticField(23)
        sage: eps = 5*a-24        # a fundamental unit
        sage: eps.multiplicative_order()
        +Infinity
        sage: eta = eps^100
        sage: discrete_log(eta,eps,bounds=(0,1000))
        100

    In this case we cannot detect negative powers::

        sage: eta = eps^(-3)
        sage: discrete_log(eta,eps,bounds=(0,100))
        Traceback (most recent call last):
        ...
        ValueError: No discrete log of -11515*a - 55224 found to base 5*a - 24

    But we can invert the base (and negate the result) instead::

        sage: - discrete_log(eta^-1,eps,bounds=(0,100))
        -3

    An additive example: elliptic curve DLOG::

        sage: F=GF(37^2,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: F.<a>=GF(37^2,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: P=E(25*a + 16 , 15*a + 7 )
        sage: P.order()
        672
        sage: Q=39*P; Q
        (36*a + 32 : 5*a + 12 : 1)
        sage: discrete_log(Q,P,P.order(),operation='+')
        39

    An example of big smooth group::

        sage: F.<a>=GF(2^63)
        sage: g=F.gen()
        sage: u=g**123456789
        sage: discrete_log(u,g)
        123456789

    AUTHORS:

    - William Stein and David Joyner (2005-01-05)
    - John Cremona (2008-02-29) rewrite using ``dict()`` and make generic

    """
    if ord is None:
        if operation in multiplication_names:
            try:
                ord = base.multiplicative_order()
            except Exception:
                ord = base.order()
        elif operation in addition_names:
            try:
                ord = base.additive_order()
            except Exception:
                ord = base.order()
        else:
            try:
                ord = base.order()
            except Exception:
                raise ValueError("ord must be specified")
    try:
        from sage.rings.infinity import Infinity
        if ord==+Infinity:
            return bsgs(base,a,bounds, operation=operation)
        if ord==1 and a!=base:
            raise ValueError
        f=ord.factor()
        l=[0]*len(f)
        for i,(pi,ri) in enumerate(f):
            for j in range(ri):
                if operation in multiplication_names:
                    c=bsgs(base**(ord//pi),(a/base**l[i])**(ord//pi**(j+1)),(0,pi),operation=operation)
                    l[i] += c*(pi**j)
                elif operation in addition_names:
                    c=bsgs(base*(ord//pi),(a-base*l[i])*(ord//pi**(j+1)),(0,pi),operation=operation)
                    l[i] += c*(pi**j)
        from sage.arith.all import CRT_list
        return  CRT_list(l,[pi**ri for pi,ri in f])
    except ValueError:
        raise ValueError("No discrete log of %s found to base %s"%(a,base))

def discrete_log_generic(a, base, ord=None, bounds=None, operation='*', identity=None, inverse=None, op=None):
    """
    Alias for ``discrete_log``.
    """
    return discrete_log(a, base, ord=None, bounds=None, operation='*', identity=None, inverse=None, op=None)

def discrete_log_lambda(a, base, bounds, operation='*', hash_function=hash):
    """
    Pollard Lambda algorithm for computing discrete logarithms. It uses
    only a logarithmic amount of memory. It's useful if you have
    bounds on the logarithm. If you are computing logarithms in a
    whole finite group, you should use Pollard Rho algorithm.

    INPUT:

    - a -- a group element
    - base -- a group element
    - bounds -- a couple (lb,ub) representing the range where we look for a logarithm
    - operation -- string: '+', '*' or 'other'
    - hash_function -- having an efficient hash function is critical for this algorithm

    OUTPUT: Returns an integer `n` such that `a=base^n` (or `a=n*base`)

    ALGORITHM: Pollard Lambda, if bounds are (lb,ub) it has time complexity
        O(sqrt(ub-lb)) and space complexity O(log(ub-lb))

    EXAMPLES::

        sage: F.<a> = GF(2^63)
        sage: discrete_log_lambda(a^1234567, a, (1200000,1250000))
        1234567

        sage: F.<a> = GF(37^5)
        sage: E = EllipticCurve(F, [1,1])
        sage: P = E.lift_x(a); P
        (a : 28*a^4 + 15*a^3 + 14*a^2 + 7 : 1)

    This will return a multiple of the order of P::

        sage: discrete_log_lambda(P.parent()(0), P, Hasse_bounds(F.order()), operation='+')
        69327408

        sage: K.<a> = GF(89**5)
        sage: hs = lambda x: hash(x) + 15
        sage: discrete_log_lambda(a**(89**3 - 3), a, (89**2, 89**4), operation = '*', hash_function = hs)  # long time (10s on sage.math, 2011)
        704966

    AUTHOR:

        -- Yann Laigle-Chapuy (2009-01-25)

    """
    from six.moves import range
    from sage.rings.integer import Integer
    from operator import mul, add, pow

    if operation in addition_names:
        mult=add
        power=mul
    elif operation in multiplication_names:
        mult=mul
        power=pow
    else:
        raise ValueError("unknown operation")

    lb,ub = bounds
    if lb<0 or ub<lb:
        raise ValueError("discrete_log_lambda() requires 0<=lb<=ub")

    # check for mutability
    mut = hasattr(base,'set_immutable')

    width = Integer(ub-lb)
    N = width.isqrt()+1

    M = dict()
    for s in range(10): #to avoid infinite loops
        #random walk function setup
        k = 0
        while (2**k<N):
            r = sage.misc.prandom.randrange(1,N)
            M[k] = (r , power(base,r))
            k += 1
        #first random walk
        H = power(base,ub)
        c = ub
        for i in range(N):
            if mut: H.set_immutable()
            r,e = M[hash_function(H)%k]
            H = mult(H,e)
            c += r
        if mut: H.set_immutable()
        mem=set([H])
        #second random walk
        H = a
        d=0
        while c-d >= lb:
            if mut: H.set_immutable()
            if ub > c-d and H in mem:
                return c-d
            r,e = M[hash_function(H)%k]
            H = mult(H,e)
            d += r

    raise ValueError("Pollard Lambda failed to find a log")


################################################################
#
# Generic linear relation finder
#
################################################################

def linear_relation(P, Q, operation='+', identity=None, inverse=None, op=None):
    r"""
    Function which solves the equation ``a*P=m*Q`` or ``P^a=Q^m``.

    Additive version: returns `(a,m)` with minimal `m>0` such that
    `aP=mQ`.  Special case: if `\left<P\right>` and `\left<Q\right>`
    intersect only in `\{0\}` then `(a,m)=(0,n)` where `n` is
    ``Q.additive_order()``.

    Multiplicative version: returns `(a,m)` with minimal `m>0` such
    that `P^a=Q^m`.  Special case: if `\left<P\right>` and
    `\left<Q\right>` intersect only in `\{1\}` then `(a,m)=(0,n)`
    where `n` is ``Q.multiplicative_order()``.

    ALGORITHM:

    Uses the generic ``bsgs()`` function, and so works in general
    finite abelian groups.

    EXAMPLES:

    An additive example (in an elliptic curve group)::

        sage: F.<a>=GF(3^6,'a')
        sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
        sage: P=E(a^5 + a^4 + a^3 + a^2 + a + 2 , 0)
        sage: Q=E(2*a^3 + 2*a^2 + 2*a , a^3 + 2*a^2 + 1)
        sage: linear_relation(P,Q,'+')
        (1, 2)
        sage: P == 2*Q
        True

    A multiplicative example (in a finite field's multiplicative group)::

        sage: F.<a>=GF(3^6,'a')
        sage: a.multiplicative_order().factor()
        2^3 * 7 * 13
        sage: b=a^7
        sage: c=a^13
        sage: linear_relation(b,c,'*')
        (13, 7)
        sage: b^13==c^7
        True
    """

    from operator import mul, add
    Z = integer_ring.ZZ

    if operation in multiplication_names:
        op = mul
        try:
            n = P.multiplicative_order()
            m = Q.multiplicative_order()
        except Exception:
            n = P.order()
            m = Q.order()
    elif operation in addition_names:
        op = add
        try:
            n = P.additive_order()
            m = Q.additive_order()
        except Exception:
            n = P.order()
            m = Q.order()
    else:
        if op is None:
            raise ValueError("operation must be specified")
        n = P.order()
        m = Q.order()

    g = sage.arith.all.gcd(n,m)
    if g==1: return (m,Z(0))
    n1 = n//g
    m1 = m//g
    P1 = multiple(P,n1,operation=operation)  # has exact order g
    Q1 = multiple(Q,m1,operation=operation)  # has exact order g

    # now see if Q1 is a multiple of P1; the only multiples we
    # need check are h*Q1 where h divides g
    for h in g.divisors(): # positive divisors!
        try:
            Q2 = multiple(Q1,h,operation=operation)
            return (n1 * bsgs(P1,Q2,(0,g-1),operation=operation),
                    m1 * h)
        except ValueError:
            pass # to next h
    raise ValueError("No solution found in linear_relation!")

################################################################
#
# Generic functions to find orders of elements
#
# 1. order_from_multiple: finds the order given a multiple of the order
#
# 2. order_from_bounds: finds the order given an interval containing a
# multiple of the order
#
################################################################

def order_from_multiple(P, m, plist=None, factorization=None, check=True,
                        operation='+'):
    r"""
    Generic function to find order of a group element given a multiple
    of its order.

    INPUT:

    - ``P`` - a Sage object which is a group element;
    - ``m`` - a Sage integer which is a multiple of the order of ``P``,
      i.e. we require that ``m*P=0`` (or ``P**m=1``);
    - ``check`` - a Boolean (default:True), indicating whether we check if ``m``
      really is a multiple of the order;
    - ``factorization`` - the factorization of ``m``, or ``None`` in which
      case this function will need to factor ``m``;
    - ``plist`` - a list of the prime factors of ``m``, or ``None`` - kept for compatibility only,
      prefer the use of ``factorization``;
    - ``operation`` - string: '+' (default) or '*'.

    .. note::

       It is more efficient for the caller to factor ``m`` and cache
       the factors for subsequent calls.

    EXAMPLES::

        sage: k.<a> = GF(5^5)
        sage: b = a^4
        sage: order_from_multiple(b,5^5-1,operation='*')
        781
        sage: E = EllipticCurve(k,[2,4])
        sage: P = E(3*a^4 + 3*a , 2*a + 1 )
        sage: M = E.cardinality(); M
        3227
        sage: F = M.factor()
        sage: order_from_multiple(P, M, factorization=F, operation='+')
        3227
        sage: Q = E(0,2)
        sage: order_from_multiple(Q, M, factorization=F, operation='+')
        7

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^50
        sage: order_from_multiple(w,230,operation='*')
        23

        sage: F=GF(2^1279,'a')
        sage: n=F.cardinality()-1 # Mersenne prime
        sage: order_from_multiple(F.random_element(),n,factorization=[(n,1)],operation='*')==n
        True

        sage: K.<a> = GF(3^60)
        sage: order_from_multiple(a, 3^60-1, operation='*', check=False)
        42391158275216203514294433200
    """
    from operator import mul, add
    Z = integer_ring.ZZ

    if operation in multiplication_names:
        identity = P.parent()(1)
    elif operation in addition_names:
        identity = P.parent()(0)
    else:
        raise ValueError("unknown group operation")

    if P == identity:
        return Z(1)

    M=Z(m)
    if check:
        assert multiple(P,M,operation=operation) == identity

    if factorization:
        F = factorization
    elif plist:
        F = [(p,M.valuation(p)) for p in plist]
    else:
        F = M.factor()

    if len(F) == 1 and list(F) == [(M,1)]:
        return M

    # Efficiency improvement (2009-10-27, implemented by Yann Laigle-Chapuy):
    # we use an internal recursive function to avoid unnecessary computations.
    def _order_from_multiple_helper(Q, L, S):
        """
        internal use, to minimize the number of group operations.
        """
        l = len(L)
        if l == 1:
            # we determine the power of p dividing the order,

            # Efficiency improvement (2009-04-01, suggested by Ryan Hinton,
            # implemented by John Cremona): avoid the last multiplication by p.
            # For example, if M itself is prime the code used to compute M*P
            # twice (unless P=0), now it does it once.
            p,e = L[0]
            e0 = 0
            while (Q != identity) and (e0<e-1):
                Q = multiple(Q,p,operation=operation)
                e0 += 1
            if (Q != identity):
                e0 += 1
            return p**e0
        else:
            # try to split the list wisely
            sum_left = 0
            i = 0
            for k in range(l):
                p,e = L[k]
                # multiplying by p**e require roughly 'e log_2(p) / 2' additions
                v = e * sage.functions.log.log(float(p))
                if abs(sum_left + v - (S / 2)) > abs(sum_left - (S / 2)):
                    break
                sum_left += v
            L1 = L[:k]
            L2 = L[k:]
            # recursive calls
            o1 = _order_from_multiple_helper(
                multiple(Q, prod([p**e for p,e in L2]), operation),
                L1,
                sum_left)
            o2 = _order_from_multiple_helper(
                multiple(Q, o1                                       , operation),
                L2,
                S-sum_left)
            return o1*o2

    return _order_from_multiple_helper(P, F, sage.functions.log.log(float(M)) )




def order_from_bounds(P, bounds, d=None, operation='+',
                         identity=None, inverse=None, op=None):
    r"""
    Generic function to find order of a group element, given only
    upper and lower bounds for a multiple of the order (e.g. bounds on
    the order of the group of which it is an element)

    INPUT:

    - ``P``      - a Sage object which is a group element

    - ``bounds`` - a 2-tuple ``(lb,ub)`` such that ``m*P=0`` (or
      ``P**m=1``) for some ``m`` with ``lb<=m<=ub``.

    - ``d`` - (optional) a positive integer; only ``m`` which are
      multiples of this will be considered.

    - ``operation`` - string: '+' (default ) or '*' or other.
      If other, the following must be supplied:

      - ``identity``: the identity element for the group;
      - ``inverse()``: a function of one argument giving the inverse
        of a group element;
      - ``op()``: a function of 2 arguments defining the group binary
        operation.


    .. note::

       Typically ``lb`` and ``ub`` will be bounds on the group order,
       and from previous calculation we know that the group order is
       divisible by ``d``.

    EXAMPLES::

        sage: k.<a> = GF(5^5)
        sage: b = a^4
        sage: order_from_bounds(b,(5^4,5^5),operation='*')
        781
        sage: E = EllipticCurve(k,[2,4])
        sage: P = E(3*a^4 + 3*a , 2*a + 1 )
        sage: bounds = Hasse_bounds(5^5)
        sage: Q = E(0,2)
        sage: order_from_bounds(Q, bounds, operation='+')
        7
        sage: order_from_bounds(P, bounds, 7, operation='+')
        3227

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^50
        sage: order_from_bounds(w,(200,250),operation='*')
        23

    """
    from operator import mul, add
    Z = integer_ring.ZZ

    if operation in multiplication_names:
        op = mul
        identity = P.parent()(1)
    elif operation in addition_names:
        op = add
        identity = P.parent()(0)
    else:
        if op is None:
            raise ValueError("operation and identity must be specified")

    Q = P
    if d is None: d = 1
    if d > 1:
        Q = multiple(P,d,operation=operation)
        lb, ub = bounds
        bounds = ( sage.arith.all.integer_ceil(lb/d),
                   sage.arith.all.integer_floor(ub/d) )

    # Use generic bsgs to find  n=d*m with lb<=n<=ub and n*P=0

    m = d * bsgs(Q, identity, bounds, operation=operation)

    # Now use the order_from_multiple() function to finish the job:

    return order_from_multiple(P, m, operation=operation, check=False)

def merge_points(P1,P2, operation='+',
                         identity=None, inverse=None, op=None, check=True):
    r"""
    Returns a group element whose order is the lcm of the given elements.

    INPUT:

    - ``P1`` -- a pair `(g_1,n_1)` where `g_1` is a group element of order `n_1`
    - ``P2`` -- a pair `(g_2,n_2)` where `g_2` is a group element of order `n_2`
    - ``operation`` -- string: '+' (default ) or '*' or other. If
      other, the following must be supplied:

        - ``identity``: the identity element for the group;
        - ``inverse()``: a function of one argument giving the inverse
          of a group element;
        - ``op()``: a function of 2 arguments defining the group
           binary operation.


    OUTPUT:

    A pair `(g_3,n_3)` where `g_3` has order `n_3=\hbox{lcm}(n_1,n_2)`.

    EXAMPLES::

        sage: F.<a>=GF(3^6,'a')
        sage: b = a^7
        sage: c = a^13
        sage: ob = (3^6-1)//7
        sage: oc = (3^6-1)//13
        sage: merge_points((b,ob),(c,oc),operation='*')
        (a^4 + 2*a^3 + 2*a^2, 728)
        sage: d,od = merge_points((b,ob),(c,oc),operation='*')
        sage: od == d.multiplicative_order()
        True
        sage: od == lcm(ob,oc)
        True

        sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
        sage: P=E(2*a^5 + 2*a^4 + a^3 + 2 , a^4 + a^3 + a^2 + 2*a + 2)
        sage: P.order()
        7
        sage: Q=E(2*a^5 + 2*a^4 + 1 , a^5 + 2*a^3 + 2*a + 2 )
        sage: Q.order()
        4
        sage: R,m = merge_points((P,7),(Q,4), operation='+')
        sage: R.order() == m
        True
        sage: m == lcm(7,4)
        True
    """
    from operator import mul, add
    Z = integer_ring.ZZ

    g1, n1 = P1
    g2, n2 = P2

    if operation in multiplication_names:
        op = mul
        identity = g1.parent()(1)
    elif operation in addition_names:
        op = add
        identity = g1.parent()(0)
    else:
        if op is None:
            raise ValueError("operation and identity must be specified")

    if check:
        assert multiple(g1,n1,operation=operation) == identity
        assert multiple(g2,n2,operation=operation) == identity

    # trivial cases
    if n1.divides(n2):
        return (g2,n2)
    if n2.divides(n1):
        return (g1,n1)

    m,k1,k2 = sage.arith.all.xlcm(n1,n2)
    m1 = n1//k1
    m2 = n2//k2
    g1 = multiple(g1,m1,operation=operation)
    g2 = multiple(g2,m2,operation=operation)
    return (op(g1,g2), m)

def structure_description(G, latex=False):
    r"""
    Return a string that tries to describe the structure of ``G``.

    This methods wraps GAP's ``StructureDescription`` method.

    Requires the *optional* ``database_gap`` package.

    For full details, including the form of the returned string and the
    algorithm to build it, see `GAP's documentation
    <http://www.gap-system.org/Manuals/doc/ref/chap39.html>`_.

    INPUT:

    - ``latex`` -- a boolean (default: ``False``). If ``True`` return a
      LaTeX formatted string.

    OUTPUT:

    - string

    .. WARNING::

        From GAP's documentation: The string returned by
        ``StructureDescription`` is **not** an isomorphism invariant:
        non-isomorphic groups can have the same string value, and two
        isomorphic groups in different representations can produce different
        strings.

    EXAMPLES::

        sage: G = CyclicPermutationGroup(6)
        sage: G.structure_description()             # optional - database_gap
        'C6'
        sage: G.structure_description(latex=True)   # optional - database_gap
        'C_{6}'
        sage: G2 = G.direct_product(G, maps=False)
        sage: LatexExpr(G2.structure_description(latex=True))   # optional - database_gap
        C_{6} \times C_{6}

    This method is mainly intended for small groups or groups with few
    normal subgroups. Even then there are some surprises::

        sage: D3 = DihedralGroup(3)
        sage: D3.structure_description()    # optional - database_gap
        'S3'

    We use the Sage notation for the degree of dihedral groups::

        sage: D4 = DihedralGroup(4)
        sage: D4.structure_description()    # optional - database_gap
        'D4'

    Works for finitely presented groups (:trac:`17573`)::

        sage: F.<x, y> = FreeGroup()
        sage: G=F / [x^2*y^-1, x^3*y^2, x*y*x^-1*y^-1]
        sage: G.structure_description()     # optional - database_gap
        'C7'

    And matrix groups (:trac:`17573`)::

        sage: groups.matrix.GL(4,2).structure_description() # optional - database_gap
        'A8'
    """
    import re
    def correct_dihedral_degree(match):
        return "%sD%d" % (match.group(1), int(match.group(2))/2)

    try:
        description = str(G._gap_().StructureDescription())
    except RuntimeError:
        from sage.misc.feature_test import SmallGroupsLibrary
        SmallGroupsLibrary().require()
        # the Small Groups Library is installed. The command failed for a different reason
        raise

    description = re.sub(r"(\A|\W)D(\d+)", correct_dihedral_degree, description)
    if not latex:
        return description
    description = description.replace("x", r"\times").replace(":", r"\rtimes")
    description = re.sub(r"([A-Za-z]+)([0-9]+)", r"\g<1>_{\g<2>}", description)
    description = re.sub(r"O([+-])", r"O^{\g<1>}", description)

    return description

