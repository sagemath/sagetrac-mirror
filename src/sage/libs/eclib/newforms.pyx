"""
Modular symbols using eclib newforms
"""

#*****************************************************************************
#       Copyright (C) 2008 Tom Boothby <boothby@u.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "cysignals/signals.pxi"

from ..eclib cimport *
from sage.libs.gmp.mpq cimport mpq_numref
from sage.libs.ntl.convert cimport mpz_to_ZZ
from sage.rings.rational_field import QQ
from sage.rings.rational cimport Rational
from sage.modular.all import Cusp


cdef class ECModularSymbol:
    """
    Modular symbol associated with an elliptic curve,  using John Cremona's newforms class.

    EXAMPLES::

        sage: from sage.libs.eclib.newforms import ECModularSymbol
        sage: E = EllipticCurve('11a')
        sage: M = ECModularSymbol(E,1); M
        Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

    By default, symbols are based at the cusp $\infty$, i.e. we evaluate $\{\infty,r\}$::

        sage: [M(1/i) for i in range(1,11)]
        [-2/5, 8/5, 3/5, -7/5, -12/5, -12/5, -7/5, 3/5, 8/5, -2/5]

    We can also switch the base point to the cusp $0$::

        sage: [M(1/i, base_at_infinity=False) for i in range(1,11)]
        [0, 2, 1, -1, -2, -2, -1, 1, 2, 0]

    For the minus symbols this makes no difference since
    $\{0,\infty\}$ is in the plus space::

        sage: M = ECModularSymbol(E,-1); M
        Modular symbol with sign -1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: [M(1/i) for i in range(1,11)]
        [0, 0, 1, 1, 0, 0, -1, -1, 0, 0]
        sage: [M(1/i, base_at_infinity=False) for i in range(1,11)]
        [0, 0, 1, 1, 0, 0, -1, -1, 0, 0]

    If the ECModularSymbol is created with sign 0 then we can ask for
    either + or - symbols, or both, though it is more work to create the full
    modular symbol space::

        sage: E = EllipticCurve('11a1')
        sage: M = ECModularSymbol(E,0); M
        Modular symbol with sign 0 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: [M(1/i) for i in range(2,11)]
        [[8/5, 0],
         [3/5, 1],
         [-7/5, 1],
         [-12/5, 0],
         [-12/5, 0],
         [-7/5, -1],
         [3/5, -1],
         [8/5, 0],
         [-2/5, 0]]

    The curve is automatically converted to its minimal model::

        sage: E = EllipticCurve([0,0,0,0,1/4])
        sage: ECModularSymbol(E)
        Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 over Rational Field
    """
    def __init__(self, E, sign=1):
        """
        Construct the modular symbol object from an elliptic curve.

        INPUT:

        - ``E``- an elliptic curve defined over Q
        - ``sign`` (int) -- 0, +1 or -1.  If +1 or -1, only modular
           symbols of this sign are availiable.  If 0, modular symbols
           of both signs are available but the construction is more expensive.


        EXAMPLES::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E, +1)
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E, -1)
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E, 0)

        TESTS:

        This one is from :trac:`8042`::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('858k2')
            sage: ECModularSymbol(E)
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + x*y = x^3 + 16353089*x - 335543012233 over Rational Field

        We allow a-invariants which are larger than 64 bits
        (:trac:`16977`)::

            sage: E = EllipticCurve([-25194941007454971, -1539281792450963687794218])  # non-minimal model of 21758k3
            sage: ECModularSymbol(E)  # long time
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + x*y = x^3 - 19440540900814*x - 32992152521343165084 over Rational Field

        The curve needs to be defined over the rationals::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: ECModularSymbol(E)
            Traceback (most recent call last):
            ...
            TypeError: the elliptic curve must have coefficients in the rational numbers
        """
        cdef ZZ_c a1, a2, a3, a4, a6, N
        cdef Curve *C
        cdef Curvedata *CD
        cdef CurveRed *CR
        cdef int n

        if E.base_field() is not QQ:
            raise TypeError("the elliptic curve must have coefficients in the rational numbers")

        E = E.minimal_model()
        self._E = E

        # The a invariants are rational numbers with denominator 1
        mpz_to_ZZ(&a1, mpq_numref((<Rational>(E.a1())).value))
        mpz_to_ZZ(&a2, mpq_numref((<Rational>(E.a2())).value))
        mpz_to_ZZ(&a3, mpq_numref((<Rational>(E.a3())).value))
        mpz_to_ZZ(&a4, mpq_numref((<Rational>(E.a4())).value))
        mpz_to_ZZ(&a6, mpq_numref((<Rational>(E.a6())).value))

        sig_on()
        C = new Curve(a1,a2,a3,a4,a6)
        CD = new Curvedata(C[0],0)
        CR = new CurveRed(CD[0])
        N = getconductor(CR[0])
        n = I2int(N)
        self.n = n
        self.sign = sign

        self.nfs = new newforms(n,0)
        self.nfs.createfromcurve(sign,CR[0])
        sig_off()

    def __repr__(self):
        """
        TESTS::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E, -1); M
            Modular symbol with sign -1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E, 0); M
            Modular symbol with sign 0 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return "Modular symbol with sign %s over Rational Field attached to %s"%(self.sign, self._E)

    def real_period(self):
        cdef double x
        #self.nfs.get_real_period(0, x, 0)
        return x

    def __call__(self, r, sign=None, base_at_infinity=True):
        """
        Computes the value of self on {0,r} or {oo,r} for rational r.

        INPUT:

        - ``r`` (rational) - a rational number
        - ``sign`` (int) - either +1, -1 or 0; must agree with the sign
           of the space unless that sign is 0 in which case either +1 or -1
           is allowed, or 0 in which case both the + and - symbols are returned.
           Default: self.sign, or +1 when self.sign=0.
        - ``base_at_infinity`` (bool) - if True, evaluates
          {oo,r}. otherwise (default) evaluates {0,r}.

        EXAMPLES::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,11)]
            [-2/5, 8/5, 3/5, -7/5, -12/5, -12/5, -7/5, 3/5, 8/5, -2/5]
            sage: [M(1/i, base_at_infinity=False) for i in range(1,11)]
            [0, 2, 1, -1, -2, -2, -1, 1, 2, 0]
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,11)]
            [0, 0, 0, 0, 1, 0, 1, 1, 0, 0]
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,11)]
            [0, 0, 0, 0, 2, 0, 1, 0, -1, 0]

        When the class is created with sign 0 we can ask for +1 or -1 symbols or (by default) both::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a1')
            sage: M = ECModularSymbol(E,0)
            sage: [M(1/i, 1) for i in range(2,11)]
            [8/5, 3/5, -7/5, -12/5, -12/5, -7/5, 3/5, 8/5, -2/5]
            sage: [M(1/i, -1) for i in range(2,11)]
            [0, 1, 1, 0, 0, -1, -1, 0, 0]
            sage: [M(1/i) for i in range(2,11)]
            [[8/5, 0],
             [3/5, 1],
             [-7/5, 1],
             [-12/5, 0],
             [-12/5, 0],
             [-7/5, -1],
             [3/5, -1],
             [8/5, 0],
             [-2/5, 0]]


        TESTS (see :trac:`11211`)::

            sage: from sage.libs.eclib.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: M(oo)
            0
            sage: M(oo, base_at_infinity=False)
            2/5
            sage: M(7/5)
            13/5
            sage: M("garbage")
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'garbage' to a cusp
            sage: M(7/5)
            13/5
        """
        cdef rational _r
        cdef rational _sp, _sm
        cdef pair[rational,rational] _spm
        cdef long n, d

        r = Cusp(r)
        d = r.denominator()
        n = r.numerator()
        if d != 0:
            n = n % d
        sig_on()
        _r = rational(n,d)
        if sign==None:
           sign = self.sign
        if sign==+1:
            if self.sign == -1:
                raise ValueError("impossible to evaluate a minus symbol on a plus space")
            else:
                _sp = self.nfs.plus_modular_symbol(_r, 0, int(base_at_infinity))
        elif sign==-1:
            if self.sign == +1:
                raise ValueError("impossible to evaluate a plus symbol on a minus space")
            else:
                _sm = self.nfs.minus_modular_symbol(_r, 0, int(base_at_infinity))
        elif sign==0:
            if self.sign != 0:
                raise ValueError("impossible to evaluate both symbols on a plus or minus space")
            else:
                _spm = self.nfs.full_modular_symbol(_r, 0, int(base_at_infinity))
                _sp = _spm.first
                _sm = _spm.second

        else:
            raise ValueError("sign must be +1, -1 or 0")
        sig_off()

        if sign==+1:
            return Rational((rational_num(_sp), rational_den(_sp)))
        elif sign==-1:
            return Rational((rational_num(_sm), rational_den(_sm)))
        else:
            return [Rational((rational_num(_sp), rational_den(_sp))),
                    Rational((rational_num(_sm), rational_den(_sm)))]
