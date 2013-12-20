r"""
Univariate polynomials over number fields

AUTHOR:

- Luis Felipe Tabera Alonso (2013-03): initial version, gcd computation.

EXAMPLES:

Define a polynomial over an absolute number field and make basic operations::

    sage: N.<a> = NumberField(x^2-2)
    sage: K.<x> = N[]
    sage: f = x - a
    sage: g = x^3 - 2*a + 1
    sage: f*(x + a)
    x^2 - 2
    sage: f + g
    x^3 + x - 3*a + 1
    sage: g // f
    x^2 + a*x + 2
    sage: g % f
    1
    sage: factor(x^3 - 2*a*x^2 - 2*x + 4*a)
    (x - 2*a) * (x - a) * (x + a)
    sage: gcd(f, x - a)
    x - a

Polynomials are aware of embeddings of the underlying field::

    sage: x = var('x')
    sage: Q7 = Qp(7)
    sage: r1 = Q7(3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8 +\
             6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16 + 7^17 +\
             4*7^18 + 6*7^19)
    sage: N.<b> = NumberField(x^2-2, embedding = r1)
    sage: K.<t> = N[]
    sage: f = t^3-2*t+1
    sage: f(r1)
    1 + O(7^20)

We can also construct polynomials over relative number fields::

    sage: N.<i, s2> = QQ[I, sqrt(2)]
    sage: K.<x> = N[]
    sage: f = x - s2
    sage: g = x^3 - 2*i*x^2 + s2*x
    sage: f*(x + s2)
    x^2 - 2
    sage: f + g
    x^3 - 2*I*x^2 + (sqrt2 + 1)*x - sqrt2
    sage: g // f
    x^2 + (-2*I + sqrt2)*x - 2*sqrt2*I + sqrt2 + 2
    sage: g % f
    -4*I + 2*sqrt2 + 2
    sage: factor(i*x^4 - 2*i*x^2 + 9*i)
    (I) * (x - I + sqrt2) * (x + I - sqrt2) * (x - I - sqrt2) * (x + I + sqrt2)
    sage: gcd(f, x-i)
    1
"""

#*****************************************************************************
#       Copyright (C) 2013 Luis Felipe Tabera Alonso <taberalf@unican.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from polynomial_element_generic import Polynomial_generic_dense_field
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.libs.ntl.ntl_ZZX import ntl_ZZX
from sage.libs.ntl.ntl_ZZ_pEX import ntl_ZZ_pEX
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext
from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pE import ntl_ZZ_pE
from sage.libs.ntl.ntl_ZZ_pEContext import ntl_ZZ_pEContext
from sage.structure.element import coerce_binop


class Polynomial_absolute_number_field_dense(Polynomial_generic_dense_field):
    """
    Class of dense univariate polynomials over an absolute number field.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new polynomial of the polynomial ring ``parent``.

        INPUT:

             - ``parent`` -- the underlying Polynomial ring.

             - ``x`` -- (default: None) An object representing the polynomial.
               e.g. a list of coefficients. See
               :meth:`sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field.__init__` for more details.

             - ``check`` -- boolean (default: True) if True make sure that the
               coefficients of the polynomial are the ring of coefficients.

             - ``is_gen`` -- boolean (defaul: False) A boolean, True if `x` is
               the disinguished generator of the polynomial ring.

             - ``construct`` -- (default: False) A boolean, unused.

        EXAMPLES::

            sage: f = QQ[I][x].random_element()
            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_number_field.Polynomial_absolute_number_field_dense'>
            sage: a = QQ[I][x](x)
            sage: a.is_gen()
            True
        """
        Polynomial_generic_dense_field.__init__(self, parent, x, check, is_gen, construct)

    @coerce_binop
    def gcd(self, other, algorithm='modular'):
        """
        Compute de monic gcd of two univariate polynomials.

        This method should not be called directly. However, it provides
        for convenience different algorithms to compute the gcd.

        INPUT:

        - ``other`` -- a polynomial with the same parent as self.
        - ``algorithm``

          - ``pari`` - use pari routines.
          - ``modular`` - modular algorithm using NTL (default).

        OUTPUT:

        - The monic gcd of ``self`` and ``other``

        EXAMPLES::

            sage: N.<a> = NumberField(x^3-1/2, 'a')
            sage: R.<r> = N['r']
            sage: f = (5/4*a^2 - 2*a + 4)*r^2 + (5*a^2 - 81/5*a - 17/2)*r + 4/5*a^2 + 24*a + 6
            sage: g = (5/4*a^2 - 2*a + 4)*r^2 + (-11*a^2 + 79/5*a - 7/2)*r - 4/5*a^2 - 24*a - 6
            sage: gcd(f, g**2)
            r - 60808/96625*a^2 - 69936/96625*a - 149212/96625
            sage: R = QQ[I]['x']
            sage: f = R.random_element(2)
            sage: g = f + 1
            sage: h = R.random_element(2).monic()
            sage: f *=h
            sage: g *=h
            sage: gcd(f, g) - h
            0
            sage: f.gcd(g, algorithm='pari') - h
            0
            sage: f.gcd(g) - h
            0
            sage: f.gcd(g, algorithm='modular') - h
            0

    TESTS:

        Test for degree one extensions::

            sage: x = var('x')
            sage: N = NumberField(x-3, 'a')
            sage: a = N.gen()
            sage: R = N[x]
            sage: f = R.random_element()
            sage: g1 = R.random_element()
            sage: g2 = g1*R.random_element() + 1
            sage: g1 *= f
            sage: g2 *= f
            sage: d = gcd(g1, g2)
            sage: f.monic() - d
            0
            sage: d.parent() is R
            True

        Test for coercion with other rings and force weird variables to test
        pari behavior::

            sage: r = var('r')
            sage: N = NumberField(r^2 - 2, 'r')
            sage: a = N.gen()
            sage: R = N['r']
            sage: r = R.gen()
            sage: f = N.random_element(4)*r + 1
            sage: g = ZZ['r']([1, 2, 3, 4, 5, 6, 7]); g
            7*r^6 + 6*r^5 + 5*r^4 + 4*r^3 + 3*r^2 + 2*r + 1
            sage: gcd(f, g) == gcd(g, f)
            True
            sage: h = f.gcd(g, algorithm='pari'); h
            1
            sage: h.parent()
            Univariate Polynomial Ring in r over Number Field in r with defining polynomial r^2 - 2
            sage: f.gcd(g, algorithm='modular')
            1

        Test keyword arguments in presence of coercion::

            sage: g.gcd(f, algorithm='pari')
            1
            sage: g.gcd(f, algorithm='modular')
            1
            sage: g.gcd(g, algorithm='pari')
            Traceback (most recent call last):
            ...
            TypeError: gcd() takes no keyword arguments
            sage: gcd([a*r+2, r^2-2])
            r + r

        ALGORITHM:

        - pari: Uses pari internal gcd routines.
        - modular (default): a combination of Langemyr-McCallum [LM89] and
          Encarnacion [E95]. Use Langemyr-McCallum algorithm but, from time to
          time, try a rational reconstruction of the gcd. Inspired from [HH09].

        REFERENCES:

        - [LM89] Langemyr, McCallum. The computation of polynomial greatest
          common divisors over an  algebraic number field. J. Symbolic Comput.
          8(5) (1989), 429--448.
        - [E95] Encarnacion. Computing gcds of polynomials over algebraic number
          fields. J. Symbolic Comput. ation, 20(3) (1995), 299-313.
        - [HH09] Hemmer, Hulse. Generic implementation of a modular gcd over
          Algebraic Extension Fields. EuroCG'09 : 25th European Workshop on
          Computational Geometry.

        """
        if self.is_zero():
            if other.is_zero():
                return self
            else:
                return other.monic()
        elif other.is_zero():
            return self.monic()
        elif self.degree() == 0:
            return self.parent().one_element()
        elif other.degree() == 0:
            return other.parent().one_element()

        #If the extension is of degree one, use the gcd from QQ[x]

        if self.base_ring().degree().is_one():
            R = self.parent()
            x = self.variable_name()
            a = QQ[x](self)
            b = QQ[x](other)
            g = a.gcd(b)
            return R(g)

        #Using pari to make the computations
        if algorithm == 'pari':
            h1 = self._pari_with_name('x')
            h2 = other._pari_with_name('x')
            g = h1.gcd(h2)
            return (self.parent()(g)).monic()
        if algorithm != 'modular':
            raise ValueError("unknown algorithm %s"%(algorithm))
        h1 = self.numerator()
        h2 = other.numerator()
        # Most cases have a small content
        #content_1 = gcd(sum([coeff._coefficients() for coeff in h1.list()]),[])
        #content_2 = gcd(sum([coeff._coefficients() for coeff in h2.list()]),[])
        #h1 = (~content_1) * h1
        #h2 = (~content_2) * h2
        N = h1.base_ring()
        R = h1.parent()
        Npol = N.polynomial().numerator()
        pol = ntl_ZZX(Npol.list())
        # leading_coefficient != 1 is not currently supported by Sage right now
        # (02-2013) but the code should work if the polynomial is not monic.
        if N.polynomial().denominator() == 1 and N.polynomial().leading_coefficient() == 1:
            # Use the denominator bound given by Langemyr, McCallum.
            # Do not assume that the leading coefficient of the gcd is an integer.
            # This is generally faster than the general bound.
            Bound = Npol.discriminant()  # D in Encarnacion paper.
            Bound = Bound * h1.leading_coefficient() * h2.leading_coefficient()
            D = ntl_ZZX(Bound.list())
            #Experimental bound IMPROVE
            #p = ZZ(3+min(2**255, (max(map(abs,Bound.list())).n()**(0.4)).floor())).next_prime(proof=False)
            p = ZZ(3+min(2**255, max(map(abs,Bound.list()))))
        else:
            # Use the denominator bound given by Encarnacion.
            Bound = Npol.discriminant()
            f= Npol.resultant(ZZ['x'](h1.leading_coefficient().polynomial()))
            g= Npol.resultant(ZZ['x'](h2.leading_coefficient().polynomial()))
            D = ntl_ZZX([Bound * f.gcd(g)])
            p = ZZ(3+min(2**255, max(map(abs,Bound.list()))))
        h1d = h1.degree()
        h2d = h2.degree()
        cd  = N.degree()
        # Save each polynomial as a list of lists for faster coercion to ntl_ZZ_pEX.
        h1ntl = [ i._coefficients() for i in h1.list()]
        h2ntl = [ i._coefficients() for i in h2.list()]
        # ss is a tuple containing: degree of the gcd, modular_gcd, modulus.
        ss = (h1d + 1,)
        # Whenever steps == nsteps, try a rational reconstruction of the gcd.
        steps = 0
        nsteps = 0
        while True:
            # We do not really need prime as long as the gcd success.
            p = p.next_prime(proof=False)
            # Recreate modular context.
            pol_p = ntl_ZZ_pX(pol, p)
            if pol_p.degree() == cd:
                c   = ntl_ZZ_pEContext(pol_p)
                h1c = ntl_ZZ_pEX(h1ntl, c)
                h2c = ntl_ZZ_pEX(h2ntl, c)
                Dc  = ntl_ZZ_pEX([ntl_ZZ_pE(D, c)])
                if h1c.degree() == h1d and h2c.degree() == h2d and Dc !=0:
                    # Compute residual gcd.
                    try:
                        gcd_pEX = h1c.gcd(h2c) * Dc
                    except (RuntimeError, ArithmeticError):
                        #RuntimeError if there is no gcd.
                        #ArithmeticError is the prime divides Dc.
                        gcd_pEX = ntl_ZZ_pEX([1],c).left_shift(h1d+3)
                    #if ss[0] < gcd_pEX.degree() discard this case.
                    if gcd_pEX.degree() == 0:
                        #The polynomials are relatively prime.
                        return R.one_element()
                    elif ss[0] > gcd_pEX.degree():
                        #All previous primes where bad primes, we start over again.
                        steps = 0
                        nsteps = 0
                        ss = gcd_pEX.degree(), gcd_pEX, p
                        cached = [[i.lift_centered() for i in j.get_as_ZZ_pX_doctest().list()] for j in gcd_pEX.list()]
                    elif ss[0] == gcd_pEX.degree():
                        #Success, apply chinese remainder to compute the
                        #residual of the gcd on a larger modulus.
                        steps +=1
                        g, c1, c2 = p.xgcd(ss[2])
                        if g<>1:
                            raise ValueError
                        m = ntl_ZZ_pContext(p*ss[2])
                        c = ntl_ZZ_pEContext(ntl_ZZ_pX(pol, m))
                        gcd_pEX = gcd_pEX.convert_to_pE(c)
                        gcd_pEX_previous = ss[1].convert_to_pE(c)
                        gcd_pEX = gcd_pEX * ntl_ZZ_pEX([c2*ss[2]],c) +\
                                  gcd_pEX_previous * ntl_ZZ_pEX([c1*p],c)
                        ss = ss[0], gcd_pEX, p * ss[2]
                        #Check if gcd_pEX has changed during the computation.
                        tentative_cached = [[i.lift_centered() for i in j.get_as_ZZ_pX_doctest().list()] for j in gcd_pEX.list()]
                        if cached != tentative_cached:
                            cached = tentative_cached
                        else:
                            # Check if we already have a gcd.
                            G = ss[1].lift_to_poly_ZZ(R).monic()
                            if (h1 % G).is_zero() and (h2 % G).is_zero():
                                return G.monic()
                        if steps >= nsteps:
                            # Try a rational reconstruction.
                            steps = 0
                            nsteps += nsteps//2 + 1
                            try:
                                # check if we already have a monic gcd.
                                gcd_pEX = gcd_pEX * ntl_ZZ_pEX([~gcd_pEX.leading_coefficient()])
                                G = gcd_pEX.lift_to_poly_QQ(R)
                                if (h1 % G).is_zero() and (h2 % G).is_zero():
                                    return G.monic()
                            except (ValueError, RuntimeError):
                                # Rational reconstruction failed.
                                # either lift_to_poly_QQ failed or
                                # ~gcd_pEX.leading_coefficient() does not exists
                                pass

class Polynomial_relative_number_field_dense(Polynomial_generic_dense_field):
    """
    Class of dense univariate polynomials over a relative number field.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new polynomial of the polynomial ring ``parent``.

        INPUT:

             - ``parent`` -- the underlying Polynomial ring.

             - ``x`` -- (default: None) An object representing the polynomial.
               e.g. a list of coefficients. See
               :meth:`sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field.__init__` for more details.

             - ``check`` -- boolean (default: True) if True make sure that the
               coefficients of the polynomial are the ring of coefficients.

             - ``is_gen`` -- boolean (defaul: False) A boolean, True if ``x`` is
               the disinguished generator of the polynomial ring.

             - ``construct`` -- (default: False) A boolean, unused.

        EXAMPLES::

            sage: f = NumberField([x^2-2, x^2-3], 'a')[x].random_element()
            sage: type(f)
            <class 'sage.rings.polynomial.polynomial_number_field.Polynomial_relative_number_field_dense'>
        """
        Polynomial_generic_dense_field.__init__(self, parent, x, check, is_gen, construct)

    @coerce_binop
    def gcd(self, other, algorithm='modular'):
        """
        Compute the monic gcd of two polynomials.

        Currently, the method checks corner cases in which one of the
        polynomials is zero or a constant. Then, computes an absolute extension
        and perform there the computations.

        INPUT:

        - ``other`` -- a polynomial with the same parent as self.
        - ``algorithm``

          - ``pari`` - use pari routines.
          - ``modular`` - modular algorithm using NTL (default).

        OUTPUT:

        - The monic gcd of ``self`` and ``other``

        See :meth:`Polynomial_absolute_number_field_dense.gcd` for more details.

        EXAMPLES::

            sage: N = QQ[sqrt(2), sqrt(3)]
            sage: s2, s3 = N.gens()
            sage: x = polygen(N)
            sage: f = x^4 - 5*x^2 +6
            sage: g = x^3 + (-2*s2 + s3)*x^2 + (-2*s3*s2 + 2)*x + 2*s3
            sage: gcd(f, g)
            x^2 + (-sqrt2 + sqrt3)*x - sqrt3*sqrt2
            sage: f.gcd(g, algorithm='pari')
            x^2 + (-sqrt2 + sqrt3)*x - sqrt3*sqrt2
            sage: f.gcd(g, algorithm='modular')
            x^2 + (-sqrt2 + sqrt3)*x - sqrt3*sqrt2


        TESTS::

            sage: x = var('x')
            sage: R = NumberField([x^2-2, x^2-3], 'a')[x]
            sage: f = R.random_element()
            sage: g1 = R.random_element()
            sage: g2 = R.random_element()*g1+1
            sage: g1 *= f
            sage: g2 *= f
            sage: f.monic() - g1.gcd(g2)
            0

        Test for degree one extensions::

            sage: R = NumberField([x-2,x+1,x-3],'a')[x]
            sage: f = R.random_element(2)
            sage: g1 = R.random_element(2)
            sage: g2 = R.random_element(2)*g1+1
            sage: g1 *= f
            sage: g2 *= f
            sage: d = gcd(g1, g2)
            sage: d - f.monic()
            0
            sage: d.parent() is R
            True

        Test for hardcoded variables::

            sage: R = N['sqrt2sqrt3']
            sage: x = R.gen()
            sage: f = x^2 - 2
            sage: g1 = x^2 - s3
            sage: g2 = x - s2
            sage: gcd(f, g1)
            1
            sage: gcd(f, g2)
            sqrt2sqrt3 - sqrt2
        """
        if self.is_zero():
            if other.is_zero():
                return self
            else:
                return other.monic()
        elif other.is_zero():
            return self.monic()
        elif self.degree() == 0:
            return self.parent().one_element()
        elif other.degree() == 0:
            return other.parent().one_element()

        L = self.parent()
        x = L.gen()
        N = self.base_ring()
        c = ''.join(map(str,N.variable_names()))
        M = N.absolute_field(c)
        M_to_N, N_to_M = M.structure()
        R = M[x]
        first = R(([N_to_M(foo) for foo in self.list()]))
        second = R(([N_to_M(foo) for foo in other.list()]))
        result = first.gcd(second, algorithm=algorithm)
        result = L(([M_to_N(foo) for foo in result.list()]))
        #the result is already monic
        return result
