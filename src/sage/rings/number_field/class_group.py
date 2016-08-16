# -*- coding: utf-8 -*-
r"""
Class Groups of Number Fields

An element of a class group is stored as a pair consisting of both an explicit
ideal in that ideal class, and a list of exponents giving that ideal class in
terms of the generators of the parent class group. These can be accessed with
the ``ideal()`` and ``exponents()`` methods respectively.

EXAMPLES::

    sage: K.<a> = NumberField(x^2 + 23)
    sage: I = K.class_group().gen(); I
    Fractional ideal class (2, 1/2*a - 1/2)
    sage: I.ideal()
    Fractional ideal (2, 1/2*a - 1/2)
    sage: I.exponents()
    (1,)

    sage: I.ideal() * I.ideal()
    Fractional ideal (4, 1/2*a + 3/2)
    sage: (I.ideal() * I.ideal()).reduce_equiv()
    Fractional ideal (2, 1/2*a + 1/2)
    sage: J = I * I; J    # class group multiplication is automatically reduced
    Fractional ideal class (2, 1/2*a + 1/2)
    sage: J.ideal()
    Fractional ideal (2, 1/2*a + 1/2)
    sage: J.exponents()
    (2,)

    sage: I * I.ideal()   # ideal classes coerce to their representative ideal
    Fractional ideal (4, 1/2*a + 3/2)

    sage: O = K.OK(); O
    Maximal Order in Number Field in a with defining polynomial x^2 + 23
    sage: O*(2, 1/2*a + 1/2)
    Fractional ideal (2, 1/2*a + 1/2)
    sage: (O*(2, 1/2*a + 1/2)).is_principal()
    False
    sage: (O*(2, 1/2*a + 1/2))^3
    Fractional ideal (1/2*a - 3/2)
"""

from sage.structure.sage_object import SageObject
from sage.groups.abelian_gps.values import AbelianGroupWithValues_class, AbelianGroupWithValuesElement
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.structure.sequence import Sequence
from sage.structure.element import MonoidElement
from sage.groups.old import Group
from sage.arith.all import LCM
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.all import ZZ
from sage.libs.pari.all import pari

def _integer_n_tuple_L1_iterator(n):
    if n == 1:
        i = 1
        while True:
            yield i
            yield -i
            i += 1
    else:
        from sage.combinat.partition import Partitions
        from sage.combinat.subset import Subsets
        N = 1
        sign_options_dict = {}
        Subsets_of_n = []
        while True:
            #print "***"
            #print N
            for k in range(1, n + 1):
                Ps = OrderedPartitions(N, k)
                for P in Ps:
                    #print "--"
                    #print P
                    try:
                        Ss = Subsets_of_n[k - 1]
                    except IndexError:
                        Ss = Subsets(range(n), k)
                        Subsets_of_n.append(Ss)
                    for S in Ss:
                        #print "="
                        #print S
                        i = [0] * n
                        for j in range(k):
                            i[S[j]] = P[j]
                        yield i
                        try:
                            sign_options = sign_options_dict[S]
                        except KeyError:
                            sign_options = Subsets(S)[1:]
                            sign_options_dict[S] = sign_options
                        for signs in sign_options:
                            ii = copy(i)
                            for index in signs:
                                ii[index] = - ii[index]
                            yield ii
            N += 1

class Modulus(SageObject):
    def __init__(self, finite, infinite=None, check=True):
        r"""
        Create a modulus of a number field.
        
        INPUT:
        
        - ``finite`` -- a non-zero fractional ideal in a number field.
        - ``infinite`` -- a list of indices corresponding to real places of the number field.
        - ``check`` (default: True) -- If ``True``, run a few checks on the input
        """
        self._finite = finite
        if infinite is None:
            infinite = []
        else:
            infinite.sort()
        self._infinite = tuple(ZZ(i) for i in infinite)
        K = self._finite.number_field()
        self._number_field = K
        if check:
            #insert various checks here
            if self._finite == 0:
                raise ValueError("Finite component of a modulus must be non-zero.")
            sgn = K.signature()[0]
            for i in self._infinite:
                if i < 0 or i >= sgn:
                    raise ValueError("Infinite component of a modulus must be a list non-negative integers less than the number of real places of K")
        return
    
    def _repr_(self):
        r"""
        EXAMPLES::
        
            sage: K.<a> = NumberField(x^2-5)
            sage: m = Modulus(K.ideal(31), [0,1]); m
            (Fractional ideal (31)) * infinity_0 * infinity_1
        """
        if len(self._infinite) == 0:
            return str(self._finite)
        str_inf = ''
        for i in self._infinite:
            str_inf += ' * infinity_%s'%(i)
        return '(' + str(self._finite) + ')' + str_inf
    
    def __eq__(self, other):
        return self._number_field == other._number_field and self._finite == other._finite and self._infinite == other._infinite
    
    def __mul__(self, other):
        r"""
        Multiply two moduli.
        
        This multiplies the two finite parts and performs an exclusive or on the real places.
        
        TESTS:
        
        ::
        
            sage: K = NumberField(x^3 - 2, 'a')
            sage: m1 = Modulus(K.ideal(2), [0])
            sage: m2 = Modulus(K.ideal(3), [0])
            sage: m1 * m2
            Fractional ideal (6)
        
        A higher degree totally real field.
        
        ::
        
            sage: K = NumberField(x^5 - x^4 - 4*x^3 + 3*x^2 + 3*x - 1, 'a')
            sage: m1 = Modulus(K.ideal(5), [2, 3])
            sage: m2 = Modulus(K.ideal(25), [0, 1, 3, 4])
            sage: m1 * m2
            (Fractional ideal (125)) * infinity_0 * infinity_1 * infinity_2 * infinity_4
            sage: _ == m2 * m1
            True
        """
        inf = list(set(self.infinite_part()).symmetric_difference(other.infinite_part()))
        return Modulus(self.finite_part() * other.finite_part(), inf, check=False)
    
    def divides(self, other):
        if not set(self.infinite_part()).issubset(other.infinite_part()):
            return False
        return self.finite_part().divides(other.finite_part())
    
    def number_field(self):
        return self._number_field
    
    def finite_part(self):
        return self._finite
    
    def infinite_part(self):
        return self._infinite
    
    def finite_factors(self):
        try:
            return self._finite_factors
        except AttributeError:
            self._finite_factors = m.finite_part().factor()
            return self._finite_factors
    
    def number_is_one_mod_star(self, a):
        K = self.number_field()
        am1 = K(a - 1)
        for P, e in self.finite_factors():
            if am1.valuation(P) < e:
                return False
        inf_places = K.places()
        for i in self.infinite_part():
            if inf_places[i](am1) <= 0:
                return False
        return True
    
    def fix_signs(self, a):
        r"""
        Given ``a`` in ``self.number_field()``, find `b` congruent to ``a`` `mod^\ast` ``self.finite_part()``
        such that `b` is positive at the infinite places dividing ``self``.
        """
        if self.is_finite() or a == 0:
            return a
        places = self.number_field().places()
        positive = []
        negative = []
        for i in self.infinite_part():
            if places[i](a) > 0:
                positive.append(i)
            else:
                negative.append(i)
        if len(negative) == 0:
            return a
        t = self.get_one_mod_star_finite_with_fixed_signs(positive, negative)
        return t * a
    
    def get_one_mod_star_finite_with_fixed_signs(self, positive, negative):
        if len(negative) == 0:
            return self.number_field().one()
        #positive = tuple(positive)
        negative = tuple(negative)
        try:
            return self._one_mod_star[negative]
        except AttributeError:
            self._one_mod_star = {}
        except KeyError:
            pass
        try:
            beta_is, Ainv = self._beta_is_Ainv
        except AttributeError:
            beta_is, Ainv = self._find_beta_is_Ainv()
        d = len(self.infinite_part())
        v = [0] * d
        for i in negative:
            v[i] = 1
        v = (GF(2)^d)(v)
        w = Ainv * v
        t = self.number_field().one()
        for i in range(d):
            if w[i] != 0:
                t *= (1 + beta_is[i])
        self._one_mod_star[negative] = t
        return t
    
    def _find_beta_is_Ainv(self):
        r"""
        Step 2 of Algorithm 4.2.20 of Cohen's Advanced...
        """
        gammas = self.finite_part().basis()
        k = len(self.infinite_part())
        beta_is = []
        Acols = []
        V = GF(2)^k
        it = _integer_n_tuple_L1_iterator(k)
        while len(beta_is) < k:
            e = it.next()
            beta = sum([e[i] * gammas[i] for i in range(k)])
            sbeta = V(self._signs(beta))
            Acols_new = Acols + [sbeta]
            #beta_is_new = beta_is + [sbeta]
            A = column_matrix(GF(2), Acols_new)
            if A.rank() == len(Acols_new):
                Acols = Acols_new
                beta_is.append(beta)
        self._beta_is_Ainv = (beta_is, ~A)
        return self._beta_is_Ainv
    
    def _signs(self, b):
        if b == 0:
            raise ValueError("Non-zero input required.")
        sigmas = K.places()
        return [(1 - sigmas[i](b).sign()).divide_knowing_divisible_by(2) for i in self.infinite_part()]
    
    def is_finite(self):
        return len(self._infinite) == 0
    
    def is_infinite(self):
        return self._finite.is_one()
    
    def _pari_(self):
        inf_mod = [0] * self._number_field.signature()[0]
        for i in self._infinite:
            inf_mod[i] = 1
        return pari([self._finite, inf_mod])
    
    def __hash__(self):
        return hash((self._finite, self._infinite))

class FractionalIdealClass(AbelianGroupWithValuesElement):
    r"""
    A fractional ideal class in a number field.

    EXAMPLES::

        sage: G = NumberField(x^2 + 23,'a').class_group(); G
        Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
        sage: I = G.0; I
        Fractional ideal class (2, 1/2*a - 1/2)
        sage: I.ideal()
        Fractional ideal (2, 1/2*a - 1/2)

        EXAMPLES::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c = C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.gens()
            (2, 1/2*w - 1/2)
    """
    def __init__(self, parent, element, ideal=None):
        """
        Returns the ideal class of this fractional ideal.

        EXAMPLE::

            sage: K.<a> = NumberField(x^2 + 23,'a'); G = K.class_group()
            sage: G(K.ideal(13, a + 4))
            Fractional ideal class (13, 1/2*a + 17/2)
        """
        if element is None:
            element = parent._ideal_log(ideal)
        AbelianGroupWithValuesElement.__init__(self, parent, element, ideal)

    def _repr_(self):
        r"""
        Return string representation of this fractional ideal class.

         EXAMPLE::

            sage: K.<a> = NumberField(x^2 + 23,'a'); G = K.class_group()
            sage: G(K.ideal(13, a + 4))._repr_()
            'Fractional ideal class (13, 1/2*a + 17/2)'
            sage: G(K.ideal(59, a+6))._repr_()
            'Trivial principal fractional ideal class'
        """
        if self.is_principal():
            return 'Trivial principal fractional ideal class'
        return 'Fractional ideal class %s'%self._value._repr_short()

    def _mul_(self, other):
        r"""
        Multiplication of two (S-)ideal classes.

        EXAMPLE::

            sage: G = NumberField(x^2 + 23,'a').class_group(); G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
            sage: I = G.0; I
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: I*I # indirect doctest
            Fractional ideal class (2, 1/2*a + 1/2)
            sage: I*I*I # indirect doctest
            Trivial principal fractional ideal class

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: G = K.ideal(3,a+1)
            sage: CS(G)*CS(G)
            Trivial S-ideal class
        """
        m = AbelianGroupElement._mul_(self, other)
        m._value = (self.ideal() * other.ideal()).reduce_equiv()
        return m

    def _div_(self, other):
        r"""
        Division of two ideal classes.

        EXAMPLE::

            sage: G = NumberField(x^2 + 23,'a').class_group(); G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
            sage: I = G.0; I
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: I*I # indirect doctest
            Fractional ideal class (2, 1/2*a + 1/2)
            sage: I*I*I # indirect doctest
            Trivial principal fractional ideal class
        """
        m = AbelianGroupElement._div_(self, other)
        m._value = (self.ideal() / other.ideal()).reduce_equiv()
        return m

    def __pow__(self, n):
        r"""
        Raise this element to the power n.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - 3*x + 8)
            sage: C=K.class_group()
            sage: c = C(2, a)
            sage: c^2
            Fractional ideal class (2, a^2 + 2*a - 1)
            sage: c^3
            Trivial principal fractional ideal class
            sage: c^1000
            Fractional ideal class (2, a)
            sage: (c^2)^2
            Fractional ideal class (2, a)
        """
        # We use MonoidElement's __pow__ routine, since that does
        # repeated squaring, and hence the ideal gets reduced as
        # we go along; actually computing self._value ** n would
        # be disastrous.
        n = n % self.order()
        return MonoidElement.__pow__(self, n)

    def inverse(self):
        r"""
        Return the multiplicative inverse of this ideal class.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - 3*x + 8); G = K.class_group()
            sage: G(2, a).inverse()
            Fractional ideal class (2, a^2 + 2*a - 1)
            sage: ~G(2, a)
            Fractional ideal class (2, a^2 + 2*a - 1)
        """
        m = AbelianGroupElement.inverse(self)
        m._value = (~self.ideal()).reduce_equiv()
        return m

    __invert__ = inverse

    def is_principal(self):
        r"""
        Returns True iff this ideal class is the trivial (principal) class

        EXAMPLES::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c=C(P2a)
            sage: c.is_principal()
            False
            sage: (c^2).is_principal()
            False
            sage: (c^3).is_principal()
            True
        """
        return self.is_one()

    def reduce(self):
        r"""
        Return representative for this ideal class that has been
        reduced using PARI's idealred.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 20072); G = k.class_group(); G
            Class group of order 76 with structure C38 x C2
            of Number Field in a with defining polynomial x^2 + 20072
            sage: I = (G.0)^11; I
            Fractional ideal class (41, 1/2*a + 5)
            sage: J = G(I.ideal()^5); J
            Fractional ideal class (115856201, 1/2*a + 40407883)
            sage: J.reduce()
            Fractional ideal class (57, 1/2*a + 44)
            sage: J == I^5
            True
        """
        return self.parent()(self.ideal().reduce_equiv())

    def ideal(self):
        r"""
        Return a representative ideal in this ideal class.

        EXAMPLE::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c=C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.ideal()
            Fractional ideal (2, 1/2*w - 1/2)
        """
        return self.value()

    def representative_prime(self, norm_bound=1000):
        r"""
        Return a prime ideal in this ideal class.

        INPUT:

        ``norm_bound`` (positive integer) -- upper bound on the norm of primes tested.

        EXAMPLE::

           sage: K.<a> = NumberField(x^2+31)
           sage: K.class_number()
           3
           sage: Cl = K.class_group()
           sage: [c.representative_prime() for c in Cl]
           [Fractional ideal (3),
           Fractional ideal (2, 1/2*a + 1/2),
           Fractional ideal (2, 1/2*a - 1/2)]

           sage: K.<a> = NumberField(x^2+223)
           sage: K.class_number()
           7
           sage: Cl = K.class_group()
           sage: [c.representative_prime() for c in Cl]
           [Fractional ideal (3),
           Fractional ideal (2, 1/2*a + 1/2),
           Fractional ideal (17, 1/2*a + 7/2),
           Fractional ideal (7, 1/2*a - 1/2),
           Fractional ideal (7, 1/2*a + 1/2),
           Fractional ideal (17, 1/2*a + 27/2),
           Fractional ideal (2, 1/2*a - 1/2)]
        """
        if self.value().is_prime():
            return self.value()
        c = self.reduce()
        if c.value().is_prime():
            return c.value()
        # otherwise we just search:
        Cl = self.parent()
        K = Cl.number_field()
        from sage.rings.all import RR
        for P in K.primes_of_bounded_norm_iter(RR(norm_bound)):
            if Cl(P)==c:
                return P
        raise RuntimeError("No prime of norm less than %s found in class %s" % (norm_bound, c))


    def gens(self):
        r"""
        Return generators for a representative ideal in this
        (S-)ideal class.

        EXAMPLES::

            sage: K.<w>=QuadraticField(-23)
            sage: OK = K.ring_of_integers()
            sage: C = OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c = C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.gens()
            (2, 1/2*w - 1/2)
       """
        return self.ideal().gens()

class RayClassGroupElement(FractionalIdealClass):
    def __init__(self, parent, element, ideal=None):
        if element is None:
            if not parent.modulus().finite_part().is_coprime(ideal):
                raise ValueError("Ideal is not coprime to the modulus.")
            element = parent._ideal_log(ideal)
        #Should treat the else case for coprime-ness as well since the code can coerce from different moduli
        FractionalIdealClass.__init__(self, parent, element, ideal)
    
    def _repr_(self):
        if self.is_one():
            return 'Trivial ray class modulo ' + str(self.parent().modulus())
        return 'Ray class of ' + self._value._repr_short() + ' modulo ' + str(self.parent().modulus())
    
    #Should be able to get rid of the operations if make the reduce function a method of the parent
    def _mul_(self, other):
        m = AbelianGroupElement._mul_(self, other)
        m._value = (self.ideal() * other.ideal()).reduce_equiv(self.parent().modulus())
        return m
    
    def _div_(self, other):
        m = AbelianGroupElement._div_(self, other)
        m._value = (self.ideal() / other.ideal()).reduce_equiv(self.parent().modulus())
        return m
    
    def inverse(self):
        m = AbelianGroupElement.inverse(self)
        m._value = (~self.ideal()).reduce_equiv(self.parent().modulus())
        return m
    
    __invert__ = inverse
    
    def reduce(self):
        return self.ideal().reduce_equiv(self.parent().modulus())

class SFractionalIdealClass(FractionalIdealClass):
    r"""
    An S-fractional ideal class in a number field for a tuple of primes S.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I)
            Trivial S-ideal class
            sage: CS(J)
            Trivial S-ideal class
            sage: CS(G)
            Fractional S-ideal class (3, a + 1)

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I).ideal()
            Fractional ideal (2, a)
            sage: CS(J).ideal()
            Fractional ideal (7, a)
            sage: CS(G).ideal()
            Fractional ideal (3, a + 1)


        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: G = K.ideal(3,a+1)
            sage: CS(G).inverse()
            Fractional S-ideal class (3, a + 2)

    TESTS::

        sage: K.<a> = QuadraticField(-14)
        sage: I = K.ideal(2,a)
        sage: S = (I,)
        sage: CS = K.S_class_group(S)
        sage: J = K.ideal(7,a)
        sage: G = K.ideal(3,a+1)
        sage: CS(I).order()
        1
        sage: CS(J).order()
        1
        sage: CS(G).order()
        2
    """

    def _repr_(self):
        r"""
        Returns a string representation of the S-ideal class of this fractional ideal.

        EXAMPLE::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: J = K.ideal(3, a + 2)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: CS(J)
            Fractional S-ideal class (3, a + 2)
            sage: CS(J^2)
            Trivial S-ideal class
        """
        if self.is_trivial():
            return 'Trivial S-ideal class'
        return 'Fractional S-ideal class %s' % self._value._repr_short()



class ClassGroup(AbelianGroupWithValues_class):
    r"""
    The class group of a number field.

    EXAMPLES::

        sage: K.<a> = NumberField(x^2 + 23)
        sage: G = K.class_group(); G
        Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
        sage: G.category()
        Category of finite commutative groups

    Note the distinction between abstract generators, their ideal, and
    exponents::

        sage: C = NumberField(x^2 + 120071, 'a').class_group(); C
        Class group of order 500 with structure C250 x C2
        of Number Field in a with defining polynomial x^2 + 120071
        sage: c = C.gen(0)
        sage: c  # random
        Fractional ideal class (5, 1/2*a + 3/2)
        sage: c.ideal()  # random
        Fractional ideal (5, 1/2*a + 3/2)
        sage: c.ideal() is c.value()   # alias
        True
        sage: c.exponents()
        (1, 0)
    """
    Element = FractionalIdealClass

    def __init__(self, gens_orders, names, number_field, gens, proof=True):
        r"""
        Create a class group.

        TESTS::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: G = K.class_group()
            sage: TestSuite(G).run()
        """
        AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
                                              values_group=number_field.ideal_monoid())
        self._proof_flag = proof
        self._number_field = number_field

    def _element_constructor_(self, *args, **kwds):
        r"""
        Create an element of this class group from the given data. This may be:
        an ideal class in this number field; an ideal class in a subfield; or
        anything from which an ideal in this number field can be constructed.

        EXAMPLES::

            sage: K.<b> = NumberField(x^2 + 389)
            sage: C = K.class_group()
            sage: C(K.ideal(b)) # indirect doctest
            Trivial principal fractional ideal class
            sage: C(K.ideal(59049, b + 35312)) # indirect doctest
            Fractional ideal class (59049, b + 35312)
            sage: C((59049, b + 35312)) # indirect doctest
            Fractional ideal class (59049, b + 35312)
            sage: C(59049, b + 35312) # indirect doctest
            Fractional ideal class (59049, b + 35312)

            sage: K.<a> = QuadraticField(-23)
            sage: L.<b> = K.extension(x^2 - 2)
            sage: CK = K.class_group()
            sage: CL = L.class_group()
            sage: [CL(I).exponents() for I in CK]
            [(0,), (2,), (4,)]
        """
        if isinstance(args[0], FractionalIdealClass):
            return self.element_class(self, None, self._number_field.ideal(args[0].ideal()))
        else:
            I = self._number_field.ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            return self.element_class(self, None, I)

    def _ideal_log(self, ideal):
        """
        Compute the exponents from the ``ideal``.

        Used by the element constructor if necessary.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23,'a')
            sage: G = K.class_group()
            sage: g = G.an_element()
            sage: G._ideal_log(g.ideal())
            (1,)
            sage: g.exponents()
            (1,)
        """
        return tuple(ZZ(order) for order in ideal.ideal_class_log(proof=self._proof_flag))

    def gens_ideals(self):
        r"""
        Return generating ideals for the (S-)class group.

        This is an alias for :meth:`gens_values`.

        OUTPUT:

        A tuple of ideals, one for each abstract Abelian group generator.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 23)
            sage: K.class_group().gens_ideals()   # random gens (platform dependent)
            (Fractional ideal (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4),)

            sage: C = NumberField(x^2 + x + 23899, 'a').class_group(); C
            Class group of order 68 with structure C34 x C2 of Number Field
            in a with defining polynomial x^2 + x + 23899
            sage: C.gens()
            (Fractional ideal class (7, a + 5), Fractional ideal class (5, a + 3))
            sage: C.gens_ideals()
            (Fractional ideal (7, a + 5), Fractional ideal (5, a + 3))
        """
        return self.gens_values()

    def __iter__(self):
        r"""
        Return an iterator of all ideal classes in this class group.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 23)
            sage: G = K.class_group()
            sage: G
            Class group of order 3 with structure C3 of Number Field
            in a with defining polynomial x^4 + 23
            sage: list(G)
            [Trivial principal fractional ideal class,
             Fractional ideal class (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4),
             Fractional ideal class (2, 1/2*a^2 + 1/2)]
            sage: G.list()
            (Trivial principal fractional ideal class,
             Fractional ideal class (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4),
             Fractional ideal class (2, 1/2*a^2 + 1/2))

        TESTS::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: G = K.class_group()
            sage: G
            Class group of order 1 of Number Field in a with defining polynomial x^2 + 1
            sage: list(G)
            [Trivial principal fractional ideal class]
            sage: G.list()
            (Trivial principal fractional ideal class,)
        """
        from sage.misc.mrange import mrange
        orders = self.gens_orders()
        T = mrange(orders)
        g = self.gens()
        for t in T:
            I = self(1)
            for i, j in enumerate(t):
                I *= g[i]**j
            yield I
        if not T:
            yield self(1)

    def _repr_(self):
        r"""
        Return string representation of self.

        EXAMPLES::

            sage: C = NumberField(x^2 + 23, 'a').class_group()
            sage: C._repr_()
            'Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23'
        """
        s = 'Class group of order %s '%self.order()
        if self.order() > 1:
            s += 'with structure %s '%self._group_notation(self.gens_orders())
        s += 'of %s'%self.number_field()
        return s

    def number_field(self):
        r"""
        Return the number field that this (S-)class group is attached to.

        EXAMPLES::

            sage: C = NumberField(x^2 + 23, 'w').class_group(); C
            Class group of order 3 with structure C3 of Number Field in w with defining polynomial x^2 + 23
            sage: C.number_field()
            Number Field in w with defining polynomial x^2 + 23

            sage: K.<a> = QuadraticField(-14)
            sage: CS = K.S_class_group(K.primes_above(2))
            sage: CS.number_field()
            Number Field in a with defining polynomial x^2 + 14
        """
        return self._number_field

class RayClassGroup(ClassGroup):
    Element = RayClassGroupElement
    
    def __init__(self, gens_orders, names, modulus, gens, proof=True, bnr=None):
        #AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
        #                                      values_group=modulus.number_field().ideal_monoid())
        ClassGroup.__init__(self, gens_orders, names, modulus.number_field(), gens, proof=proof)
        #self._proof_flag = proof
        self._modulus = modulus
        #self._number_field = modulus.number_field()
        self._bnr = bnr
    
    #Can remove this I think if make ClassGroup's version say isinstance(args[0], Element):
    def _element_constructor_(self, *args, **kwds):
        if isinstance(args[0], RayClassGroupElement):
            return self.element_class(self, None, self._number_field.ideal(args[0].ideal()))
        else:
            I = self._number_field.ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            return self.element_class(self, None, I)
    
    def _repr_(self):
        return "Ray class group of " + str(self._number_field) + " of modulus " + str(self._modulus)
    
    def _ideal_log(self, ideal):
        return tuple(ZZ(c) for c in self._bnr.bnrisprincipal(ideal, flag = 0))
    
    def modulus(self):
        return self._modulus
    
    def pari_bnr(self):
        return self._bnr
    
    def pari_gens(self):
        return self._bnr[4][2]

class SClassGroup(ClassGroup):
    r"""
    The S-class group of a number field.

    EXAMPLES::

        sage: K.<a> = QuadraticField(-14)
        sage: S = K.primes_above(2)
        sage: K.S_class_group(S).gens()   # random gens (platform dependent)
        (Fractional S-ideal class (3, a + 2),)

        sage: K.<a> = QuadraticField(-974)
        sage: CS = K.S_class_group(K.primes_above(2)); CS
        S-class group of order 18 with structure C6 x C3
        of Number Field in a with defining polynomial x^2 + 974
        sage: CS.gen(0) # random
        Fractional S-ideal class (3, a + 2)
        sage: CS.gen(1) # random
        Fractional S-ideal class (31, a + 24)
    """
    Element = SFractionalIdealClass

    def __init__(self, gens_orders, names, number_field, gens, S, proof=True):
        r"""
        Create an S-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: K.S_class_group(S)
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14
            sage: K.<a> = QuadraticField(-105)
            sage: K.S_class_group([K.ideal(13, a + 8)])
            S-class group of order 4 with structure C2 x C2 of Number Field in a with defining polynomial x^2 + 105
        """
        AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
                                              values_group=number_field.ideal_monoid())
        self._proof_flag = proof
        self._number_field = number_field
        self._S = S

    def S(self):
        r"""
        Return the set (or rather tuple) of primes used to define this class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S);CS
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14
            sage: T = tuple([])
            sage: CT = K.S_class_group(T);CT
            S-class group of order 4 with structure C4 of Number Field in a with defining polynomial x^2 + 14
            sage: CS.S()
            (Fractional ideal (2, a),)
            sage: CT.S()
            ()
        """
        return self._S

    def _ideal_log(self, ideal):
        """
        Compute the exponents from the ``ideal``.

        Used by the element constructor if necessary.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: s = CS.an_element()
            sage: CS._ideal_log(s.ideal())
            (1,)
            sage: s.exponents()
            (1,)
        """
        return tuple(ZZ(order) for order in ideal.S_ideal_class_log(self.S()))

    def _element_constructor_(self, *args, **kwds):
        r"""
        Create an element of this class group from the given data.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I)
            Trivial S-ideal class
            sage: CS(J)
            Trivial S-ideal class
            sage: CS(G)
            Fractional S-ideal class (3, a + 1)
        """
        if isinstance(args[0], FractionalIdealClass):
            return self.element_class(self, None, args[0].ideal())
        else:
            I = self.number_field().ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            return self.element_class(self, None, I)

    def _repr_(self):
        r"""
        Return string representation of this S-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: CS = K.S_class_group(K.primes_above(2))
            sage: CS._repr_()
            'S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14'
        """
        s = 'S-class group of order %s ' % self.order()
        if self.order() > 1:
            s += 'with structure %s ' % self._group_notation(self.gens_orders())
        s += 'of %s' % self.number_field()
        return s
