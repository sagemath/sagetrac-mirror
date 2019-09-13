r"""
Base class for multivariate polynomial rings
"""
from __future__ import absolute_import, print_function

import sage.misc.latex
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod

from sage.structure.element cimport parent
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport rich_to_bool, richcmp
from cpython.object cimport Py_NE

import sage.categories as categories
from sage.categories.morphism import IdentityMorphism
from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()

from sage.arith.all import binomial

from sage.combinat.integer_vector import IntegerVectors

from sage.rings.integer_ring import ZZ

from .polydict import PolyDict
from . import (multi_polynomial_ideal, polynomial_ring, multi_polynomial_element)
from .term_order import TermOrder
from .polynomial_ring_constructor import (PolynomialRing, polynomial_default_category)


def is_MPolynomialRing(x):
    return isinstance(x, MPolynomialRing_base)


cdef class MPolynomialRing_base(sage.rings.ring.CommutativeRing):
    def __init__(self, base_ring, n, names, order):
        """
        Create a polynomial ring in several variables over a commutative ring.

        EXAMPLES::

            sage: R.<x,y> = ZZ['x,y']; R
            Multivariate Polynomial Ring in x, y over Integer Ring
            sage: class CR(CommutativeRing):
            ....:     def __init__(self):
            ....:         CommutativeRing.__init__(self,self)
            ....:     def __call__(self,x):
            ....:         return None
            sage: cr = CR()
            sage: cr.is_commutative()
            True
            sage: cr['x,y']
            Multivariate Polynomial Ring in x, y over <__main__.CR_with_category object at ...>

        TESTS:

        Check that containment works correctly (:trac:`10355`)::

            sage: A1.<a> = PolynomialRing(QQ)
            sage: A2.<a,b> = PolynomialRing(QQ)
            sage: 3 in A2
            True
            sage: A1(a) in A2
            True

        Check that :trac:`26958` is fixed::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: class Foo(MPolynomialRing_libsingular):
            ....:     pass
            sage: Foo(QQ, 2, ['x','y'], 'degrevlex')
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        if base_ring not in _CommutativeRings:
            raise TypeError("The base ring %s is not a commutative ring" % base_ring)

        n = int(n)
        if n < 0:
            raise ValueError("Multivariate Polynomial Rings must " + \
                  "have more than 0 variables.")
        order = TermOrder(order, n)
        self.__ngens = n
        self.__term_order = order
        self._has_singular = False #cannot convert to Singular by default
        self._magma_cache = {}
        # Ring.__init__ already does assign the names.
        # It would be a mistake to call ParentWithGens.__init__
        # as well, assigning the names twice.
        #ParentWithGens.__init__(self, base_ring, names)
        if base_ring.is_zero():
            category = categories.rings.Rings().Finite()
        else:
            category = polynomial_default_category(base_ring.category(), n)

        # Avoid calling __init_extra__ of Algebras(...).parent_class
        self._no_generic_basering_coercion = True

        sage.rings.ring.Ring.__init__(self, base_ring, names, category=category)

    def is_integral_domain(self, proof = True):
        """
        EXAMPLES::

            sage: ZZ['x,y'].is_integral_domain()
            True
            sage: Integers(8)['x,y'].is_integral_domain()
            False
        """
        return self.base_ring().is_integral_domain(proof)

    def is_noetherian(self):
        """
        EXAMPLES::

            sage: ZZ['x,y'].is_noetherian()
            True
            sage: Integers(8)['x,y'].is_noetherian()
            True
        """
        return self.base_ring().is_noetherian()

    @cached_method
    def flattening_morphism(self):
        r"""
        Return the flattening morphism of this polynomial ring

        EXAMPLES::

            sage: QQ['a','b']['x','y'].flattening_morphism()
            Flattening morphism:
              From: Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in a, b over Rational Field
              To:   Multivariate Polynomial Ring in a, b, x, y over Rational Field

            sage: QQ['x,y'].flattening_morphism()
            Identity endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
        """
        base = self.base_ring()
        if is_MPolynomialRing(base) or polynomial_ring.is_PolynomialRing(base):
            from .flatten import FlatteningMorphism
            return FlatteningMorphism(self)
        else:
            return IdentityMorphism(self)

    def construction(self):
        """
        Returns a functor F and base ring R such that F(R) == self.

        EXAMPLES::

            sage: S = ZZ['x,y']
            sage: F, R = S.construction(); R
            Integer Ring
            sage: F
            MPoly[x,y]
            sage: F(R) == S
            True
            sage: F(R) == ZZ['x']['y']
            False

        """
        from .polynomial_ring_constructor import PolynomialRing
        from sage.categories.pushout import MultiPolynomialFunctor
        return MultiPolynomialFunctor(self.variable_names(), self.term_order()), self.base_ring()

    def irrelevant_ideal(self):
        """
        Return the irrelevant ideal of this multivariate polynomial ring.

        This is the ideal generated by all of the indeterminate
        generators of this ring.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: R.irrelevant_ideal()
            Ideal (x, y, z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.ideal(self.gens(), check=False)

    def completion(self, names, prec=20):
        """
        Return the completion of self with respect to the ideal
        generated by the variable(s) ``names``.

        INPUT:

        - ``names`` -- variable or list/tuple of variables (given either
          as elements of the polynomial ring or as strings)

        - ``prec`` -- default precision of resulting power series ring

        - ``extras`` -- deprecated and ignored

        EXAMPLES::

            sage: P.<x,y,z,w> = PolynomialRing(ZZ)
            sage: P.completion('w')
            Power Series Ring in w over Multivariate Polynomial Ring in x, y, z over Integer Ring
            sage: P.completion((w,x,y))
            Multivariate Power Series Ring in w, x, y over Univariate Polynomial Ring in z over Integer Ring
            sage: Q.<w,x,y,z> = P.completion(); Q
            Multivariate Power Series Ring in w, x, y, z over Integer Ring

            sage: H = PolynomialRing(PolynomialRing(ZZ,3,'z'),4,'f'); H
            Multivariate Polynomial Ring in f0, f1, f2, f3 over
            Multivariate Polynomial Ring in z0, z1, z2 over Integer Ring

            sage: H.completion(H.gens())
            Multivariate Power Series Ring in f0, f1, f2, f3 over
            Multivariate Polynomial Ring in z0, z1, z2 over Integer Ring

            sage: H.completion(H.gens()[2])
            Power Series Ring in f2 over
            Multivariate Polynomial Ring in f0, f1, f3 over
            Multivariate Polynomial Ring in z0, z1, z2 over Integer Ring

        TESTS::

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: P.completion([]) is P
            True
            sage: P.completion(SR.var('x'))
            Traceback (most recent call last):
            ...
            TypeError: x is not an element of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: P.completion(x + y)
            Traceback (most recent call last):
            ...
            ValueError: x + y is not a variable of Multivariate Polynomial Ring in x, y over Integer Ring
            sage: P.completion('q')
            Traceback (most recent call last):
            ...
            ValueError: q is not a variable of Multivariate Polynomial Ring in x, y over Integer Ring
        """
        if not isinstance(names, (list, tuple)):
            names = [names]  # Single variable
        elif not names:
            return self      # 0 variables => completion is self

        vars = []
        for v in names:
            # Convert variable names to str and check that they really
            # are variables of self
            if isinstance(v, str):
                pass
            elif parent(v) is self:
                v = str(v)
            else:
                raise TypeError(f"{v!r} is not an element of {self}")
            if v not in self.variable_names():
                raise ValueError(f"{v} is not a variable of {self}")
            vars.append(v)

        from sage.rings.power_series_ring import PowerSeriesRing
        new_base = self.remove_var(*vars)
        return PowerSeriesRing(new_base, names=vars, default_prec=prec)

    def remove_var(self, *var, order=None):
        """
        Remove a variable or sequence of variables from self.

        If ``order`` is not specified, then the subring inherits the
        term order of the original ring, if possible.

        EXAMPLES::

            sage: P.<x,y,z,w> = PolynomialRing(ZZ)
            sage: P.remove_var(z)
            Multivariate Polynomial Ring in x, y, w over Integer Ring
            sage: P.remove_var(z,x)
            Multivariate Polynomial Ring in y, w over Integer Ring
            sage: P.remove_var(y,z,x)
            Univariate Polynomial Ring in w over Integer Ring

        Removing all variables results in the base ring::

            sage: P.remove_var(y,z,x,w)
            Integer Ring

        If possible, the term order is kept::

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, order='deglex')
            sage: R.remove_var(y).term_order()
            Degree lexicographic term order

            sage: R.<x,y,z,w> = PolynomialRing(ZZ, order='lex')
            sage: R.remove_var(y).term_order()
            Lexicographic term order

        Be careful with block orders when removing variables::

            sage: R.<x,y,z,u,v> = PolynomialRing(ZZ, order='deglex(2),lex(3)')
            sage: R.remove_var(x,y,z)
            Traceback (most recent call last):
            ...
            ValueError: impossible to use the original term order (most likely because it was a block order). Please specify the term order for the subring
            sage: R.remove_var(x,y,z, order='degrevlex')
            Multivariate Polynomial Ring in u, v over Integer Ring

        """
        vars = list(self.variable_names())
        for v in var:
            vars.remove(str(v))
        if len(vars) == 0:
            return self.base_ring()
        from .polynomial_ring_constructor import PolynomialRing
        if order is None:
            try:
                return PolynomialRing(self.base_ring(), vars,  order=self.term_order())
            except ValueError:
                raise ValueError("impossible to use the original term order (most likely because it was a block order). Please specify the term order for the subring")
        else:
            return PolynomialRing(self.base_ring(), vars, order=order)

    def univariate_ring(self, x):
        """
        Return a univariate polynomial ring whose base ring comprises all but one variables of self.

        INPUT:

        - ``x`` -- a variable of self.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.univariate_ring(y)
            Univariate Polynomial Ring in y over Multivariate Polynomial Ring in x, z over Rational Field
        """
        return self.remove_var(x)[str(x)]

    cdef _coerce_c_impl(self, x):
        """
        Return the canonical coercion of x to this multivariate
        polynomial ring, if one is defined, or raise a TypeError.

        The rings that canonically coerce to this polynomial ring are:

        - this ring itself
        - polynomial rings in the same variables over any base ring that
          canonically coerces to the base ring of this ring
        - polynomial rings in a subset of the variables over any base ring that
          canonically coerces to the base ring of this ring
        - any ring that canonically coerces to the base ring of this polynomial
          ring.

        TESTS:

        This fairly complicated code (from Michel Vandenbergh) ends up
        implicitly calling ``_coerce_c_impl``::

            sage: z = polygen(QQ, 'z')
            sage: W.<s>=NumberField(z^2+1)
            sage: Q.<u,v,w> = W[]
            sage: W1 = FractionField (Q)
            sage: S.<x,y,z> = W1[]
            sage: u + x
            x + u
            sage: x + 1/u
            x + 1/u
        """
        try:
            P = x.parent()
            # polynomial rings in the same variable over the any base that coerces in:
            if is_MPolynomialRing(P):
                if P.variable_names() == self.variable_names():
                    if self.has_coerce_map_from(P.base_ring()):
                        return self(x)
                elif self.base_ring().has_coerce_map_from(P._mpoly_base_ring(self.variable_names())):
                    return self(x)

            elif polynomial_ring.is_PolynomialRing(P):
                if P.variable_name() in self.variable_names():
                    if self.has_coerce_map_from(P.base_ring()):
                        return self(x)

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this polynomial ring.
        return self._coerce_try(x, [self.base_ring()])

    def _extract_polydict(self, x):
        """
        Assuming other_vars is a subset of ``self.variable_names()``,
        convert the dict of ETuples with respect to other_vars to
        a dict with respect to ``self.variable_names()``.
        """
        # This is probably horribly inefficient
        from .polydict import ETuple
        other_vars = list(x.parent().variable_names())
        name_mapping = [(other_vars.index(var) if var in other_vars else -1) for var in self.variable_names()]
        K = self.base_ring()
        D = {}
        var_range = xrange(len(self.variable_names()))
        for ix, a in x.dict().iteritems():
            ix = ETuple([0 if name_mapping[t] == -1 else ix[name_mapping[t]] for t in var_range])
            D[ix] = K(a)
        return D

    def __richcmp__(left, right, int op):
        if left is right:
            return rich_to_bool(op, 0)

        if not isinstance(right, Parent) or not isinstance(left, Parent):
            # One is not a parent -- not equal and not ordered
            return op == Py_NE

        if not is_MPolynomialRing(right):
            return op == Py_NE

        lft = <MPolynomialRing_base>left
        other = <MPolynomialRing_base>right

        lx = (lft.base_ring(), lft.__ngens,
              lft.variable_names(),
              lft.__term_order)
        rx = (other.base_ring(), other.__ngens,
              other.variable_names(),
              other.__term_order)
        return richcmp(lx, rx, op)

    def _repr_(self):
        """
        Return string representation of this object.

        EXAMPLES::

            sage: PolynomialRing(QQ, names=[])
            Multivariate Polynomial Ring in no variables over Rational Field

            sage: PolynomialRing(QQ, names=['x', 'y'])
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        if self.ngens() == 0:
            generators_rep = "no variables"
        else:
            generators_rep = ", ".join(self.variable_names())
        return "Multivariate Polynomial Ring in %s over %s"%(generators_rep, self.base_ring())

    def repr_long(self):
        """
        Return structured string representation of self.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,order=TermOrder('degrevlex',1)+TermOrder('lex',2))
            sage: print(P.repr_long())
            Polynomial Ring
             Base Ring : Rational Field
                  Size : 3 Variables
              Block  0 : Ordering : degrevlex
                         Names    : x
              Block  1 : Ordering : lex
                         Names    : y, z
        """
        from sage.rings.polynomial.term_order import inv_singular_name_mapping
        n = self.ngens()
        k = self.base_ring()
        names = self.variable_names()
        T = self.term_order()
        _repr =  "Polynomial Ring\n"
        _repr += "  Base Ring : %s\n"%(k,)
        _repr += "       Size : %d Variables\n"%(n,)
        offset = 0
        i = 0
        for order in T.blocks():
            _repr += "   Block % 2d : Ordering : %s\n"%(i,inv_singular_name_mapping.get(order.singular_str(), order.singular_str()))
            _repr += "              Names    : %s\n"%(", ".join(names[offset:offset + len(order)]))
            offset += len(order)
            i+=1
        return _repr

    def _latex_(self):
        vars = ', '.join(self.latex_variable_names())
        return "%s[%s]"%(sage.misc.latex.latex(self.base_ring()), vars)

    def _ideal_class_(self, n=0):
        return multi_polynomial_ideal.MPolynomialIdeal

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: T.<t> = ZZ[]
            sage: K.<i> = NumberField(t^2 + 1)
            sage: R.<x,y> = K[]
            sage: Q5 = Qp(5); i5 = Q5(-1).sqrt()
            sage: R._is_valid_homomorphism_(Q5, [Q5.teichmuller(2), Q5(6).log()]) # no coercion
            False
            sage: R._is_valid_homomorphism_(Q5, [Q5.teichmuller(2), Q5(6).log()], base_map=K.hom([i5]))
            True
        """
        if base_map is None:
            # all that is needed is that elements of the base ring
            # of the polynomial ring canonically coerce into codomain.
            # Since poly rings are free, any image of the gen
            # determines a homomorphism
            return codomain.has_coerce_map_from(self._base)
        return True

    def _magma_init_(self, magma):
        """
        Used in converting this ring to the corresponding ring in Magma.

        EXAMPLES::

            sage: R.<a,b,c,d,e,f,g,h,i,j> = PolynomialRing(GF(127),10)
            sage: R._magma_init_(magma)                      # optional - magma
            'SageCreateWithNames(PolynomialRing(_sage_ref...,10,"grevlex"),["a","b","c","d","e","f","g","h","i","j"])'
            sage: R.<y,z,w> = PolynomialRing(QQ,3)
            sage: magma(R)                                   # optional - magma
            Polynomial ring of rank 3 over Rational Field
            Order: Graded Reverse Lexicographical
            Variables: y, z, w

        A complicated nested example::

            sage: R.<a,b,c> = PolynomialRing(GF(9,'a')); S.<T,W> = R[]; S
            Multivariate Polynomial Ring in T, W over Multivariate Polynomial Ring in a, b, c over Finite Field in a of size 3^2
            sage: magma(S)                                   # optional - magma
            Polynomial ring of rank 2 over Polynomial ring of rank 3 over GF(3^2)
            Order: Graded Reverse Lexicographical
            Variables: T, W


            sage: magma(PolynomialRing(GF(7),4, 'x'))        # optional - magma
            Polynomial ring of rank 4 over GF(7)
            Order: Graded Reverse Lexicographical
            Variables: x0, x1, x2, x3

            sage: magma(PolynomialRing(GF(49,'a'),10, 'x'))  # optional - magma
            Polynomial ring of rank 10 over GF(7^2)
            Order: Graded Reverse Lexicographical
            Variables: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9

            sage: magma(PolynomialRing(ZZ['a,b,c'],3, 'x'))  # optional - magma
            Polynomial ring of rank 3 over Polynomial ring of rank 3 over Integer Ring
            Order: Graded Reverse Lexicographical
            Variables: x0, x1, x2
        """
        B = magma(self.base_ring())
        Bref = B._ref()
        s = 'PolynomialRing(%s,%s,%s)'%(Bref, self.ngens(), self.term_order().magma_str())
        return magma._with_names(s, self.variable_names())

    def _gap_init_(self, gap=None):
        """
        Return a string that yields a representation of ``self`` in GAP.

        INPUT:

        ``gap`` -- (optional GAP instance) Interface to which the
                   string is addressed.

        NOTE:

        - If the optional argument ``gap`` is provided, the base ring
          of ``self`` will be represented as ``gap(self.base_ring()).name()``.
        - The result of applying the GAP interface to ``self`` is cached.

        EXAMPLES::

            sage: F = CyclotomicField(8)
            sage: P.<x,y> = F[]
            sage: gap(P)     # indirect doctest
            PolynomialRing( CF(8), ["x", "y"] )
            sage: libgap(P)
            <field in characteristic 0>[x,y]
        """
        L = ['"%s"'%t for t in self.variable_names()]
        if gap is not None:
            return 'PolynomialRing(%s,[%s])'%(gap(self.base_ring()).name(),','.join(L))
        return 'PolynomialRing(%s,[%s])'%(self.base_ring()._gap_init_(),','.join(L))

    def is_field(self, proof = True):
        """
        Test whether this multivariate polynomial ring is a field.

        A polynomial ring is a field when there are no variable and the base
        ring is a field.

        EXAMPLES::

            sage: PolynomialRing(QQ, 'x', 2).is_field()
            False
            sage: PolynomialRing(QQ, 'x', 0).is_field()
            True
            sage: PolynomialRing(ZZ, 'x', 0).is_field()
            False
            sage: PolynomialRing(Zmod(1), names=['x','y']).is_finite()
            True
        """
        if not self.ngens():
            return self.base_ring().is_field(proof)
        return False

    def term_order(self):
        return self.__term_order

    def characteristic(self):
        """
        Return the characteristic of this polynomial ring.

        EXAMPLES::

            sage: R = PolynomialRing(QQ, 'x', 3)
            sage: R.characteristic()
            0
            sage: R = PolynomialRing(GF(7),'x', 20)
            sage: R.characteristic()
            7
        """
        return self.base_ring().characteristic()

    def gen(self, n=0):
        if n < 0 or n >= self.__ngens:
            raise ValueError("Generator not defined.")
        return self._gens[int(n)]

    def variable_names_recursive(self, depth=sage.rings.infinity.infinity):
        r"""
        Returns the list of variable names of this and its base rings, as if
        it were a single multi-variate polynomial.

        EXAMPLES::

            sage: R = QQ['x,y']['z,w']
            sage: R.variable_names_recursive()
            ('x', 'y', 'z', 'w')
            sage: R.variable_names_recursive(3)
            ('y', 'z', 'w')

        """
        if depth <= 0:
            all = ()
        elif depth == 1:
            all = self.variable_names()
        else:
            my_vars = self.variable_names()
            try:
               all = self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
            except AttributeError:
                all = my_vars
        if len(all) > depth:
            all = all[-depth:]
        return all

    def _mpoly_base_ring(self, vars=None):
        """
        Returns the base ring if this is viewed as a polynomial ring over vars.
        See also MPolynomial._mpoly_dict_recursive.
        """
        if vars is None:
            vars = self.variable_names_recursive()
        vars = list(vars)
        my_vars = list(self.variable_names())
        if vars == list(my_vars):
            return self.base_ring()
        elif not my_vars[-1] in vars:
            return self
        elif not set(my_vars).issubset(set(vars)):
            while my_vars[-1] in vars:
                my_vars.pop()
            from .polynomial_ring_constructor import PolynomialRing
            return PolynomialRing(self.base_ring(), my_vars)
        else:
            try:
                return self.base_ring()._mpoly_base_ring(vars[:vars.index(my_vars[0])])
            except AttributeError:
                return self.base_ring()


    def krull_dimension(self):
        return self.base_ring().krull_dimension() + self.ngens()

    def ngens(self):
        return self.__ngens

    def _monomial_order_function(self):
        raise NotImplementedError

    def __reduce__(self):
        """
        """

        base_ring = self.base_ring()
        n = self.ngens()
        names = self.variable_names()
        order = self.term_order()

        return unpickle_MPolynomialRing_generic_v1,(base_ring, n, names, order)

    def _precomp_counts(self,n, d):
        """
        Given a number of variable n and a degree d return a tuple (C,t)
        such that C is a list of the cardinalities of the sets of
        monomials up to degree d (including) in n variables and t is the
        sum of these cardinalities.

        INPUT:

        - ``n`` -- number of variables
        - ``d`` -- degree

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: C,t = P._precomp_counts(10,2)
            sage: C[0]
            1
            sage: C[1]
            10
            sage: C[2]
            55
            sage: t
            66

        TESTS::

            sage: P.<x,y> = PolynomialRing(ZZ)
            sage: C,t = P._precomp_counts(1,2)
            sage: C[0]
            1
            sage: C[1]
            1
            sage: C[2]
            1
            sage: t
            3

        """
        C = [1]  #d = 0
        for dbar in xrange(1, d+1):
            C.append(binomial(n+dbar-1, dbar))
        t = sum(C)
        return C, t

    def _to_monomial(self,i, n, d):
        """
        Given an index i, a number of variables n and a degree d return
        the i-th monomial of degree d in n variables.

        INPUT:

        - ``i`` -- index: 0 <= i < binom(n+d-1,n-1)
        - ``n`` -- number of variables
        - ``d`` -- degree

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ)
            sage: P._to_monomial(0,10,2)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 2)
            sage: P._to_monomial(8,10,2)
            (0, 0, 0, 0, 0, 0, 1, 1, 0, 0)
            sage: P._to_monomial(54,10,2)
            (2, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        .. NOTE::

            We do not check if the provided index/rank is within the allowed
            range. If it is not an infinite loop will occur.
        """
        from sage.combinat import combination
        comb = combination.from_rank(i, n+d-1, n-1)
        if comb == []:
            return (d,)
        monomial = [ comb[0] ]
        res = []
        for j in range(n-2):
            res.append(comb[j+1]-comb[j]-1)
        monomial += res
        monomial.append( n+d-1-comb[-1]-1 )
        return tuple(monomial)

    def _random_monomial_upto_degree_class(self,n, degree, counts=None, total=None):
        """
        Choose a random exponent tuple for `n` variables with a random
        degree `d`, i.e. choose the degree uniformly at random first
        before choosing a random monomial.

        INPUT:

        - ``n`` -- number of variables
        - ``degree`` -- degree of monomials
        - ``counts`` -- ignored
        - ``total`` -- ignored

        EXAMPLES::

            sage: K.<x,y,z,w> = QQ[]
            sage: K._random_monomial_upto_degree_class(5, 7)
            (0, 0, 0, 3, 0)
            """
        # bug: doesn't handle n=1
        #Select random degree
        d = ZZ.random_element(0,degree+1)
        total = binomial(n+d-1, d)

        #Select random monomial of degree d
        random_index = ZZ.random_element(0, total-1)
        #Generate the corresponding monomial
        return self._to_monomial(random_index, n, d)

    def _random_monomial_upto_degree_uniform(self,n, degree, counts=None, total=None):
        """
        Choose a random exponent tuple for `n` variables with a random
        degree up to `d`, i.e. choose a random monomial uniformly random
        from all monomials up to degree `d`. This discriminates against
        smaller degrees because there are more monomials of bigger
        degrees.

        INPUT:

        - ``n`` -- number of variables
        - ``degree`` -- degree of monomials
        - ``counts`` -- ignored
        - ``total`` -- ignored

        EXAMPLES::

            sage: K.<x,y,z,w> = QQ[]
            sage: K._random_monomial_upto_degree_uniform(4, 3)
            (1, 0, 0, 1)
            """
        if counts is None or total is None:
            counts, total = self._precomp_counts(n, degree)

        #Select a random one
        random_index = ZZ.random_element(0, total-1)
        #Figure out which degree it corresponds to
        d = 0
        while random_index >= counts[d]:
            random_index -= counts[d]
            d += 1
        #Generate the corresponding monomial
        return self._to_monomial(random_index, n, d)

    def random_element(self, degree=2, terms=None, choose_degree=False,*args, **kwargs):
        """
        Return a random polynomial of at most degree `d` and at most `t`
        terms.

        First monomials are chosen uniformly random from the set of all
        possible monomials of degree up to `d` (inclusive). This means
        that it is more likely that a monomial of degree `d` appears than
        a monomial of degree `d-1` because the former class is bigger.

        Exactly `t` *distinct* monomials are chosen this way and each one gets
        a random coefficient (possibly zero) from the base ring assigned.

        The returned polynomial is the sum of this list of terms.

        INPUT:

        - ``degree`` -- maximal degree (likely to be reached) (default: 2)

        - ``terms`` -- number of terms requested (default: 5). If more
          terms are requested than exist, then this parameter is
          silently reduced to the maximum number of available terms.

        - ``choose_degree`` -- choose degrees of monomials randomly first
          rather than monomials uniformly random.

        - ``**kwargs`` -- passed to the random element generator of the base
          ring

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: P.random_element(2, 5)
            -6/5*x^2 + 2/3*z^2 - 1

            sage: P.random_element(2, 5, choose_degree=True)
            -1/4*x*y - x - 1/14*z - 1

        Stacked rings::

            sage: R = QQ['x,y']
            sage: S = R['t,u']
            sage: S.random_element(degree=2, terms=1)
            -1/2*x^2 - 1/4*x*y - 3*y^2 + 4*y
            sage: S.random_element(degree=2, terms=1)
            (-x^2 - 2*y^2 - 1/3*x + 2*y + 9)*u^2

        Default values apply if no degree and/or number of terms is
        provided::

            sage: random_matrix(QQ['x,y,z'], 2, 2)
            [357*x^2 + 1/4*y^2 + 2*y*z + 2*z^2 + 28*x      2*x*y + 3/2*y^2 + 2*y*z - 2*z^2 - z]
            [                       x*y - y*z + 2*z^2         -x^2 - 4/3*x*z + 2*z^2 - x + 4*y]

            sage: random_matrix(QQ['x,y,z'], 2, 2, terms=1, degree=2)
            [ 1/2*y -1/4*x]
            [   1/2  1/3*x]

            sage: P.random_element(0, 1)
            1

            sage: P.random_element(2, 0)
            0

            sage: R.<x> = PolynomialRing(Integers(3), 1)
            sage: R.random_element()
            2*x^2 + x

        To produce a dense polynomial, pick ``terms=Infinity``::

            sage: P.<x,y,z> = GF(127)[]
            sage: P.random_element(degree=2, terms=Infinity)
            -55*x^2 - 51*x*y + 5*y^2 + 55*x*z - 59*y*z + 20*z^2 + 19*x - 55*y - 28*z + 17
            sage: P.random_element(degree=3, terms=Infinity)
            -54*x^3 + 15*x^2*y - x*y^2 - 15*y^3 + 61*x^2*z - 12*x*y*z + 20*y^2*z - 61*x*z^2 - 5*y*z^2 + 62*z^3 + 15*x^2 - 47*x*y + 31*y^2 - 14*x*z + 29*y*z + 13*z^2 + 61*x - 40*y - 49*z + 30
            sage: P.random_element(degree=3, terms=Infinity, choose_degree=True)
            57*x^3 - 58*x^2*y + 21*x*y^2 + 36*y^3 + 7*x^2*z - 57*x*y*z + 8*y^2*z - 11*x*z^2 + 7*y*z^2 + 6*z^3 - 38*x^2 - 18*x*y - 52*y^2 + 27*x*z + 4*y*z - 51*z^2 - 63*x + 7*y + 48*z + 14

        The number of terms is silently reduced to the maximum
        available if more terms are requested::

            sage: P.<x,y,z> = GF(127)[]
            sage: P.random_element(degree=2, terms=1000)
            5*x^2 - 10*x*y + 10*y^2 - 44*x*z + 31*y*z + 19*z^2 - 42*x - 50*y - 49*z - 60

        """
        k = self.base_ring()
        n = self.ngens()

        counts, total = self._precomp_counts(n, degree)

        # Note that 'terms' could be None while 'total' is a
        # nonnegative integer, so the comparison 'terms > total' could
        # fail in Python 3.
        if terms and terms > total:
            terms = total

        if terms is None:
            if total >= 5:
                terms = 5
            else:
                terms = total

        if terms < 0:
            raise TypeError("Cannot compute polynomial with a negative number of terms.")
        elif terms == 0:
            return self._zero_element
        if degree == 0:
            return k.random_element(**kwargs)


        from sage.combinat.integer_vector import IntegerVectors

        #total is 0. Just return
        if total == 0:
            return self._zero_element

        elif terms < total/2:
            # we choose random monomials if t < total/2 because then we
            # expect the algorithm to be faster than generating all
            # monomials and picking a random index from the list. if t ==
            # total/2 we expect every second random monomial to be a
            # double such that our runtime is doubled in the worst case.
            M = set()
            if not choose_degree:
                while terms:
                    m = self._random_monomial_upto_degree_uniform(n, degree, counts, total)
                    if not m in M:
                        M.add(m)
                        terms -= 1
            else:
                while terms:
                    m = self._random_monomial_upto_degree_class(n, degree)
                    if not m in M:
                        M.add(m)
                        terms -= 1
        elif terms <= total:
            # generate a list of all monomials and choose among them
            if not choose_degree:
                M = sum([list(IntegerVectors(_d,n)) for _d in xrange(degree+1)],[])
                for mi in xrange(total - terms): # we throw away those we don't need
                    M.pop( ZZ.random_element(0,len(M)-1) )
                M = [tuple(m) for m in M]
            else:
                M = [list(IntegerVectors(_d,n)) for _d in xrange(degree+1)]
                Mbar = []
                for mi in xrange(terms):
                    d = ZZ.random_element(0,len(M)) #choose degree at random
                    m = ZZ.random_element(0,len(M[d])) # choose monomial at random
                    Mbar.append( M[d].pop(m) ) # remove and insert
                    if len(M[d]) == 0:
                        M.pop(d) # bookkeeping
                M = [tuple(m) for m in Mbar]

        C = [k.random_element(*args,**kwargs) for _ in range(len(M))]

        return self(dict(zip(M,C)))

    def change_ring(self, base_ring=None, names=None, order=None):
        """
        Return a new multivariate polynomial ring which isomorphic to
        self, but has a different ordering given by the parameter
        'order' or names given by the parameter 'names'.

        INPUT:

        - ``base_ring`` -- a base ring
        - ``names`` -- variable names
        - ``order`` -- a term order

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(GF(127),3,order='lex')
            sage: x > y^2
            True
            sage: Q.<x,y,z> = P.change_ring(order='degrevlex')
            sage: x > y^2
            False
        """
        if base_ring is None:
            base_ring = self.base_ring()
        if names is None:
            names = self.variable_names()
        if order is None:
            order = self.term_order()

        from .polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(base_ring, self.ngens(), names, order=order)

    def monomial(self,*exponents):
        """
        Return the monomial with given exponents.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(ZZ, 3)
            sage: R.monomial(1,1,1)
            x*y*z
            sage: e=(1,2,3)
            sage: R.monomial(*e)
            x*y^2*z^3
            sage: m = R.monomial(1,2,3)
            sage: R.monomial(*m.degrees()) == m
            True
        """
        return self({exponents:self.base_ring().one()})

    def _macaulay_resultant_getS(self,mon_deg_tuple,dlist):
        r"""
        In the Macaulay resultant algorithm the list of all monomials of the total degree is partitioned into sets `S_i`.
        This function returns the index `i` for the set `S_i` for the given monomial.

        INPUT:

        - ``mon_deg_tuple`` -- a list representing a monomial of a degree `d`
        - ``dlist`` -- a list of degrees `d_i` of the polynomials in question, where
        `d =  sum(dlist) - len(dlist) + 1`

        OUTPUT:

        - the index `i` such that the input monomial is in `S_i`

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 2)
            sage: R._macaulay_resultant_getS([1,1,0],[2,1,1]) # the monomial xy where the total degree = 2
            1

            sage: R._macaulay_resultant_getS([29,21,8],[10,20,30])
            0

            sage: R._macaulay_resultant_getS(list(range(9))+[10],list(range(1,11)))
            9
        """
        for i in xrange(len(dlist)):
            if mon_deg_tuple[i] - dlist[i] >= 0:
                return i

    def _macaulay_resultant_is_reduced(self,mon_degs,dlist):
        r"""
        Helper function for the Macaulay resultant algorithm.
        A monomial in the variables `x_0,...,x_n` is called reduced with respect to the list of degrees `d_0,...,d_n`
        if the degree of `x_i` in the monomial is `>= d_i` for exactly one `i`. This function checks this property for a monomial.

        INPUT:

        - ``mon_degs`` -- a monomial represented by a vector of degrees
        - ``dlist`` -- a list of degrees with respect to which we check reducedness

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R._macaulay_resultant_is_reduced([2,3,1],[2,3,3]) # the monomial x^2*y^3*z is not reduced w.r.t. degrees vector [2,3,3]
            False

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R._macaulay_resultant_is_reduced([1,3,2],[2,3,3]) # the monomial x*y^3*z^2 is not reduced w.r.t. degrees vector [2,3,3]
            True
        """
        diff = [mon_degs[i] - dlist[i] for i in xrange(0,len(dlist))]
        return len([1 for d in diff if d >= 0]) == 1

    def _macaulay_resultant_universal_polynomials(self, dlist):
        r"""
        Given a list of degrees, this function returns a list of ``len(dlist)`` polynomials with ``len(dlist)`` variables,
        with generic coefficients. This is useful for generating polynomials for tests,
        and for getting a universal macaulay resultant for the given degrees.

        INPUT:

        - ``dlist`` -- a list of degrees.

        OUTPUT:

        - a list of polynomials of the given degrees with general coefficients.
        - a polynomial ring over ``self`` generated by the coefficients of the output polynomials.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(ZZ, 2)
            sage: R._macaulay_resultant_universal_polynomials([1,1,2])
            ([u0*x0 + u1*x1 + u2*x2, u3*x0 + u4*x1 + u5*x2, u6*x0^2 + u7*x0*x1 + u9*x1^2 + u8*x0*x2 + u10*x1*x2 + u11*x2^2], Multivariate Polynomial Ring in x0, x1, x2 over Multivariate Polynomial Ring in u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11 over Integer Ring)
        """
        n  = len(dlist) - 1
        number_of_coeffs = sum([binomial(n+di,di) for di in dlist])
        U = PolynomialRing(ZZ,'u',number_of_coeffs)
        d = sum(dlist) - len(dlist) + 1
        flist = []
        R = PolynomialRing(U,'x',n+1)
        ulist  = U.gens()
        for d in dlist:
            xlist = R.gens()
            degs = IntegerVectors(d, n+1)
            mon_d = [prod([xlist[i]**(deg[i]) for i in xrange(0,len(deg))])
                     for deg in degs]

            f = sum([mon_d[i]*ulist[i] for i in xrange(0,len(mon_d))])
            flist.append (f)
            ulist = ulist[len(mon_d):]
        return flist, R

    def macaulay_resultant(self, *args, **kwds):
        r"""
        This is an implementation of the Macaulay Resultant. It computes
        the resultant of universal polynomials as well as polynomials
        with constant coefficients. This is a project done in
        sage days 55. It's based on the implementation in Maple by
        Manfred Minimair, which in turn is based on the references listed below:
        It calculates the Macaulay resultant for a list of polynomials,
        up to sign!

        REFERENCES:

        - [CLO2005]_

        - [Can1990]_

        - [Mac1916]_

        AUTHORS:

        - Hao Chen, Solomon Vishkautsan (7-2014)

        INPUT:

        - ``args`` -- a list of `n` homogeneous polynomials in `n` variables.
                  works when ``args[0]`` is the list of polynomials,
                  or ``args`` is itself the list of polynomials

        kwds:

        - ``sparse`` -- boolean (optional - default: ``False``)
                     if ``True`` function creates sparse matrices.

        OUTPUT:

        - the macaulay resultant, an element of the base ring of ``self``

        .. TODO::
            Working with sparse matrices should usually give faster results,
            but with the current implementation it actually works slower.
            There should be a way to improve performance with regards to this.

        EXAMPLES:

        The number of polynomials has to match the number of variables::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R.macaulay_resultant([y,x+z])
            Traceback (most recent call last):
            ...
            TypeError: number of polynomials(= 2) must equal number of variables (= 3)

        The polynomials need to be all homogeneous::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R.macaulay_resultant([y, x+z, z+x^3])
            Traceback (most recent call last):
            ...
            TypeError: resultant for non-homogeneous polynomials is not supported

        All polynomials must be in the same ring::

            sage: S.<x,y> = PolynomialRing(QQ, 2)
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: S.macaulay_resultant([y, z+x])
            Traceback (most recent call last):
            ...
            TypeError: not all inputs are polynomials in the calling ring

        The following example recreates Proposition 2.10 in Ch.3 in [CLO2005]::

            sage: K.<x,y> = PolynomialRing(ZZ, 2)
            sage: flist,R = K._macaulay_resultant_universal_polynomials([1,1,2])
            sage: R.macaulay_resultant(flist)
            u2^2*u4^2*u6 - 2*u1*u2*u4*u5*u6 + u1^2*u5^2*u6 - u2^2*u3*u4*u7 + u1*u2*u3*u5*u7 + u0*u2*u4*u5*u7 - u0*u1*u5^2*u7 + u1*u2*u3*u4*u8 - u0*u2*u4^2*u8 - u1^2*u3*u5*u8 + u0*u1*u4*u5*u8 + u2^2*u3^2*u9 - 2*u0*u2*u3*u5*u9 + u0^2*u5^2*u9 - u1*u2*u3^2*u10 + u0*u2*u3*u4*u10 + u0*u1*u3*u5*u10 - u0^2*u4*u5*u10 + u1^2*u3^2*u11 - 2*u0*u1*u3*u4*u11 + u0^2*u4^2*u11

        The following example degenerates into the determinant of a `3*3` matrix::

            sage: K.<x,y> = PolynomialRing(ZZ, 2)
            sage: flist,R = K._macaulay_resultant_universal_polynomials([1,1,1])
            sage: R.macaulay_resultant(flist)
            -u2*u4*u6 + u1*u5*u6 + u2*u3*u7 - u0*u5*u7 - u1*u3*u8 + u0*u4*u8

        The following example is by Patrick Ingram (:arxiv:`1310.4114`)::

            sage: U = PolynomialRing(ZZ,'y',2); y0,y1 = U.gens()
            sage: R = PolynomialRing(U,'x',3); x0,x1,x2 = R.gens()
            sage: f0 = y0*x2^2 - x0^2 + 2*x1*x2
            sage: f1 = y1*x2^2 - x1^2 + 2*x0*x2
            sage: f2 = x0*x1 - x2^2
            sage: flist = [f0,f1,f2]
            sage: R.macaulay_resultant([f0,f1,f2])
            y0^2*y1^2 - 4*y0^3 - 4*y1^3 + 18*y0*y1 - 27

        a simple example with constant rational coefficients::

            sage: R.<x,y,z,w> = PolynomialRing(QQ,4)
            sage: R.macaulay_resultant([w,z,y,x])
            1

        an example where the resultant vanishes::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R.macaulay_resultant([x+y,y^2,x])
            0

        an example of bad reduction at a prime `p = 5`::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R.macaulay_resultant([y,x^3+25*y^2*x,5*z])
            125

        The input can given as an unpacked list of polynomials::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R.macaulay_resultant(y,x^3+25*y^2*x,5*z)
            125

        an example when the coefficients live in a finite field::

            sage: F = FiniteField(11)
            sage: R.<x,y,z,w> = PolynomialRing(F,4)
            sage: R.macaulay_resultant([z,x^3,5*y,w])
            4

        example when the denominator in the algorithm vanishes(in this case
        the resultant is the constant term of the quotient of
        char polynomials of numerator/denominator)::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: R.macaulay_resultant([y, x+z, z^2])
            -1

        when there are only 2 polynomials, macaulay resultant degenerates to the traditional resultant::

            sage: R.<x> = PolynomialRing(QQ,1)
            sage: f =  x^2+1; g = x^5+1
            sage: fh = f.homogenize()
            sage: gh = g.homogenize()
            sage: RH = fh.parent()
            sage: f.resultant(g) == RH.macaulay_resultant([fh,gh])
            True

        """
        from sage.matrix.constructor import matrix
        from sage.matrix.constructor import zero_matrix

        if len(args) == 1 and isinstance(args[0],list):
            flist = args[0]
        else:
            flist = args

        if len(flist) <= 0:
            raise TypeError('input list should contain at least 1 polynomial')
        if not all(f.is_homogeneous() for f in flist):
            raise TypeError('resultant for non-homogeneous polynomials is not supported')
        if not all(self.is_parent_of(f) for f in flist):
            raise TypeError('not all inputs are polynomials in the calling ring')

        sparse = kwds.pop('sparse', False)

        U = self.base_ring() # ring of coefficients of self
        dlist = [f.degree() for f in flist]
        xlist = self.gens()
        if len(xlist) != len(dlist):
            raise TypeError('number of polynomials(= %d) must equal number of variables (= %d)'%(len(dlist),len(xlist)))
        n = len(dlist) - 1
        d = sum(dlist) - len(dlist) + 1
        mons = IntegerVectors(d, n+1).list()  # list of exponent-vectors(/lists) of monomials of degree d
        mons_idx = { str(mon) : idx for idx,mon in enumerate(mons)} # a reverse index lookup for monomials
        mons_num = len(mons)
        mons_to_keep = []
        newflist = []
        flist=[[f.exponents(),f.coefficients()] for f in flist] # strip coefficients of the input polynomials
        numer_matrix = zero_matrix(self.base_ring(),mons_num, sparse=sparse)

        for j,mon in enumerate(mons):
            # if monomial is not reduced, then we keep it in the denominator matrix:
            if not self._macaulay_resultant_is_reduced(mon,dlist):
                mons_to_keep.append(j)
            si_mon = self._macaulay_resultant_getS(mon, dlist)
            # Monomial is in S_i under the partition, now we reduce the i'th degree of the monomial
            new_mon = list(mon)
            new_mon[si_mon] -= dlist[si_mon]
            new_f = [[[g[k] + new_mon[k] for k in range(n+1)] for g in flist[si_mon][0]], flist[si_mon][1]]

            for i,mon in enumerate(new_f[0]):
                k = mons_idx[str(mon)]
                numer_matrix[j,k]=new_f[1][i]

        denom_matrix = numer_matrix.matrix_from_rows_and_columns(mons_to_keep,mons_to_keep)
        if denom_matrix.dimensions()[0] == 0: # here we choose the determinant of an empty matrix to be 1
            return U(numer_matrix.det())
        denom_det = denom_matrix.det()
        if denom_det != 0:
            return U(numer_matrix.det()/denom_det)
        # if we get to this point, the determinant of the denominator was 0, and we get the resultant
        # by taking the free coefficient of the quotient of two characteristic polynomials
        poly_num = numer_matrix.characteristic_polynomial('T')
        poly_denom = denom_matrix.characteristic_polynomial('T')
        poly_quo = poly_num.quo_rem(poly_denom)[0]
        return U(poly_quo(0))

    def weyl_algebra(self):
        """
        Return the Weyl algebra generated from ``self``.

        EXAMPLES::

            sage: R = QQ['x,y,z']
            sage: W = R.weyl_algebra(); W
            Differential Weyl algebra of polynomials in x, y, z over Rational Field
            sage: W.polynomial_ring() == R
            True
        """
        from sage.algebras.weyl_algebra import DifferentialWeylAlgebra
        return DifferentialWeylAlgebra(self)

####################
# Leave *all* old versions!

def unpickle_MPolynomialRing_generic_v1(base_ring, n, names, order):
    from .polynomial_ring_constructor import PolynomialRing
    return PolynomialRing(base_ring, n, names=names, order=order)


def unpickle_MPolynomialRing_generic(base_ring, n, names, order):
    from .polynomial_ring_constructor import PolynomialRing

    return PolynomialRing(base_ring, n, names=names, order=order)
