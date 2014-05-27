r"""
Base class for multivariate polynomial rings
"""

include 'sage/ext/stdsage.pxi'


from sage.structure.parent_gens cimport ParentWithGens
import sage.misc.latex
import multi_polynomial_ideal
from term_order import TermOrder
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polydict import PolyDict
from sage.misc.latex import latex_variable_name
import multi_polynomial_element
import polynomial_ring
from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()
from sage.rings.polynomial.polynomial_ring_constructor import polynomial_default_category

def is_MPolynomialRing(x):
    return bool(PY_TYPE_CHECK(x, MPolynomialRing_generic))

cdef class MPolynomialRing_generic(sage.rings.ring.CommutativeRing):
    def __init__(self, base_ring, n, names, order):
        """
        Create a polynomial ring in several variables over a commutative ring.

        EXAMPLES::

            sage: R.<x,y> = ZZ['x,y']; R
            Multivariate Polynomial Ring in x, y over Integer Ring
            sage: class CR(CommutativeRing):
            ...       def __init__(self):
            ...           CommutativeRing.__init__(self,self)
            ...       def __call__(self,x):
            ...           return None
            sage: cr = CR()
            sage: cr.is_commutative()
            True
            sage: cr['x,y']
            Multivariate Polynomial Ring in x, y over <class '....CR_with_category'>

        TESTS:

        Check that containment works correctly (ticket #10355)::

            sage: A1.<a> = PolynomialRing(QQ)
            sage: A2.<a,b> = PolynomialRing(QQ)
            sage: 3 in A2
            True
            sage: A1(a) in A2
            True
        """
        if base_ring not in _CommutativeRings:
            raise TypeError, "The base ring %s is not a commutative ring"%base_ring

        n = int(n)
        if n < 0:
            raise ValueError, "Multivariate Polynomial Rings must " + \
                  "have more than 0 variables."
        order = TermOrder(order,n)
        self.__ngens = n
        self.__term_order = order
        self._has_singular = False #cannot convert to Singular by default
        self._magma_cache = {}
        # Ring.__init__ already does assign the names.
        # It would be a mistake to call ParentWithGens.__init__
        # as well, assigning the names twice.
        #ParentWithGens.__init__(self, base_ring, names)
        sage.rings.ring.Ring.__init__(self, base_ring, names,
                                      category=polynomial_default_category(base_ring.category(),n>1))

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
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.categories.pushout import MultiPolynomialFunctor
        return MultiPolynomialFunctor(self.variable_names(), self.term_order()), self.base_ring()

    def irrelevant_ideal(self):
        """
        Return the irrelevant ideal of this multivariate polynomial
        ring, which is the ideal generated by all of the indeterminate
        generators of this ring.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: R.irrelevant_ideal()
            Ideal (x, y, z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.ideal(self.gens(), check=False)

    def completion(self, p, prec=20, extras=None):
        """
        Return the completion of self with respect to the ideal generated
        by the variable(s) ``p``.

        INPUT:

        - ``p`` -- variable or tuple of variables

        - ``prec`` -- default precision of resulting power series ring

        - ``extras`` -- ignored; present for backward compatibility

        EXAMPLES::

            sage: P.<x,y,z,w> = PolynomialRing(ZZ)
            sage: P.completion((w,x,y))
            Multivariate Power Series Ring in w, x, y over Univariate Polynomial Ring in z over Integer Ring
            sage: P.completion((w,x,y,z))
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

        """
        if p in self or isinstance(p,str) and set(p).issubset(set([str(g) for g in self.gens()])):
            p = tuple([p])
        elif isinstance(p,(list,tuple)):
            p = tuple(p)
        else:
            raise TypeError("input %s is not of type list or tuple, and is not a variable of %s" % (p,self))

        try:
            from sage.rings.power_series_ring import PowerSeriesRing
            if len(p) > 1:
                return PowerSeriesRing(self.remove_var(*p), names=p, default_prec=prec)
            elif len(p) == 1:
                p = p[0]
                return PowerSeriesRing(self.remove_var(p), name=str(p), default_prec=prec)
            else:
                return self

            # p is one variable; possibly as a string
            # check to make sure it is one of the variables of self
            #if p in self or type(p) == str and set(p).issubset(set([str(g) for g in self.gens()])):

        except ValueError:
            raise TypeError, "Cannot complete %s with respect to %s" % (self, p)

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
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
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

        EXAMPLE::

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
        from polydict import ETuple
        other_vars = list(x.parent().variable_names())
        name_mapping = [(other_vars.index(var) if var in other_vars else -1) for var in self.variable_names()]
        K = self.base_ring()
        D = {}
        var_range = range(len(self.variable_names()))
        for ix, a in x.dict().iteritems():
            ix = ETuple([0 if name_mapping[t] == -1 else ix[name_mapping[t]] for t in var_range])
            D[ix] = K(a)
        return D

    def __richcmp__(left, right, int op):
        return (<ParentWithGens>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        if not is_MPolynomialRing(right):
            return cmp(type(left),type(right))
        else:
            return cmp((left.base_ring(), left.__ngens, left.variable_names(), left.__term_order),
                       (right.base_ring(), (<MPolynomialRing_generic>right).__ngens, right.variable_names(), (<MPolynomialRing_generic>right).__term_order))

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
            sage: print P.repr_long()
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

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            # all that is needed is that elements of the base ring
            # of the polynomial ring canonically coerce into codomain.
            # Since poly rings are free, any image of the gen
            # determines a homomorphism
            codomain._coerce_(self.base_ring()(1))
        except TypeError:
            return False
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

    def is_finite(self):
        if self.ngens() == 0:
            return self.base_ring().is_finite()
        return False

    def is_field(self, proof = True):
        """
        Return True if this multivariate polynomial ring is a field, i.e.,
        it is a ring in 0 generators over a field.
        """
        if self.ngens() == 0:
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
            raise ValueError, "Generator not defined."
        return self._gens[int(n)]

    #def gens(self):
        #return self._gens

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
            from polynomial_ring_constructor import PolynomialRing
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
        from sage.rings.arith import binomial
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

        .. note::

            We do not check if the provided index/rank is within the allowed
            range. If it is not an infinite loop will occur.
        """
        from sage.combinat import choose_nk
        comb = choose_nk.from_rank(i, n+d-1, n-1)
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
        from sage.rings.arith import binomial
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
        - ``terms`` -- number of terms requested (default: 5)
        - ``choose_degree`` -- choose degrees of monomials randomly first
          rather than monomials uniformly random.
        - ``**kwargs`` -- passed to the random element generator of the base
          ring

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: P.random_element(2, 5)
            -6/5*x^2 + 2/3*z^2 - 1

            sage: P.random_element(2, 5, choose_degree=True)
            -1/4*x*y - 1/5*x*z - 1/14*y*z - z^2

        Stacked rings::

            sage: R = QQ['x,y']
            sage: S = R['t,u']
            sage: S.random_element(degree=2, terms=1)
            -3*x*y + 5/2*y^2 - 1/2*x - 1/4*y + 4
            sage: S.random_element(degree=2, terms=1)
            (-1/2*x^2 - x*y - 2/7*y^2 + 3/2*x - y)*t*u

        Default values apply if no degree and/or number of terms is
        provided::

            sage: random_matrix(QQ['x,y,z'], 2, 2)
            [        2*y^2 - 2/27*y*z - z^2 + 2*z        1/2*x*y - 1/2*y^2 + 2*x - 2*y]
            [-1/27*x^2 + 2/5*y^2 - 1/10*z^2 - 2*z              -13*y^2 + 2/3*z^2 + 2*y]

            sage: random_matrix(QQ['x,y,z'], 2, 2, terms=1, degree=2)
            [-1/4*x    1/2]
            [ 1/3*x    x*y]

            sage: P.random_element(0, 1)
            -1

            sage: P.random_element(2, 0)
            0

            sage: R.<x> = PolynomialRing(Integers(3), 1)
            sage: R.random_element()
            x + 1
        """
        k = self.base_ring()
        n = self.ngens()

        counts, total = self._precomp_counts(n, degree)

        if terms is None:
            if total >= 5:
                terms = 5
            else:
                terms = total

        if terms < 0:
            raise TypeError, "Cannot compute polynomial with a negative number of terms."
        elif terms == 0:
            return self._zero_element
        if degree == 0:
            if terms != 1:
                raise TypeError, "Cannot compute polynomial with more terms than exist."
            return k.random_element(**kwargs)


        from sage.combinat.integer_vector import IntegerVectors
        from sage.rings.arith import binomial

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
                M = map(tuple, M)
            else:
                M = [list(IntegerVectors(_d,n)) for _d in xrange(degree+1)]
                Mbar = []
                for mi in xrange(terms):
                    d = ZZ.random_element(0,len(M)) #choose degree at random
                    m = ZZ.random_element(0,len(M[d])) # choose monomial at random
                    Mbar.append( M[degree].pop(m) ) # remove and insert
                    if len(M[degree]) == 0:
                        M.pop(degree) # bookkeeping
                M = map(tuple, Mbar)

        else:
            raise TypeError, "Cannot compute polynomial with more terms than exist."

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

        from polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(base_ring, self.ngens(), names, order=order)


####################
# Leave *all* old versions!

def unpickle_MPolynomialRing_generic_v1(base_ring, n, names, order):
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    return PolynomialRing(base_ring, n, names=names, order=order)


def unpickle_MPolynomialRing_generic(base_ring, n, names, order):
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    return PolynomialRing(base_ring, n, names=names, order=order)
