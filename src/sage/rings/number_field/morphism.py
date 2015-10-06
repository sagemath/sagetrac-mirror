r"""
Morphisms between number fields

This module provides classes to represent ring homomorphisms between number
fields (i.e. field embeddings).
"""

from sage.misc.cachefunc import cached_method

from sage.rings.homset import RingHomset_generic
from sage.rings.morphism import RingHomomorphism_im_gens, RingHomomorphism
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.structure.sequence import Sequence

from sage.rings.real_mpfr import RealField_class, RealField
from sage.rings.complex_field import ComplexField_class
from sage.rings.real_mpfi import RealIntervalField_class, RealIntervalField
from sage.rings.complex_interval_field import ComplexIntervalField_class, ComplexIntervalField


class NumberFieldHomset(RingHomset_generic):
    """
    Set of homomorphisms with domain a given number field.

    TESTS::

        sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))
        sage: TestSuite(H).run()
          Failure in _test_category:
        ...
        The following tests failed: _test_elements
    """
    def __init__(self, *args, **kwds):
        """
        Initialize this homset. The ``element_class`` is set, depending
        on the codomain.

        EXAMPLES::

            sage: K = QuadraticField(-1)
            sage: Hom(K, K).element_class
            <class 'sage.rings.number_field.morphism.NumberFieldHomomorphism_im_gens'>
            sage: Hom(K, RR).element_class
            <class 'sage.rings.number_field.morphism.NumberFieldRealEmbedding'>
            sage: Hom(K, CC).element_class
            <class 'sage.rings.number_field.morphism.NumberFieldComplexEmbedding'>
        """
        super(NumberFieldHomset, self).__init__(*args, **kwds)

        if isinstance(self.codomain(),
                (RealField_class, RealIntervalField_class)):
            self.element_class = NumberFieldRealEmbedding
        elif isinstance(self.codomain(),
                (ComplexField_class, ComplexIntervalField_class)):
            self.element_class = NumberFieldComplexEmbedding
        else:
            self.element_class = NumberFieldHomomorphism_im_gens

    def __call__(self, im_gens, check=True):
        """
        Create the homomorphism sending the generators to ``im_gens``.

        EXAMPLES::

            sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))
            sage: phi = H([H.domain().gen()]); phi # indirect doctest
            Ring morphism:
            From: Number Field in a with defining polynomial x^2 + 1
            To:   Number Field in b with defining polynomial x^2 + 1
            Defn: a |--> b
        """
        if isinstance(im_gens, RingHomomorphism_im_gens):
            return self._coerce_impl(im_gens)
        try:
            return self.element_class(self, im_gens, check=check)
        except (NotImplementedError, ValueError):
            raise TypeError("images do not define a valid homomorphism")

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into this homset. The only things that
        coerce canonically into self are elements of self and of homsets equal
        to self.

        EXAMPLES::

            sage: H1 = End(QuadraticField(-1, 'a'))
            sage: H1.coerce(loads(dumps(H1[1]))) # indirect doctest
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 1
              Defn: a |--> -a
        """
        if not isinstance(x, RingHomomorphism_im_gens):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return NumberFieldHomomorphism_im_gens(self, x.im_gens(), check=False)
        raise TypeError

    def _an_element_(self):
        r"""
        Return an element of this set of embeddings.

        EXAMPLES::

            sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))
            sage: H.an_element() # indirect doctest
            Ring morphism:
            From: Number Field in a with defining polynomial x^2 + 1
            To:   Number Field in b with defining polynomial x^2 + 1
            Defn: a |--> b

            sage: H = Hom(QuadraticField(-1, 'a'), QuadraticField(-2, 'b'))
            sage: H.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError: There is no morphism from Number Field in a with
            defining polynomial x^2 + 1 to Number Field in b with defining
            polynomial x^2 + 2
        """
        L = self.list()
        if len(L) != 0:
            return L[0]
        else:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("There is no morphism from {} to {}".format(
                                              self.domain(), self.codomain()))

    def _repr_(self):
        r"""
        String representation of this homset.

        EXAMPLES::

            sage: repr(Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'b'))) # indirect doctest
            'Set of field embeddings from Number Field in a with defining polynomial x^2 + 1 to Number Field in b with defining polynomial x^2 + 1'
            sage: repr(Hom(QuadraticField(-1, 'a'), QuadraticField(-1, 'a'))) # indirect doctest
            'Automorphism group of Number Field in a with defining polynomial x^2 + 1'
        """
        D = self.domain()
        C = self.codomain()
        if C == D:
            return "Automorphism group of {}".format(D)
        else:
            return "Set of field embeddings from {} to {}".format(D, C)

    def order(self):
        """
        Return the order of this set of field homomorphism.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: End(k)
            Automorphism group of Number Field in a with defining polynomial x^2 + 1
            sage: End(k).order()
            2
            sage: k.<a> = NumberField(x^3 + 2)
            sage: End(k).order()
            1

            sage: K.<a> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: End(K).order()
            6
        """
        return Integer(len(self.list()))

    cardinality = order

    @cached_method
    def list(self):
        """
        Return a list of all the elements of self.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 3*x + 1)
            sage: End(K).list()
            [
            Ring endomorphism of Number Field in a with defining polynomial x^3 - 3*x + 1
              Defn: a |--> a,
            Ring endomorphism of Number Field in a with defining polynomial x^3 - 3*x + 1
              Defn: a |--> a^2 - 2,
            Ring endomorphism of Number Field in a with defining polynomial x^3 - 3*x + 1
              Defn: a |--> -a^2 - a + 2
            ]
            sage: Hom(K, CyclotomicField(9))[0] # indirect doctest
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 3*x + 1
              To:   Cyclotomic Field of order 9 and degree 6
              Defn: a |--> -zeta9^4 + zeta9^2 - zeta9

        An example where the codomain is a relative extension::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<b> = K.extension(x^2 + 3)
            sage: Hom(K, L).list()
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in b with defining polynomial x^2 + 3 over its base field
              Defn: a |--> a,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in b with defining polynomial x^2 + 3 over its base field
              Defn: a |--> -1/2*a*b - 1/2*a,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in b with defining polynomial x^2 + 3 over its base field
              Defn: a |--> 1/2*a*b - 1/2*a
            ]
        """
        D = self.domain()
        C = self.codomain()
        if D.degree().divides(C.absolute_degree()):
            roots = D.polynomial().roots(ring=C, multiplicities=False)
            v = [D.hom([r], codomain=C, check=False) for r in roots]
        else:
            v = []
        return Sequence(v, universe=self, check=False, immutable=True, cr=v!=[])

    def __getitem__(self, n):
        r"""
        Return the ``n``th element of ``self.list()``.

        EXAMPLES::

            sage: End(CyclotomicField(37))[3] # indirect doctest
            Ring endomorphism of Cyclotomic Field of order 37 and degree 36
              Defn: zeta37 |--> zeta37^4
        """
        return self.list()[n]


class NumberFieldHomomorphism_im_gens(RingHomomorphism_im_gens):
    def __invert__(self):
        r"""
        Return the inverse of an isomorphism of absolute number fields

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 5)
            sage: tau1, tau2 = K.automorphisms(); tau1, tau2
            (Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
              Defn: a |--> a,
             Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
              Defn: a |--> -a)
            sage: ~tau1
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
             Defn: a |--> a
            sage: ~tau2
            Ring endomorphism of Number Field in a with defining polynomial x^2 + 5
             Defn: a |--> -a

            sage: L.<z> = CyclotomicField(5)
            sage: tau1, tau2, tau3, tau4 = L.automorphisms()
            sage: (tau1, ~tau1)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z)
            sage: (tau2, ~tau2)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^2,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^3)
            sage: (tau4, ~tau4)
            (Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^3,
             Ring endomorphism of Cyclotomic Field of order 5 and degree 4
              Defn: z |--> z^2)

             sage: M.<w> = NumberField(x^4 - 5*x + 5)
             sage: phi = M.hom([z - z^2]); phi
             Ring morphism:
               From: Number Field in w with defining polynomial x^4 - 5*x + 5
               To:   Cyclotomic Field of order 5 and degree 4
               Defn: w |--> -z^2 + z
             sage: phi^-1
             Ring morphism:
               From: Cyclotomic Field of order 5 and degree 4
               To:   Number Field in w with defining polynomial x^4 - 5*x + 5
               Defn: z |--> 3/11*w^3 + 4/11*w^2 + 9/11*w - 14/11
        """
        K = self.domain()
        L = self.codomain()
        if K.degree() != L.degree():
            raise TypeError("Can only invert isomorphisms")
        V, V_into_K, _ = K.vector_space()
        _, _, L_into_W = L.vector_space()
        linear_inverse = ~V.hom([(L_into_W*self*V_into_K)(_) for _ in V.basis()])
        return L.hom([(V_into_K*linear_inverse*L_into_W)(_) for _ in [L.gen()]])

    def preimage(self, y):
        r"""
        Computes a preimage of `y` in the domain, provided one exists.
        Raises a ValueError if `y` has no preimage.

        INPUT:

        - `y` -- an element of the codomain of self.

        OUTPUT:

        Returns the preimage of `y` in the domain, if one exists.
        Raises a ValueError if `y` has no preimage.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 7)
            sage: L.<b> = NumberField(x^4 - 7)
            sage: f = K.embeddings(L)[0]
            sage: f.preimage(3*b^2 - 12/7)
            3*a - 12/7
            sage: f.preimage(b)
            Traceback (most recent call last):
            ...
            ValueError: Element 'b' is not in the image of this homomorphism.

        ::

            sage: F.<b> = QuadraticField(23)
            sage: G.<a> = F.extension(x^3+5)
            sage: f = F.embeddings(G)[0]
            sage: f.preimage(a^3+2*b+3)
            2*b - 2
        """
        # Throughout this method I am using the convention that self is a homomorphism from the number field K to the number field L
        # Therefore, I use the names K and L in place of domain and codomain

        # try to get the cached transformation matrix and vector space isomorphisms if they exist
        try:
            M,LtoV,VtoK = self._transformation_data
        except Exception:
            # get the identifications of K and L with vector spaces over Q
            V,VtoL,LtoV = self.codomain().absolute_vector_space()
            V,VtoK,KtoV = self.domain().absolute_vector_space()
            # construct the transformation matrix from K to L by making the columns be the image of the basis of V_K in V_L using the homomorphism
            from sage.matrix.constructor import matrix
            from sage.rings.all import QQ
            M = matrix(QQ, [LtoV(self(VtoK(e))) for e in V.basis()]).transpose()
            self._transformation_data = (M,LtoV,VtoK)

        # get the coordinate vector of y, solve the linear system, pass to domain
        yvec = LtoV(y)                  # pass from a point in L to its vector space representation
        try:
            xvec = M.solve_right(yvec)      # solve the linear system, throws an exception if there is no solution
        except ValueError:
            raise ValueError("Element '{}' is not in the image of this homomorphism.".format(y))
        return VtoK(xvec)               # pass from the vector space representation of K back to a point in K

# For old pickles
CyclotomicFieldHomomorphism_im_gens = NumberFieldHomomorphism_im_gens


class NumberFieldComplexEmbedding(RingHomomorphism_im_gens):
    """
    An embedding of a number field into a floating-point complex field
    (either a :func:`ComplexField` or :func:`ComplexIntervalField`).
    The image may or may not be actually real.

    .. SEEALSO::

        :class:`NumberFieldRealEmbedding` for embeddings into real
        fields.

    EXAMPLES::

        sage: K.<zeta7> = CyclotomicField(7)
        sage: phi = K.embeddings(CC)[0]
        sage: phi
        Number field embedding morphism:
          From: Cyclotomic Field of order 7 and degree 6
          To:   Complex Field with 53 bits of precision
          Defn: zeta7 |--> -0.900968867902419 - 0.433883739117558*I
        sage: type(phi)
        <class 'sage.rings.number_field.morphism.NumberFieldComplexEmbedding'>
        sage: phi(zeta7)
        -0.900968867902419 - 0.433883739117558*I
        sage: parent(phi(zeta7))
        Complex Field with 53 bits of precision

    Into the complex interval field::

        sage: phi = K.embeddings(CIF)[0]
        sage: phi
        Number field embedding morphism:
          From: Cyclotomic Field of order 7 and degree 6
          To:   Complex Interval Field with 53 bits of precision
          Defn: zeta7 |--> -0.9009688679024191? - 0.4338837391175581?*I
        sage: type(phi)
        <class 'sage.rings.number_field.morphism.NumberFieldComplexEmbedding'>
        sage: phi(zeta7)
        -0.9009688679024191? - 0.4338837391175581?*I
        sage: parent(phi(zeta7))
        Complex Interval Field with 53 bits of precision

    In the following example, the image is real::

        sage: K.<a> = NumberField(x^3 - 3*x + 1)
        sage: phi = K.embeddings(CC)[2]
        sage: phi(a)
        1.53208888623796
        sage: parent(phi(a))
        Complex Field with 53 bits of precision

    AUTHORS:

    - Jeroen Demeyer (2015-10-02): see :trac:`19288`
    """
    def __init__(self, parent, im_gens, check=None):
        """
        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: homset = Hom(K, CC)
            sage: from sage.rings.number_field.morphism import NumberFieldComplexEmbedding
            sage: NumberFieldComplexEmbedding(homset, CC(-3).sqrt())
            Number field embedding morphism:
              From: Number Field in a with defining polynomial x^2 + 3
              To:   Complex Field with 53 bits of precision
              Defn: a |--> 1.73205080756888*I

        It still works if the given root is approximate::

            sage: NumberFieldComplexEmbedding(homset, CC(0,2))
            Number field embedding morphism:
              From: Number Field in a with defining polynomial x^2 + 3
              To:   Complex Field with 53 bits of precision
              Defn: a |--> 1.73205080756888*I

        If we give a completely bogus root, we get an error::

            sage: phi = NumberFieldComplexEmbedding(homset, CC(1))
            Traceback (most recent call last):
            ...
            ArithmeticError: cannot refine polynomial root (not enough steps?)
        """
        # Compute an algebraic version of the given image g0
        try:
            n = len(im_gens)
        except TypeError:
            # No len() => assume it's a single element
            g0 = im_gens
        else:
            if n != 1:
                raise ValueError("number of images must equal number of generators")
            g0 = im_gens[0]

        pol = parent.domain().polynomial()
        from sage.rings.all import AA, QQbar
        if not g0.imag():
            ga = AA.polynomial_root(pol, g0)
        else:
            ga = QQbar.polynomial_root(pol, g0)
        self._im_gen_algebraic = ga

        # We now initialize the RingHomomorphism_im_gens with our
        # refined root (in the correct codomain).  This root isn't
        # actually used to compute images, but it is used for display
        # and other purposes.  That's why it is still useful that this
        # class inherits from RingHomomorphism_im_gens.
        #
        # We always use check=False since the image is almost certainly
        # not an exact root.
        g = parent.codomain()(ga)
        super(NumberFieldComplexEmbedding, self).__init__(parent, (g,), check=False)

    interval_field = staticmethod(ComplexIntervalField)

    def _repr_type(self):
        """
        Return a string used to display this morphism.

        EXAMPLES::

            sage: K = CyclotomicField(11)
            sage: K.embeddings(ComplexField(20))[1]
            Number field embedding morphism:
              From: Cyclotomic Field of order 11 and degree 10
              To:   Complex Field with 20 bits of precision
              Defn: zeta11 |--> -0.95949 + 0.28173*I
        """
        return "Number field embedding"

    def _call_(self, x):
        """
        Evaluate this embedding at ``x``.

        ALGORITHM: to avoid catastrophic cancellation, we compute the
        embedding at higher precision, increasing it if needed until
        the result is precise enough.

        EXAMPLES:

            sage: K.<a> = NumberField(x^2 - x - 104)
            sage: phi = K.embeddings(RR)[0]
            sage: phi(0)
            0.000000000000000
            sage: phi(a)
            -9.71028892833107
            sage: phi(10^100)
            1.00000000000000e100

        We check that :trac:`19288` is fixed::

            sage: D = -4754362903238080086532729679384769785488*a - 46166237460580684378417136025570042177152
            sage: phi(D)
            -4.32015718504111e-35

        Compare with the naive computation which is totally wrong::

            sage: D.polynomial()(phi.im_gens()[0])
            -9.67140655691703e24
        """
        if not x:  # map 0 to 0
            return self.codomain().zero()

        # We evaluate this using interval arithmetic to avoid
        # catastrophic cancellation
        xpol = x.polynomial()
        C = self.codomain()
        prec = C.precision()
        diameter_target = RealField(prec).one() >> (prec-1)
        prec += 7  # arbitrary extra bits
        prec1 = prec

        while True:
            g = self._im_gen_algebraic.interval_fast(self.interval_field(prec))
            value = xpol(g)  # Result as interval
            if value.contains_zero():
                # This is an *embedding*, so the result can only be zero
                # if we start from zero. We already ensured that x isn't
                # zero, so we just keep increasing precision...
                prec *= 2
            elif value.diameter() >= diameter_target:
                # If the interval does not contain zero, we probably
                # have the magnitude correct. So we just need some
                # extra bits of precision, doubling is overkill.
                prec += prec1
            else:
                return self.codomain()(value)

    def im_gens_algebraic(self):
        """
        Return the image of the generator in an algebraic field
        (``AA`` or ``QQbar``) as a 1-tuple.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: phis = K.embeddings(CC)
            sage: [phi.im_gens_algebraic() for phi in phis]
            [(-0.6299605249474365? - 1.091123635971722?*I,),
             (-0.6299605249474365? + 1.091123635971722?*I,),
             (1.259921049894873?,)]
            sage: [parent(phi.im_gens_algebraic()[0]) for phi in phis]
            [Algebraic Field, Algebraic Field, Algebraic Real Field]
        """
        return (self._im_gen_algebraic,)


class NumberFieldRealEmbedding(NumberFieldComplexEmbedding):
    """
    An embedding of a number field into a floating-point real field
    (either a :func:`RealField` or :func:`RealIntervalField`).

    .. SEEALSO::

        :class:`NumberFieldComplexEmbedding` which contains the
        implementation for real and complex embeddings.

    EXAMPLES::

        sage: K.<a> = QuadraticField(3)
        sage: phi = K.embeddings(RR)[0]
        sage: phi
        Number field embedding morphism:
          From: Number Field in a with defining polynomial x^2 - 3
          To:   Real Field with 53 bits of precision
          Defn: a |--> -1.73205080756888
        sage: type(phi)
        <class 'sage.rings.number_field.morphism.NumberFieldRealEmbedding'>
        sage: phi(a)
        -1.73205080756888
        sage: parent(phi(a))
        Real Field with 53 bits of precision

    Into the complex interval field::

        sage: phi = K.embeddings(RIF)[0]
        sage: phi
        Number field embedding morphism:
          From: Number Field in a with defining polynomial x^2 - 3
          To:   Real Interval Field with 53 bits of precision
          Defn: a |--> -1.732050807568878?
        sage: type(phi)
        <class 'sage.rings.number_field.morphism.NumberFieldRealEmbedding'>
        sage: phi(a)
        -1.732050807568878?
        sage: parent(phi(a))
        Real Interval Field with 53 bits of precision
    """
    interval_field = staticmethod(RealIntervalField)


class RelativeNumberFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given relative number field.

    EXAMPLES:

    We construct a homomorphism from a relative field by giving
    the image of a generator::

        sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2)
        sage: phi = L.hom([cuberoot2 * zeta3]); phi
        Relative number field endomorphism of Number Field in cuberoot2 with defining polynomial x^3 - 2 over its base field
          Defn: cuberoot2 |--> zeta3*cuberoot2
                zeta3 |--> zeta3
        sage: phi(cuberoot2 + zeta3)
        zeta3*cuberoot2 + zeta3

    In fact, this phi is a generator for the Kummer Galois group of this
    cyclic extension::

        sage: phi(phi(cuberoot2 + zeta3))
        (-zeta3 - 1)*cuberoot2 + zeta3
        sage: phi(phi(phi(cuberoot2 + zeta3)))
        cuberoot2 + zeta3
    """
    def __call__(self, im_gen, base_hom=None, check=True):
        r"""
        Create a homomorphism in this homset from the given data, which can be:

        - A homomorphism from this number field.
        - A homomorphism from the absolute number field corresponding to this
          relative number field.
        - An element (specifying the image of the generator) of a ring into
          which the base ring coerces.
        - A pair consisting of an element of a ring R and a homomorphism from
          the base ring to R.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: L.<b> = K.extension(x^4 - 2)
            sage: E = End(L)
            sage: E(E[0]) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> b
                    a |--> a
            sage: E(L.absolute_field('c').hom(b+a, L)) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> b
                    a |--> -a
            sage: E(-b*a) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> -a*b
                    a |--> a
            sage: E(-a*b, K.hom([-a])) # indirect doctest
            Relative number field endomorphism of Number Field in b with defining polynomial x^4 - 2 over its base field
              Defn: b |--> -a*b
                    a |--> -a

        It is possible to construct homomorphisms into fields such as
        ``CC``, where calculations are only approximate::

            sage: K.<a> = QuadraticField(-7)
            sage: f = K.hom([CC(sqrt(-7))])
            sage: x = polygen(K)
            sage: L.<b> = K.extension(x^2 - a - 5)
            sage: L.Hom(CC)(f(a + 5).sqrt(), f)
            Relative number field morphism:
              From: Number Field in b with defining polynomial x^2 - a - 5 over its base field
              To:   Complex Field with 53 bits of precision
              Defn: b |--> 2.30833860703888 + 0.573085617291335*I
                    a |--> -3.46944695195361e-18 + 2.64575131106459*I
        """
        if isinstance(im_gen, RingHomomorphism_im_gens):
            # Then it must be a homomorphism from the corresponding
            # absolute number field
            abs_hom = im_gen
            K = abs_hom.domain()
            if K != self.domain().absolute_field(K.variable_name()):
                raise TypeError("domain of morphism must be absolute field of domain.")
            from_K, to_K = K.structure()
            if abs_hom.domain() != K:
                raise ValueError("domain of absolute homomorphism must be absolute field of domain.")
            if abs_hom.codomain() != self.codomain():
                raise ValueError("codomain of absolute homomorphism must be codomain of this homset.")
            return RelativeNumberFieldHomomorphism_from_abs(self, abs_hom)
        if isinstance(im_gen, RelativeNumberFieldHomomorphism_from_abs):
            return self._coerce_impl(im_gen)
        if base_hom is None:
            base_hom = self.default_base_hom()
        if isinstance(im_gen, (list, tuple)) and len(im_gen) == 1:
            im_gen = im_gen[0]
        if check:
            im_gen = self.codomain()(im_gen)
        return self._from_im(im_gen, base_hom, check=check)

    def _coerce_impl(self, x):
        r"""
        Canonically coerce ``x`` into this homset. This will only work if ``x``
        is already in the homset.

        EXAMPLES::

            sage: L.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
            sage: E = End(L)
            sage: E.coerce(loads(dumps(E[0])))  # indirect doctest
            Relative number field endomorphism of Number Field in a with defining polynomial x^3 - x + 1 over its base field
              Defn: a |--> a
                    b |--> b
        """
        if not isinstance(x, RelativeNumberFieldHomomorphism_from_abs):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return RelativeNumberFieldHomomorphism_from_abs(self, x.abs_hom())
        raise TypeError

    def _from_im(self, im_gen, base_hom, check=True):
        """
        Return the homomorphism that acts on the base as given and
        sends the generator of the domain to im_gen.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: L.<b> = K.extension(x^3 - x + 1)
            sage: End(L)._from_im( -3/23*a*b^2 + (-9/46*a - 1/2)*b + 2/23*a, K.hom([-a], K))
            Relative number field endomorphism of Number Field in b with defining polynomial x^3 - x + 1 over its base field
              Defn: b |--> -3/23*a*b^2 + (-9/46*a - 1/2)*b + 2/23*a
                    a |--> -a
        """
        K = self.domain().absolute_field('a')
        from_K, to_K = K.structure()
        a = from_K(K.gen())
        # We just have to figure out where a goes to
        # under the morphism defined by im_gen and base_hom.
        L = self.codomain()
        R = L['x']
        f = R([base_hom(x) for x in a.list()])
        b = f(im_gen)
        abs_hom = K.hom([b], check=check)
        return RelativeNumberFieldHomomorphism_from_abs(self, abs_hom)

    def default_base_hom(self):
        r"""
        Pick an embedding of the base field of self into the codomain of this
        homset. This is done in an essentially arbitrary way.

        EXAMPLES::

            sage: L.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
            sage: M.<c> = NumberField(x^4 + 80*x^2 + 36)
            sage: Hom(L, M).default_base_hom()
            Ring morphism:
              From: Number Field in b with defining polynomial x^2 + 23
              To:   Number Field in c with defining polynomial x^4 + 80*x^2 + 36
              Defn: b |--> 1/12*c^3 + 43/6*c
        """
        try:
            return self.__default_base_hom
        except AttributeError:
            pass
        v = self.domain().base_field().embeddings(self.codomain())
        if len(v) == 0:
            raise ValueError("no way to map base field to codomain.")
        self.__default_base_hom = v[0]
        return v[0]

    @cached_method
    def list(self):
        """
        Return a list of all the elements of self (for which the domain
        is a relative number field).

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^2 + x + 1, x^3 + 2])
            sage: End(K).list()
            [
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: a |--> a
                    b |--> b,
            ...
            Relative number field endomorphism of Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: a |--> a
                    b |--> -b*a - b
            ]

        An example with an absolute codomain::

            sage: K.<a, b> = NumberField([x^2 - 3, x^2 + 2])
            sage: Hom(K, CyclotomicField(24, 'z')).list()
            [
            Relative number field morphism:
              From: Number Field in a with defining polynomial x^2 - 3 over its base field
              To:   Cyclotomic Field of order 24 and degree 8
              Defn: a |--> z^6 - 2*z^2
                    b |--> -z^5 - z^3 + z,
            ...
            Relative number field morphism:
              From: Number Field in a with defining polynomial x^2 - 3 over its base field
              To:   Cyclotomic Field of order 24 and degree 8
              Defn: a |--> -z^6 + 2*z^2
                    b |--> z^5 + z^3 - z
            ]
        """
        D = self.domain()
        C = self.codomain()
        D_abs = D.absolute_field('a')
        v = [self(f, check=False) for f in D_abs.Hom(C).list()]
        v = Sequence(v, universe=self, check=False, immutable=True, cr=v!=[])
        return v


class RelativeNumberFieldHomomorphism_from_abs(RingHomomorphism):
    r"""
    A homomorphism from a relative number field to some other ring, stored as a
    homomorphism from the corresponding absolute field.
    """

    def __init__(self, parent, abs_hom):
        r"""
        EXAMPLES::

            sage: K.<a, b> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: f = K.hom(-a*b - a, K); f
            Relative number field endomorphism of Number Field in a with defining polynomial x^3 + 2 over its base field
              Defn: a |--> (-b - 1)*a
                    b |--> b
            sage: type(f)
            <class 'sage.rings.number_field.morphism.RelativeNumberFieldHomomorphism_from_abs'>
        """
        RingHomomorphism.__init__(self, parent)
        self.__abs_hom = abs_hom
        K = abs_hom.domain()
        from_K, to_K = K.structure()
        self.__K = K
        self.__from_K = from_K
        self.__to_K = to_K

    def abs_hom(self):
        r"""
        Return the corresponding homomorphism from the absolute number field.

        EXAMPLES::

            sage: K.<a, b> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: K.hom(a, K).abs_hom()
            Ring morphism:
              From: Number Field in a with defining polynomial x^6 - 3*x^5 + 6*x^4 - 3*x^3 - 9*x + 9
              To:   Number Field in a with defining polynomial x^3 + 2 over its base field
              Defn: a |--> a - b
        """
        return self.__abs_hom

    def _repr_type(self):
        r"""
        A short string to identify the type of this homomorphism.

        EXAMPLES::

            sage: K.<a, b> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: K.hom(a, K)._repr_type()
            'Relative number field'
        """
        return "Relative number field"

    def im_gens(self):
        r"""
        Return the images of the generators under this map.

        EXAMPLES::

            sage: K.<a, b> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: K.hom(a, K).im_gens()
            [a, b]
        """
        try:
            return self.__im_gens
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        v = Sequence([self(x) for x in D.gens()], universe=C, check=False, immutable=True)
        self.__im_gens = v
        return v

    def _cmp_(self, other):
        """
        Compare

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^2 - 2, x^2 - 3])
            sage: e, u, v, w = End(K)
            sage: all([u^2 == e, u*v == w, u != e])
            True
        """
        return cmp(self.abs_hom(), other.abs_hom())

    __cmp__ = _cmp_

    def _repr_defn(self):
        r"""
        Return a string describing the images of the generators under this map.

        EXAMPLES::

            sage: K.<a, b> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: K.hom(a, K)._repr_defn()
            'a |--> a\nb |--> b'
        """
        D = self.domain()
        ig = self.im_gens()
        return '\n'.join(['%s |--> %s'%(D.gen(i), ig[i]) for\
                       i in range(D.ngens())])

    def _call_(self, x):
        r"""
        Evaluate this map at the element ``x``. This is done by first
        converting ``x`` to an element of the absolute field and then
        evaluating ``self.abs_hom()`` on it.

        EXAMPLES::

            sage: K.<a, b> = NumberField( [x^3 + 2, x^2 + x + 1] )
            sage: K.hom(a*b, K)(17 + 3*a + 2*b) # indirect doctest
            3*b*a + 2*b + 17
        """
        return self.__abs_hom(self.__to_K(x))


class CyclotomicFieldHomset(NumberFieldHomset):
    """
    Set of homomorphisms with domain a given cyclotomic field.

    EXAMPLES::

        sage: End(CyclotomicField(16))
        Automorphism group of Cyclotomic Field of order 16 and degree 8

    TESTS::

        sage: K.<z> = CyclotomicField(16)
        sage: E = End(K)
        sage: E(E[0]) # indirect doctest
        Ring endomorphism of Cyclotomic Field of order 16 and degree 8
          Defn: z |--> z
        sage: E(z^5) # indirect doctest
        Ring endomorphism of Cyclotomic Field of order 16 and degree 8
          Defn: z |--> z^5
        sage: E(z^6) # indirect doctest
        Traceback (most recent call last):
        ...
        TypeError: images do not define a valid homomorphism

    ::

        sage: E = End(CyclotomicField(16))
        sage: E.coerce(E[0]) # indirect doctest
        Ring endomorphism of Cyclotomic Field of order 16 and degree 8
          Defn: zeta16 |--> zeta16
        sage: E.coerce(17) # indirect doctest
        Traceback (most recent call last):
        ...
        TypeError: no canonical coercion from Integer Ring to Automorphism group of Cyclotomic Field of order 16 and degree 8
    """
    @cached_method
    def list(self):
        """
        Return a list of all the elements of self (for which the domain
        is a cyclotomic field).

        EXAMPLES::

            sage: K.<z> = CyclotomicField(12)
            sage: G = End(K); G
            Automorphism group of Cyclotomic Field of order 12 and degree 4
            sage: [g(z) for g in G]
            [z, z^3 - z, -z, -z^3 + z]
            sage: L.<a, b> = NumberField([x^2 + x + 1, x^4 + 1])
            sage: L
            Number Field in a with defining polynomial x^2 + x + 1 over its base field
            sage: Hom(CyclotomicField(12), L)[3]
            Ring morphism:
              From: Cyclotomic Field of order 12 and degree 4
              To:   Number Field in a with defining polynomial x^2 + x + 1 over its base field
              Defn: zeta12 |--> -b^2*a
            sage: list(Hom(CyclotomicField(5), K))
            []
            sage: Hom(CyclotomicField(11), L).list()
            []
        """
        D = self.domain()
        C = self.codomain()
        z = D.gen()
        n = z.multiplicative_order()
        if not n.divides(C.zeta_order()):
            v =[]
        else:
            if D == C:
                w = z
            else:
                w = C.zeta(n)
            v = [self([w**k], check=False) for k in Zmod(n) if k.is_unit()]
        v = Sequence(v, universe=self, check=False, immutable=True, cr=v!=[])
        return v
