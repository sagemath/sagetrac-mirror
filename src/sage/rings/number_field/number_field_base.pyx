"""
Base class for all number fields


TESTS::

    sage: k = NumberField(x^2 + 1, 'i'); k == loads(dumps(k))
    True
"""

def is_NumberField(x):
    """
    Return True if x is of number field type.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_base import is_NumberField
        sage: is_NumberField(NumberField(x^2+1,'a'))
        True
        sage: is_NumberField(QuadraticField(-97,'theta'))
        True
        sage: is_NumberField(CyclotomicField(97))
        True

    Note that the rational numbers QQ are a number field.::

        sage: is_NumberField(QQ)
        True
        sage: is_NumberField(ZZ)
        False
    """
    return isinstance(x, NumberField)

from sage.rings.ring cimport Field

cdef class NumberField(Field):
    r"""
    Base class for all number fields.
    """
    # This token docstring is mostly there to prevent Sphinx from pasting in
    # the docstring of the __init__ method inherited from IntegralDomain, which
    # is rather confusing.

    def ring_of_integers(self, *args, **kwds):
        r"""
        Synomym for ``self.maximal_order(...)``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.ring_of_integers()
            Gaussian Integers in Number Field in a with defining polynomial x^2 + 1
        """
        return self.maximal_order(*args, **kwds)

    def OK(self, *args, **kwds):
        r"""
        Synomym for ``self.maximal_order(...)``.

        EXAMPLES::

            sage: NumberField(x^3 - 2,'a').OK()
            Maximal Order in Number Field in a with defining polynomial x^3 - 2
        """
        return self.maximal_order(*args, **kwds)

    def maximal_order(self):
        """
        Return the maximal order, i.e., the ring of integers of this
        number field.

        EXAMPLES::

            sage: NumberField(x^3 - 2,'b').maximal_order()
            Maximal Order in Number Field in b with defining polynomial x^3 - 2
        """
        raise NotImplementedError

    def is_finite(self):
        """
        Return False since number fields are not finite.

        EXAMPLES::

            sage: z = polygen(QQ)
            sage: K.<theta, beta> = NumberField([z^3 - 3, z^2 + 1])
            sage: K.is_finite()
            False
            sage: K.order()
            +Infinity
        """
        return False

    def is_absolute(self):
        """
        Return True if self is viewed as a single extension over Q.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3+2)
            sage: K.is_absolute()
            True
            sage: y = polygen(K)
            sage: L.<b> = NumberField(y^2+1)
            sage: L.is_absolute()
            False
            sage: QQ.is_absolute()
            True
        """
        raise NotImplementedError

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES::

            sage: NumberField(x^3 - 2, 'a').signature()
            (1, 1)
        """
        raise NotImplementedError

    def degree(self):
        """
        Return the degree of this number field.

        EXAMPLES::

            sage: NumberField(x^3 + 9, 'a').degree()
            3
        """
        raise NotImplementedError

    def discriminant(self):
        """
        Return the discriminant of this number field.

        EXAMPLES::

            sage: NumberField(x^3 + 9, 'a').discriminant()
            -243
        """
        raise NotImplementedError

    def minkowski_bound(self):
        r"""
        Return the Minkowski bound associated to this number field,
        which is a bound B so that every integral ideal is equivalent
        modulo principal fractional ideals to an integral ideal of
        norm at most B.

        .. SEEALSO::

            :meth:`~bach_bound`

        OUTPUT:

        symbolic expression or Rational

        EXAMPLES:

        The Minkowski bound for `\QQ[i]` tells us that the class
        number is 1::

            sage: K = QQ[I]
            sage: B = K.minkowski_bound(); B
            4/pi
            sage: B.n()
            1.27323954473516

        We compute the Minkowski bound for `\QQ[\sqrt[3]{2}]`::

            sage: K = QQ[2^(1/3)]
            sage: B = K.minkowski_bound(); B
            16/3*sqrt(3)/pi
            sage: B.n()
            2.94042077558289
            sage: int(B)
            2

        We compute the Minkowski bound for `\QQ[\sqrt{10}]`, which has class
        number 2::

            sage: K = QQ[sqrt(10)]
            sage: B = K.minkowski_bound(); B
            sqrt(10)
            sage: int(B)
            3
            sage: K.class_number()
            2

        We compute the Minkowski bound for `\QQ[\sqrt{2}+\sqrt{3}]`::

            sage: K.<y,z> = NumberField([x^2-2, x^2-3])
            sage: L.<w> = QQ[sqrt(2) + sqrt(3)]
            sage: B = K.minkowski_bound(); B
            9/2
            sage: int(B)
            4
            sage: B == L.minkowski_bound()
            True
            sage: K.class_number()
            1

        The bound of course also works for the rational numbers::

            sage: QQ.minkowski_bound()
            1
        """
        _, s = self.signature()
        n = self.absolute_degree()
        d = self.absolute_discriminant().abs().sqrt()
        from sage.symbolic.constants import pi
        if s > 0:
            return d * (4/pi)**s * n.factorial() / (n**n)
        else:
            return d * n.factorial() / (n**n)

    def bach_bound(self):
        r"""
        Return the Bach bound associated to this number field.
        Assuming the General Riemann Hypothesis, this is a bound B so
        that every integral ideal is equivalent modulo principal
        fractional ideals to an integral ideal of norm at most B.

        .. SEEALSO::

            :meth:`~minkowski_bound`

        OUTPUT:

        symbolic expression or the Integer 1

        EXAMPLES:

        We compute both the Minkowski and Bach bounds for a quadratic
        field, where the Minkowski bound is much better::

            sage: K = QQ[sqrt(5)]
            sage: K.minkowski_bound()
            1/2*sqrt(5)
            sage: K.minkowski_bound().n()
            1.11803398874989
            sage: K.bach_bound()
            12*log(5)^2
            sage: K.bach_bound().n()
            31.0834847277628

        We compute both the Minkowski and Bach bounds for a bigger
        degree field, where the Bach bound is much better::

            sage: K = CyclotomicField(37)
            sage: K.minkowski_bound().n()
            7.50857335698544e14
            sage: K.bach_bound().n()
            191669.304126267

        The bound of course also works for the rational numbers:
            sage: QQ.minkowski_bound()
            1
        """
        ans = 12 * abs(self.discriminant()).log()**2
        if ans == 0: # rational numbers
            from sage.rings.integer import Integer
            return Integer(1)
        return ans


    # Approximate embeddings for comparisons with respect to the order of RR or
    # CC

    def _init_embedding_approx(self):
        r"""
        Initialize the approximation of embeddings.

        This should be called only once.

        TESTS::

            sage: K.<a> = NumberField(x^3 - x^2 - x - 1, embedding=1)
            sage: K._get_embedding_approx(0)   # indirect doctest
            1.839286755214161?
        """

        if self._gen_approx is not None or self._embedding is None:
            return

        from sage.rings.qqbar import AA
        from sage.rings.real_lazy import RLF
        codomain = self._embedding.codomain()
        if codomain is AA or codomain is RLF:
            self._gen_approx = []
            self._embedded_real = 1

    cpdef _get_embedding_approx(self, size_t i):
        r"""
        Return an interval approximation of the generator of this number field.

        OUTPUT:

        A real interval element with precision `53 \times 2^i`.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: p = x^5 - 3*x + 1
            sage: a_AA = AA.polynomial_root(p, RIF(0,1))
            sage: K.<a> = NumberField(p, embedding=a_AA)
            sage: K._get_embedding_approx(2)
            0.3347341419433526870750989624732833071257517550374285560578335863?
            sage: K._get_embedding_approx(1)
            0.33473414194335268707509896247329?
            sage: K._get_embedding_approx(1).str(style='brackets')
            '[0.334734141943352687075098962473280 .. 0.334734141943352687075098962473287]'


            sage: K._get_embedding_approx(2).prec()
            212
            sage: K._get_embedding_approx(1).prec()
            106
            sage: K._get_embedding_approx(0).prec()
            53

        If a real embedding is not specified, this method will result in an error::

            sage: N.<g> = NumberField(x^3+2)
            sage: N._get_embedding_approx(1)
            Traceback (most recent call last):
            ...
            ValueError: No embedding set. You need to specify a a real embedding.


        .. SEEALSO::

            :class:` RealIntervalField_class <sage.rings.real_mpfi.RealIntervalField_class>`
        """
        if self._embedded_real and i < len(self._gen_approx):
            return self._gen_approx[i]

        cdef size_t j
        if self._embedded_real:
            j = len(self._gen_approx)
            from sage.rings.real_mpfi import RealIntervalField
            gen = self._embedding.gen_image()
            while j <= i:
                self._gen_approx.append(RealIntervalField(53 << j)(gen))
                j += 1
            return self._gen_approx[i]
        else:
            raise ValueError("No embedding set. You need to specify a a real embedding.")

