r"""
Extension of rings

See :func:`RingExtension` for documentation.

AUTHOR:

- Xavier Caruso (2016)
"""

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.pushout import pushout
from sage.categories.algebras import Algebras
from sage.rings.integer_ring import IntegerRing
from sage.rings.ring import CommutativeRing, CommutativeAlgebra
from sage.categories.action import BaseActionOnRing


class AlgebraFromMorphism(CommutativeAlgebra, UniqueRepresentation):
    r"""
    Create a ring extension

    Do not call this function directly
    Use :func:`RingExtension` instead

    INPUT:

    - ``defining_morphism`` -- a ring homomorphism

    - ``coerce`` -- boolean (specify whether defining_morphism
    is a coercion map or not)

    OUTPUT:

    The extension defined by ``defining_morphism``

    EXAMPLES::

        sage: K = GF(5^2)
        sage: L = GF(5^4)
        sage: E = L/K; E
        Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

        sage: from sage.rings.algebra_from_morphism import AlgebraFromMorphism
        sage: isinstance(E, AlgebraFromMorphism)
        True

    See :func:`RingExtension` for more documentation

    AUTHOR:

    - Xavier Caruso (2016)
    """
    def __init__(self, defining_morphism, coerce):
        r"""
        TESTS::

        The attribute _coerce indicates if the defining morphism is
        a coercion map::

            sage: K = GF(5^2)
            sage: E = RingExtension(K, K, K.frobenius_endomorphism())
            sage: E._coerce
            False

            sage: L = GF(5^4)
            sage: E1 = L/K
            sage: E1._coerce
            True

        Creating an extension of extensions is not allowed::

            sage: M = GF(5^8)
            sage: E2 = M/K
            sage: E2/E1
            Traceback (most recent call last):
            ...
            TypeError: Creating an extension of extensions is not allowed
        """
        base = defining_morphism.domain()
        ring = defining_morphism.codomain()

        # To avoid problems...
        if isinstance(base, AlgebraFromMorphism) or isinstance(ring, AlgebraFromMorphism):
            raise TypeError("Creating an extension of extensions is not allowed")

        # We do not want to have a coercion map base -> self
        # So we first set the base to ring.base() (which indeed coerces to self)
        CommutativeAlgebra.__init__(self, ring.base(), category=Algebras(base)) 
        # and then update the base
        self._base = base
        self._ring = ring
        self._defining_morphism = defining_morphism
        self._coerce = coerce

        # Set left action of base on self
        # and a coercion map self -> ring
        from sage.rings.morphism import AlgebraToRing_coercion
        self._unset_coercions_used()
        self._populate_coercion_lists_(
            action_list = [BaseActionOnRing(self)],
            embedding = AlgebraToRing_coercion(self.Hom(ring), check=False))

        from sage.structure.element import AlgebraFMElement
        self.element_class = AlgebraFMElement

    def _element_constructor_(self, x, *args, **kwargs):
        r"""
        Convert ``x`` into an element of this parent

        EXAMPLES::

            sage: K = GF(5^4)
            sage: L = GF(5^8)
            sage: E = L/K

        Conversion of an element of the ring::

            sage: a = L.random_element()
            sage: a.parent()
            Finite Field in z8 of size 5^8
            sage: a.parent() is L
            True

            sage: aE = E(a)
            sage: aE.parent()
            Finite Field in z8 of size 5^8 viewed as an algebra over Finite Field in z4 of size 5^4
            sage: aE.parent() is E
            True

        Conversion from another extension::

            sage: k = GF(5^2)
            sage: F = K/k
            sage: b = K.gen(); b
            z4
            sage: bF = F(b); bF
            z4

            sage: bE = E(bF); bE
            4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4
            sage: bE.parent()
            Finite Field in z8 of size 5^8 viewed as an algebra over Finite Field in z4 of size 5^4
            sage: bE.parent() is E
            True
        """
        from sage.structure.element import AlgebraFMElement
        if isinstance(x, AlgebraFMElement):
            x = x.element_in_ring()
        elt = self._ring(x, *args, **kwargs)
        return self.element_class(self, elt)

    def from_base_ring(self, r):
        r"""
        Return the canonical embedding of ``r`` into ``self``.

        INPUT:

        - ``r`` -- an element of the base of the ring of this extension

        IMPORTANT REMARK:

        The base in which ``r`` is supposed to lie is not the base of this 
        extension (``self``) but the base of its ring (``self.ring()``) which
        is canonically attached to ``self.ring()`` and in particular in not
        affected by the defining morphism of ``self``.

        Note in addition that ``r`` is mapped into ``self`` using the
        canonical embedding. Again the defining morphism does not interfere
        at this stage.

        EXAMPLES::

            sage: k = GF(5)
            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = L/K; E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

            sage: E.base()
            Finite Field in z2 of size 5^2
            sage: E.ring().base()
            Finite Field of size 5

            sage: x = E.from_base_ring(k(2)); x
            2
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

            sage: x = E.from_base_ring(z2);
            Traceback (most recent call last):
            ...
            TypeError: z2 is not an element of the base of Finite Field in z4 of size 5^4 (= Finite Field of size 5)

        TESTS::

            sage: R.<t> = K[]
            sage: F = RingExtension(R, K, K.frobenius_endomorphism())
            sage: F.from_base_ring(z2)
            z2
            sage: F.defining_morphism()(z2)
            4*z2 + 1
        """
        ring = self._ring
        if r not in ring.base():
            raise TypeError("%s is not an element of the base of %s (= %s)" % (r, ring, ring.base()))
        elt = ring(r)
        return self.element_class(self, elt)

    def _pushout_(self, other):
        r"""
        Construct the pushout of this extension and ``other``
        This is called by :func:`sage.categories.pushout.pushout`.

        ALGORITHM:

        If ``other`` is not an extension, we return the pushout of
        ``self.ring()`` and ``other``

        Otherwise, if the defining morphism of ``self`` or ``other`` 
        is not a coercion map, we return the pushout of ``self.ring()`` 
        and ``other.ring()``

        Otherwise, we compute the pushout of ``self.ring()`` and 
        ``other.ring()`` and call it ``ring``. We then check whether
        there exists a coercion map between ``self.base()`` and 
        ``other.base()``. If no coercion map is found, we return 
        ``ring``

        Otherwise, we call ``base`` the domain of the coercion
        map and return the extension ``ring/base``

        .. TODO::

            Handle the case where defining morphisms are not
            coercion maps

            If no coercion map is discovered between ``self.base()`` 
            and ``other.base()``, compute the pullover ``base`` of 
            these two rings (i.e. the greatest ring coercing to both) 
            and return the extension ``ring/base``

        TESTS::

            sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = L1/K1
            sage: K2 = GF(3^2); L2 = GF(3^16); E2 = L2/K2
            sage: E2.has_coerce_map_from(E1)
            True
            sage: E2._pushout_(E1) is E2
            True
            sage: E1._pushout_(E2) is E2
            True

            sage: E1._pushout_(L2) is L2
            True

        When the defining morphism is not a coercion map, the pushout
        is never an extension, always a ring::

            sage: E1p = RingExtension(L1, K1, K1.frobenius_endomorphism())
            sage: E1p._pushout_(E2) is L2
            True

            sage: E2p = RingExtension(L2, K2, K2.frobenius_endomorphism())
            sage: E1p._pushout_(E2p) is L2
            True
        """
        base = None
        if isinstance(other,AlgebraFromMorphism):
            # TODO: implement pushout when:
            # - defining_morphism are not coercion maps
            #   (this implies to check the equality of two morphisms)
            # - there is no direct coercion map between the bases
            #   but there exists a third ring which coerces to eash base
            #   Question: how can one find the greatest such third ring?
            if self._coerce and other._coerce:
                sbase = self._base
                obase = other._base
                if sbase.has_coerce_map_from(obase):
                    base = obase
                elif obase.has_coerce_map_from(sbase):
                    base = sbase
            ring = pushout(self._ring, other.ring())
        else:
            ring = pushout(self._ring, other)
        if base is None:
            return ring
        else:
            from sage.rings.algebra_from_morphism_constructor import RingExtension
            return RingExtension(ring,base)

    def _repr_(self):
        r"""
        Return a string representation of this extension.

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = L/K; 
            sage: E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

            sage: E._repr_()
            'Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2'
        """
        return "%s viewed as an algebra over %s" % (self._ring, self._base)

    def _latex_(self):
        r"""
        Return a latex representation of this extension.

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = L/K; 
            sage: latex(E)
            \Bold{F}_{5^{4}}/\Bold{F}_{5^{2}}
        """
        from sage.misc.latex import latex
        return "%s/%s" % (latex(self._ring), latex(self._base))

    def _coerce_map_from_(self, other):
        r"""
        An extension L1/K1 coerces to another extension L2/K2 when L1 coerces 
        to L2 and K2 coerces to K1 (be careful at the direction) and the 
        defining morphisms are coercion maps.

        .. TODO::

            Relax the condition on defining morphisms (and only require
            that they agree on K2)

        TESTS::

            sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = L1/K1
            sage: K2 = GF(3^2); L2 = GF(3^16); E2 = L2/K2
            sage: E2.has_coerce_map_from(E1)  # indirect doctest
            True

        In the next example, the defining morphisms are not coercion maps.
        So there is no coercion between the extensions (through it probably
        should).

            sage: E1p = RingExtension(L1, K1, K1.frobenius_endomorphism())
            sage: E2p = RingExtension(L2, K2, K2.frobenius_endomorphism())
            sage: E2p.has_coerce_map_from(E1p)  # indirect doctest
            False
        """
        if isinstance(other, AlgebraFromMorphism):
            if self._coerce and other._coerce:
                return other.base().has_coerce_map_from(self._base) and self._ring.has_coerce_map_from(other.ring())                

    def defining_morphism(self):
        r"""
        Return the defining morphism of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = L/K; E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z4 of size 5^4
              Defn: z2 |--> z4^3 + z4^2 + z4 + 3

            sage: E2 = RingExtension(L, K, K.frobenius_endomorphism()); E2
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E2.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z4 of size 5^4
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
        """
        return self._defining_morphism

    def base(self):
        r"""
        Return the base of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = L/K; E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E.base()
            Finite Field in z2 of size 5^2
            sage: E.base() is K
            True

        Note that the base of the extension is generally different
        from the base of the ring of the extension::

            sage: E.ring().base()
            Finite Field of size 5
        """
        return self._base

    def ring(self):
        r"""
        Return the ring of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = L/K; E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E.ring()
            Finite Field in z4 of size 5^4
            sage: E.ring() is L
            True
        """
        return self._ring

    def _an_element_(self):
        r"""
        Return an element of this extension

        This element is constructed by converting an element of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = L/K; E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

            sage: x = E.an_element()  # indirect doctest
            sage: x
            0
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

            sage: a = E.ring().an_element()
            sage: a == x
            True
        """
        elt = self._ring.an_element()
        return self.element_class(self, elt)

    def random_element(self):
        r"""
        Return a random element of this extension

        This element is constructed by converting a random element of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = L/K; E
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2

            sage: x = E.random_element(); x  # random
            0
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over Finite Field in z2 of size 5^2
        """
        elt = self._ring.random_element()
        return self.element_class(self, elt)

    def scalar_restriction(self, newbase):
        r"""
        Return the scalar restriction of this extension to 
        the new base ``newbase``.

        INPUT:

        - ``newbase`` -- it can be either
          * a ring which coerces to the base of ``self`` or
          * an extension whose ring coerces to the base of ``self`` or
          * a morphism whose codomain coerces to the base of ``self``

        OUTPUT:

        The scalar restriction of ``self`` to the ``newbase``.

        It is defined as follows.
        In each above three cases, we get a ring homorphism ``f`` 
        from a domain ``domain`` to the base of ``self``.
        The scalar restriction of ``self`` to ``newbase`` is the
        extension ``self.ring()/domain`` with defining morphism the 
        compositum of ``f`` with the defining morphism of ``self``

        EXAMPLES::

            sage: k = GF(5^2)
            sage: K = GF(5^4)
            sage: L = GF(5^8)
            sage: E = L/K; E
            Finite Field in z8 of size 5^8 viewed as an algebra over Finite Field in z4 of size 5^4

            sage: E1 = E.scalar_restriction(k); E1
            Finite Field in z8 of size 5^8 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E1.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn: z2 |--> z8^7 + z8^6 + z8^5 + z8^4 + 3*z8^3 + 3*z8^2 + z8 + 3

            sage: Frob = k.frobenius_endomorphism()
            sage: E2 = E.scalar_restriction(Frob); E2
            Finite Field in z8 of size 5^8 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E2.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
                    then
                      Ring morphism:
                      From: Finite Field in z4 of size 5^4
                      To:   Finite Field in z8 of size 5^8
                      Defn: z4 |--> 4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4

            sage: F = RingExtension(K, k, Frob)
            sage: E3 = E.scalar_restriction(F); E3
            Finite Field in z8 of size 5^8 viewed as an algebra over Finite Field in z2 of size 5^2
            sage: E3.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
                    then
                      Ring morphism:
                      From: Finite Field in z4 of size 5^4
                      To:   Finite Field in z8 of size 5^8
                      Defn: z4 |--> 4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4

            sage: E2 is E3
            True
        """
        coerce = self._coerce
        if isinstance(newbase, AlgebraFromMorphism):
            if coerce:
                coerce = newbase._coerce
            newbase = newbase.defining_morphism()
        elif isinstance(newbase, CommutativeRing):
            if self._base.has_coerce_map_from(newbase):
                newbase = self._base.coerce_map_from(newbase)
            else:
                raise TypeError("No coercion map from %s to %s" % (newbase, self._base))
        else:
            domain = newbase.domain()
            codomain = newbase.codomain()
            if coerce:
                coerce = (newbase == codomain.coerce_map_from(domain))
            if self._base.has_coerce_map_from(codomain):
                newbase = newbase.post_compose(self._base.coerce_map_from(codomain))
            else:
                raise TypeError("No coercion map from %s to %s" % (codomain, self._base))
        defining_morphism = self._defining_morphism.pre_compose(newbase)
        return AlgebraFromMorphism(defining_morphism, coerce)
