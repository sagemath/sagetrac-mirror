"""
Galois Groups of Number Fields

AUTHORS:

- William Stein (2004, 2005): initial version
- David Loeffler (2009): rewrite to give explicit homomorphism groups

TESTS:

Standard test of pickleability::

    sage: G = NumberField(x^3 + 2, 'alpha').galois_group(type="pari"); G
    Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in alpha with defining polynomial x^3 + 2
    sage: G == loads(dumps(G))
    True

    sage: G = NumberField(x^3 + 2, 'alpha').galois_group(names='beta'); G
    Galois group of Galois closure in beta of Number Field in alpha with defining polynomial x^3 + 2
    sage: G == loads(dumps(G))
    True
"""

from sage.structure.sage_object import SageObject
from sage.groups.perm_gps.permgroup import PermutationGroup_generic, standardize_generator
from sage.structure.category_object import normalize_names

from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.libs.pari.all import pari
from sage.rings.infinity import infinity
from sage.rings.number_field.number_field import refine_embedding
from sage.rings.number_field.morphism import NumberFieldHomomorphism_im_gens
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet


class GaloisGroup_v1(SageObject):
    r"""
    A wrapper around a class representing an abstract transitive group.

    This is just a fairly minimal object at present.  To get the underlying
    group, do ``G.group()``, and to get the corresponding number field do
    ``G.number_field()``. For a more sophisticated interface use the
    ``type=None`` option.

    EXAMPLES::

        sage: K = QQ[2^(1/3)]
        sage: G = K.galois_group(type="pari"); G
        Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in a with defining polynomial x^3 - 2 with a = 1.259921049894873?
        sage: G.order()
        6
        sage: G.group()
        PARI group [6, -1, 2, "S3"] of degree 3
        sage: G.number_field()
        Number Field in a with defining polynomial x^3 - 2 with a = 1.259921049894873?
    """

    def __init__(self, group, number_field):
        """
        Create a Galois group.

        EXAMPLES::

            sage: NumberField([x^2 + 1, x^2 + 2],'a').galois_group(type="pari")
            Galois group PARI group [4, 1, 2, "E(4) = 2[x]2"] of degree 4 of the Number Field in a0 with defining polynomial x^2 + 1 over its base field
        """
        self.__group = group
        self.__number_field = number_field

    def __eq__(self, other):
        """
        Compare two number field Galois groups.

        First the number fields are compared, then the Galois groups
        if the number fields are equal.  (Of course, if the number
        fields are the same, the Galois groups are automatically
        equal.)

        EXAMPLES::

            sage: G = NumberField(x^3 + 2, 'alpha').galois_group(type="pari")
            sage: H = QQ[sqrt(2)].galois_group(type="pari")
            sage: H == G
            False
            sage: H == H
            True
            sage: G == G
            True
        """
        if not isinstance(other, GaloisGroup_v1):
            return False
        if self.__number_field == other.__number_field:
            return True
        if self.__group == other.__group:
            return True
        return False

    def __ne__(self, other):
        """
        Test for unequality.

        EXAMPLES::

            sage: G = NumberField(x^3 + 2, 'alpha').galois_group(type="pari")
            sage: H = QQ[sqrt(2)].galois_group(type="pari")
            sage: H != G
            True
            sage: H != H
            False
            sage: G != G
            False
        """
        return not (self == other)

    def __repr__(self):
        """
        Display print representation of a Galois group.

        EXAMPLES::

            sage: G = NumberField(x^4 + 2*x + 2, 'a').galois_group(type="pari")
            sage: G.__repr__()
            'Galois group PARI group [24, -1, 5, "S4"] of degree 4 of the Number Field in a with defining polynomial x^4 + 2*x + 2'
        """
        return "Galois group %s of the %s" % (self.__group,
                                              self.__number_field)

    def group(self):
        """
        Return the underlying abstract group.

        EXAMPLES::

            sage: G = NumberField(x^3 + 2*x + 2, 'theta').galois_group(type="pari")
            sage: H = G.group(); H
            PARI group [6, -1, 2, "S3"] of degree 3
            sage: P = H.permutation_group(); P
            Transitive group number 2 of degree 3
            sage: sorted(P)
            [(), (2,3), (1,2), (1,2,3), (1,3,2), (1,3)]
        """
        return self.__group

    def order(self):
        """
        Return the order of this Galois group.

        EXAMPLES::

            sage: G = NumberField(x^5 + 2, 'theta_1').galois_group(type="pari"); G
            Galois group PARI group [20, -1, 3, "F(5) = 5:4"] of degree 5 of the Number Field in theta_1 with defining polynomial x^5 + 2
            sage: G.order()
            20
        """
        return self.__group.order()

    def number_field(self):
        """
        Return the number field of which this is the Galois group.

        EXAMPLES::

            sage: G = NumberField(x^6 + 2, 't').galois_group(type="pari"); G
            Galois group PARI group [12, -1, 3, "D(6) = S(3)[x]2"] of degree 6 of the Number Field in t with defining polynomial x^6 + 2
            sage: G.number_field()
            Number Field in t with defining polynomial x^6 + 2
        """
        return self.__number_field


class GaloisGroup_v2(PermutationGroup_generic):
    r"""
    The Galois group of an (absolute) number field.

    .. NOTE::

        We define the Galois group of a non-normal field K to be the
        Galois group of its Galois closure L, and elements are stored as
        permutations of the roots of the defining polynomial of L, *not* as
        permutations of the roots (in L) of the defining polynomial of K. The
        latter would probably be preferable, but is harder to implement. Thus
        the permutation group that is returned is always simply-transitive.

        The 'arithmetical' features (decomposition and ramification groups,
        Artin symbols etc) are only available for Galois fields.
    """

    def __init__(self, number_field, algorithm='pari', names=None, gc_numbering=None):
        r"""
        Create a Galois group.

        EXAMPLES::

            sage: QuadraticField(-23,'a').galois_group()
            Galois group of Number Field in a with defining polynomial x^2 + 23 with a = 4.795831523312720?*I
            sage: NumberField(x^3 - 2, 'b').galois_group()
            Traceback (most recent call last):
            ...
            TypeError: You must specify the name of the generator.
            sage: NumberField(x^3 - 2, 'b').galois_group(names="c")
            Galois group of Galois closure in c of Number Field in b with defining polynomial x^3 - 2

        TESTS::

            sage: F.<z> = CyclotomicField(7)
            sage: G = F.galois_group()

        We test that a method inherited from PermutationGroup_generic returns
        the right type of element (see :trac:`133`)::

            sage: phi = G.random_element()
            sage: type(phi) is G.element_class
            True
            sage: phi(z) # random
            z^3
        """
        self._number_field = number_field
        self._default_algorithm = algorithm
        if gc_numbering is None:
            gc_numbering = False if algorithm == 'magma' else True
        self._gc_numbering = gc_numbering
        if names is None:
            # add a c for Galois closure
            names = number_field.variable_name() + 'c'
        self._gc_names = normalize_names(1, names)
        # We do only the parts of the initialization of PermutationGroup_generic
        # that don't depend on _gens
        from sage.categories.permutation_groups import PermutationGroups
        category = PermutationGroups().FinitelyGenerated().Finite()
        # Note that we DON'T call the __init__ method for PermutationGroup_generic
        # Instead, the relevant attributes are computed lazily
        super(PermutationGroup_generic, self).__init__(category=category)

    def _get_algorithm(self, algorithm):
        return self._default_algorithm if algorithm is None else algorithm

    @cached_method(key=_get_algorithm)
    def _pol_galgp(self, algorithm=None):
        algorithm = self._get_algorithm(algorithm)
        f = self._number_field.absolute_polynomial()
        return f.galois_group(pari_group=True, algorithm=algorithm)

    def order(self, algorithm=None, recompute=False):
        algorithm = self._get_algorithm(algorithm)
        # We cache manually since we're computing the same quantity using different backends
        if not recompute and '_size' in self.__dict__:
            return self._size # _order is a method on permgroup
        K = self._number_field
        if K.absolute_degree() < 12 or algorithm != "pari":
            self._size = self._pol_galgp(algorithm=algorithm).order()
        else:
            self._size = self._galois_closure.degree()
        return self._size

    def easy_order(self, algorithm=None):
        algorithm = self._get_algorithm(algorithm)
        if '_size' in self.__dict__:
            return self._size
        K = self._number_field
        if K.absolute_degree() < 12 or algorithm != "pari":
            self._size = self._pol_galgp(algorithm=algorithm).order()
            return self._size

    def transitive_number(self, algorithm=None, recompute=False):
        """
        Regardless of the value of ``gc_numbering``, this gives the transitive number
        for the action on the roots of the defining polynomial of the original number field,
        not the Galois closure.
        """
        algorithm = self._get_algorithm(algorithm)
        # We cache manually since we're computing the same quantity using different backends
        if not recompute and '_t' in self.__dict__:
            return self._t
        K = self._number_field
        if K.absolute_degree() < 12 or algorithm != "pari":
            self._t = self._pol_galgp(algorithm=algorithm).transitive_number()
        elif not self._gc_numbering:
            self._t = ZZ(self.gap().TransitiveIdentification())
        else:
            self._t = ZZ(PermutationGroup(self._short_gens(self._gens)).TransitiveIdentification())
        return self._t

    def transitive_label(self):
        return "%sT%s" % (self._number_field.degree(), self.transitive_number())

    def pari_label(self):
        return self._pol_galgp().label()

    @cached_method
    def signature(self):
        """
        Returns 1 if contained in the alternating group, -1 otherwise.
        """
        if self._number_field.degree() < 12:
            return self._pol_galgp().signature()
        elif self._number_field.absolute_polynomial().discriminant().is_square():
            return ZZ(1)
        else:
            return ZZ(-1)

    # We compute various attributes lazily so that we can support quick lookup
    # of some that are more easily computed.  This allows us to emulate
    # having initialized as a permutation group.
    @lazy_attribute
    def _gcdata(self):
        K = self._number_field
        if self.is_galois():
            return K, K.hom(K.gen(), K)
        else:
            return K.galois_closure(names=self._gc_names, map=True)

    @lazy_attribute
    def _galois_closure(self):
        return self._gcdata[0]

    @lazy_attribute
    def _gc_map(self):
        return self._gcdata[1]

    @lazy_attribute
    def _pari_data(self):
        return self._galois_closure.__pari__().galoisinit()

    @lazy_attribute
    def _elts(self):
        # PARI computes all the elements of self anyway, so we might as well store them
        return sorted([self(x, check=False) for x in self._pari_data[5]])

    @lazy_attribute
    def _deg(self):
        """
        The number of moved points in the permutation representation.

        This will be the degree of the original number field if `_gc_numbering``
        is ``False``, or the degree of the Galois closure otherwise.

        EXAMPES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._deg
            5
            sage: G = K.galois_group(gc_numbering=True); G._deg
            20
        """
        if self._gc_numbering:
            return self.order()
        else:
            return self._number_field.degree()

    @lazy_attribute
    def _domain(self):
        """
        The integers labeling the roots on which this Galois group acts.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._domain
            {1, 2, 3, 4, 5}
            sage: G = K.galois_group(gc_numbering=True); G._domain
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}
        """
        return FiniteEnumeratedSet(range(1, self._deg+1))

    @lazy_attribute
    def _domain_to_gap(self):
        """
        Dictionary implementing the identity (used by PermutationGroup_generic).

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False)
            sage: G._domain_to_gap[5]
            5
        """
        return dict((key, i+1) for i, key in enumerate(self._domain))

    @lazy_attribute
    def _domain_from_gap(self):
        """
        Dictionary implementing the identity (used by PermutationGroup_generic).

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=True)
            sage: G._domain_from_gap[20]
            20
        """
        return dict((i+1, key) for i, key in enumerate(self._domain))

    @lazy_attribute
    def _gens(self):
        """
        Computes the generators as permutations.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5-2)
            sage: G = K.galois_group(gc_numbering=False); G
            Galois group 5T3 (5:4) with order 20 of x^5 - 2
            sage: G._gens
            [(1,2,3,5), (1,4,3,2,5)]
            sage: G = K.galois_group(gc_numbering=True)
            sage: G._gens
            [(1,2,15,3)(4,19,11,8)(5,20,13,7)(6,9,10,16)(12,17,18,14),
             (1,7,17,11,6)(2,8,5,9,18)(3,12,16,13,19)(4,14,20,15,10)]
        """
        if self._gc_numbering:
            gens = [standardize_generator(x, as_cycles=True) for x in self._pari_data[6]]
            if not gens:
                gens = [()]
            gens = [self.element_class(x, self, check=False) for x in gens]
            return sorted(set(gens))
        else:
            self._gc_numbered = G = self._number_field.galois_group(algorithm=self._default_algorithm, names=self._gc_names, gc_numbering=True)
            self._galois_closure = L = G._galois_closure
            gens = [g.as_hom() for g in G._gens]
            if gens:
                # We add None so that we're 1-indexed
                roots = [None] + self._number_field.absolute_polynomial().roots(L, multiplicities=False)
                new_gens = []
                for g in gens:
                    seen = set()
                    cycles = []
                    for start in range(1, len(roots)):
                        if start in seen:
                            continue
                        cycle = [start]
                        r = roots[start]
                        while True:
                            r = g(r)
                            i = roots.index(r)
                            seen.add(i)
                            if i == start:
                                break
                            cycle.append(i)
                        cycles.append(tuple(cycle))
                    new_gens.append(cycles)
            else:
                new_gens = [()]
            # Want order to match G's, so don't sort
            return [self.element_class(x, self, check=False) for x in new_gens]

    def __call__(self, x, check=True):
        r"""
        Create an element of this Galois group.

        INPUT:

        - ``x`` -- should be one of:

          - the integer 1, denoting the identity of G
          - an element of G
          - a permutation of the right length which defines an element of G, or anything that
            coerces into a permutation of the right length
          - an abstract automorphism of the underlying number field.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-23)
            sage: G = K.galois_group()
            sage: G(1)
            ()
            sage: G(G.gens()[0])
            (1,2)
            sage: G([(1,2)])
            (1,2)
            sage: G(K.hom(-a, K))
            (1,2)
        """
        if x == 1:
            return self.identity()

        from sage.rings.number_field.morphism import NumberFieldHomomorphism_im_gens
        if isinstance(x, NumberFieldHomomorphism_im_gens) and x.parent() == self.number_field().Hom(self.number_field()):
            l = [g for g in self if g.as_hom() == x]
            if len(l) != 1:
                raise ArithmeticError
            return l[0]
        return self.element_class(x, self, check=check)

    def is_galois(self):
        r"""
        Whether the underlying number field is Galois

        EXAMPLES::

            sage: NumberField(x^3 - x + 1,'a').galois_group(names='b').is_galois()
            False
            sage: NumberField(x^2 - x + 1,'a').galois_group().is_galois()
            True
        """
        K = self._number_field
        if K.degree() < 12:
            return self._pol_galgp().order() == K.degree()
        else:
            return len(K.automorphisms()) == K.degree()

    def ngens(self):
        r"""
        Number of generators of this Galois group

        EXAMPLES::

            sage: QuadraticField(-23, 'a').galois_group().ngens()
            1
        """
        return len(self._gens)

    def _repr_(self):
        r"""
        String representation of this Galois group

        EXAMPLES::

            sage: G = QuadraticField(-23, 'a').galois_group()
            sage: G._repr_()
            'Galois group of x^2 + 23'
            sage: G = NumberField(x^3 - 2, 'a').galois_group(names='b')
            sage: G._repr_()
            'Galois group of Galois closure in b of Number Field in a with defining polynomial x^3 - 2'
        """
        K = self.number_field()
        f = K.defining_polynomial()
        if K.degree() < 12:
            plabel = self.pari_label().split('=')[-1].strip()
            tlabel = "%sT%s (%s) with order %s " % (K.degree(), self.transitive_number(), plabel, self.order())
        else:
            tlabel = ""
        if K.degree() < 12 or self.is_galois():
            return "Galois group %sof %s" % (tlabel, f)
        else:
            return "Galois group %sof (non-Galois) %s" % (tlabel, f)

    def number_field(self):
        r"""
        The ambient number field.

        EXAMPLES::

            sage: K = NumberField(x^3 - x + 1, 'a')
            sage: K.galois_group(names='b').number_field() is K
            True
        """
        return self._number_field

    def splitting_field(self):
        r"""
        The Galois closure of the ambient number field.

        EXAMPLES::

            sage: K = NumberField(x^3 - x + 1, 'a')
            sage: K.galois_group(names='b').splitting_field()
            Number Field in b with defining polynomial x^6 - 6*x^4 + 9*x^2 + 23
            sage: L = QuadraticField(-23, 'c'); L.galois_group().splitting_field() is L
            True
        """
        return self._galois_closure

    def list(self):
        r"""
        List of the elements of self.

        EXAMPLES::

            sage: NumberField(x^3 - 3*x + 1,'a').galois_group().list()
            [(), (1,2,3), (1,3,2)]
        """
        return self._elts

    def unrank(self, i):
        r"""
        Return the ``i``-th element of ``self``.

        INPUT:

        - ``i`` -- integer between ``0`` and ``n-1`` where
          ``n`` is the cardinality of this set

        EXAMPLES::

            sage: G = NumberField(x^3 - 3*x + 1,'a').galois_group()
            sage: [G.unrank(i) for i in range(G.cardinality())]
            [(), (1,2,3), (1,3,2)]

        TESTS::

            sage: G = NumberField(x^3 - 3*x + 1,'a').galois_group()
            sage: L = [G.unrank(i) for i in range(G.cardinality())]
            sage: L == G.list()
            True
        """
        return self._elts[i]

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: G = NumberField(x^3 - 3*x + 1,'a').galois_group()
            sage: list(G) == G.list()
            True
        """
        return iter(self._elts)

    def subgroup(self, elts):
        r"""
        Return the subgroup of self with the given elements. Mostly for internal use.

        EXAMPLES::

            sage: G = NumberField(x^3 - x - 1, 'a').galois_closure('b').galois_group()
            sage: G.subgroup([ G(1), G([(1,2,3),(4,5,6)]), G([(1,3,2),(4,6,5)]) ])
            Subgroup [(), (1,2,3)(4,5,6), (1,3,2)(4,6,5)] of Galois group of Number Field in b with defining polynomial x^6 - 6*x^4 + 9*x^2 + 23

        Subgroups can be specified using generators (:trac:`26816`)::

            sage: K.<a> = NumberField(x^6 - 6*x^4 + 9*x^2 + 23)
            sage: G = K.galois_group()
            sage: list(G)
            [(),
             (1,2,3)(4,5,6),
             (1,3,2)(4,6,5),
             (1,4)(2,6)(3,5),
             (1,5)(2,4)(3,6),
             (1,6)(2,5)(3,4)]
            sage: g = G[1]
            sage: h = G[3]
            sage: list(G.subgroup([]))
            [()]
            sage: list(G.subgroup([g]))
            [(), (1,2,3)(4,5,6), (1,3,2)(4,6,5)]
            sage: list(G.subgroup([h]))
            [(), (1,4)(2,6)(3,5)]
            sage: list(G.subgroup([g,h])) == list(G)
            True
        """
        return GaloisGroup_subgroup(self, elts)

    # Proper number theory starts here. All the functions below make no sense
    # unless the field is Galois.

    @cached_method
    def _ramgroups(self, P):
        """
        Compute ramification data using Pari.

        INPUT:

        - ``P`` -- a prime ideal.

        OUTPUT:

        A pari vector holding the decomposition group, inertia groups, and higher ramification groups.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 2*x^2 + 2,'b').galois_closure()
            sage: P = K.ideal([17, a^2])
            sage: G = K.galois_group()
            sage: G._ramgroups(P)
            [[[Vecsmall([8, 7, 6, 5, 4, 3, 2, 1])], Vecsmall([2])]]
        """
        K = self.number_field()
        P = K.ideal_monoid()(P).pari_prime()
        return pari(K).idealramgroups(self._pari_data, P)

    def decomposition_group(self, P):
        r"""
        Decomposition group of a prime ideal P, i.e. the subgroup of elements
        that map P to itself. This is the same as the Galois group of the
        extension of local fields obtained by completing at P.

        This function will raise an error if P is not prime or the given number
        field is not Galois.

        P can also be an infinite prime, i.e. an embedding into `\RR` or `\CC`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 2*x^2 + 2,'b').galois_closure()
            sage: P = K.ideal([17, a^2])
            sage: G = K.galois_group()
            sage: G.decomposition_group(P)
            Subgroup [(), (1,8)(2,7)(3,6)(4,5)] of Galois group of Number Field in a with defining polynomial x^8 - 20*x^6 + 104*x^4 - 40*x^2 + 1156
            sage: G.decomposition_group(P^2)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (...) is not a prime ideal
            sage: G.decomposition_group(17)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (17) is not a prime ideal

        An example with an infinite place::

            sage: L.<b> = NumberField(x^3 - 2,'a').galois_closure(); G=L.galois_group()
            sage: x = L.places()[0]
            sage: G.decomposition_group(x).order()
            2
        """
        if not self.is_galois():
            raise TypeError("Decomposition groups only defined for Galois extensions")

        if isinstance(P, NumberFieldHomomorphism_im_gens):
            if self.number_field().is_totally_real():
                return self.subgroup([self.identity()])
            else:
                return self.subgroup([self.identity(), self.complex_conjugation(P)])
        else:
            return self.ramification_group(P, -1)

    def complex_conjugation(self, P=None):
        """
        Return the unique element of self corresponding to complex conjugation,
        for a specified embedding P into the complex numbers. If P is not
        specified, use the "standard" embedding, whenever that is well-defined.

        EXAMPLES::

            sage: L.<z> = CyclotomicField(7)
            sage: G = L.galois_group()
            sage: conj = G.complex_conjugation(); conj
            (1,4)(2,5)(3,6)
            sage: conj(z)
            -z^5 - z^4 - z^3 - z^2 - z - 1

        An example where the field is not CM, so complex conjugation really
        depends on the choice of embedding::

            sage: L = NumberField(x^6 + 40*x^3 + 1372,'a')
            sage: G = L.galois_group()
            sage: [G.complex_conjugation(x) for x in L.places()]
            [(1,3)(2,6)(4,5), (1,5)(2,4)(3,6), (1,2)(3,4)(5,6)]
        """
        if P is None:
            Q = self.number_field().specified_complex_embedding()
            if Q is None:
                raise ValueError("No default complex embedding specified")
            P = Q

        P = refine_embedding(P, infinity)

        if not self.number_field().is_galois():
            raise TypeError("Extension is not Galois")
        if self.number_field().is_totally_real():
            raise TypeError("No complex conjugation (field is real)")

        g = self.number_field().gen()
        gconj = P(g).conjugate()
        elts = [s for s in self if P(s(g)) == gconj]
        if len(elts) != 1:
            raise ArithmeticError("Something has gone very wrong here")
        return elts[0]

    def ramification_group(self, P, v):
        """
        Return the vth ramification group of self for the prime P, i.e. the set
        of elements s of self such that s acts trivially modulo P^(v+1). This
        is only defined for Galois fields.

        EXAMPLES::

            sage: K.<b> = NumberField(x^3 - 3,'a').galois_closure()
            sage: G=K.galois_group()
            sage: P = K.primes_above(3)[0]
            sage: G.ramification_group(P, 3)
            Subgroup [(), (1,2,4)(3,5,6), (1,4,2)(3,6,5)] of Galois group of Number Field in b with defining polynomial x^6 + 243
            sage: G.ramification_group(P, 5)
            Subgroup [()] of Galois group of Number Field in b with defining polynomial x^6 + 243
        """
        if not self.is_galois():
            raise TypeError("Ramification groups only defined for Galois extensions")
        ramdata = self._ramgroups(P)
        if v < -1:
            raise ValueError("v must be at least -1")
        elif v + 1 >= len(ramdata):
            return self.subgroup([])
        else:
            return self.subgroup(ramdata[v + 1][0])

    def inertia_group(self, P):
        """
        Return the inertia group of the prime P, i.e. the group of elements acting
        trivially modulo P. This is just the 0th ramification group of P.

        EXAMPLES::

            sage: K.<b> = NumberField(x^2 - 3,'a')
            sage: G = K.galois_group()
            sage: G.inertia_group(K.primes_above(2)[0])
            Subgroup [(), (1,2)] of Galois group of Number Field in b with defining polynomial x^2 - 3
            sage: G.inertia_group(K.primes_above(5)[0])
            Subgroup [()] of Galois group of Number Field in b with defining polynomial x^2 - 3
        """
        if not self.is_galois():
            raise TypeError("Inertia groups only defined for Galois extensions")
        return self.ramification_group(P, 0)

    def ramification_breaks(self, P):
        r"""
        Return the set of ramification breaks of the prime ideal P, i.e. the
        set of indices i such that the ramification group `G_{i+1} \ne G_{i}`.
        This is only defined for Galois fields.

        EXAMPLES::

            sage: K.<b> = NumberField(x^8 - 20*x^6 + 104*x^4 - 40*x^2 + 1156)
            sage: G = K.galois_group()
            sage: P = K.primes_above(2)[0]
            sage: G.ramification_breaks(P)
            {1, 3, 5}
            sage: min( [ G.ramification_group(P, i).order() / G.ramification_group(P, i+1).order() for i in G.ramification_breaks(P)] )
            2
        """
        if not self.is_galois():
            raise TypeError("Ramification breaks only defined for Galois extensions")
        ramdata = self._ramgroups(P)
        from sage.sets.set import Set
        return Set([i - 1 for (i, (v, w)) in enumerate(zip(ramdata[:-1], ramdata[1:]))
                    if v[1] != w[1]] + [len(ramdata) - 2])

    def artin_symbol(self, P):
        r"""
        Return the Artin symbol `\left(\frac{K /
        \QQ}{\mathfrak{P}}\right)`, where K is the number field of self,
        and `\mathfrak{P}` is an unramified prime ideal. This is the unique
        element s of the decomposition group of `\mathfrak{P}` such that `s(x) = x^p \bmod
        \mathfrak{P}`, where p is the residue characteristic of `\mathfrak{P}`.

        EXAMPLES::

            sage: K.<b> = NumberField(x^4 - 2*x^2 + 2, 'a').galois_closure()
            sage: G = K.galois_group()
            sage: [G.artin_symbol(P) for P in K.primes_above(7)]
            [(1,5)(2,6)(3,7)(4,8), (1,5)(2,6)(3,7)(4,8), (1,4)(2,3)(5,8)(6,7), (1,4)(2,3)(5,8)(6,7)]
            sage: G.artin_symbol(17)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (17) is not prime
            sage: QuadraticField(-7,'c').galois_group().artin_symbol(13)
            (1,2)
            sage: G.artin_symbol(K.primes_above(2)[0])
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (...) is ramified
        """
        if not self.is_galois():
            raise TypeError("Artin symbols only defined for Galois extensions")

        P = self.number_field().ideal_monoid()(P)
        if not P.is_prime():
            raise ValueError("%s is not prime" % P)
        p = P.smallest_integer()
        t = []
        gens = self.number_field().ring_of_integers().ring_generators()
        for s in self.decomposition_group(P):
            w = [(s(g) - g**p).valuation(P) for g in gens]
            if min(w) >= 1:
                t.append(s)
        if len(t) > 1:
            raise ValueError("%s is ramified" % P)
        return t[0]


class GaloisGroup_subgroup(GaloisGroup_v2):
    r"""
    A subgroup of a Galois group, as returned by functions such as ``decomposition_group``.
    """

    def __init__(self, ambient, elts):
        r"""
        Return the subgroup of this Galois group generated by the
        given elements.

        It is generally better to use the :meth:`subgroup` method of
        the parent group.

        EXAMPLES::

            sage: from sage.rings.number_field.galois_group import GaloisGroup_subgroup
            sage: G = NumberField(x^3 - x - 1, 'a').galois_closure('b').galois_group()
            sage: GaloisGroup_subgroup( G, [ G(1), G([(1,2,3),(4,5,6)]), G([(1,3,2),(4,6,5)])])
            Subgroup [(), (1,2,3)(4,5,6), (1,3,2)(4,6,5)] of Galois group of Number Field in b with defining polynomial x^6 - 6*x^4 + 9*x^2 + 23

        TESTS:

        Check that :trac:`17664` is fixed::

            sage: L.<c> = QuadraticField(-1)
            sage: P = L.primes_above(5)[0]
            sage: G = L.galois_group()
            sage: H = G.decomposition_group(P)
            sage: H.domain()
            {1, 2}
            sage: G.artin_symbol(P)
            ()
        """
        # XXX This should be fixed so that this can use GaloisGroup_v2.__init__
        PermutationGroup_generic.__init__(self, elts, canonicalize=True,
                                          domain=ambient.domain())
        self._ambient = ambient
        self._number_field = ambient.number_field()
        self._galois_closure = ambient._galois_closure
        self._pari_data = ambient._pari_data
        self._gc_map = ambient._gc_map
        self._elts = sorted(self.iteration())

    def fixed_field(self):
        r"""
        Return the fixed field of this subgroup (as a subfield of the Galois
        closure of the number field associated to the ambient Galois group).

        EXAMPLES::

            sage: L.<a> = NumberField(x^4 + 1)
            sage: G = L.galois_group()
            sage: H = G.decomposition_group(L.primes_above(3)[0])
            sage: H.fixed_field()
            (Number Field in a0 with defining polynomial x^2 + 2 with a0 = a^3 + a,
             Ring morphism:
               From: Number Field in a0 with defining polynomial x^2 + 2 with a0 = a^3 + a
               To:   Number Field in a with defining polynomial x^4 + 1
               Defn: a0 |--> a^3 + a)

        """
        if self.order() == 1:
            return self._galois_closure  # work around a silly error

        vecs = [pari(g.domain()).Vecsmall() for g in self._elts]
        v = self._ambient._pari_data.galoisfixedfield(vecs)
        x = self._galois_closure(v[1])
        return self._galois_closure.subfield(x)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: G = NumberField(x^3 - x - 1, 'a').galois_closure('b').galois_group()
            sage: H = G.subgroup([ G(1), G([(1,2,3),(4,5,6)]), G([(1,3,2),(4,6,5)])])
            sage: H._repr_()
            'Subgroup [(), (1,2,3)(4,5,6), (1,3,2)(4,6,5)] of Galois group of Number Field in b with defining polynomial x^6 - 6*x^4 + 9*x^2 + 23'
        """
        return "Subgroup %s of %s" % (self._elts, self._ambient)


class GaloisGroupElement(PermutationGroupElement):
    r"""
    An element of a Galois group. This is stored as a permutation, but may also
    be made to act on elements of the field (generally returning elements of
    its Galois closure).

    EXAMPLES::

        sage: K.<w> = QuadraticField(-7); G = K.galois_group()
        sage: G[1]
        (1,2)
        sage: G[1](w + 2)
        -w + 2

        sage: L.<v> = NumberField(x^3 - 2); G = L.galois_group(names='y')
        sage: G[4]
        (1,5)(2,4)(3,6)
        sage: G[4](v)
        1/18*y^4
        sage: G[4](G[4](v))
        -1/36*y^4 - 1/2*y
        sage: G[4](G[4](G[4](v)))
        1/18*y^4
    """
    @cached_method
    def as_hom(self):
        r"""
        Return the homomorphism L -> L corresponding to self, where L is the
        Galois closure of the ambient number field.

        EXAMPLES::

            sage: G = QuadraticField(-7,'w').galois_group()
            sage: G[1].as_hom()
            Ring endomorphism of Number Field in w with defining polynomial x^2 + 7 with w = 2.645751311064591?*I
              Defn: w |--> -w

        TESTS:

        Number fields defined by non-monic and non-integral
        polynomials are supported (:trac:`252`)::

            sage: R.<x> = QQ[]
            sage: f = 7/9*x^3 + 7/3*x^2 - 56*x + 123
            sage: K.<a> = NumberField(f)
            sage: G = K.galois_group()
            sage: G[1].as_hom()
            Ring endomorphism of Number Field in a with defining polynomial 7/9*x^3 + 7/3*x^2 - 56*x + 123
              Defn: a |--> -7/15*a^2 - 18/5*a + 96/5
            sage: prod(x - sigma(a) for sigma in G) == f.monic()
            True
        """
        G = self.parent()
        L = G.splitting_field()
        # First compute the image of the standard generator of the
        # PARI number field.
        a = G._pari_data.galoispermtopol(pari(self.domain()).Vecsmall())
        # Now convert this to a conjugate of the standard generator of
        # the Sage number field.
        P = L._pari_absolute_structure()[1].lift()
        a = L(P(a.Mod(L.pari_polynomial('y'))))
        return L.hom(a, L)

    def __call__(self, x):
        r"""
        Return the action of self on an element x in the number field of self
        (or its Galois closure).

        EXAMPLES::

            sage: K.<w> = QuadraticField(-7)
            sage: f = K.galois_group()[1]
            sage: f(w)
            -w
        """
        if x.parent() == self.parent().splitting_field():
            return self.as_hom()(x)
        else:
            return self.as_hom()(self.parent()._gc_map(x))

    def ramification_degree(self, P):
        """
        Return the greatest value of v such that s acts trivially modulo P^v.
        Should only be used if P is prime and s is in the decomposition group of P.

        EXAMPLES::

            sage: K.<b> = NumberField(x^3 - 3, 'a').galois_closure()
            sage: G = K.galois_group()
            sage: P = K.primes_above(3)[0]
            sage: s = hom(K, K, 1/18*b^4 - 1/2*b)
            sage: G(s).ramification_degree(P)
            4
        """
        if not self.parent().is_galois():
            raise TypeError("Ramification degree only defined for Galois extensions")
        gens = self.parent().number_field().ring_of_integers().ring_generators()
        w = [(self(g) - g).valuation(P) for g in gens]
        return min(w)

GaloisGroup_v2.Element = GaloisGroupElement

# For unpickling purposes we rebind GaloisGroup as GaloisGroup_v1.

GaloisGroup = GaloisGroup_v1
