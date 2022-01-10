r"""
Finitely presented graded modules

Let `R` be a connected graded algebra.  A finitely presented module over `R`
is isomorphic to the cokernel of an `R`-linear homomorphism `f: F_1 \to F_0`
of finitely generated free modules: The generators of `F_0` corresponds to the
generators of the module, and the generators of `F_1` corresponds to its
relations, via the map `f`.

The class constructor of this module class is given a set of generators and
relations, and uses them to construct a presentation, using the class
:class:`sage.modules.fp_graded.free_morphism.FreeGradedModuleMorphism`.

This package was designed with homological algebra in mind, and its API
focuses on maps rather than objects.  A good example of this is the kernel
function :meth:`sage.modules.fp_graded.morphism.FPModuleMorphism.kernel_inclusion`
which computes the kernel of a homomorphism `f: M\to N`.  Its return value
is not an instance of the module class, but rather an injective homomorphism
`i: K\to M` with the property that `\operatorname{im}(i) = \ker(f)`.

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Hom
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import PlusInfinity
from sage.categories.graded_modules import GradedModules
from sage.modules.module import Module
from sage.structure.indexed_generators import IndexedGenerators
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.element import parent

from sage.modules.fp_graded.free_module import FreeGradedModule
from sage.modules.fp_graded.free_element import FreeGradedModuleElement
from sage.modules.fp_graded.element import FPElement

# These are not free modules over the algebra, but they are free as
# vector spaces. They have a distinguished set of generators over the
# algebra, and as long as the algebra has a vector space basis
# implemented in Sage, the modules will have a vector space basis as well.
class FPModule(UniqueRepresentation, IndexedGenerators, Module):
    r"""
    Create a finitely presented module over a connected graded algebra.

    INPUT:

    - ``algebra`` -- the graded connected algebra over which the module is
      defined; this algebra must be equipped with a graded basis

    - ``generator_degrees`` -- tuple of integer degrees

    - ``relations`` -- tuple of relations; a relation is a tuple of
      coefficients `(c_1, \ldots, c_n)`, ordered so that they
      correspond to the module generators, that is, such a tuple
      corresponds to the relation

      .. MATH::

          c_1 g_1 + \ldots + c_n g_n = 0

      if the generators are `(g_1, \ldots, g_n)`

    EXAMPLES::

        sage: from sage.modules.fp_graded.module import FPModule

        sage: E.<x,y> = ExteriorAlgebra(QQ)
        sage: M = FPModule(E, [0, 1], [[x, 1]])
        sage: a, b = M.generators()
        sage: x*a + b == 0
        True
        sage: (x*a + b).normalize()
        0

        sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
        sage: M = FPModule(A3, [0, 1], [[Sq(2), Sq(1)]])
        sage: M.generators()
        (G[0], G[1])
        sage: M.relations()
        (Sq(2)*G[0] + Sq(1)*G[1],)
        sage: M.is_trivial()
        False

        sage: Z = FPModule(A3, [])
        sage: Z.generators()
        ()
        sage: Z.relations()
        ()
        sage: Z.is_trivial()
        True
    """
    @staticmethod
    def __classcall_private__(cls, algebra, generator_degrees, relations=(), names=None):
        r"""
        Normalize input to ensure a unique representation.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M1.<m0,m1> = FPModule(A3, [0, 1], [[Sq(2), Sq(1)]])
            sage: M2 = FPModule(A3, (0, 1), [[Sq(2), Sq(1)]], names='m0,m1')
            sage: M1 is M2
            True
        """
        if names is not None:
            from sage.structure.category_object import normalize_names
            names = normalize_names(-1, names)
        return super(FPModule, cls).__classcall__(cls,
            algebra=algebra,
            generator_degrees=tuple(generator_degrees),
            relations=tuple([tuple([algebra(x) for x in r]) for r in relations]),
            names=names)

    def __init__(self, algebra, generator_degrees, relations=(), names=None):
        r"""
        Create a finitely presented module over a connected graded algebra.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPModule(A3, [0, 1], [[Sq(2), Sq(1)]])
            sage: TestSuite(M).run()
        """
        self._relations = relations
        self._generator_degrees = generator_degrees

        # The free module on the generators of the module.
        generator_module = FreeGradedModule(algebra, generator_degrees, names=names)
        # Use the coefficients given for the relations and make module elements
        # from them.  Filter out the zero elements, as they are redundant.
        rels = [v for v in [generator_module(r) for r in relations] if not v.is_zero()]

        # The free module for the relations of the module.
        relations_module = FreeGradedModule(algebra,
                                            tuple([r.degree() for r in rels]))

        # The module we want to model is the cokernel of the
        # following morphism.
        self.j = Hom(relations_module, generator_module)(rels)

        # Call the base class constructors.
        keys = generator_module.basis().keys()
        cat = GradedModules(algebra).WithBasis().FinitelyPresented()
        IndexedGenerators.__init__(self, keys)
        Module.__init__(self, algebra, category=cat)

        from sage.combinat.family import Family
        self._spanning_set = Family(self._indices, self.monomial)

    Element = FPElement

    def _free_module(self):
        """
        Return the free module of which ``self`` is a quotient.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra()
            sage: M.<x,y> = FPModule(A, [0, 1], [[Sq(2), Sq(1)]])
            sage: M.generators()
            (x, y)
            sage: F = M._free_module()
            sage: F.generators()
            (x, y)
        """
        return self.j.codomain()


    @classmethod
    def from_free_module(cls, free_module):
        r"""
        Initialize from a finitely generated free module.

        INPUT:

        - ``free_module`` -- a finitely generated free module

        OUTPUT: the finitely presented module having same set of generators
        as ``free_module``, and no relations.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = FreeGradedModule(A, (-2,2,4))
            sage: FPModule.from_free_module(F)
            Finitely presented left module on 3 generators and 0 relations over
             mod 2 Steenrod algebra, milnor basis
        """
        return cls(algebra=free_module.base_ring(),
                   generator_degrees=free_module.generator_degrees(),
                   relations=())


    @classmethod
    def from_free_module_morphism(cls, morphism):
        r"""
        Create a finitely presented module from a morphism of finitely
        generated free modules.

        INPUT:

        - ``morphism`` -- a morphism between finitely generated free modules

        OUTPUT:

        The finitely presented module having presentation equal to the
        homomorphism ``morphism``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F1 = FreeGradedModule(A, (2,))
            sage: F2 = FreeGradedModule(A, (0,))
            sage: v = F2([Sq(2)])
            sage: pres = Hom(F1, F2)([v])
            sage: M = FPModule.from_free_module_morphism(pres); M
            Finitely presented left module on 1 generator and 1 relation over
             mod 2 Steenrod algebra, milnor basis
            sage: M.generator_degrees()
            (0,)
            sage: M.relations()
            (Sq(2)*G[0],)
        """
        return cls(algebra=morphism.base_ring(),
                   generator_degrees=morphism.codomain().generator_degrees(),
                   relations=tuple([r.coefficients() for r in morphism.values()]))


    def change_ring(self, algebra):
        r"""
        Change the base ring of ``self``.

        INPUT:

        - ``algebra`` -- a connected graded algebra

        OUTPUT:

        The finitely presented module over ``algebra`` defined with the
        exact same number of generators of the same degrees and relations
        as ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = M.change_ring(A2); N
            Finitely presented left module on 2 generators and 1 relation over
             sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]

        Changing back yields the original module::

            sage: N.change_ring(A) is M
            True
        """
        # self.relations() consists of module elements. We need to extra the coefficients.
        relations = tuple(r.coefficients() for r in self.relations())
        return FPModule(algebra, self.generator_degrees(), relations)


    def __contains__(self, x):  # Probably can be removed
        """
        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: x = M([Sq(1), 1])
            sage: x in M
            True
            sage: N = FPModule(SteenrodAlgebra(2), [0], [[Sq(2)]])
            sage: y = Sq(2)*N.generator(0)
            sage: y in M
            False
        """
        return parent(x) is self


    def _from_dict(self, d, coerce=False, remove_zeros=True):
        r"""
        Construct an element of ``self`` from an ``{index: coefficient}``
        dictionary.

        INPUT:

        - ``d`` -- a dictionary

        This code is taken from the method of the same name for
        ``sage.combinat.free_module.FreeModule``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: x = M.an_element()
            sage: M._from_dict(x.monomial_coefficients()) == x
            True
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            d = {key: R(coeff) for key, coeff in d.items()}
        if remove_zeros:
            d = {key: coeff for key, coeff in d.items() if coeff}
        return self.element_class(self, d)


    def _monomial(self, index):
        """
        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: M._monomial(0)
            G[0]
        """
        return self._from_dict({index: self.base_ring().one()}, remove_zeros=False)

    @lazy_attribute
    def monomial(self):
        """
        Return the basis element indexed by ``i``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: M.monomial(0)
            G[0]
            sage: M.monomial(1)
            G[1]
        """
        # Should use a real Map, as soon as combinatorial_classes are enumerated sets, and therefore parents
        from sage.categories.poor_man_map import PoorManMap
        return PoorManMap(self._monomial, domain=self._indices, codomain=self, name="Term map")


    @cached_method
    def zero(self):
        r"""
        Return the zero element.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: M = FPModule(SteenrodAlgebra(2), [0,1], [[Sq(4), Sq(3)]])
            sage: M.zero()
            0
        """
        return self.element_class(self, {})


    def _element_constructor_(self, x):
        r"""
        Construct any element of ``self``.

        INPUT:

        - ``x`` -- a tuple of coefficients, an element of FPModule, or the
          zero integer constant

        OUTPUT:

        An instance of the element class with coefficients from ``x``, the
        element ``x`` if it already was an element, or the zero element.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M.<a0,a2,a4> = FPModule(A, [0,2,4], [[Sq(4), Sq(2), 0]])

            sage: # Creating an element from coefficients:
            sage: e = M((Sq(6), 0, Sq(2))); e
            Sq(6)*a0 + Sq(2)*a4
            sage: e in M
            True

            sage: # Special syntax for creating the zero element:
            sage: z = M(0); z
            0
            sage: z.is_zero()
            True

            sage: # Creating an element from another element returns a reference to itself:
            sage: M(e)
            Sq(6)*a0 + Sq(2)*a4
            sage: e is M(e)
            True
        """
        if isinstance(x, self.element_class):
            return x
        if not self._generator_degrees: # the trivial module
            return self.zero()
        if not x:
            return self.zero()
        ngens = len(self._generator_degrees)
        B = self._spanning_set
        if isinstance(x, FreeGradedModuleElement):
            if x.parent() is self._free_module():
                # x.parent() should have the same generator list as self.
                coeffs = x.monomial_coefficients()
                return sum(coeffs[idx] * B[idx] for idx in coeffs)
            raise ValueError("element is not in this module")
        return self._from_dict({b: c for (c,b) in zip(x, self._indices) if c},
                               remove_zeros=False)


    def _repr_(self):
        r"""
        Construct a string representation of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,2,4], [[Sq(4),Sq(2),0]]); M
            Finitely presented left module on 3 generators and 1 relation over
             mod 2 Steenrod algebra, milnor basis
            sage: N = FPModule(A, [0,1], [[Sq(2),Sq(1)], [Sq(2)*Sq(1),Sq(2)]]); N
            Finitely presented left module on 2 generators and 2 relations over
             mod 2 Steenrod algebra, milnor basis
            sage: F = FPModule(A, [2]); F
            Finitely presented left module on 1 generator and 0 relations over
             mod 2 Steenrod algebra, milnor basis
        """
        return "Finitely presented left module on %s generator%s and %s relation%s over %s"\
            %(len(self._free_module().generator_degrees()),
              "" if len(self._free_module().generator_degrees()) == 1 else "s",
              len(self.j.values()), "" if len(self.j.values()) == 1 else "s",
              self.base_ring())


    def _repr_term(self, m):
        """
        Return a string representing the generator indexed by ``m``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,2,4], [[Sq(4),Sq(2),0]])
            sage: M._repr_term(4)
            'G[4]'
        """
        return self._free_module()._repr_term(m)


    def _latex_term(self, m):
        """
        Return a LaTeX representing the generator indexed by ``m``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,2,4], [[Sq(4),Sq(2),0]])
            sage: M._latex_term(4)
            'G_{4}'
        """
        return self._free_module()._latex_term(m)


    def connectivity(self):
        r"""
        The connectivity of ``self``.

        Since a finitely presented module over a connected algebra is in
        particular bounded below, the connectivity is an integer when the
        module is non-trivial, and `+\infty` when the module is trivial.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)

            sage: M = FPModule(A, [0,2,4], [[0, Sq(5), Sq(3)], [Sq(7), 0, Sq(2)*Sq(1)]])
            sage: M.connectivity()
            0

            sage: G = FPModule(A, [0,2], [[1,0]])
            sage: G.connectivity()
            2

        TESTS::

            sage: C = FPModule(SteenrodAlgebra(2, profile=(3,2,1)), [0], relations=[[Sq(1)], [0]])
            sage: C.connectivity()
            0

            sage: F = FPModule(A, [-1])
            sage: F.connectivity()
            -1

            sage: F = FPModule(A, [])
            sage: F.connectivity()
            +Infinity

            sage: F = FPModule(A, [0], [[1]])
            sage: F.connectivity()
            +Infinity
        """
        # In case there are no relations, the connectivity is the equal to
        # the connectivity of the free module on the generators.
        if self.j._degree is None:
            return self._free_module().connectivity()

        # We must check that the generator(s) in the free generator module are
        # not hit by relations, since we are not guaranteed that the
        # presentation we have is minimal.
        X = sorted(x for x in self.generator_degrees())

        previous = None
        for k in X:
            if previous is not None and k == previous:
                continue
            if not self.j.vector_presentation(k - self.j._degree).is_surjective():
                return k
            previous = k

        return PlusInfinity()


    def is_trivial(self):
        r"""
        Return ``True`` if ``self`` is isomorphic to the trivial module
        and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: M = FPModule(A2, [])
            sage: M.is_trivial()
            True

            sage: N = FPModule(A, [1,2])
            sage: N.is_trivial()
            False

            sage: P = FPModule(A, [1,2], [[1,0], [0,1]])
            sage: P.is_trivial()
            True

        TESTS::

            sage: C = FPModule(SteenrodAlgebra(2, profile=(3,2,1)), [0], [[Sq(1)], [0]])
            sage: C.is_trivial()
            False

            sage: C = FPModule(SteenrodAlgebra(2), [0], [[Sq(1)], [1]])
            sage: C.is_trivial()
            True
        """
        return self.connectivity() == PlusInfinity()


    def has_relations(self):
        r"""
        Return ``True`` if no relations are defined, and ``False``
        otherwise.

        .. NOTE::

            This module is free if this function returns ``False``, but a free
            module can have (redundant) relations.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: F = FPModule(A2, [1,2])
            sage: F.has_relations()
            False

            sage: M = FPModule(A2, [1,2], [[Sq(2), Sq(1)]])
            sage: M.has_relations()
            True

            sage: # A free module constructed with a redundant
            ....: # generator and relation.
            sage: N = FPModule(A2, [0,0], [[0, 1]])
            sage: N.has_relations()
            True
            sage: # Computing a minimal presentation reveals an
            ....: # isomorphic module with no relations.
            sage: N_min = N.min_presentation().domain()
            sage: N_min.has_relations()
            False
        """
        return not self.j.is_zero()


    def an_element(self, n=None):
        r"""
        An element of this module.

        This function chooses deterministically an element, i.e. the output
        depends only on the module and its input ``n``.

        INPUT:

        - ``n`` --  (optional) the degree of the element to construct

        OUTPUT:

        A module element of the given degree.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FPModule(A2, [0,2,4], [[0, Sq(5), Sq(3)], [Sq(7), 0, Sq(2)*Sq(1)]])

            sage: [M.an_element(i) for i in range(10)]
            [G[0],
             Sq(1)*G[0],
             Sq(2)*G[0] + G[2],
             Sq(0,1)*G[0] + Sq(1)*G[2],
             Sq(1,1)*G[0] + Sq(2)*G[2] + G[4],
             Sq(2,1)*G[0] + Sq(0,1)*G[2] + Sq(1)*G[4],
             Sq(0,2)*G[0] + Sq(1,1)*G[2] + Sq(2)*G[4],
             Sq(0,0,1)*G[0] + Sq(2,1)*G[2] + Sq(0,1)*G[4],
             Sq(1,0,1)*G[0] + Sq(6)*G[2] + Sq(1,1)*G[4],
             Sq(2,0,1)*G[0] + Sq(4,1)*G[2] + Sq(2,1)*G[4]]
        """
        a_free_element = self._free_module().an_element(n)
        return self(a_free_element)


    @cached_method
    def basis_elements(self, n, verbose=False):
        r"""
        A basis for the vector space of degree ``n`` module elements.
        Note that this returns a basis as a vector space over the
        base field.

        INPUT:

        - ``n`` -- an integer
        - ``verbose`` -- (default: ``False``) a boolean to control if log
          messages should be emitted

        OUTPUT:

        A list of homogeneous module elements of degree ``n`` which is
        a basis for the vector space of all degree ``n`` module elements.

        .. SEEALSO::

            :meth:`vector_presentation`, :meth:`element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M.<m0,m2> = FPModule(A2, [0,2], [[Sq(4), Sq(2)], [0, Sq(6)]])

            sage: M.basis_elements(4)
            (Sq(1,1)*m0, Sq(4)*m0)

            sage: M.basis_elements(5)
            (Sq(2,1)*m0, Sq(5)*m0, Sq(0,1)*m2)

            sage: M.basis_elements(25)
            ()

            sage: M.basis_elements(0)
            (m0,)

            sage: M.basis_elements(2)
            (Sq(2)*m0, m2)

        TESTS::

            sage: Z0 = FPModule(A2, [])
            sage: Z0.basis_elements(n=10)
            ()

            sage: Z1 = FPModule(A2, [1], [[1]])
            sage: Z1.basis_elements(n=10)
            ()
        """
        return tuple([self.element_from_coordinates(x, n) for
                      x in self.vector_presentation(n, verbose).basis()])


    @cached_method
    def element_from_coordinates(self, coordinates, n):
        r"""
        The module element in degree ``n`` having the given coordinates with
        respect to the basis returned by :meth:`basis_elements`.

        This function is inverse to
        :meth:`sage.modules.fp_graded.element.FPElement.vector_presentation`.

        INPUT:

        - ``coordinates`` -- a vector of coordinates
        - ``n`` -- the degree of the element to construct

        OUTPUT:

        A module element of degree ``n`` having the given coordinates
        with respect to the basis returned by :meth:`basis_elements`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0], [[Sq(4)], [Sq(7)], [Sq(4)*Sq(9)]])

            sage: M.vector_presentation(12).dimension()
            3
            sage: x = M.element_from_coordinates((0,1,0), 12); x
            Sq(0,4)*G[0]

        Applying the inverse function brings us back to the coordinate form::

            sage: x.vector_presentation()
            (0, 1, 0)

        TESTS::

            sage: M.element_from_coordinates((0,1,0,0), 12)
            Traceback (most recent call last):
             ...
            ValueError: the given coordinate vector has incorrect length (4);
             it should have length 3

        .. SEEALSO::

            :meth:`sage.modules.fp_graded.module.FPModule.vector_presentation`
        """
        M_n = self.vector_presentation(n)

        if len(coordinates) != M_n.dimension():
            raise ValueError('the given coordinate vector has incorrect length (%d); '
                             'it should have length %d' % (len(coordinates), M_n.dimension()))

        free_element = self._free_module().element_from_coordinates(
            M_n.lift(coordinates), n)

        return self(free_element.coefficients())


    def __getitem__(self, n):
        r"""
        A basis for the vector space of degree ``n`` module elements.

        INPUT:

        - ``n`` -- an integer

        OUTPUT:

        A list of homogeneous module elements of degree ``n`` which is
        a basis for the vector space of all degree ``n`` module elements.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,2,4], [[Sq(4),Sq(2),0]])

            sage: M[4]
            (Sq(1,1)*G[0], Sq(4)*G[0], G[4])

        .. SEEALSO::

            :meth:`basis_elements`
        """
        return self.basis_elements(n)


    @cached_method
    def vector_presentation(self, n, verbose=False):
        r"""
        A vector space isomorphic to the vector space of module elements of
        degree ``n``.

        INPUT:

        - ``n`` -- the degree of the presentation

        OUTPUT:

        A vector space.

        .. SEEALSO::

            :meth:`basis_elements`, :meth:`element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,2,4], [[Sq(4),Sq(2),0]])

            sage: V = M.vector_presentation(4)
            sage: V.dimension()
            3

            sage: len(M.basis_elements(4))
            3
        """
        # Get the vector space presentation of the free module on the
        # module generators.
        F_n = self._free_module().vector_presentation(n)

        # Compute the sub vector space generated by the relations.
        spanning_set = []

        if verbose:
            num_total_iterations = 0
            for relation in self.j.values():
                if relation.is_zero():
                    continue

                num_total_iterations += len(self.base_ring().basis(n - relation.degree()))

            progress = 0
            iteration_count = 0

        for relation in self.j.values():
            if relation.is_zero():
                continue

            for a in self.base_ring().basis(n - relation.degree()):
                if verbose:
                    iteration_count += 1
                    prog = int(100*iteration_count/num_total_iterations)
                    if prog > progress:
                        progress = prog
                        print('Progress: %d/100' % prog)

                # assert: isinstance(FreeElement, relation)
                v = (a*relation).vector_presentation()
                if v is not None:
                    spanning_set.append(v)

        R_n = F_n.subspace(spanning_set)

        # Return the quotient of the free part by the relations.
        return F_n / R_n


    def _Hom_(self, Y, category):
        r"""
        The internal hook used by the free function
        :meth:`sage.categories.homset.hom.Hom` to create homsets
        involving ``self``.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = FPModule(A, [1,3]);
            sage: L = FPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]]);

            sage: homset = Hom(F, L); homset
            Set of Morphisms from Finitely presented left module on 2 generators ...
        """
        from .homspace import FPModuleHomspace
        if not isinstance(Y, self.__class__):
            raise ValueError('cannot create homspace between incompatible types:\n%s  ->\n%s' % (self.__class__, type(Y)))
        if Y.base_ring() != self.base_ring():
            raise ValueError('the modules are not defined over the same base ring')

        return FPModuleHomspace(self, Y, category)


    def generator_degrees(self):
        r"""
        The degrees of the generators for ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A4 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: N = FPModule(A4, [0, 1], [[Sq(2), Sq(1)]])

            sage: N.generator_degrees()
            (0, 1)
        """
        return self._generator_degrees


    def generators(self):
        r"""
        The generators of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A4 = SteenrodAlgebra(2, profile=(4,3,2,1))

            sage: M = FPModule(A4, [0,0,2,3])
            sage: M.generators()
            (g_{0,0}, g_{0,1}, g_{2,0}, g_{3,0})

            sage: N = FPModule(A4, [0, 1], [[Sq(2), Sq(1)]], names='h')
            sage: N.generators()
            (h[0], h[1])

            sage: Z = FPModule(A4, [])
            sage: Z.generators()
            ()
        """
        return tuple([self.generator(i) for i in range(len(self.generator_degrees()))])

    gens = generators


    def generator(self, index):
        r"""
        Return the module generator with the given index.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A4 = SteenrodAlgebra(2, profile=(4,3,2,1))

            sage: M = FPModule(A4, [0,2,3])
            sage: M.generator(0)
            G[0]

            sage: N = FPModule(A4, [0, 1], [[Sq(2), Sq(1)]], names='h')
            sage: N.generator(1)
            h[1]
        """
        return self(self._free_module().generator(index))

    gen = generator


    def relations(self):
        r"""
        Return the relations of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A4 = SteenrodAlgebra(2, profile=(4,3,2,1))

            sage: M = FPModule(A4, [0,2,3])
            sage: M.relations()
            ()

            sage: N = FPModule(A4, [0, 1], [[Sq(2), Sq(1)]])
            sage: N.relations()
            (Sq(2)*G[0] + Sq(1)*G[1],)

            sage: Z = FPModule(A4, [])
            sage: Z.relations()
            ()
        """
        return self.j.values()


    def relation(self, index):
        r"""
        Return the module relation of the given index.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A4 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: N = FPModule(A4, [0, 1], [[Sq(2), Sq(1)]])
            sage: N.relation(0)
            Sq(2)*G[0] + Sq(1)*G[1]
        """
        return self.j.values()[index]


    def min_presentation(self, top_dim=None, verbose=False):
        r"""
        A minimal presentation of ``self``.

        OUTPUT:

        An isomorphism `M \to S`, where `M` has minimal presentation
        and `S` is ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: M = FPModule(A2, [0,1], [[Sq(2),Sq(1)],[0,Sq(2)],[Sq(3),0]])
            sage: i = M.min_presentation()
            sage: M_min = i.domain()

        ``i`` is an isomorphism between ``M_min`` and ``M``::

            sage: i.codomain() is M
            True
            sage: i.is_injective()
            True
            sage: i.is_surjective()
            True

        There are more relations in ``M`` than in ``M_min``::

            sage: M.relations()
            (Sq(2)*G[0] + Sq(1)*G[1], Sq(2)*G[1], Sq(3)*G[0])
            sage: M_min.relations()
            (Sq(2)*G[0] + Sq(1)*G[1], Sq(2)*G[1])

        TESTS::

            sage: T = FPModule(A2, [0], [[1]])
            sage: T_min = T.min_presentation().domain()
            sage: T_min.is_trivial()
            True
            sage: T_min
            Finitely presented left module on 0 generators and 0 relations over ...
        """
        return Hom(self, self).identity().image(top_dim, verbose)


    def suspension(self, t):
        r"""
        The suspension of ``self`` by degree ``t``.

        INPUT:

        - ``t`` -- an integer degree by which the module is suspended

        OUTPUT:

        A module which is identical to this module by a shift of
        degrees by the integer ``t``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: Y = FPModule(A2, [0], [[Sq(1)]])
            sage: X = Y.suspension(4)
            sage: X.generator_degrees()
            (4,)
            sage: X.relations()
            (Sq(1)*G[4],)

            sage: M = FPModule(A, [2,3], [[Sq(2), Sq(1)], [0, Sq(2)]])
            sage: Q = M.suspension(1)
            sage: Q.generator_degrees()
            (3, 4)
            sage: Q.relations()
            (Sq(2)*G[3] + Sq(1)*G[4], Sq(2)*G[4])
            sage: Q = M.suspension(-3)
            sage: Q.generator_degrees()
            (-1, 0)
            sage: Q = M.suspension(0)
            sage: Q.generator_degrees()
            (2, 3)
        """
        return FPModule(algebra=self.base_ring(),
                        generator_degrees=tuple([g + t for g in self._generator_degrees]),
                        relations=self._relations)


    def submodule_inclusion(self, spanning_elements):
        r"""
        Return the inclusion morphism of the submodule of ``self`` spanned
        by the given elements.

        INPUT:

        -  ``spanning_elements``  -- an iterable of elements

        OUTPUT:

        The inclusion of the submodule into this module.

        Because a submodule of a finitely presented module need not be
        finitely presented, this method will only work if the
        underlying algebra is finite-dimensional. Indeed, the current
        implementation only works if the algebra has a ``top_class``
        method, which gets used in
        :meth:`sage.modules.fp_graded.morphism._resolve_kernel`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: M = FPModule(A2, [0,1], [[Sq(2),Sq(1)]])
            sage: i = M.submodule_inclusion([M.generator(0)])
            sage: i.codomain() is M
            True
            sage: i.is_injective()
            True
            sage: i.domain().generator_degrees()
            (0,)
            sage: i.domain().relations()
            (Sq(3)*G[0],)
        """
        # Create the free graded module on the set of spanning elements.
        degs = [x.degree() for x in spanning_elements]
        F = FPModule(self.base_ring(), tuple(degs))

        # The submodule is the module generated by the spanning elements.
        return Hom(F, self)(spanning_elements).image()


    def resolution(self, k, top_dim=None, verbose=False):
        r"""
        A resolution of this module of length ``k``.

        INPUT:

        - ``k`` -- an non-negative integer
        - ``verbose`` -- (default: ``False``) a boolean to control if
          log messages should be emitted

        OUTPUT:

        A list of homomorphisms `[\epsilon, f_1, \ldots, f_k]` such that

        .. MATH::

            f_i: F_i \to F_{i-1} \text{ for } 1 < i \leq k,
            \qquad
            \epsilon: F_0 \to M,

        where each `F_i` is a finitely generated free module, and the
        sequence

        .. MATH::

            F_k \xrightarrow{\mathit{f_k}} F_{k-1}
            \xrightarrow{\mathit{f_{k-1}}} \ldots \rightarrow F_0
            \xrightarrow{\epsilon} M \rightarrow 0

        is exact.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FPModule(A2, [0,1], [[Sq(2), Sq(1)]])
            sage: M.resolution(0)
            [Module homomorphism of degree 0 defined by sending the generators
               [g_{0}, g_{1}]
             to
               (g_{0}, g_{1})]
            sage: res = M.resolution(4, verbose=True)
            Computing f_1 (1/4)
            Computing f_2 (2/4)
            Resolving the kernel in the range of dimensions [2, 25]:
             2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25.
            Computing f_3 (3/4)
            Resolving the kernel in the range of dimensions [8, 31]:
             8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31.
            Computing f_4 (4/4)
            Resolving the kernel in the range of dimensions [9, 33]:
             9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33.
            sage: len(res)
            5
            sage: res
            [Module homomorphism of degree 0 defined by sending the generators
               [G[0], G[1]]
             to
               (G[0], G[1]),
             Module homomorphism of degree 0 defined by sending the generators
               [g_{2}]
             to
               (Sq(2)*G[0] + Sq(1)*G[1],),
             Module homomorphism of degree 0 defined by sending the generators
               [g_{8}]
             to
               (Sq(3,1)*G[2],),
             Module homomorphism of degree 0 defined by sending the generators
               [g_{9}, g_{10}]
             to
               (Sq(1)*g_{8}, Sq(2)*g_{8}),
             Module homomorphism of degree 0 defined by sending the generators
               [g_{10}, g_{12}]
             to
               (Sq(1)*g_{9}, Sq(0,1)*g_{9} + Sq(2)*g_{10})]
            sage: for i in range(len(res)-1):
            ....:     assert (res[i]*res[i+1]).is_zero(), 'the result is not a complex'
        """
        def _print_progress(i, k):
            if verbose:
                print ('Computing f_%d (%d/%d)' % (i, i, k))

        if k < 0:
            raise ValueError('the length of the resolution must be non-negative')

        ret_complex = []

        # Epsilon: F_0 -> M
        F_0 = FPModule.from_free_module(self._free_module())
        epsilon = Hom(F_0, self)(self.generators())
        ret_complex.append(epsilon)

        if k == 0:
            return ret_complex

        # f_1: F_1 -> F_0
        _print_progress(1, k)
        F_1 = FPModule.from_free_module(self.j.domain())
        pres = Hom(F_1, F_0)(tuple([ F_0(x.coefficients()) for x in self.j.values() ]))

        ret_complex.append(pres)

        from .morphism import FPModuleMorphism

        # f_i: F_i -> F_i-1, for i > 1
        for i in range(2, k+1):
            _print_progress(i, k)

            f = ret_complex[i-1]
            ret_complex.append(FPModuleMorphism._resolve_kernel(f, top_dim=top_dim,
                                                                verbose=verbose))

        return ret_complex

