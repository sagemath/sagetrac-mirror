r"""
Homomorphisms of finitely generated free graded modules

This class implements construction and basic manipulation of
elements of the Sage parent
:class:`sage.modules.fp_over_steenrod_algebra.free_homspace.FreeModuleHomspace`, which models
homomorphisms of free graded modules over connected algebras.

.. NOTE:: This class is intended for private use by
    :class:`sage.modules.fp_over_steenrod_algebra.fp_morphism.FP_ModuleMorphism`.

For an overview of the free module API, see :doc:`free_module`.

AUTHORS:

    - Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
    - Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
      original software to Sage version 8.9.
    - Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
      new documentation and tests.

"""

#*****************************************************************************
#       Copyright (C) 2019 Robert R. Bruner <rrb@math.wayne.edu>
#                     and  Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from inspect import isfunction

from sage.categories.homset import Hom
from sage.categories.morphism import Morphism as SageMorphism
from sage.misc.cachefunc import cached_method

from .free_homspace import is_FreeModuleHomspace


class FreeModuleMorphism(SageMorphism):

    def __init__(self, parent, values):
        r"""
        Create a homomorphism between finitely generated free graded modules.

        INPUT::

        - ``parent`` -- A homspace in the category of finitely generated free
            modules.

        - ``values`` -- A list of elements in the codomain.  Each element
            corresponds (by their ordering) to a module generator in the domain.

        OUTPUT:: A module homomorphism defined by sending each generator to its
        corresponding value.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: A = SteenrodAlgebra(2)
            sage: F1 = FreeModule((4,5), A)
            sage: F2 = FreeModule((3,4), A)
            sage: F3 = FreeModule((2,3), A)
            sage: H1 = Hom(F1, F2)
            sage: H2 = Hom(F2, F3)
            sage: f = H1( ( F2((Sq(4), 0)), F2((0, Sq(4))) ) )
            sage: g = H2( ( F3((Sq(2), 0)), F3((Sq(3), Sq(2))) ) )
            sage: g*f
            Module homomorphism of degree 4 defined by sending the generators
              [<1, 0>, <0, 1>]
            to
              [<Sq(0,2) + Sq(3,1) + Sq(6), 0>, <Sq(1,2) + Sq(7), Sq(0,2) + Sq(3,1) + Sq(6)>]

        """

        if not is_FreeModuleHomspace(parent):
            raise TypeError('the parent (%s) must be a f.p. free module homset' % parent)

        # Get the values.
        C = parent.codomain()
        D = parent.domain()
        if isfunction(values):
            _values = [C(values(g)) for g in D.generators()]
        elif values == 0:
            _values = len(D.generator_degrees())*[C(0)]
        else:
            _values = [C(a) for a in values]

        # Check the homomorphism is well defined.
        if len(D.generator_degrees()) != len(_values):
            raise ValueError('the number of values must equal the number of '\
                'generators in the domain.  Invalid argument: %s' % _values)

        if all(v.is_zero() for v in _values):
            # The zero homomorphism does not get a degree.
            self._degree = None
        else:
            # Find the first non-zero value, and and record the shift
            # of degrees imposed by this homomorphism.
            for i, value in enumerate(_values):
                if not value.is_zero():
                    x = value.degree()
                    xx = D.generator_degrees()[i]
                    self._degree = x-xx
                    break

            # Check that all generators are shifted by the same degree.
            if not all(not v.degree() or self._degree == (v.degree() - g) \
                       for g, v in zip(D.generator_degrees(), _values)):
                errorMessage = "Ill defined homomorphism (degrees do not match)\n"
                gen_index = 0
                for g, v in zip(D.generator_degrees(), _values):
                    errorMessage += "  Generator #%d (degree %d) -> %s (degree %d)"\
                        " shifts degrees by %d\n" % (
                        gen_index, g, v, v.degree(), v.degree() - g)
                    gen_index += 1
                raise ValueError(errorMessage)

        self._values = _values

        SageMorphism.__init__(self, parent)


    def degree(self):
        r"""
        The degree of this homomorphism.

        OUTPUT:: The integer degree of this homomorphism, or ``None`` if this is
        the trivial homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeModule((0,1), A), FreeModule((0,), A))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.degree()
            5

        The zero homomorphism has no degree::

            sage: homspace.zero().degree() is None
            True

        """
        return self._degree


    def values(self):
        r"""
        The values under this homomorphism corresponding to the generators of
        the domain module.

        OUTPUT:: A sequence of elements of the codomain module.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeModule((0,1), A), FreeModule((2,), A))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.values()
            [<Sq(5)>, <Sq(3,1)>]
            sage: homspace.zero().values()
            [<0>, <0>]

        """
        return self._values


    def _richcmp_(self, other, op):
        r"""
        Compare this homomorphism to the given homomorphism.

        INPUT::

        - ``other`` -- An instance of this class.

        - ``op`` -- An integer specifying the comparison operation to be
          carried out: If ``op`` == 2, then return ``True`` if and only if the
          homomorphisms are equal.  If ``op`` == 3, then return ``True `` if
          and only if the homomorphisms are not equal.  Otherwise,
          return ``False``.

        OUTPUT:: A boolean.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeModule((0,1), A), FreeModule((2,), A))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f._richcmp_(f, op=2)
            True
            sage: f._richcmp_(f, op=3)
            False

        """

        try:
            same = (self - other).is_zero()
        except ValueError:
            return False

        # Equality
        if op == 2:
            return same

        # Non-equality
        if op == 3:
            return not same

        return False


    def __add__(self, g):
        r"""
        The pointwise sum of this and the given homomorphism.

        Pointwise addition of two homomorphisms `f` and `g` with the same domain
        and codomain is given by the formula `(f+g)(x) = f(x) + g(x)` for
        every `x` in the domain of `f`.

        INPUT::

        - ``g`` -- A homomorphism with the same domain and codomain as this
          homomorphism.

        OUTPUT:: The pointwise sum homomorphism of this and the given
        homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeModule((0,1), A), FreeModule((2,), A))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: ff = f.__add__(f)
            sage: ff.is_zero()
            True
            sage: ff.__add__(f) == f
            True

        """

        if self.domain() != g.domain():
            raise ValueError('morphisms do not have the same domain')
        elif self.codomain() != g.codomain():
            raise ValueError('morphisms do not have the same codomain')
        elif self._degree and g.degree() and self._degree != g.degree():
            raise ValueError('morphisms do not have the same degree')

        v = [self(x) + g(x) for x in self.domain().generators()]

        return self.parent()(v)


    def __neg__(self):
        r"""
        The additive inverse of this homomorphism with respect to the group
        structure given by pointwise sum.

        OUTPUT:: An instance of this class.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeModule((0,1), A), FreeModule((2,), A))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f_inverse = f.__neg__(); f_inverse
            Module homomorphism of degree 7 defined by sending the generators
              [<1, 0>, <0, 1>]
            to
              [<Sq(5)>, <Sq(3,1)>]
            sage: (f + f_inverse).is_zero()
            True

        """

        return self.parent()([-x for x in self._values])


    def __sub__(self, g):
        r"""
        The pointwise difference between this and the given homomorphism.

        OUTPUT:: An instance of this class.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeModule((0,1), A), FreeModule((2,), A))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: values2 = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: g = homspace(values2)
            sage: f.__sub__(g)
            The trivial homomorphism.

        """

        return self.__add__(g.__neg__())


    def __mul__(self, g):
        r"""
        The composition of the given homomorphism ``g``, followed by this
        homomorphisms.

        OUTPUT:: A homomorphism from the domain of this homomorphism, into the
        codomain of the homomorphism ``g``.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: N = FreeModule((2,), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: values2 = [Sq(2)*M.generator(0)]
            sage: g = Hom(N, M)(values2)
            sage: fg = f.__mul__(g); fg
            Module homomorphism of degree 7 defined by sending the generators
              [<1>]
            to
              [<Sq(4,1) + Sq(7)>]
            sage: fg.is_endomorphism()
            True

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.__mul__(f)
            Traceback (most recent call last):
            ...
            ValueError: morphisms are not composable

        """

        if self.parent().domain() != g.parent().codomain():
            raise ValueError('morphisms are not composable')
        homset = Hom(g.parent().domain(), self.parent().codomain())
        return homset([self(g(x)) for x in g.domain().generators()])


    def is_zero(self):
        r"""
        Decide if this homomomorphism is trivial.

        OUTPUT:: The boolean value ``True`` if this homomorphism is trivial, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: N = FreeModule((2,), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.is_zero()
            False
            sage: (f-f).is_zero()
            True

        """
        return self._degree == None


    def is_identity(self):
        r"""
        Decide if this homomomorphism is the identity endomorphism.

        OUTPUT:: The boolean value ``True`` if this homomorphism is the
        identity, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: N = FreeModule((2,), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.is_identity()
            False
            sage: id = Hom(M, M)(M.generators()); id
            The identity homomorphism.
            sage: id.is_identity()
            True
        """

        if self.parent().is_endomorphism_set():
            return self.parent().identity() == self
        else:
            return False


    def __call__(self, x):
        r"""
        Evaluate the homomorphism at the given domain element ``x``.

        INPUT::

        - ``x`` -- An element of the domain of this morphism.

        OUTPUT:: The module element of the codomain which is the value of ``x``
        under this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: N = FreeModule((2,), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.__call__(M.generator(0))
            <Sq(5)>
            sage: f.__call__(M.generator(1))
            <Sq(3,1)>

        """

        if x.parent() != self.domain():
            raise ValueError('cannot evaluate morphism on element not in the domain')

        value = sum([c*v for c, v in zip(
            x.coefficients(), self._values)], self.codomain()(0))

        return value


    def _repr_(self):
        r"""
        A string representation of this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: N = FreeModule((2,), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]

            sage: Hom(M, N)(values)._repr_()
            'Module homomorphism of degree 7 defined by sending the generators\n  [<1, 0>, <0, 1>]\nto\n  [<Sq(5)>, <Sq(3,1)>]'

            sage: Hom(M, N).zero()._repr_()
            'The trivial homomorphism.'

            sage: Hom(M, M).identity()._repr_()
            'The identity homomorphism.'

        """

        if self.is_zero():
            return "The trivial homomorphism."
        elif self.is_identity():
            return "The identity homomorphism."
        else:
            r = "Module homomorphism of degree {} defined by sending the generators\n  {}\nto\n  {}"
            return r.format(self._degree, self.domain().generators(), self._values)


    @cached_method
    def vector_presentation(self, n):
        r"""
        The restriction of this homomorphism to the domain module elements of
        degree ``n``.

        The restriction of a non-zero module homomorphism to the vectorspace of
        module elements of degree `n` is a linear function into the vectorspace
        of elements of degree `n+d` belonging to the codomain.  Here `d` is the
        degree of this homomorphism.

        When this homomorphism is zero, it has no well defined degree so the
        function cannot be presented since we do not know the degree of its
        codomain.  In this case, the return value is ``None``.

        INPUT::

        - ``n`` -- An integer degree.

        OUTPUT:: A linear function of finite dimensional vectorspaces over the
        ground field of the algebra for this module.  The domain is isomorphic
        to the vectorspace of domain elements of degree ``n`` of this free
        module, via the choice of basis given by
        :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.basis_elements`.
        If the morphism is zero, the value ``None`` is returned.

        .. SEEALSO::

            :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.vector_presentation`,
            :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.basis_elements`.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: N = FreeModule((2,), A)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.vector_presentation(0)
            Vector space morphism represented by the matrix:
            [0 1]
            Domain: Vector space of dimension 1 over Finite Field of size 2
            Codomain: Vector space of dimension 2 over Finite Field of size 2
            sage: f.vector_presentation(1)
            Vector space morphism represented by the matrix:
            [0 0 0]
            [0 1 0]
            Domain: Vector space of dimension 2 over Finite Field of size 2
            Codomain: Vector space of dimension 3 over Finite Field of size 2
            sage: f.vector_presentation(2)
            Vector space morphism represented by the matrix:
            [0 0 1 1]
            [0 0 0 0]
            Domain: Vector space of dimension 2 over Finite Field of size 2
            Codomain: Vector space of dimension 4 over Finite Field of size 2

        TESTS:

            sage: F = FreeModule((0,), A)
            sage: z = Hom(F, F)([0])
            sage: z.is_zero()
            True
            sage: z.vector_presentation(0) is None
            True

        """

        # The trivial map has no degree, so we can not create the codomain
        # of the linear transformation.
        if self._degree is None:
            return None

        D_n = self.domain().vector_presentation(n)
        C_n = self.codomain().vector_presentation(n + self._degree)

        values = [self(e) for e in self.domain().basis_elements(n)]
        return Hom(D_n, C_n)([
            C_n.zero() if e.is_zero() else e.vector_presentation() for e in values])






