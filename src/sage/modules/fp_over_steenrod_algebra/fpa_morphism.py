r"""
Homomorphisms of finitely presented modules over the Steenrod algebra

This class implements construction and basic manipulation of homomorphisms
between finitely presented graded modules over the `\operatorname{mod} p`
Steenrod algebra.

For an overview of the API, see :doc:`fpa_module`.

AUTHORS:

    - Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
    - Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
      original software to Sage version 8.9.
    - Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
      new documentation and tests.

"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu>
#             and          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from sage.modules.fp_over_steenrod_algebra.fp_morphism import FP_ModuleMorphism
from sage.modules.fp_over_steenrod_algebra.profile import enveloping_profile_elements


class FPA_ModuleMorphism(FP_ModuleMorphism):
    r"""
    Create a homomorphism between finitely presented graded modules over
    the `\operatorname{mod} p` Steenrod algebra.

    INPUT::

    - ``parent`` -- A homspace object.

    - ``values`` -- A list of elements in the codomain.  Each element
      corresponds to a module generator in the domain.

    OUTPUT:: A module homomorphism defined by sending the generator with
    index `i` to the corresponding element in ``values``.

    .. NOTE:: Never use this constructor explicitly, but rather the parent's
        call method, or this class' __call__ method.  The reason for this
        is that the dynamic type of the element class changes as a
        consequence of the category system.

    TESTS:

        sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
        sage: # Trying to map the generators of a non-free module into a
        sage: # free module:
        sage: A = SteenrodAlgebra(2)
        sage: F = FPA_Module([2,3], A)
        sage: Q = FPA_Module([2,3], A, relations=[[Sq(6), Sq(5)]])
        sage: m = Hom(F, Q)( (F((Sq(1), 0)), F((0, 1))) )
        Traceback (most recent call last):
         ...
        ValueError: Ill defined homomorphism (degrees do not match)
              Generator #0 (degree 2) -> <Sq(1), 0> (degree 3) shifts degrees by 1
              Generator #1 (degree 3) -> <0, 1> (degree 3) shifts degrees by 0

        sage: # Trying to map the generators of a non-free module into a
        sage: # free module:
        sage: w = Hom(Q, F)( (F((1, 0)), F((0, 1))) )
        Traceback (most recent call last):
         ...
        ValueError: relation <Sq(6), Sq(5)> is not sent to zero

    """

    def __init__(self, parent, values):
        r"""
        Create a homomorphism between finitely presented graded modules over
        the `\operatorname{mod} p` Steenrod algebra.

        """
        # Call the base class constructor.
        FP_ModuleMorphism.__init__(self, parent, values)


    def profile(self):
        r"""
        A finite profile over which this homomorphism can be defined.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module([0,1], A, [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: id = Hom(M,M).identity()
            sage: id.profile()
            (2, 1)
            sage: zero = Hom(M,M).zero()
            sage: zero.profile()
            (2, 1)
            sage: A_fin = SteenrodAlgebra(2, profile=(2,1))
            sage: M_fin = M.change_ring(A_fin)

        Change the ring of the module M::

            sage: M_fin.change_ring(A) is M
            True

        We can change rings to the finite sub-Hopf algebra defined by
        the profile we just computed::

            sage: id_fin = id.change_ring(A_fin); id_fin
            The identity homomorphism.
            sage: id_fin.domain()
            Finitely presented module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]

        And if we change back to the full Steenrod algebra, we are back were
        we started::

            sage: id_fin.change_ring(A) == id
            True

        """
        def _flatten(f):
            return [coeffifient for value in f.values()\
                for coeffifient in value.coefficients()]

        elements = _flatten(self.domain().j) +\
            _flatten(self.codomain().j) +\
            _flatten(self)


        profile = enveloping_profile_elements(elements)

        # Avoid returning the zero profile because it triggers a corner case
        # in FP_Module.resolution().
        #
        # XXX: Fix FP_Module.resolution().
        #
        return (1,) if profile == (0,) else profile


    def is_injective(self, verbose=False):
        r"""
        Determine if this homomorphism is injective.

        INPUT::

        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT:: The boolean value ``True`` if this homomorphism has a trivial
        kernel, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module([0,1], A, [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = FPA_Module([0], A, [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: g = Hom(S, M).zero()
            sage: g.is_injective()
            False
            sage: z = Hom(FPA_Module([], A), M).zero()
            sage: z.is_injective()
            True
            sage: z.is_zero()
            True

        """
        algebra = self.base_ring()

        finite_algebra = algebra.__class__(algebra.prime(), profile=self.profile())

        return FP_ModuleMorphism.is_injective(
            self.change_ring(finite_algebra),
            verbose=verbose)


    def kernel(self, top_dim=None, verbose=False):
        r"""
        The kernel of this homomorphism.

        INPUT::

        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT:: An injective homomorphism into the domain ``self`` which is
        onto the kernel of this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module([0,1], A, [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = FPA_Module([0], A, [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: k = f.kernel(); k
            The trivial homomorphism.

        Since k is both trivial and injective, its domain should
        be the zero module::

            sage: k.domain().is_trivial()
            True

            sage: g = Hom(S, M)([M([Sq(3),Sq(2)])])
            sage: h = g.kernel(); h
            The identity homomorphism.
            sage: ker = h.domain();
            sage: ker is S
            True

        So `g` had to be trivial::

            sage: g.is_zero()
            True

        """
        return self._action(FP_ModuleMorphism.kernel, verbose)


    def image(self, verbose=False):
        r"""
        Compute the image of this homomorphism.

        INPUT::

        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT:: An injective homomorphism into the codomain of ``self`` which is
        onto the image of this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: M = FPA_Module([0,1], A, [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = FPA_Module([0], A, [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: i = f.image(); i
            Module homomorphism of degree 0 defined by sending the generators
              [<1>]
            to
              [<0, 1>]
            sage: i.codomain() is M
            True

        Lift the map ``f`` over the inclusion ``i``::

            sage: f_ = f.lift(i)
            sage: f_.is_injective()
            True
            sage: f_.is_surjective()
            True

            sage: g = Hom(S, M)([M([Sq(3),Sq(2)])])
            sage: j = g.image(); j
            The trivial homomorphism.

        So ``g`` had to be trivial::

            sage: g.is_zero()
            True

        """
        return self._action(FP_ModuleMorphism.image, verbose)


    def _resolve_kernel(self, top_dim=None, verbose=False):
        r"""
        Resolve the kernel of this homomorphism by a free module.

        INPUT::

        - ``verbose`` -- A boolean to enable progress messages. (optional,
          default: ``False``)

        OUTPUT:: A homomorphism `j: F \rightarrow D` where `D` is the domain of
        this homomorphism, `F` is free and such that `\ker(self) = \operatorname{im}(j)`.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: F = FPA_Module([0,0], A)
            sage: L = FPA_Module([0,0], A, [[Sq(3),Sq(0,1)], [0,Sq(2)]])
            sage: f = Hom(F, L)([L([Sq(2),0]), L([0, Sq(2)])])
            sage: f._resolve_kernel()
            Module homomorphism of degree 0 defined by sending the generators
              [<1, 0, 0>, <0, 1, 0>, <0, 0, 1>]
            to
              [<0, 1>, <Sq(0,1), 0>, <Sq(3), 0>]

        """
        return self._action(FP_ModuleMorphism._resolve_kernel, verbose)


    def _resolve_image(self, top_dim=None, verbose=False):
        r"""
        Resolve the image of this homomorphism by a free module.

        INPUT::

        - ``verbose`` -- A boolean to enable progress messages. (optional,
          default: ``False``)

        OUTPUT:: A homomorphism `j: F \rightarrow C` where `C` is the codomain
        of this homomorphism, `F` is free, and
        `\operatorname{im}(self) = \operatorname{im}(j)`.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: A = SteenrodAlgebra(2)
            sage: F = FPA_Module([0,0], A)
            sage: L = FPA_Module([0,0], A, [[Sq(3),Sq(0,1)], [0,Sq(2)]])
            sage: f = Hom(F, L)([L([Sq(2),0]), L([0, Sq(2)])])
            sage: f._resolve_image()
            Module homomorphism of degree 0 defined by sending the generators
              [<1>]
            to
              [<Sq(2), 0>]

        """
        return self._action(FP_ModuleMorphism._resolve_image, verbose)


    def _action(self, method, verbose=False):
        r"""
        Changes the ground ring to a finite algebra, acts by the given method
        and changes back into the original ground ring before returning.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: from sage.modules.fp_over_steenrod_algebra.fp_morphism import FP_ModuleMorphism
            sage: A = SteenrodAlgebra(2)
            sage: F = FPA_Module([0], A)
            sage: L = FPA_Module([0], A, [[Sq(3)]])
            sage: f = Hom(F, L)([L([Sq(2)])])
            sage: f._action(FP_ModuleMorphism._resolve_image, verbose=True)
            Computing the kernel using the profile:
            (2, 1)
            Resolving the image in the range of dimensions [0, 8]: 0 1 2 3 4 5 6 7 8.
            Module homomorphism of degree 0 defined by sending the generators
              [<1>]
            to
              [<Sq(2)>]

        """
        small_profile = self.profile()

        if verbose:
            print('Computing the kernel using the profile:')
            print(small_profile)

        algebra = self.base_ring()

        # Choose a finite sub Hopf-algebra of the original algebra.
        finite_algebra = algebra.__class__(algebra.prime(), profile=small_profile)

        # Perform the chosen action on the module after having changed rings
        # to the finite algebra.
        fp_result = method(
            self.change_ring(finite_algebra),
            verbose=verbose)

        # Change back to the original algebra and return the result.
        return fp_result.change_ring(self.base_ring())



