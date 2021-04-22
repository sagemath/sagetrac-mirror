"""
Category of chain complexes
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw
#                2009 Mike Hansen
#                2013 Volker Braun
#                2013, 2015 Travis Scrimshaw
#                2021 Michael Jung
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .category_types import Category_module
from .modules import Modules
from .functor import Functor
from .morphism import SetMorphism
from sage.misc.abstract_method import abstract_method

#############################################################
# ChainComplex
#############################################################
class ChainComplexes(Category_module):
    """
    The category of all chain complexes over a base ring.

    EXAMPLES::

        sage: ChainComplexes(RationalField())
        Category of chain complexes over Rational Field
        sage: ChainComplexes(Integers(9))
        Category of chain complexes over Ring of integers modulo 9

     TESTS::

        sage: TestSuite(ChainComplexes(RationalField())).run()

    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: ChainComplexes(Integers(9)).super_categories()
            [Category of modules over Ring of integers modulo 9]

        """
        from sage.categories.all import Fields, Modules, VectorSpaces
        base_ring = self.base_ring()
        if base_ring in Fields():
            return [VectorSpaces(base_ring)]
        return [Modules(base_ring)]

    class ParentMethods:
        @abstract_method(optional=True)
        def homology(self, n=None):
            r"""
            Return the homology of the chain complex.

            INPUT:

            - ``n`` -- (default: ``None``) degree of the homology; if none is
              provided, the direct sum homology will be used

            EXAMPLES::

                sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
                sage: C.homology(0)
                Z x Z
                sage: C.homology(1)
                Z x C3
                sage: C.homology(2)
                0

            ::

                sage: A.<x,y,z> = GradedCommutativeAlgebra(QQ, degrees=(2, 2, 3))
                sage: C = A.cdg_algebra({z: x*y})
                sage: C.homology(0)
                Free module generated by {[1]} over Rational Field
                sage: C.homology(1)
                Free module generated by {} over Rational Field
                sage: C.homology(2)
                Free module generated by {[x], [y]} over Rational Field
                sage: C.homology(3)
                Free module generated by {} over Rational Field
                sage: C.homology(4)
                Free module generated by {[x^2], [y^2]} over Rational Field

            """

        @abstract_method
        def differential(self, *args, **kwargs):
            r"""
            Return the differentials (or boundary maps) of the chain complex.

            EXAMPLES::

                sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
                sage: C.differential(0)
                [3 0 0]
                [0 0 0]

            ::

                sage: A.<x,y,z> = GradedCommutativeAlgebra(QQ, degrees=(2, 2, 3))
                sage: C = A.cdg_algebra({z: x*y})
                sage: C.differential()
                Differential of Commutative Differential Graded Algebra with generators ('x', 'y', 'z') in degrees (2, 2, 3) over Rational Field
                  Defn: x --> 0
                        y --> 0
                        z --> x*y

            """

        def lift_from_homology(self, x):
            r"""
            Lift the homology element ``x`` to the corresponding module.

            EXAMPLES::

                ...

            """
            try:
                # is homology already a quotient?
                x.lift()
            except AttributeError:
                # if not, this methods needs to be overwritten by parent
                raise NotImplementedError

        def reduce_to_homology(self, x, n=None):
            r"""
            Reduce a cycle to the corresponding quotient in homology.

            INPUT:

            - ``x`` -- a cycle
            - ``n`` -- (default: ``None``) degree of the homology; if none is
              provided, the direct sum homology will be used


            EXAMPLES::

                ...

            """
            try:
                # try coercion
                self.homology(n)(x)
            except TypeError:
                # if not, this methods needs to be overwritten by parent
                raise NotImplementedError

class HomologyFunctor(Functor):
    r"""
    Homology functor.

    INPUT:

    - ``domain`` -- must be a category of chain complexes
    - ``n`` -- (default: ``None``) degree of the homology; if none is provided,
      the direct sum homology will be used

    EXAMPLES::

        sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
        sage: H = HomologyFunctor(ChainComplexes(ZZ), 1)
        sage: H(C)
        Z x C3

    ::

        sage: A.<x,y,z> = GradedCommutativeAlgebra(QQ, degrees=(2, 2, 3))
        sage: C = A.cdg_algebra({z: x*y})
        sage: H = HomologyFunctor(ChainComplexes(QQ), 2)
        sage: H(C)
        Free module generated by {[x], [y]} over Rational Field

    Applying to a chain map::

        ...

    """
    def __init__(self, domain, n=None):
        r"""
        Construct the homology functor.

        TESTS::

            sage: H = HomologyFunctor(ChainComplexes(QQ), 1); H
            Functor from Category of chain complexes over Rational Field to
             Category of vector spaces over Rational Field

        """
        if not isinstance(domain, ChainComplexes):
            raise TypeError(f'{domain} must be a category of chain complexes')
        codomain = Modules(domain.base_ring())
        super().__init__(domain, codomain)
        self.__n = n

    def _apply_functor(self, x):
        r"""
        Apply ``self`` to a chain complex.

        TESTS::

            sage: C = ChainComplex({0: matrix(ZZ, 2, 3, [3, 0, 0, 0, 0, 0])})
            sage: H = HomologyFunctor(ChainComplexes(ZZ), 1)
            sage: H._apply_functor(C)
            Z x C3

        """
        return x.homology(self.__n)

    def _apply_functor_to_morphism(self, f):
        r"""
        Apply ``self`` to a chain map.

        TESTS:

        ...

        """
        lift = self.domain.lift_from_homology
        reduce = self.codomain.reduce_to_homology
        apply_f_star = lambda x: reduce(f(lift(x)), self.__n)
        return SetMorphism(self.domain.homology(self.__n),
                           self.codomain.homology(self.__n),
                           apply_f_star)
