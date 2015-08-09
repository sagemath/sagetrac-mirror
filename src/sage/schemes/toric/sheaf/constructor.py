r"""
Construct sheaves on toric varieties.

A toric vector bundle (on a toric variety) is a vector bundle that is
equivariant with respect to the algebraic torus action.
"""


#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.toric.variety import is_ToricVariety
from sage.modules.filtered_vector_space import FilteredVectorSpace


def TangentBundle(X):
    r"""
    Construct the tangent bundle of a toric variety.

    INPUT:

    - ``X`` -- a toric variety. The base space of the bundle.

    OUTPUT:

    The tangent bundle as a Klyachko bundle.

    EXAMPLES::

        sage: dP7 = toric_varieties.dP7()
        sage: from sage.schemes.toric.sheaf.constructor import TangentBundle
        sage: TangentBundle(dP7)
        Rank 2 bundle on 2-d CPR-Fano toric variety covered by 5 affine patches.
    """
    if not is_ToricVariety(X):
        raise ValueError('not a toric variety')
    base_ring = X.base_ring()
    fan = X.fan()
    filtrations = dict()
    from sage.modules.filtered_vector_space import FilteredVectorSpace
    for i, ray in enumerate(fan.rays()):
        F = FilteredVectorSpace(fan.rays(), {0:range(fan.nrays()), 1:[i]})
        filtrations[ray] = F
    import klyachko
    return klyachko.Bundle(X, filtrations, check=True)


def CotangentBundle(X):
    r"""
    Construct the cotangent bundle of a toric variety.

    INPUT:

    - ``X`` -- a toric variety. The base space of the bundle.

    OUTPUT:

    The cotangent bundle as a Klyachko bundle.

    EXAMPLES::

        sage: dP7 = toric_varieties.dP7()
        sage: from sage.schemes.toric.sheaf.constructor import CotangentBundle
        sage: CotangentBundle(dP7)
        Rank 2 bundle on 2-d CPR-Fano toric variety covered by 5 affine patches.
    """
    return TangentBundle(X).dual()


def TrivialBundle(X, rank=1):
    r"""
    Return the trivial bundle of rank ``r``.

    INPUT:

    - ``X`` -- a toric variety. The base space of the bundle.

    - ``rank`` -- the rank of the bundle.

    OUTPUT:

    The trivial bundle as a Klyachko bundle.

    EXAMPLES::

        sage: P2 = toric_varieties.P2()
        sage: from sage.schemes.toric.sheaf.constructor import TrivialBundle
        sage: I3 = TrivialBundle(P2, 3);  I3
        Rank 3 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
        sage: I3.cohomology(weight=(0,0), dim=True)
        (3, 0, 0)
    """
    if not is_ToricVariety(X):
        raise ValueError('not a toric variety')
    from sage.modules.free_module import VectorSpace
    base_ring = X.base_ring()
    filtrations = dict([ray, FilteredVectorSpace(rank, 0, base_ring=base_ring)]
                       for ray in X.fan().rays())
    import klyachko
    return klyachko.Bundle(X, filtrations, check=True)


def LineBundle(X, D):
    """
    Construct the rank-1 bundle `O(D)`.

    INPUT:

    - ``X`` -- a toric variety. The base space of the bundle.

    - ``D`` -- a toric divisor.

    OUTPUT:

    The line bundle `O(D)` as a Klyachko bundle of rank 1.

    EXAMPLES::

        sage: X = toric_varieties.dP8()
        sage: D = X.divisor(0)
        sage: from sage.schemes.toric.sheaf.constructor import LineBundle
        sage: O_D = LineBundle(X, D)
        sage: O_D.cohomology(dim=True, weight=(0,0))
        (1, 0, 0)
    """
    if not is_ToricVariety(X):
        raise ValueError('not a toric variety')
    from sage.modules.free_module import VectorSpace
    base_ring = X.base_ring()
    filtrations = dict([X.fan().ray(i),
                        FilteredVectorSpace(1, D.function_value(i), base_ring=base_ring)]
                       for i in range(X.fan().nrays()))
    import klyachko
    return klyachko.Bundle(X, filtrations, check=True)



class SheafLibrary(object):

    def __init__(self, toric_variety):
        """
        Utility object to construct sheaves on toric varieties.

        .. warning::

            You should never constuct instances manually. Can be
            accessed from a toric variety via the
            :attr:`sage.schemes.toric.variety.ToricVariety_field.sheaves`
            attribute.

        EXAMPLES::

            sage: type(toric_varieties.P2().sheaves)
            <class 'sage.schemes.toric.sheaf.constructor.SheafLibrary'>
        """
        self._variety = toric_variety

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: toric_varieties.P2().sheaves      # indirect doctest
            Sheaf constructor on 2-d CPR-Fano toric variety covered by 3 affine patches
        """
        return 'Sheaf constructor on ' + repr(self._variety)

    def trivial_bundle(self, rank=1):
        r"""
        Return the trivial bundle of rank ``r``.

        INPUT:

        - ``rank`` -- integer (optional; default: `1`). The rank of
          the bundle.

        OUTPUT:

        The trivial bundle as a Klyachko bundle.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: I3 = P2.sheaves.trivial_bundle(3);  I3
            Rank 3 bundle on 2-d CPR-Fano toric variety covered by 3 affine patches.
            sage: I3.cohomology(weight=(0,0), dim=True)
            (3, 0, 0)
        """
        return TrivialBundle(self._variety, rank)

    def line_bundle(self, divisor):
        """
        Construct the rank-1 bundle `O(D)`.

        INPUT:

        - ``divisor`` -- a toric divisor.

        OUTPUT:

        The line bundle `O(D)` for the given divisor as a Klyachko
        bundle of rank 1.

        EXAMPLES::

            sage: X = toric_varieties.dP8()
            sage: D = X.divisor(0)
            sage: O_D = X.sheaves.line_bundle(D)
            sage: O_D.cohomology(dim=True, weight=(0,0))
            (1, 0, 0)
        """
        return LineBundle(self._variety, divisor)

    def tangent_bundle(self):
        r"""
        Return the tangent bundle of the toric variety.

        OUTPUT:

        The tangent bundle as a Klyachko bundle.

        EXAMPLES::

            sage: toric_varieties.dP6().sheaves.tangent_bundle()
            Rank 2 bundle on 2-d CPR-Fano toric variety covered by 6 affine patches.
        """
        return TangentBundle(self._variety)

    def cotangent_bundle(self):
        r"""
        Return the cotangent bundle of the toric variety.

        OUTPUT:

        The cotangent bundle as a Klyachko bundle.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: TX = dP6.sheaves.tangent_bundle()
            sage: TXdual = dP6.sheaves.cotangent_bundle()
            sage: TXdual == TX.dual()
            True
        """
        return CotangentBundle(self._variety)

    def Klyachko(self, multi_filtration):
        """
        Construct a Klyachko bundle (sheaf) from filtration data.

        INPUT:

        - ``multi_filtration`` -- a multi-filtered vectors space with
          multiple filtrations being indexed by the rays of the
          fan. Either an instance of
          :func:`~sage.modules.multi_filtered_vector_space.MultiFilteredVectorSpace`
          or something (like a dictionary of ordinary filtered vector
          spaces).

        OUTPUT:

        The Klyachko bundle defined by the filtrations, one for each
        ray, of a vector space.

        EXAMPLES::

            sage: P1 = toric_varieties.P1()
            sage: v1, v2, v3 = [(1,0,0),(0,1,0),(0,0,1)]
            sage: F1 = FilteredVectorSpace({1:[v1, v2, v3], 3:[v1]})
            sage: F2 = FilteredVectorSpace({0:[v1, v2, v3], 2:[v2, v3]})
            sage: P1 = toric_varieties.P1()
            sage: r1, r2 = P1.fan().rays()
            sage: F = MultiFilteredVectorSpace({r1:F1, r2:F2});  F
            Filtrations
                N(-1): QQ^3 >= QQ^2 >= QQ^2 >=  0   >= 0
                 N(1): QQ^3 >= QQ^3 >= QQ^1 >= QQ^1 >= 0
            sage: P1.sheaves.Klyachko(F)
            Rank 3 bundle on 1-d CPR-Fano toric variety covered by 2 affine patches.
        """
        from klyachko import Bundle
        return Bundle(self._variety, multi_filtration, check=True)

    def divisor(self, *args, **kwds):
        """
        Return a toric divisor.

        INPUT:

        This is just an alias for
        :meth:`sage.schemes.toric.variety.ToricVariety_field.divisor`,
        see there for details.

        By abuse of notation, you can usually use the divisor `D`
        interchangeably with the line bundle `O(D)`.

        OUTPUT:

        A toric divisor.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: dP6.inject_variables()
            Defining x, u, y, v, z, w
            sage: D = dP6.sheaves.divisor(x*u^3);  D
            V(x) + 3*V(u)
            sage: D == dP6.divisor(x*u^3)
            True
        """
        return self._variety.divisor(*args, **kwds)
