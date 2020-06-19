r"""
Lie Conformal Algebras.

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.

.. include:: ../../../lie_conformal_algebras/lie_conformal_algebra_desc.rst

EXAMPLES:

The base class for all Lie conformal algebras is
:class:`LieConformalAlgebra`.
All subclasses are called through its method ``__classcall_private__``.
We provide some convenience functions to define named Lie conformal
algebras like :meth:`VirasoroLieConformalAlgebra` and
:meth:`AffineLieConformalAlgebra`

- We construct the Virasoro Lie conformal algebra, its universal
  enveloping vertex algebra and lift some elements::

        sage: Vir = VirasoroLieConformalAlgebra(QQ)
        sage: Vir.inject_variables()
        Defining L, C
        sage: L.bracket(L)
        {0: TL, 1: 2*L, 3: 1/2*C}
        sage: cp = Family({C:1/2})
        sage: V = Vir.universal_enveloping_algebra(cp)
        sage: V
        The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
        sage: L.lift()
        L_-2|0>
        sage: L*L
        L_-2L_-2|0>
        sage: sorted(L.bracket(L*L).items())
        [(0, 2*L_-3L_-2|0> + L_-5|0>),
         (1, 4*L_-2L_-2|0>),
         (2, 3*L_-3|0>),
         (3, 17/2*L_-2|0>),
         (5, 3/2*|0>)]

- We construct the Current algebra for `\mathfrak{sl}_2`::

        sage: V = AffineLieConformalAlgebra(QQ, 'A1')
        sage: V.gens()
        (alpha[1], alphacheck[1], -alpha[1], K)
        sage: e = V.0; f = V.2; e.nproduct(f,1)
        K

- We construct the `\beta-\gamma` system by directly giving the
  `\lambda`-brackets of the generators::

        sage: betagamma_dict = {('b','a'):{0:{('K',0):1}}}
        sage: V = LieConformalAlgebra(QQ, betagamma_dict, names=('a','b'), weights=(1,0), central_elements=('K',))
        sage: V.category()
        Category of finitely generated H-graded Lie conformal algebras with basis over Rational Field
        sage: V.inject_variables()
        Defining a, b, K
        sage: a.bracket(b)
        {0: -K}
"""


#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.sets.family import Family
from sage.categories.commutative_rings import CommutativeRings
from sage.structure.parent import Parent

class LieConformalAlgebra(Parent, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, **kwds):

            sage: betagamma_dict = {('b','a'):{0:{('K',0):1}}}
            sage: V = LieConformalAlgebra(QQ, betagamma_dict, names=('a','b'), weights=(1,0), central_elements=('K',))
            sage: type(V)
            <class 'sage.algebras.lie_conformal_algebras.graded_lie_conformal_algebra.GradedLieConformalAlgebra_with_category'>
        """
        if not R in CommutativeRings():
            raise ValueError("arg0 must be a commutative ring got {}".format(R))

        #This is the only exposed class so we clean keywords here
        known_keywords = ['names', 'central_elements', 'index_set',
                          'weights','graded', 'category', 'prefix',
                          'bracket', 'latex_bracket', 'string_quotes',
                          'sorting_key','parity','super']
        for key in kwds:
            if key not in known_keywords:
                raise ValueError("LieConformalAlgebra(): got an unexpected " +
                                "keyword argument '%s'"%key)

        if isinstance(arg0,dict) and arg0:
            graded=kwds.get("graded", False)
            weights=kwds.get("weights",None)
            if weights is not None or graded:
                #if not arg0:
                    #return AbelianLieConformalAlgebra(R, names, index_set,
                    #    weights, graded=True)
                from .graded_lie_conformal_algebra import \
                                                    GradedLieConformalAlgebra
                return GradedLieConformalAlgebra(R, Family(arg0),**kwds)
            else:
                #if not arg0:
                #    return AbelianLieConformalAlgebra(R, names, index_set)
                from .lie_conformal_algebra_with_structure_coefs import \
                        LieConformalAlgebraWithStructureCoefficients
                return LieConformalAlgebraWithStructureCoefficients(R,
                                                    Family(arg0), **kwds)
        return NotImplementedError("Not implemented")

    def __init__(self, R, names=None, category=None):
        r"""
        Lie conformal algebras base class and factory.

        INPUT:

        - ``R`` -- a commutative ring (default: ``None``); the base
          ring of this Lie conformal algebra. Behaviour is undefined
          if it is not a field of characteristic zero.

        - ``arg0`` -- dictionary (default: ``None``);
          a dictionary containing the `\lambda` brackets of the
          generators of this Lie conformal algebra. The keys of this
          dictionary are pairs of either names or indices of the
          generators and the values are themselves dictionaries. For a
          pair of generators `a` and `b`, the value of
          ``arg0[('a','b')]`` is a dictionary whose keys are positive
          integer numbers and the corresponding value for the
          key `j` is a dictionary itself representing the j-th product
          `a_{(j)}b`. Thus, for a positive integer number `j`, the
          value of ``arg0[('a','b')][j]`` is a dictionary whose entries
          are pairs ``('c',n)`` where ``'c'`` is the name of a generator
          and `n` is a positive number. The value for this key is the
          coefficient of `\frac{T^{n}}{n!} c` in `a_{(j)}b`. For
          example the ``arg0`` for the *Virasoro* Lie conformal algebra
          is::

                {('L','L'):{0:{('L',1):1}, 1:{('L',0):2}, 3:{('C',0):1/2}}}


          Do not include central elements in this dictionary. Also, if
          the key ``('a','b')`` is present, there is no need to include
          ``('b','a')`` as it is defined by skew-symmetry. Any missing
          pair (besides the ones defined by skew-symmetry) is assumed
          to have vanishing `\lambda`-bracket.

        - ``names`` -- tuple of ``str`` (default: ``None``); the list of
          names for generators of this Lie conformal algebra. Do not
          include central elements in this list.

        - ``central_elements`` -- tuple of ``str`` (default: ``None``);
          A list of names for central elements of this Lie conformal
          algebra.

        - ``index_set`` -- enumerated set (default: ``None``); an
          indexing set for the generators of this Lie conformal algebra.
          Do not include central elements in this list.

        - ``weights`` -- tuple of non-negative rational numbers
           (default: ``None``); a list of degrees for this Lie
           conformal algebra.
           The returned Lie conformal algebra is H-Graded. This tuple
           needs to have the same cardinality as ``index_set`` or
           ``names``. Central elements are assumed to have weight `0`.

        - ``parity`` -- tuple of `0` or `1` (default: tuple of `0`);
          if this is a super Lie conformal algebra, this tuple
          specifies the parity of each of the non-central generators of
          this Lie conformal algebra. Central elements are assumed to
          be even. Notice that if this tuple is present, the category
          of this Lie conformal algebra is set to be a subcategory of
          ``LieConformalAlgebras(R).Super()``, even if all generators
          are even.

        - ``category`` The category that this Lie conformal algebra
          belongs to.

        In addition we accept the following keywords:

        - ``graded`` -- a boolean (default: ``False``);
          if `True`, the returned algebra is H-Graded.
          If ``weights`` is not specified, all non-central generators
          are assigned degree `1`. This keyword is unnecessary if
          ``weights`` is specified

        - ``super`` -- a boolean (default: ``False``);
          if `True`, the returned algebra is a super
          Lie conformal algebra. Even if all generators have even
          parity. If ``parity`` is not specified, all generators are
          assigned even parity. This keyword is unnecessary if
          ``parity`` is specified.

        EXAMPLES:

        We construct the `\beta-\gamma` system or *Weyl* Lie conformal
        algebra::

            sage: betagamma_dict = {('b','a'):{0:{('K',0):1}}}
            sage: V = LieConformalAlgebra(QQbar, betagamma_dict, names=('a','b'), weights=(1,0), central_elements=('K',))
            sage: V.category()
            Category of finitely generated H-graded Lie conformal algebras with basis over Algebraic Field
            sage: V.inject_variables()
            Defining a, b, K
            sage: a.bracket(b)
            {0: -K}

        We construct the current algebra for `\mathfrak{sl}_2`::

            sage: sl2dict = {('e','f'):{0:{('h',0):1}, 1:{('K',0):1}}, ('e','h'):{0:{('e',0):-2}}, ('f','h'):{0:{('f',0):2}}, ('h', 'h'):{1:{('K', 0):2}}}
            sage: V = LieConformalAlgebra(QQ, sl2dict, names=('e', 'h', 'f'), central_elements=('K',), graded=True)
            sage: V.inject_variables()
            Defining e, h, f, K
            sage: e.bracket(f)
            {0: h, 1: K}
            sage: h.bracket(e)
            {0: 2*e}
            sage: e.bracket(f.T())
            {0: Th, 1: h, 2: 2*K}
            sage: V.category()
            Category of finitely generated H-graded Lie conformal algebras with basis over Rational Field
            sage: e.degree()
            1

        .. SEEALSO::

            :mod:`sage.algebras.lie_conformal_algebras.graded_lie_conformal_algebra`

        """
        category = LieConformalAlgebras(R).or_subcategory(category)
        super(LieConformalAlgebra,self).__init__(base=R, names=names,
                                                 category=category)

    def _element_constructor_(self,x):
        """
        Construct an element of this Lie conformal algebra.

        EXAMPLES::

            sage: R = VirasoroLieConformalAlgebra(QQbar)
            sage: R(0)
            0
            sage: R(1)
            Traceback (most recent call last):
            ...
            ValueError: can only convert the scalar 0 into a Lie conformal algebra element
        """
        if x in self.base_ring():
            if not self.base_ring()(x).is_zero():
                raise ValueError("can only convert the scalar 0 into a"\
                                 " Lie conformal algebra element")
            return self.zero()
        try:
            x = self.module()(x)
        except (TypeError, ValueError):
            raise ValueError("Don't know how to convert {0} into an element"\
                             " of {1}".format(x,self))
        return self.element_class(self,x)

    def set_lift(self, liftmorphism):
        """
        Register ``liftmorphsm`` as a coercion between this Lie
        conformal algebra and its universal enveloping vertex algebra.
        """
        self.lift = liftmorphism
        try:
            self.lift.register_as_coercion()
        except AssertionError:
            #we already constructed this morphisms and its fine
            pass
