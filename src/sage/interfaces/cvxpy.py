r"""
CVXPY --> Sage conversion
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.real_set import RealSet
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import minus_infinity, infinity
from sage.rings.real_mpfr import RR
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.geometry.cone_catalog import nonnegative_orthant
from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices
from sage.manifolds.manifold import Manifold
from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
from sage.manifolds.subsets.pullback import ManifoldSubsetPullback

# Leaf subclasses

def _cvxpy_Variable_sage_(self):
    r"""
    Return an equivalent Sage object.

    A CVXPY ``Variable`` is converted as follows:

    - the actual domain of the variable is converted to a Riemannian manifold
      (or a subset thereof) with a default chart corresponding to the entries
      of the variable;

    - return a chart function, vector of chart functions, or matrix of chart
      functions.

    EXAMPLES::

        sage: from sage.interfaces.cvxpy import cvxpy_init
        sage: cvxpy_init()
        sage: import cvxpy as cp

    A scalar variable::

        sage: a = cp.Variable()
        sage: s_a = a._sage_(); s_a
        1-dimensional Euclidean space dom_var0

    The result is cached::

        sage: s_a is a._sage_()
        True

    A vector variable with shape ``(5,)``::

        sage: x = cp.Variable(5, name='x')
        sage: x._sage_()
        5-dimensional Euclidean space dom_x

    A matrix variable with shape ``(5, 1)``::

        sage: y = cp.Variable((5, 1), name='y')
        sage: y._sage_()
        5-dimensional Euclidean space dom_y

    A positive semidefinite matrix variable::

        sage: X = cp.Variable((3, 3), PSD=True, name='X')
        sage: X._sage_()
        Subset dom_X of the 9-dimensional Euclidean space amb_dom_X

    A 10-vector constrained to have boolean valued entries::

        sage: b = cp.Variable(10, boolean=True, name='b')
        sage: b._sage_()
        Subset dom_b of the 10-dimensional Euclidean space amb_dom_b

    A 5 by 7 matrix constrained to have integer valued entries::

        sage: Z = cp.Variable((5, 7), integer=True, name='Z')
        sage: Z._sage_()
        Subset dom_Z of the 35-dimensional Euclidean space amb_dom_Z
    """
    try:
        return self._sage_object
    except AttributeError:
        pass

    is_complex = any(self.attributes[key]
                     for key in ('complex', 'hermitian', 'imag'))
    if is_complex:
        raise NotImplementedError

    name = self.name()
    if self.ndim == 0:
        names = [name]
        ambient_model_space = RealSet(minus_infinity, infinity)
    elif self.ndim == 1:
        names = [f'{name}_{i}'
                 for i in range(self.shape[0])]
        ambient_model_space = VectorSpace(RR, self.shape[0])
    elif self.ndim == 2:
        names = [f'{name}_{i}_{j}'
                 for i in range(self.shape[0]) for j in range(self.shape[1])]
        ambient_model_space = MatrixSpace(RR, *self.shape)
    else:
        # As of CVXPY 1.1.13, higher shapes are not implemented
        raise NotImplementedError

    # "A leaf may carry attributes that constrain the set values permissible for
    # it. Leafs can have no more than one attribute, with the exception that a
    # leaf may be both nonpos and nonneg or both boolean in some indices and
    # integer in others."
    is_ambient = not any(attribute
                         for attribute in self.attributes.values())
    is_integer = any(self.attributes[key]
                     for key in ('boolean', 'integer'))
    is_convex = not is_integer
    is_full_dimensional_convex = (is_convex
                                  and not any(self.attributes[key]
                                             for key in ('symmetric', 'diag', 'sparsity',
                                                         'PSD', 'NSD')))
    is_relatively_open_convex = (is_convex
                                 and not any(self.attributes[key]
                                             for key in ('nonneg', 'nonpos',
                                                         'PSD', 'NSD')))

    # In convex analysis, the word "domain" refers to the closure.
    # We use the term "actual domain" - which can be closed ('nonneg') or open ('pos')
    dom_name = f'dom_{name}'
    if is_ambient:
        ambient_name = dom_name
    else:
        ambient_name = f'amb_dom_{name}'
    if is_full_dimensional_convex:
        aff_name = ambient_name
    else:
        aff_name = f'aff_dom_{name}'
    if is_relatively_open_convex and not is_ambient:
        actual_dom_name = f'ri_dom_{name}'
    else:
        actual_dom_name = dom_name

    ambient_space = EuclideanSpace(name=ambient_name, names=names,
                                   start_index=0, unique_tag=self)
    if is_ambient:
        actual_dom = ambient_space
    else:
        if is_full_dimensional_convex or is_integer:
            aff_space = ambient_space
        else:
            if any(self.attributes[key] for key in ('diag', 'sparsity')):
                raise NotImplementedError
            elif any(self.attributes[key] for key in ('symmetric', 'PSD', 'NSD')):
                n = ZZ(self.shape[0])
                dimension = n * (n + 1) // 2
                aff_space = Manifold(dimension, name=aff_name, ambient=ambient_space)

        if any(self.attributes[key] for key in ('nonneg', 'pos', 'nonpos', 'neg')):
            actual_model_dom = cones.nonnegative_orthant(dim)
        elif any(self.attributes[key] for key in ('PSD', 'NSD')):
            actual_model_dom = PositiveSemidefiniteMatrices(ambient_model_space)
        else:
            actual_model_dom = ambient_model_space
        if any(self.attributes[key] for key in ('nonpos', 'neg', 'NSD')):
            actual_model_dom = -actual_model_dom
        if is_relatively_open_convex:
            actual_model_dom = actual_model_dom.relative_interior()

        actual_dom = ManifoldSubsetPullback(ambient_space.default_chart(), None,
                                            actual_model_dom, name=actual_dom_name)

    self._sage_object = actual_dom
    return self._sage_object

#

from sage.repl.ipython_extension import run_once

@run_once
def cvxpy_init():
    r"""
    Add ``_sage_()`` methods to CVXPY classes where needed.
    """

    from cvxpy.expressions.variable import Variable
    try:
        if Variable._sage_ == _cvxpy_Variable_sage_:
            return
    except AttributeError:
        pass
    Variable._sage_ = _cvxpy_Variable_sage_
