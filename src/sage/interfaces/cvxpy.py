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
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.geometry.cone_catalog import nonnegative_orthant
from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices
from sage.manifolds.manifold import Manifold
from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
from sage.functions.log import log

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

        sage: a = cp.Variable(); a
        Variable(())
        sage: s_a = a._sage_(); s_a
        var0
        sage: s_a.parent()
        Ring of chart functions on Chart (dom_var0, (var0,))

    The result is cached::

        sage: s_a is a._sage_()
        True

    A vector variable with shape ``(5,)``::

        sage: x = cp.Variable(5, name='x'); x
        Variable((5,))
        sage: x._sage_()
        Coordinate functions (x_0, x_1, x_2, x_3, x_4)
         on the Chart (dom_x, (x_0, x_1, x_2, x_3, x_4))
        sage: x._sage_().parent()
        Ambient free module of rank 5 over
         Ring of chart functions
          on Chart (dom_x, (x_0, x_1, x_2, x_3, x_4))

    A matrix variable with shape ``(5, 1)``::

        sage: y = cp.Variable((5, 1), name='y'); y
        Variable((5, 1))
        sage: y._sage_()
        [y_0_0]
        [y_1_0]
        [y_2_0]
        [y_3_0]
        [y_4_0]
        sage: y._sage_().parent()
        Full MatrixSpace of 5 by 1 dense matrices over
         Ring of chart functions
          on Chart (dom_y, (y_0_0, y_1_0, y_2_0, y_3_0, y_4_0))

    A positive semidefinite matrix variable::

        sage: X = cp.Variable((3, 3), PSD=True, name='X'); X
        Variable((3, 3), PSD=True)
        sage: X._sage_()
        [X_0_0 X_0_1 X_0_2]
        [X_1_0 X_1_1 X_1_2]
        [X_2_0 X_2_1 X_2_2]
        sage: X._sage_().parent()
        Full MatrixSpace of 3 by 3 dense matrices over
         Ring of chart functions
          on Chart (amb_dom_X, (X_0_0, X_0_1, X_0_2,
                                X_1_0, X_1_1, X_1_2,
                                X_2_0, X_2_1, X_2_2))

    A 10-vector constrained to have boolean valued entries::

        sage: b = cp.Variable(10, boolean=True, name='b'); b
        Variable((10,), boolean=True)
        sage: b._sage_()
        Coordinate functions (b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)
         on the Chart (amb_dom_b, (b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9))
        sage: b._sage_().parent()
        Ambient free module of rank 10
         over Ring of chart functions
          on Chart (amb_dom_b, (b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9))

    A 3 by 4 matrix constrained to have integer valued entries::

        sage: Z = cp.Variable((3, 4), integer=True, name='Z'); Z
        Variable((3, 4), integer=True)
        sage: Z._sage_().parent()
        Full MatrixSpace of 3 by 4 dense matrices
         over Ring of chart functions
          on Chart (amb_dom_Z, (Z_0_0, Z_0_1, Z_0_2, Z_0_3,
                                Z_1_0, Z_1_1, Z_1_2, Z_1_3,
                                Z_2_0, Z_2_1, Z_2_2, Z_2_3))
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
    ambient_chart = ambient_space.default_chart()

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
            else:
                raise NotImplementedError

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

        actual_dom = ManifoldSubsetPullback(ambient_chart, None,
                                            actual_model_dom, name=actual_dom_name)

    if self.ndim == 0:
        sage_object = ambient_chart.function(ambient_chart[0])
    elif self.ndim == 1:
        # like a multifunction, but better for matrix-vector multiplication
        sage_object = vector(ambient_chart.function(variable)
                             for variable in ambient_chart)
    elif self.ndim == 2:
        sage_object = matrix(ambient_chart.function_ring(),
                             self.shape[0], self.shape[1],
                             [ambient_chart.function(variable)
                              for variable in ambient_chart])

    self._sage_actual_dom = actual_dom
    self._sage_object = sage_object
    return self._sage_object


def _cvxpy_Constant_sage_(self):
    r"""
    Return an equivalent Sage object.

    A CVXPY ``Constant`` becomes a constant expression, vector, or matrix.

    Note that all constant values are inexact.  If exact constants are needed, use
    parameters instead and substitute the exact value of the constant.

    EXAMPLES::

        sage: from sage.interfaces.cvxpy import cvxpy_init
        sage: cvxpy_init()
        sage: import cvxpy as cp

    A scalar constant::

        sage: print("possible deprecation warning"); a = cp.Constant(28/5); a
        possible deprecation warning...
        Constant(CONSTANT, NONNEGATIVE, ())
        sage: s_a = a._sage_(); s_a
        5.60000000000000
        sage: s_a.parent()
        Real Field with 53 bits of precision

    A vector constant with shape ``(5,)``::

        sage: x = cp.Constant([1/7, 2/7, 3/7, 4/7, 5/7]); x
        Constant(CONSTANT, NONNEGATIVE, (5,))
        sage: x._sage_()
        (0.142857142857143, 0.285714285714286, 0.428571428571429, 0.571428571428571, 0.714285714285714)
        sage: x._sage_().parent()
        Vector space of dimension 5 over Real Field with 53 bits of precision

    A matrix variable with shape ``(2, 3)``::

        sage: A = cp.Constant([[1, 2, 3], [4, 5, 6]]); A
        Constant(CONSTANT, NONNEGATIVE, (3, 2))
        sage: A._sage_()
        [1.00000000000000 4.00000000000000]
        [2.00000000000000 5.00000000000000]
        [3.00000000000000 6.00000000000000]
        sage: A._sage_().parent()
        Full MatrixSpace of 3 by 2 dense matrices over Real Field with 53 bits of precision
    """
    try:
        return self._sage_object
    except AttributeError:
        pass

    is_complex = any(self.attributes[key]
                     for key in ('complex', 'hermitian', 'imag'))
    if is_complex:
        raise NotImplementedError

    if self.ndim == 0:
        self._sage_object = RR(self.value)
    elif self.ndim == 1:
        self._sage_object = vector(RR, self.value)
    elif self.ndim == 2:
        self._sage_object = matrix(RR, self.value)
    else:
        # As of CVXPY 1.1.13, higher shapes are not implemented
        raise NotImplementedError

    return self._sage_object


# Binary operators

def _cvxpy_MulExpression_sage_(self):
    r"""
    Return an equivalent Sage object.

    EXAMPLES::

        sage: from sage.interfaces.cvxpy import cvxpy_init
        sage: cvxpy_init()
        sage: import cvxpy as cp

    A scalar variable::

        sage: a = cp.Variable(name='a'); a
        Variable(())
        sage: s_a = a._sage_(); s_a
        a
        sage: a2 = a * a; a2
        Expression(UNKNOWN, UNKNOWN, ())
        sage: s_a2 = a2._sage_(); s_a2
        a^2
        sage: s_a2 is s_a2._sage_()
        True
        sage: s_a2 == s_a * s_a
        True

    A vector variable with shape ``(5,)``::

        sage: x = cp.Variable(5, name='x'); x
        Variable((5,))
        sage: s_x = x._sage_(); s_x
        Coordinate functions (x_0, x_1, x_2, x_3, x_4) on the Chart (dom_x, (x_0, x_1, x_2, x_3, x_4))
        sage: x_dot_x = x @ x; x_dot_x
        Expression(UNKNOWN, UNKNOWN, ())
        sage: s_x_dot_x = x_dot_x._sage_(); s_x_dot_x
        x_0^2 + x_1^2 + x_2^2 + x_3^2 + x_4^2
        sage: s_x_dot_x == s_x * s_x
        True
    """
    try:
        return self._sage_object
    except AttributeError:
        pass

    left = self.args[0]._sage_()
    right = self.args[1]._sage_()

    self._sage_object = left * right
    return self._sage_object


# Affine

def _cvxpy_AddExpression_sage_(self):
    r"""
    Return an equivalent Sage object.

    EXAMPLES::

        sage: from sage.interfaces.cvxpy import cvxpy_init
        sage: cvxpy_init()
        sage: import cvxpy as cp

    A scalar variable::

        sage: a = cp.Variable(name='a'); a
        Variable(())
        sage: s_a = a._sage_(); s_a
        a
        sage: _2a = a + a; _2a
        Expression(AFFINE, UNKNOWN, ())
        sage: s_2a = _2a._sage_(); s_2a
        2*a
        sage: s_2a is s_2a._sage_()
        True
        sage: s_2a == s_a + s_a
        True
    """
    try:
        return self._sage_object
    except AttributeError:
        pass

    summands = [arg._sage_() for arg in self.args]

    self._sage_object = sum(summands)
    return self._sage_object


# Elementwise

def _cvxpy_log1p_sage_(self):
    r"""
    Return an equivalent Sage object.

    EXAMPLES::

        sage: from sage.interfaces.cvxpy import cvxpy_init
        sage: cvxpy_init()
        sage: import cvxpy as cp

    A vector variable with shape ``(5,)``::

        sage: x = cp.Variable(5, name='x'); x
        Variable((5,))
        sage: s_x = x._sage_(); s_x
        Coordinate functions (x_0, x_1, x_2, x_3, x_4)
         on the Chart (dom_x, (x_0, x_1, x_2, x_3, x_4))

        sage: log1p_x = cp.log1p(x); log1p_x
        Expression(CONCAVE, UNKNOWN, (5,))
        sage: s_log1p_x = log1p_x._sage_(); s_log1p_x
        Coordinate functions (log(x_0 + 1), log(x_1 + 1), log(x_2 + 1),
                              log(x_3 + 1), log(x_4 + 1))
         on the Chart (dom_x, (x_0, x_1, x_2, x_3, x_4))
    """
    try:
        return self._sage_object
    except AttributeError:
        pass

    arg = self.args[0]._sage_()

    def m(x):
        return log(x+1)

    try:
        self._sage_object = arg.apply_map(m)
    except AttributeError:
        self._sage_object = m(arg)
    return self._sage_object


# Other atoms

def _cvxpy_log_det_sage_(self):
    r"""
    Return an equivalent Sage object.

    EXAMPLES::

        sage: from sage.interfaces.cvxpy import cvxpy_init
        sage: cvxpy_init()
        sage: import cvxpy as cp

    A positive semidefinite matrix variable::

        sage: X = cp.Variable((3, 3), PSD=True, name='X'); X
        Variable((3, 3), PSD=True)
        sage: X._sage_()
        [X_0_0 X_0_1 X_0_2]
        [X_1_0 X_1_1 X_1_2]
        [X_2_0 X_2_1 X_2_2]

        sage: log_det_X = cp.log_det(X); log_det_X
        Expression(CONCAVE, NONNEGATIVE, ())
        sage: s_log_det_X = log_det_X._sage_(); s_log_det_X
        log(-(X_0_2*X_1_1 - X_0_1*X_1_2)*X_2_0
            + (X_0_2*X_1_0 - X_0_0*X_1_2)*X_2_1
            - (X_0_1*X_1_0 - X_0_0*X_1_1)*X_2_2)
    """
    try:
        return self._sage_object
    except AttributeError:
        pass

    arg = self.args[0]._sage_()

    self._sage_object = log(arg.det())
    return self._sage_object


# Monkey patching

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

    from cvxpy.expressions.constants import Constant
    Constant._sage_ = _cvxpy_Constant_sage_

    from cvxpy.atoms.affine.binary_operators import MulExpression
    MulExpression._sage_ = _cvxpy_MulExpression_sage_

    from cvxpy.atoms.affine.add_expr import AddExpression
    AddExpression._sage_ = _cvxpy_AddExpression_sage_

    from cvxpy.atoms.elementwise.log1p import log1p
    log1p._sage_ = _cvxpy_log1p_sage_

    from cvxpy.atoms.log_det import log_det
    log_det._sage_ = _cvxpy_log_det_sage_
