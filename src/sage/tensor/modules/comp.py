r"""
Components as indexed sets of ring elements

The class :class:`Components` is a technical class to take in charge the
storage and manipulation of **indexed elements of a commutative ring** that
represent the components of some "mathematical entity" with respect to some
"frame". Examples of *entity/frame* are *vector/vector-space basis* or
*vector field/vector frame on some manifold*. More generally, the components
can be those of a tensor on a free module or those of a tensor field on a
manifold. They can also be non-tensorial quantities, like connection
coefficients or structure coefficients of a vector frame.

The individual components are assumed to belong to a given commutative ring
and are labelled by *indices*, which are *tuples of integers*.
The following operations are implemented on components with respect
to a given frame:

* arithmetics (addition, subtraction, multiplication by a ring element)
* handling of symmetries or antisymmetries on the indices
* symmetrization and antisymmetrization
* tensor product
* contraction

Various subclasses of class :class:`Components` are

* :class:`CompWithSym` for components with symmetries or antisymmetries w.r.t.
  index permutations

  * :class:`CompFullySym` for fully symmetric components w.r.t. index
    permutations

    * :class:`KroneckerDelta` for the Kronecker delta symbol

  * :class:`CompFullyAntiSym` for fully antisymmetric components w.r.t. index
    permutations

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Joris Vankerschaver (2010): for the idea of storing only the non-zero
  components as dictionaries, whose keys are the component indices (implemented
  in the old class ``DifferentialForm``; see :trac:`24444`)
- Marco Mancini (2015) : parallelization of some computations
- Michael Jung (2021): refactor to parent model

EXAMPLES:

Set of components with 2 indices on a 3-dimensional vector space, the frame
being some basis of the vector space::

    sage: from sage.tensor.modules.comp import Components
    sage: V = VectorSpace(QQ,3)
    sage: basis = V.basis() ; basis
    [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]
    sage: c = Components(QQ, basis, 2) ; c
    2-index components w.r.t. [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]

Actually, the frame can be any object that has some length, i.e. on which
the function :func:`len()` can be called::

    sage: basis1 = V.gens() ; basis1
    ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    sage: c1 = Components(QQ, basis1, 2) ; c1
    2-index components w.r.t. ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    sage: basis2 = ['a', 'b' , 'c']
    sage: c2 = Components(QQ, basis2, 2) ; c2
    2-index components w.r.t. ['a', 'b', 'c']

A just created set of components is initialized to zero::

    sage: c.is_zero()
    True
    sage: c == 0
    True

This can also be checked on the list of components, which is returned by
the operator ``[:]``::

    sage: c[:]
    [0 0 0]
    [0 0 0]
    [0 0 0]

Individual components are accessed by providing their indices inside
square brackets::

    sage: c[1,2] = -3
    sage: c[:]
    [ 0  0  0]
    [ 0  0 -3]
    [ 0  0  0]
    sage: v = Components(QQ, basis, 1)
    sage: v[:]
    [0, 0, 0]
    sage: v[0]
    0
    sage: v[:] = (-1,3,2)
    sage: v[:]
    [-1, 3, 2]
    sage: v[0]
    -1

Sets of components with 2 indices can be converted into a matrix::

    sage: matrix(c)
    [ 0  0  0]
    [ 0  0 -3]
    [ 0  0  0]
    sage: matrix(c).parent()
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

By default, the indices range from `0` to `n-1`, where `n` is the length
of the frame. This can be changed via the argument ``start_index`` in
the :class:`Components` constructor::

    sage: v1 = Components(QQ, basis, 1, start_index=1)
    sage: v1[:]
    [0, 0, 0]
    sage: v1[0]
    Traceback (most recent call last):
    ...
    IndexError: index out of range: 0 not in range(1, 4)
    sage: v1[1]
    0
    sage: v1[:] = v[:]  # list copy of all components
    sage: v1[:]
    [-1, 3, 2]
    sage: v1[1], v1[2], v1[3]
    (-1, 3, 2)
    sage: v[0], v[1], v[2]
    (-1, 3, 2)

If some formatter function or unbound method is provided via the argument
``output_formatter`` in the :class:`Components` constructor, it is used to
change the output of the access operator ``[...]``::

    sage: a = Components(QQ, basis, 2, output_formatter=Rational.numerical_approx)
    sage: a[1,2] = 1/3
    sage: a[1,2]
    0.333333333333333

The format can be passed to the formatter as the last argument of the
access operator ``[...]``::

    sage: a[1,2,10] # here the format is 10, for 10 bits of precision
    0.33
    sage: a[1,2,100]
    0.33333333333333333333333333333

The raw (unformatted) components are then accessed by the double bracket
operator::

    sage: a[[1,2]]
    1/3

For sets of components declared without any output formatter, there is no
difference between ``[...]`` and ``[[...]]``::

    sage: c[1,2] = 1/3
    sage: c[1,2], c[[1,2]]
    (1/3, 1/3)

The formatter is also used for the complete list of components::

    sage: a[:]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    [0.000000000000000 0.000000000000000 0.333333333333333]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    sage: a[:,10] # with a format different from the default one (53 bits)
    [0.00 0.00 0.00]
    [0.00 0.00 0.33]
    [0.00 0.00 0.00]

The complete list of components in raw form can be recovered by the double
bracket operator, replacing ``:`` by ``slice(None)`` (since ``a[[:]]``
generates a Python syntax error)::

    sage: a[[slice(None)]]
    [  0   0   0]
    [  0   0 1/3]
    [  0   0   0]

Another example of formatter: the Python built-in function :func:`str`
to generate string outputs::

    sage: b = Components(QQ, V.basis(), 1, output_formatter=str)
    sage: b[:] = (1, 0, -4)
    sage: b[:]
    ['1', '0', '-4']

For such a formatter, 2-indices components are no longer displayed as a
matrix::

    sage: b = Components(QQ, basis, 2, output_formatter=str)
    sage: b[0,1] = 1/3
    sage: b[:]
    [['0', '1/3', '0'], ['0', '0', '0'], ['0', '0', '0']]

But unformatted outputs still are::

    sage: b[[slice(None)]]
    [  0 1/3   0]
    [  0   0   0]
    [  0   0   0]

Internally, the components are stored as a dictionary (:attr:`_comp`) whose
keys are the indices; only the non-zero components are stored::

    sage: a[:]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    [0.000000000000000 0.000000000000000 0.333333333333333]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    sage: a._comp
    {(1, 2): 1/3}
    sage: v[:] = (-1, 0, 3)
    sage: v._comp  # random output order of the component dictionary
    {(0,): -1, (2,): 3}

In case of symmetries, only non-redundant components are stored::

    sage: from sage.tensor.modules.comp import CompFullyAntiSym
    sage: c = CompFullyAntiSym(QQ, basis, 2)
    sage: c[0,1] = 3
    sage: c[:]
    [ 0  3  0]
    [-3  0  0]
    [ 0  0  0]
    sage: c._comp
    {(0, 1): 3}

"""

# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#       Copyright (C) 2020 Michael Jung <m.jung@vu.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.comp_parent import CompParent, CompParentWithSym, \
                                    CompParentFullySym, CompParentFullyAntiSym
from sage.tensor.modules.comp_element_dict import ComponentsKroneckerDelta_dict
from sage.misc.cachefunc import cached_function

def Components(ring, frame, nb_indices, start_index=0, output_formatter=None,
               sym=None, antisym=None):
    r"""

    """
    if not sym and not antisym:
        parent = CompParent(ring, nb_indices)
    else:
        parent = CompParentWithSym(ring, nb_indices,
                                   sym=sym, antisym=antisym)
    return parent(frame, start_index=start_index, output_formatter=output_formatter)

CompWithSym = Components

def CompFullySym(ring, frame, nb_indices, start_index=0, output_formatter=None):
    r"""

    """
    parent = CompParentFullySym(ring, nb_indices)
    return parent(frame, start_index=start_index, output_formatter=output_formatter)

def CompFullyAntiSym(ring, frame, nb_indices, start_index=0,
                     output_formatter=None):
    r"""

    """
    parent = CompParentFullyAntiSym(ring, nb_indices)
    return parent(frame, start_index=start_index, output_formatter=output_formatter)

@cached_function
def KroneckerDelta(ring, frame, start_index=0, output_formatter=None):
    r"""

    """
    parent = CompParentFullySym(ring, 2)
    return ComponentsKroneckerDelta_dict(parent, frame, start_index=start_index,
                                         output_formatter=output_formatter)
