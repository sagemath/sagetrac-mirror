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

* :class:`ComponentsWithSym` for components with symmetries or antisymmetries w.r.t.
  index permutations

  * :class:`ComponentsFullySym` for fully symmetric components w.r.t. index
    permutations

    * :class:`ComponentsKroneckerDelta` for the Kronecker delta symbol

  * :class:`ComponentsFullyAntiSym` for fully antisymmetric components w.r.t. index
    permutations

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Joris Vankerschaver (2010): for the idea of storing only the non-zero
  components as dictionaries, whose keys are the component indices (implemented
  in the old class ``DifferentialForm``; see :trac:`24444`)
- Marco Mancini (2015) : parallelization of some computations

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

For such a formatter, 2-index components are no longer displayed as a
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
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.element import Element
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.parallel.decorate import parallel
from sage.parallel.parallelism import Parallelism
from operator import itemgetter

class Components_base(Element):
    r"""
    Indexed set of ring elements forming some components with respect
    to a given "frame".

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities,
    such as connection coefficients or structure coefficients. The symmetries
    over some indices are dealt by subclasses of the class :class:`Components`.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have a method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of integer indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an element of ``ring`` and
      the second one, if any, some format specification.

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

    By default, the indices range from `0` to `n-1`, where `n` is the length
    of the frame. This can be changed via the argument ``start_index``::

        sage: c1 = Components(QQ, basis, 2, start_index=1)
        sage: c1[0,1]
        Traceback (most recent call last):
        ...
        IndexError: index out of range: 0 not in range(1, 4)
        sage: c[0,1]  # for c, the index 0 is OK
        0
        sage: c[0,1] = -3
        sage: c1[:] = c[:] # list copy of all components
        sage: c1[1,2]  # (1,2) = (0,1) shifted by 1
        -3

    If some formatter function or unbound method is provided via the argument
    ``output_formatter``, it is used to change the output of the access
    operator ``[...]``::

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

    For such a formatter, 2-index components are no longer displayed as a
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
        sage: v = Components(QQ, basis, 1)
        sage: v[:] = (-1, 0, 3)
        sage: v._comp  # random output order of the component dictionary
        {(0,): -1, (2,): 3}


    .. RUBRIC:: ARITHMETIC EXAMPLES:

    Unary plus operator::

        sage: a = Components(QQ, basis, 1)
        sage: a[:] = (-1, 0, 3)
        sage: s = +a ; s[:]
        [-1, 0, 3]
        sage: +a == a
        True

    Unary minus operator::

        sage: s = -a ; s[:]
        [1, 0, -3]

    Addition::

        sage: b = Components(QQ, basis, 1)
        sage: b[:] = (2, 1, 4)
        sage: s = a + b ; s[:]
        [1, 1, 7]
        sage: a + b == b + a
        True
        sage: a + (-a) == 0
        True

    Subtraction::

        sage: s = a - b ; s[:]
        [-3, -1, -1]
        sage: s + b == a
        True
        sage: a - b == - (b - a)
        True

    Multiplication by a scalar::

        sage: s = 2*a ; s[:]
        [-2, 0, 6]

    Division by a scalar::

        sage: s = a/2 ; s[:]
        [-1/2, 0, 3/2]
        sage: 2*(a/2) == a
        True

    Tensor product (by means of the operator ``*``)::

        sage: c = a*b ; c
        2-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:]
        ([-1, 0, 3], [2, 1, 4])
        sage: c[:]
        [-2 -1 -4]
        [ 0  0  0]
        [ 6  3 12]
        sage: d = c*a ; d
        3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: d[:]
        [[[2, 0, -6], [1, 0, -3], [4, 0, -12]],
         [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
         [[-6, 0, 18], [-3, 0, 9], [-12, 0, 36]]]
        sage: d[0,1,2] == a[0]*b[1]*a[2]
        True

    """
    def __init__(self, parent, ring, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import Components
            sage: Components(ZZ, [1,2,3], 2)
            2-index components w.r.t. [1, 2, 3]

        """
        # For efficiency, no test is performed regarding the type and range of
        # the arguments:
        dim = len(frame)
        # parent
        Element.__init__(self, parent)
        self._dim = dim
        self._nid = parent._nid
        self._sindex = start_index
        # element
        self._ring = ring
        self._frame = frame
        r = range(start_index, dim + start_index)
        self._ranges = tuple(r for x in range(self._nid))
        self._output_formatter = output_formatter

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._repr_()
            '2-index components w.r.t. [1, 2, 3]'

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: CompWithSym(ZZ, [1,2,3], 4, sym=(0,1))
            4-index components w.r.t. [1, 2, 3],
             with symmetry on the index positions (0, 1)
            sage: CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            4-index components w.r.t. [1, 2, 3],
             with symmetry on the index positions (0, 1),
             with antisymmetry on the index positions (2, 3)

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: CompFullySym(ZZ, (1,2,3), 4)
            Fully symmetric 4-index components w.r.t. (1, 2, 3)

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: CompFullyAntiSym(ZZ, (1,2,3), 4)
            Fully antisymmetric 4-index components w.r.t. (1, 2, 3)
        """
        prefix, suffix = self.parent()._repr_symmetry()
        return prefix + f"{self._nid}-index components w.r.t. {self._frame}" + suffix

    def _new_instance(self, parent=None):
        r"""
        Creates a :class:`Components` instance of the same number of indices
        and w.r.t. the same frame, but possibly with different symmetries. If
        ``parent`` is ``None``, the parent of ``self`` is used and the same
        symmetries are imposed.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._new_instance()
            2-index components w.r.t. [1, 2, 3]

        """
        if parent is None:
            parent = self.parent()
        return parent(self._ring, self._frame, start_index=self._sindex,
                      output_formatter=self._output_formatter)

    def _check_indices(self, indices):
        r"""
        Check the validity of a list of indices and returns a tuple from it

        INPUT:

        - ``indices`` -- list of indices (possibly a single integer if
          self is a 1-index object)

        OUTPUT:

        - a tuple containing valid indices

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._check_indices((0,1))
            (0, 1)
            sage: c._check_indices([0,1])
            (0, 1)
            sage: c._check_indices([2,1])
            (2, 1)
            sage: c._check_indices([2,3])
            Traceback (most recent call last):
            ...
            IndexError: index out of range: 3 not in range(0, 3)
            sage: c._check_indices(1)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 1 are provided
            sage: c._check_indices([1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 3 are provided

        """
        parent = self.parent()
        if isinstance(indices, (int, Integer)):
            ind = (indices,)
        else:
            ind = tuple(indices)
        return parent._check_indices(ind, self._ranges)

    def _get_list(self, ind_slice, no_format=True, format_type=None):
        r"""
        Return the list of components (as nested list or matrix).

        INPUT:

        - ``ind_slice`` --  a slice object. Unless the dimension is 1,
          this must be ``[:]``.
        - ``no_format`` -- (default: ``True``) determines whether some
          formatting of the components is to be performed
        - ``format_type`` -- (default: ``None``) argument to be passed
          to the formatting function ``self._output_formatter``, as the
          second (optional) argument

        OUTPUT:

        - general case: the nested list of components in the form
          ``T[i][j]...`` for the components `T_{ij...}`.

        - in the 1-dim case, a slice of that list if
          ``ind_slice = [a:b]``.

        - in the 2-dim case, a matrix (over the base ring of the components or
          of the formatted components if ``no_format`` is ``False``) is
          returned instead, except if the formatted components do not belong
          to any ring (for instance if they are strings).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c[0,1], c[1,2] = 5, -4
            sage: c._get_list(slice(None))
            [ 0  5  0]
            [ 0  0 -4]
            [ 0  0  0]
            sage: v = Components(ZZ, [1,2,3], 1)
            sage: v[:] = 4, 5, 6
            sage: v._get_list(slice(None))
            [4, 5, 6]
            sage: v._get_list(slice(0,1))
            [4]
            sage: v._get_list(slice(0,2))
            [4, 5]
            sage: v._get_list(slice(2,3))
            [6]
        """
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 1:
            if ind_slice.start is None:
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None:
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("function [start:stop:step] not implemented")
            if no_format:
                return [self[[i]] for i in range(start, stop)]
            else:
                return [self[i, format_type] for i in range(start, stop)]
        if ind_slice.start is not None or ind_slice.stop is not None:
            raise NotImplementedError("function [start:stop] not " +
                      "implemented for components with {} indices".format(self._nid))
        resu = [self._gen_list([i], no_format, format_type)
                for i in range(si, nsi)]
        if self._nid == 2:
            # 2-dim case: convert to matrix for a nicer output
            from sage.matrix.constructor import matrix
            from sage.structure.element import parent
            from sage.categories.rings import Rings
            if parent(resu[0][0]) in Rings():
                return matrix(resu)
        return resu

    def _gen_list(self, ind, no_format=True, format_type=None):
        r"""
        Recursive function to generate the list of values.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c[0,1], c[1,2] = 5, -4
            sage: c._gen_list([])
            [[0, 5, 0], [0, 0, -4], [0, 0, 0]]
            sage: c._gen_list([0])
            [0, 5, 0]
            sage: c._gen_list([1])
            [0, 0, -4]
            sage: c._gen_list([2])
            [0, 0, 0]
            sage: c._gen_list([0,1])
            5

        """
        if len(ind) == self._nid:
            if no_format:
                return self[ind]
            else:
                args = tuple(ind + [format_type])
                return self[args]
        else:
            si = self._sindex
            nsi = si + self._dim
            return [self._gen_list(ind + [i], no_format, format_type)
                    for i in range(si, nsi)]

    def _set_list(self, ind_slice, format_type, values):
        r"""
        Set the components from a list.

        INPUT:

        - ``ind_slice`` --  a slice object
        - ``format_type`` -- format possibly used to construct a ring element
        - ``values`` -- list of values for the components : the full list if
          ``ind_slice = [:]``, in the form ``T[i][j]...`` for the
          component `T_{ij...}`; in the 1-D case, ``ind_slice`` can be
          a slice of the full list, in the form  ``[a:b]``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._set_list(slice(None), None, [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: c[:]
            [0 1 2]
            [3 4 5]
            [6 7 8]

        """
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 1:
            if ind_slice.start is None:
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None:
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("function [start:stop:step] not implemented")
            for i in range(start, stop):
                self[i, format_type] = values[i-start]
        else:
            if ind_slice.start is not None or ind_slice.stop is not None:
                raise NotImplementedError("function [start:stop] not " +
                      "implemented for components with {} indices".format(self._nid))
            for i in range(si, nsi):
                self._set_value_list([i], format_type, values[i-si])

    def _set_value_list(self, ind, format_type, val):
        r"""
        Recursive function to set a list of values to ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._set_value_list([], None, [[1,2,3], [4,5,6], [7,8,9]])
            sage: c[:]
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: c._set_value_list([0], None, [-1,-2,-3])
            sage: c[:]
            [-1 -2 -3]
            [ 4  5  6]
            [ 7  8  9]
            sage: c._set_value_list([2,1], None, -8)
            sage: c[:]
            [-1 -2 -3]
            [ 4  5  6]
            [ 7 -8  9]

        """
        if len(ind) == self._nid:
            if format_type is not None:
                ind = tuple(ind + [format_type])
            self[ind] = val
        else:
            si = self._sindex
            nsi = si + self._dim
            for i in range(si, nsi):
                self._set_value_list(ind + [i], format_type, val[i-si])

    def display(self, symbol, latex_symbol=None, index_positions=None,
                index_labels=None, index_latex_labels=None,
                format_spec=None, only_nonzero=True, only_nonredundant=False):
        r"""
        Display all the components, one per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``symbol`` -- string (typically a single letter) specifying the
          symbol for the components
        - ``latex_symbol`` -- (default: ``None``) string specifying the LaTeX
          symbol for the components; if ``None``, ``symbol`` is used
        - ``index_positions`` -- (default: ``None``) string of length the
          number of indices of the components and composed of characters 'd'
          (for "down") or 'u' (for "up") to specify the position of each index:
          'd' corresponds to a subscript and 'u' to a superscript. If
          ``index_positions`` is ``None``, all indices are printed as
          subscripts
        - ``index_labels`` -- (default: ``None``) list of strings representing
          the labels of each of the individual indices within the index range
          defined at the construction of the object; if ``None``, integer
          labels are used
        - ``index_latex_labels`` -- (default: ``None``) list of strings
          representing the LaTeX labels of each of the individual indices
          within the index range defined at the construction of the object; if
          ``None``, integers labels are used
        - ``format_spec`` -- (default: ``None``) format specification passed
          to the output formatter declared at the construction of the object
        - ``only_nonzero`` -- (default: ``True``) boolean; if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- (default: ``False``) boolean; if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of 3-index components w.r.t. to the canonical basis of the
        free module `\ZZ^2` over the integer ring::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, (ZZ^2).basis(), 3)
            sage: c[0,1,0], c[1,0,1], c[1,1,1] = -2, 5, 3
            sage: c.display('c')
            c_010 = -2
            c_101 = 5
            c_111 = 3

        By default, only nonzero components are shown; to display all the
        components, it suffices to set the parameter ``only_nonzero`` to
        ``False``::

            sage: c.display('c', only_nonzero=False)
            c_000 = 0
            c_001 = 0
            c_010 = -2
            c_011 = 0
            c_100 = 0
            c_101 = 5
            c_110 = 0
            c_111 = 3

        By default, all indices are printed as subscripts, but any index
        position can be specified::

            sage: c.display('c', index_positions='udd')
            c^0_10 = -2
            c^1_01 = 5
            c^1_11 = 3
            sage: c.display('c', index_positions='udu')
            c^0_1^0 = -2
            c^1_0^1 = 5
            c^1_1^1 = 3
            sage: c.display('c', index_positions='ddu')
            c_01^0 = -2
            c_10^1 = 5
            c_11^1 = 3

        The LaTeX output is performed as an array, with the symbol adjustable
        if it differs from the text symbol::

            sage: latex(c.display('c', latex_symbol=r'\Gamma', index_positions='udd'))
            \begin{array}{lcl}
             \Gamma_{\phantom{\, 0}\,1\,0}^{\,0\phantom{\, 1}\phantom{\, 0}} & = & -2 \\
             \Gamma_{\phantom{\, 1}\,0\,1}^{\,1\phantom{\, 0}\phantom{\, 1}} & = & 5 \\
             \Gamma_{\phantom{\, 1}\,1\,1}^{\,1\phantom{\, 1}\phantom{\, 1}} & = & 3
            \end{array}

        The index labels can differ from integers::

            sage: c.display('c', index_labels=['x','y'])
            c_xyx = -2
            c_yxy = 5
            c_yyy = 3

        If the index labels are longer than a single character, they are
        separated by a comma::

            sage: c.display('c', index_labels=['r', 'th'])
            c_r,th,r = -2
            c_th,r,th = 5
            c_th,th,th = 3

        The LaTeX labels for the indices can be specified if they differ
        from the text ones::

            sage: c.display('c', index_labels=['r', 'th'],
            ....:           index_latex_labels=['r', r'\theta'])
            c_r,th,r = -2
            c_th,r,th = 5
            c_th,th,th = 3

        The display of components with symmetries is governed by the parameter
        ``only_nonredundant``::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, (ZZ^2).basis(), 3, sym=(1,2)) ; c
            3-index components w.r.t. [
            (1, 0),
            (0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: c[0,0,1] = 2
            sage: c.display('c')
            c_001 = 2
            c_010 = 2
            sage: c.display('c', only_nonredundant=True)
            c_001 = 2

        If some nontrivial output formatter has been set, the format can be
        specified by means of the argument ``format_spec``::

            sage: c = Components(QQ, (QQ^3).basis(), 2,
            ....:                output_formatter=Rational.numerical_approx)
            sage: c[0,1] = 1/3
            sage: c[2,1] = 2/7
            sage: c.display('C')  # default format (53 bits of precision)
            C_01 = 0.333333333333333
            C_21 = 0.285714285714286
            sage: c.display('C', format_spec=10)  # 10 bits of precision
            C_01 = 0.33
            C_21 = 0.29

        Check that the bug reported in :trac:`22520` is fixed::

            sage: c = Components(SR, [1, 2], 1)
            sage: c[0] = SR.var('t', domain='real')
            sage: c.display('c')
            c_0 = t

        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        si = self._sindex
        nsi = si + self._dim
        if latex_symbol is None:
            latex_symbol = symbol
        if index_positions is None:
            index_positions = self._nid * 'd'
        elif len(index_positions) != self._nid:
            raise ValueError("the argument 'index_positions' must contain " +
                             "{} characters".format(self._nid))
        if index_labels is None:
            index_labels = [str(i) for i in range(si, nsi)]
        elif len(index_labels) != self._dim:
            raise ValueError("the argument 'index_labels' must contain " +
                             "{} items".format(self._dim))
        # Index separator:
        max_len_symbols = max(len(s) for s in index_labels)
        if max_len_symbols == 1:
            sep = ''
        else:
            sep = ','
        if index_latex_labels is None:
            index_latex_labels = index_labels
        elif len(index_latex_labels) != self._dim:
            raise ValueError("the argument 'index_latex_labels' must " +
                             "contain {} items".format(self._dim))
        if only_nonredundant:
            generator = self.non_redundant_index_generator()
        else:
            generator = self.index_generator()
        rtxt = ''
        rlatex = r'\begin{array}{lcl}'
        for ind in generator:
            ind_arg = ind + (format_spec,)
            val = self[ind_arg]
            # Check whether the value is zero, preferably via the
            # fast method is_trivial_zero():
            if hasattr(val, 'is_trivial_zero'):
                zero_value = val.is_trivial_zero()
            else:
                zero_value = val == 0
            if not zero_value or not only_nonzero:
                indices = ''  # text indices
                d_indices = '' # LaTeX down indices
                u_indices = '' # LaTeX up indices
                previous = None  # position of previous index
                for k in range(self._nid):
                    i = ind[k] - si
                    if index_positions[k] == 'd':
                        if previous == 'd':
                            indices += sep + index_labels[i]
                        else:
                            indices += '_' + index_labels[i]
                        d_indices += r'\,' + index_latex_labels[i]
                        u_indices += r'\phantom{{\, {}}}'.format(index_latex_labels[i])
                        previous = 'd'
                    else:
                        if previous == 'u':
                            indices += sep + index_labels[i]
                        else:
                            indices += '^' + index_labels[i]
                        d_indices += r'\phantom{{\, {}}}'.format(index_latex_labels[i])
                        u_indices += r'\,' + index_latex_labels[i]
                        previous = 'u'
                rtxt += symbol + indices + ' = {} \n'.format(val)
                rlatex += (latex_symbol + r'_{' + d_indices + r'}^{'
                           + u_indices + r'} & = & ' + latex(val) + r'\\')
        if rtxt == '':
            # no component has been displayed
            rlatex = ''
        else:
            # closing the display
            rtxt = rtxt[:-1]  # remove the last new line
            rlatex = rlatex[:-2] + r'\end{array}'
        return FormattedExpansion(rtxt, rlatex)

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a set of components or 0

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c.__eq__(0)  # uninitialized components are zero
            True
            sage: c[0,1], c[1,2] = 5, -4
            sage: c.__eq__(0)
            False
            sage: c1 = Components(ZZ, [1,2,3], 2)
            sage: c1[0,1] = 5
            sage: c.__eq__(c1)
            False
            sage: c1[1,2] = -4
            sage: c.__eq__(c1)
            True
            sage: v = Components(ZZ, [1,2,3], 1)
            sage: c.__eq__(v)
            False

        """
        if isinstance(other, (int, Integer)): # other is 0
            if other == 0:
                return self.is_zero()
            else:
                raise TypeError("cannot compare a set of components to a number")
        elif not isinstance(other, Components_base):
            return False
        else:
            # other is another Components
            if other._frame != self._frame:
                return False
            if other._nid != self._nid:
                return False
            if other._sindex != self._sindex:
                return False
            if other._output_formatter != self._output_formatter:
                return False
            return (self - other).is_zero()

    def __ne__(self, other):
        r"""
        Non-equality operator.

        INPUT:

        - ``other`` -- a set of components or 0

        OUTPUT:

        - True if ``self`` is different from ``other``,  or False otherwise

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 1)
            sage: c.__ne__(0)  # uninitialized components are zero
            False
            sage: c1 = Components(ZZ, [1,2,3], 1)
            sage: c.__ne__(c1)  # c and c1 are both zero
            False
            sage: c[0] = 4
            sage: c.__ne__(0)
            True
            sage: c.__ne__(c1)
            True

        """
        return not self == other

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 1)
            sage: c[:] = 5, 0, -4
            sage: a = c.__pos__() ; a
            1-index components w.r.t. [1, 2, 3]
            sage: a[:]
            [5, 0, -4]
            sage: a == +c
            True
            sage: a == c
            True

        """
        return self.copy()

    def __radd__(self, other):
        r"""
        Reflected addition (addition on the right to `other``)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__radd__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [5, 5, 3]
            sage: s == a+b
            True
            sage: s = 0 + a ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s == a
            True

        """
        return self + other

    def __sub__(self, other):
        r"""
        Component subtraction.

        INPUT:

        - ``other`` -- components, of the same type as ``self``

        OUTPUT:

        - components resulting from the subtraction of ``other`` from ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__sub__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [-3, -5, -9]
            sage: s == a - b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__sub__(b) ; s_par
            1-index components w.r.t. [1, 2, 3]
            sage: s_par[:]
            [-3, -5, -9]
            sage: s_par == s
            True
            sage: b.__sub__(a) == -s # test of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        return self + (-other)  #!# correct, deals properly with
                                # symmetries, but is probably not optimal

    def __rsub__(self, other):
        r"""
        Reflected subtraction (subtraction from ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__rsub__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [3, 5, 9]
            sage: s == b - a
            True
            sage: s = 0 - a ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [-1, 0, 3]
            sage: s == -a
            True

        """
        return (-self) + other


    def contract(self, *args):
        r"""
        Contraction on one or many indices with another instance of
        :class:`Components`.

        INPUT:

        - ``pos1`` -- positions of the indices in ``self`` involved in the
          contraction; ``pos1`` must be a sequence of integers, with 0 standing
          for the first index position, 1 for the second one, etc. If ``pos1``
          is not provided, a single contraction on the last index position of
          ``self`` is assumed
        - ``other`` -- the set of components to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``. If ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT:

        - set of components resulting from the contraction

        EXAMPLES:

        Contraction of a 1-index set of components with a 2-index one::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = (-1, 2, 3)
            sage: b = Components(QQ, V.basis(), 2)
            sage: b[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s0 = a.contract(0, b, 0) ; s0
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s0[:]
            [28, 32, 36]
            sage: s0[:] == [sum(a[j]*b[j,i] for j in range(3)) for i in range(3)]  # check
            True
            sage: s1 = a.contract(0, b, 1) ; s1[:]
            [12, 24, 36]
            sage: s1[:] == [sum(a[j]*b[i,j] for j in range(3)) for i in range(3)]  # check
            True

        Parallel computations (see
        :class:`~sage.parallel.parallelism.Parallelism`)::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s0_par = a.contract(0, b, 0) ; s0_par
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s0_par[:]
            [28, 32, 36]
            sage: s0_par == s0
            True
            sage: s1_par = a.contract(0, b, 1) ; s1_par[:]
            [12, 24, 36]
            sage: s1_par == s1
            True
            sage: Parallelism().set('tensor', nproc = 1)  # switch off parallelization

        Contraction on 2 indices::

            sage: c = a*b ; c
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s = c.contract(1,2, b, 0,1) ; s
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [-285, 570, 855]
            sage: [sum(sum(c[i,j,k]*b[j,k] for k in range(3)) # check
            ....:      for j in range(3)) for i in range(3)]
            [-285, 570, 855]

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: c_par = a*b ; c_par
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c_par == c
            True
            sage: s_par = c_par.contract(1,2, b, 0,1) ; s_par
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s_par[:]
            [-285, 570, 855]
            sage: s_par == s
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Consistency check with :meth:`trace`::

            sage: b = a*a ; b   # the tensor product of a with itself
            Fully symmetric 2-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: b[:]
            [ 1 -2 -3]
            [-2  4  6]
            [-3  6  9]
            sage: b.trace(0,1)
            14
            sage: a.contract(0, a, 0) == b.trace(0,1)
            True

        """
        #
        # Treatment of the input
        #
        nargs = len(args)
        for i, arg in enumerate(args):
            if isinstance(arg, Components_base):
                other = arg
                it = i
                break
        else:
            raise ValueError("a set of components must be provided in the " +
                             "argument list")
        if it == 0:
            pos1 = (self._nid - 1,)
        else:
            pos1 = args[:it]
        if it == nargs-1:
            pos2 = (0,)
        else:
            pos2 = args[it+1:]
        ncontr = len(pos1) # number of contractions
        if len(pos2) != ncontr:
            raise TypeError("Different number of indices for the contraction.")
        if other._frame != self._frame:
            raise TypeError("The two sets of components are not defined on " +
                            "the same frame.")
        if other._sindex != self._sindex:
            raise TypeError("The two sets of components do not have the " +
                            "same starting index.")
        contractions = [(pos1[i], pos2[i]) for i in range(ncontr)]
        res_nid = self._nid + other._nid - 2*ncontr
        #
        # Special case of a scalar result
        #
        if res_nid == 0:
            # To generate the indices tuples (of size ncontr) involved in the
            # the contraction, we create an empty instance of Components with
            # ncontr indices and call the method index_generator() on it:
            from .comp import Components
            comp_for_contr = Components(self._ring, self._frame, ncontr,
                                        start_index=self._sindex)
            res = 0


            if Parallelism().get('tensor') != 1:
                # parallel contraction to scalar

                # parallel multiplication
                @parallel(p_iter='multiprocessing',ncpus=Parallelism().get('tensor'))
                def compprod(a,b):
                    return a*b

                # parallel list of inputs
                partial = list(compprod([(other[[ind]],self[[ind]]) for ind in
                                     comp_for_contr.index_generator()
                    ]))
                res = sum(map(itemgetter(1),partial))
            else:
                # sequential
                res = 0
                for ind in comp_for_contr.index_generator():
                    res += self[[ind]] * other[[ind]]

            return res


        #
        # Positions of self and other indices in the result
        #  (None = the position is involved in a contraction and therefore
        #   does not appear in the final result)
        #
        pos_s = [None for i in range(self._nid)]  # initialization
        pos_o = [None for i in range(other._nid)] # initialization
        shift = 0
        for pos in range(self._nid):
            for contract_pair in contractions:
                if pos == contract_pair[0]:
                    shift += 1
                    break
            else:
                pos_s[pos] = pos - shift
        for pos in range(other._nid):
            for contract_pair in contractions:
                if pos == contract_pair[1]:
                    shift += 1
                    break
            else:
                pos_o[pos] = self._nid + pos - shift
        rev_s = [pos_s.index(i) for i in range(self._nid-ncontr)]
        rev_o = [pos_o.index(i) for i in range(self._nid-ncontr, res_nid)]
        #
        # Determination of the symmetries of the result
        #
        max_len_sym = 0 # maximum length of symmetries in the result
        max_len_antisym = 0 # maximum length of antisymmetries in the result
        if res_nid <= 1:
            # no need to search for symmetries if res_nid == 1
            res_sym = ()
            res_antisym = ()
        else:
            from .comp_element_dict import ComponentsWithSym_dict, ComponentsWithSym_dict
            if isinstance(self, ComponentsWithSym_dict):
                s_sym = self.parent()._sym
                s_antisym = self.parent()._antisym
            else:
                s_sym = []
                s_antisym = []
            if isinstance(other, ComponentsWithSym_dict):
                o_sym = other.parent()._sym
                o_antisym = other.parent()._antisym
            else:
                o_sym = []
                o_antisym = []
            res_sym = []
            res_antisym = []
            for isym in s_sym:
                r_isym = []
                for pos in isym:
                    if pos_s[pos] is not None:
                        r_isym.append(pos_s[pos])
                if len(r_isym) > 1:
                    res_sym.append(r_isym)
                    max_len_sym = max(max_len_sym, len(r_isym))
            for isym in s_antisym:
                r_isym = []
                for pos in isym:
                    if pos_s[pos] is not None:
                        r_isym.append(pos_s[pos])
                if len(r_isym) > 1:
                    res_antisym.append(r_isym)
                    max_len_antisym = max(max_len_antisym, len(r_isym))
            for isym in o_sym:
                r_isym = []
                for pos in isym:
                    if pos_o[pos] is not None:
                        r_isym.append(pos_o[pos])
                if len(r_isym) > 1:
                    res_sym.append(r_isym)
                    max_len_sym = max(max_len_sym, len(r_isym))
            for isym in o_antisym:
                r_isym = []
                for pos in isym:
                    if pos_o[pos] is not None:
                        r_isym.append(pos_o[pos])
                if len(r_isym) > 1:
                    res_antisym.append(r_isym)
                    max_len_antisym = max(max_len_antisym, len(r_isym))
        #
        # Construction of the result object in view of the remaining symmetries:
        #
        from .comp import Components
        res = Components(self._ring, self._frame, res_nid,
                         start_index=self._sindex,
                         output_formatter=self._output_formatter,
                         sym=res_sym, antisym=res_antisym)
        #
        # Performing the contraction
        #
        # To generate the indices tuples (of size ncontr) involved in the
        # the contraction, we create an empty instance of Components with
        # ncontr indices and call the method index_generator() on it:
        comp_for_contr = Components(self._ring, self._frame, ncontr,
                                    start_index=self._sindex)
        shift_o = self._nid - ncontr

        if Parallelism().get('tensor') != 1:
            # parallel computation
            nproc = Parallelism().get('tensor')
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = [ind for ind in res.non_redundant_index_generator()]
            ind_step = max(1,int(len(ind_list)/nproc/2))
            local_list = lol(ind_list,ind_step)

            listParalInput = []
            for ind_part in local_list:
                listParalInput.append((self,other,ind_part,rev_s,rev_o,shift_o,contractions,comp_for_contr))

            # definition of the parallel function
            @parallel(p_iter='multiprocessing',ncpus=nproc)
            def make_Contraction(this,other,local_list,rev_s,rev_o,shift_o,contractions,comp_for_contr):
                local_res = []
                for ind in local_list:
                    ind_s = [None for i in range(this._nid)]  # initialization
                    ind_o = [None for i in range(other._nid)] # initialization
                    for i, pos in enumerate(rev_s):
                        ind_s[pos] = ind[i]
                    for i, pos in enumerate(rev_o):
                        ind_o[pos] = ind[shift_o+i]
                    sm = 0
                    for ind_c in comp_for_contr.index_generator():
                        ic = 0
                        for pos_s, pos_o in contractions:
                            k = ind_c[ic]
                            ind_s[pos_s] = k
                            ind_o[pos_o] = k
                            ic += 1
                        sm += this[[ind_s]] * other[[ind_o]]
                    local_res.append([ind,sm])
                return local_res

            for ii, val in make_Contraction(listParalInput):
                for jj in val:
                      res[[jj[0]]] = jj[1]
        else:
            # sequential
            for ind in res.non_redundant_index_generator():
                ind_s = [None for i in range(self._nid)]  # initialization
                ind_o = [None for i in range(other._nid)] # initialization
                for i, pos in enumerate(rev_s):
                    ind_s[pos] = ind[i]
                for i, pos in enumerate(rev_o):
                    ind_o[pos] = ind[shift_o+i]
                sm = 0
                for ind_c in comp_for_contr.index_generator():
                    ic = 0
                    for pos_s, pos_o in contractions:
                        k = ind_c[ic]
                        ind_s[pos_s] = k
                        ind_o[pos_o] = k
                        ic += 1
                    sm += self[[ind_s]] * other[[ind_o]]
                res[[ind]] = sm

        return res

    def index_generator(self):
        r"""
        Generator of indices.

        OUTPUT:

        - an iterable index

        EXAMPLES:

        Indices on a 3-dimensional vector space::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: c = Components(QQ, V.basis(), 1)
            sage: list(c.index_generator())
            [(0,), (1,), (2,)]
            sage: c = Components(QQ, V.basis(), 1, start_index=1)
            sage: list(c.index_generator())
            [(1,), (2,), (3,)]
            sage: c = Components(QQ, V.basis(), 2)
            sage: list(c.index_generator())
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]

        """
        parent = self.parent()
        yield from parent.index_generator(self._ranges)

    def non_redundant_index_generator(self):
        r"""
        Generator of non redundant indices.

        In the absence of declared symmetries, all possible indices are
        generated. So this method is equivalent to :meth:`index_generator`.
        Only versions for derived classes with symmetries or antisymmetries
        are not trivial.

        OUTPUT:

        - an iterable index

        EXAMPLES:

        Indices on a 3-dimensional vector space::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]
            sage: c = Components(QQ, V.basis(), 2, start_index=1)
            sage: list(c.non_redundant_index_generator())
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1),
             (3, 2), (3, 3)]

        """
        parent = self.parent()
        yield from parent.non_redundant_index_generator(self._ranges)

    def symmetrize(self, *pos):
        r"""
        Symmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the
          symmetrization (with the convention position=0 for the first slot);
          if none, the symmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`ComponentsWithSym` describing the symmetrized
          components

        EXAMPLES:

        Symmetrization of 2-index components::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s = c.symmetrize() ; s
            Fully symmetric 2-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            (
            [1 2 3]  [1 3 5]
            [4 5 6]  [3 5 7]
            [7 8 9], [5 7 9]
            )
            sage: c.symmetrize() == c.symmetrize(0,1)
            True

        Full symmetrization of 3-index components::

            sage: c = Components(QQ, V.basis(), 3)
            sage: c[:] = [[[1,2,3], [4,5,6], [7,8,9]], [[10,11,12], [13,14,15], [16,17,18]], [[19,20,21], [22,23,24], [25,26,27]]]
            sage: s = c.symmetrize() ; s
            Fully symmetric 3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 16/3, 29/3], [16/3, 29/3, 14], [29/3, 14, 55/3]],
              [[16/3, 29/3, 14], [29/3, 14, 55/3], [14, 55/3, 68/3]],
              [[29/3, 14, 55/3], [14, 55/3, 68/3], [55/3, 68/3, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[i,k,j]+c[j,k,i]+c[j,i,k]+c[k,i,j]+c[k,j,i])/6  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: c.symmetrize() == c.symmetrize(0,1,2)
            True

        Partial symmetrization of 3-index components::

            sage: s = c.symmetrize(0,1) ; s   # symmetrization on the first two indices
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
              [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
              [[13, 14, 15], [19, 20, 21], [25, 26, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[j,i,k])/2   # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.symmetrize(1,2) ; s   # symmetrization on the last two indices
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 3, 5], [3, 5, 7], [5, 7, 9]],
              [[10, 12, 14], [12, 14, 16], [14, 16, 18]],
              [[19, 21, 23], [21, 23, 25], [23, 25, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[i,k,j])/2   # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.symmetrize(0,2) ; s   # symmetrization on the first and last indices
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 2)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 6, 11], [4, 9, 14], [7, 12, 17]],
              [[6, 11, 16], [9, 14, 19], [12, 17, 22]],
              [[11, 16, 21], [14, 19, 24], [17, 22, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[k,j,i])/2   # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True

        """
        new_parent = self.parent().symmetrize(*pos)
        result = self._new_instance(parent=new_parent)
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        n_sym = len(pos)  # number of indices involved in the antisymmetry
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            for perm in sym_group.list():
                # action of the permutation on [0,1,...,n_sym-1]:
                perm_action = [x - 1 for x in perm.domain()]
                ind_perm = list(ind)
                for k in range(n_sym):
                    ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                sum += self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def antisymmetrize(self, *pos):
        r"""
        Antisymmetrization over the given index positions

        INPUT:

        - ``pos`` -- list of index positions involved in the antisymmetrization
          (with the convention position=0 for the first slot); if none, the
          antisymmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`ComponentsWithSym` describing the antisymmetrized
          components.

        EXAMPLES:

        Antisymmetrization of 2-index components::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 2-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            (
            [1 2 3]  [ 0 -1 -2]
            [4 5 6]  [ 1  0 -1]
            [7 8 9], [ 2  1  0]
            )
            sage: c.antisymmetrize() == c.antisymmetrize(0,1)
            True

        Full antisymmetrization of 3-index components::

            sage: c = Components(QQ, V.basis(), 3)
            sage: c[:] = [[[-1,-2,3], [4,-5,4], [-7,8,9]], [[10,10,12], [13,-14,15], [-16,17,19]], [[-19,20,21], [1,2,3], [-25,26,27]]]
            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, 0, 0], [0, 0, -13/6], [0, 13/6, 0]],
              [[0, 0, 13/6], [0, 0, 0], [-13/6, 0, 0]],
              [[0, -13/6, 0], [13/6, 0, 0], [0, 0, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[i,k,j]+c[j,k,i]-c[j,i,k]+c[k,i,j]-c[k,j,i])/6  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: c.symmetrize() == c.symmetrize(0,1,2)
            True

        Partial antisymmetrization of 3-index components::

            sage: s = c.antisymmetrize(0,1) ; s  # antisymmetrization on the first two indices
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, 0, 0], [-3, -15/2, -4], [6, -6, -6]],
              [[3, 15/2, 4], [0, 0, 0], [-17/2, 15/2, 8]],
              [[-6, 6, 6], [17/2, -15/2, -8], [0, 0, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[j,i,k])/2  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.antisymmetrize(1,2) ; s  # antisymmetrization on the last two indices
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, -3, 5], [3, 0, -2], [-5, 2, 0]],
              [[0, -3/2, 14], [3/2, 0, -1], [-14, 1, 0]],
              [[0, 19/2, 23], [-19/2, 0, -23/2], [-23, 23/2, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[i,k,j])/2  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.antisymmetrize(0,2) ; s  # antisymmetrization on the first and last indices
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 2)
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, -6, 11], [0, -9, 3/2], [0, 12, 17]],
              [[6, 0, -4], [9, 0, 13/2], [-12, 0, -7/2]],
              [[-11, 4, 0], [-3/2, -13/2, 0], [-17, 7/2, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[k,j,i])/2  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True

        The order of index positions in the argument does not matter::

            sage: c.antisymmetrize(1,0) == c.antisymmetrize(0,1)
            True
            sage: c.antisymmetrize(2,1) == c.antisymmetrize(1,2)
            True
            sage: c.antisymmetrize(2,0) == c.antisymmetrize(0,2)
            True

        """
        new_parent = self.parent().antisymmetrize(*pos)
        result = self._new_instance(parent=new_parent)
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        n_sym = len(pos)  # number of indices involved in the antisymmetry
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            for perm in sym_group.list():
                # action of the permutation on [0,1,...,n_sym-1]:
                perm_action = [x - 1 for x in perm.domain()]
                ind_perm = list(ind)
                for k in range(n_sym):
                    ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                if perm.sign() == 1:
                    sum += self[[ind_perm]]
                else:
                    sum -= self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def _matrix_(self):
        r"""
        Convert a set of ring components with 2 indices into a matrix.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2, start_index=1)
            sage: c[:] = [[-1,2,3], [4,-5,6], [7,8,-9]]
            sage: c._matrix_()
            [-1  2  3]
            [ 4 -5  6]
            [ 7  8 -9]

            sage: matrix(c) == c._matrix_()
            True

        """
        from sage.matrix.constructor import matrix
        if self._nid != 2:
            raise ValueError("the set of components must have 2 indices")
        si = self._sindex
        nsi = self._dim + si
        tab = [[self[[i,j]] for j in range(si, nsi)] for i in range(si, nsi)]
        return matrix(tab)
