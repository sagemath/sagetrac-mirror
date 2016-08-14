"""
Matrix/Vector-Valued Linear Functions: Elements

Here is an example of a linear function tensored with a vector space::

    sage: mip.<x> = MixedIntegerLinearProgram('ppl')   # base ring is QQ
    sage: lt = x[0] * vector([3,4]) + 1;   lt
    (1, 1) + (3, 4)*x_0
    sage: type(lt)
    <type 'sage.numerical.linear_tensor_element.LinearTensor'>
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport *

from sage.misc.fast_methods cimport hash_by_id
from sage.structure.element cimport ModuleElement, RingElement
from sage.numerical.linear_functions cimport LinearFunction, is_LinearFunction


#*****************************************************************************
#
# Elements of linear functions tensored with a free module
#
#*****************************************************************************

cdef class LinearTensor(ModuleElement):
    r"""
    A linear function tensored with a free module

    .. warning::

        You should never instantiate :class:`LinearTensor`
        manually. Use the element constructor in the parent
        instead.

    EXAMPLES::

        sage: parent = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
        sage: parent({0: [1,2], 3: [-7,-8]})
        (1.0, 2.0)*x_0 + (-7.0, -8.0)*x_3
    """

    def __init__(self, parent, f):
        r"""
        Constructor taking a dictionary as its argument.

        INPUT:
        
        - ``parent`` -- the parent
          :class:`~sage.numerical.linear_tensor.LinearTensorParent_class`.
        
        - ``f`` -- A linear function tensored by a free module is
          represented as a dictionary. The values are the coefficient
          (free module elements) of the variable represented by the
          keys. The key ``-1`` corresponds to the constant term.

        EXAMPLES:

        With a dictionary::

            sage: LT = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: LT({0: [1,2], 3: [-7,-8]})
            (1.0, 2.0)*x_0 + (-7.0, -8.0)*x_3
        
            sage: TestSuite(LT).run(skip=['_test_an_element', '_test_elements_eq_reflexive',
            ....:     '_test_elements_eq_symmetric', '_test_elements_eq_transitive',
            ....:     '_test_elements_neq', '_test_additive_associativity',
            ....:     '_test_elements', '_test_pickling', '_test_zero'])
        """
        ModuleElement.__init__(self, parent)
        assert isinstance(f, dict)
        self._f = f

    def __getitem__(self, indices):
        """
        Return the linear function component with given tensor indices.

        INPUT:

        - ``indices`` -- one or more integers. The basis indices of
          the free module. E.g. a single integer for vectors, two for
          matrices.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: lt = p({0:[1,2], 3:[4,5]});  lt
            (1.0, 2.0)*x_0 + (4.0, 5.0)*x_3
            sage: lt[0]
            x_0 + 4*x_3
            sage: lt[1]
            2*x_0 + 5*x_3
        """
        f = dict([key, value[indices]] for key, value in self._f.items())
        LF = self.parent().linear_functions()
        return LF(f)

    def dict(self):
        r"""
        Return the dictionary corresponding to the tensor product.

        OUTPUT:

        The linear function tensor product is represented as a
        dictionary. The value are the coefficient (free module
        elements) of the variable represented by the keys (which are
        integers). The key ``-1`` corresponds to the constant term.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: lt = p({0:[1,2], 3:[4,5]})
            sage: lt.dict()
            {0: (1.0, 2.0), 3: (4.0, 5.0)}
        """
        return dict(self._f)

    def coefficient(self, x):
        r"""
        Return one of the the coefficients.

        INPUT:

        - ``x`` -- a linear variable or an integer. If an integer `i`
          is passed, then `x_i` is used as linear variable. Pass
          ``-1`` for the constant term.

        OUTPUT:

        A constant, that is, an element of the free module factor. The
        coefficient of ``x`` in the linear function.

        EXAMPLE::

            sage: mip.<b> = MixedIntegerLinearProgram()
            sage: lt = vector([1,2]) * b[3] + vector([4,5]) * b[0] - 5;  lt
            (-5.0, -5.0) + (1.0, 2.0)*x_0 + (4.0, 5.0)*x_1
            sage: lt.coefficient(b[3])
            (1.0, 2.0)
            sage: lt.coefficient(0)      # x_0 is b[3]
            (1.0, 2.0)
            sage: lt.coefficient(4)
            (0.0, 0.0)
            sage: lt.coefficient(-1)
            (-5.0, -5.0)

        TESTS::

            sage: lt.coefficient(b[3] + b[4])
            Traceback (most recent call last):
            ...
            ValueError: x is a sum, must be a single variable
            sage: lt.coefficient(2*b[3])
            Traceback (most recent call last):
            ...
            ValueError: x must have a unit coefficient
            sage: mip.<q> = MixedIntegerLinearProgram(solver='ppl')
            sage: lt.coefficient(q[0])
            Traceback (most recent call last):
            ...
            ValueError: x is from a different linear functions module
        """
        if is_LinearFunction(x):
            if self.parent().linear_functions() != x.parent():
                raise ValueError('x is from a different linear functions module')
            if len((<LinearFunction>x)._f) != 1:
                raise ValueError('x is a sum, must be a single variable')
            i = (<LinearFunction>x)._f.keys()[0]
            if (<LinearFunction>x)._f[i] != 1:
                raise ValueError('x must have a unit coefficient')
        else:
            i = int(x)
        try:
            return self._f[i]
        except KeyError:
            return self.parent().free_module().zero()

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: R.<s,t> = RDF[]
            sage: LT = LinearFunctionsParent(RDF).tensor(R)
            sage: LT.an_element()  # indirect doctest
            (s) + (5.0*s)*x_2 + (7.0*s)*x_5

            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: LT.an_element()  # indirect doctest
            (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
        """
        if self.parent().is_matrix_space():
            return self._repr_matrix()
        terms = []
        for key in sorted(self._f.keys()):
            coeff = self._f[key]
            if coeff._is_atomic():
                if key == -1:
                    term = '({1})'.format(key, coeff)
                else:
                    term = '({1})*x_{0}'.format(key, coeff)
            else:
                if key == -1:
                    term = '{1}'.format(key, coeff)
                else:
                    term = '{1}*x_{0}'.format(key, coeff)
            terms.append(term)
        return ' + '.join(terms)
            
    def _repr_matrix(self):
        """
        Return a matrix-like string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^(2,2))
            sage: LT.an_element()  # indirect doctest
            [1 + 5*x_2 + 7*x_5 0]
            [0                 0]
        """
        MS = self.parent().free_module()
        assert self.parent().is_matrix_space()
        col_lengths = []
        columns = []
        for c in range(MS.ncols()):
            column = []
            for r in range(MS.nrows()):
                cell = repr(self[r, c])
                column.append(cell)
            columns.append(column)
            col_lengths.append(max(map(len, column)))
        s = ''
        for r in range(MS.nrows()):
            if r > 0:
                s += '\n'
            s += '['
            for c in range(MS.ncols()):
                if c > 0:
                    s += ' '
                s += columns[c][r].ljust(col_lengths[c])
            s += ']'
        return s

    cpdef _add_(self, b):
        r"""
        Return sum.

        INPUT:
        
        - ``b`` -- a :class:`LinearTensor`.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: LT({0: [1,2], 3: [-7,-8]}) + LT({2: [5,6], 3: [2,-2]}) + 16
            (16.0, 16.0) + (1.0, 2.0)*x_0 + (5.0, 6.0)*x_2 + (-5.0, -10.0)*x_3
        """
        result = dict(self._f)
        for key, coeff in b.dict().iteritems():
            result[key] = self._f.get(key, 0) + coeff
        return self.parent()(result)

    cpdef _neg_(self):
        r"""
        Return the negative.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: -LT({0: [1,2], 3: [-7,-8]})
            (-1.0, -2.0)*x_0 + (7.0, 8.0)*x_3
        """
        result = dict()
        for key, coeff in self._f.items():
            result[key] = -coeff
        return self.parent()(result)

    cpdef _sub_(self, b):
        r"""
        Return difference.

        INPUT:
        
        - ``b`` -- a :class:`LinearTensor`.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: LT({0: [1,2], 3: [-7,-8]}) - LT({1: [1,2]})
            (1.0, 2.0)*x_0 + (-1.0, -2.0)*x_1 + (-7.0, -8.0)*x_3
            sage: LT({0: [1,2], 3: [-7,-8]}) - 16
            (-16.0, -16.0) + (1.0, 2.0)*x_0 + (-7.0, -8.0)*x_3
        """
        result = dict(self._f)
        for key, coeff in b.dict().iteritems():
            result[key] = self._f.get(key, 0) - coeff
        return self.parent()(result)

    cpdef _rmul_(self, RingElement b):
        r"""
        Return right multiplication by scalar.

        INPUT:

        - ``b`` -- base ring element. The scalar to multiply by.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: 10 * LT({0: [1,2], 3: [-7,-8]})
            (10.0, 20.0)*x_0 + (-70.0, -80.0)*x_3
        """
        result = dict()
        for key, coeff in self._f.items():
            result[key] = b * coeff
        return self.parent()(result)

    def __richcmp__(left, right, int op):
        """
        Create an inequality or equality object.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: lt0 = x[0] * vector([1,2])
            sage: lt1 = x[1] * vector([2,3])
            sage: lt0.__le__(lt1)    # indirect doctest
            (1.0, 2.0)*x_0 <= (2.0, 3.0)*x_1

        ::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: from sage.numerical.linear_functions import LinearFunction
            sage: x[0] * vector([1,2]) <= x[1] * vector([2,3])
            (1.0, 2.0)*x_0 <= (2.0, 3.0)*x_1

            sage: x[0] * vector([1,2]) >= x[1] * vector([2,3])
            (2.0, 3.0)*x_1 <= (1.0, 2.0)*x_0

            sage: x[0] * vector([1,2]) == x[1] * vector([2,3])
            (1.0, 2.0)*x_0 == (2.0, 3.0)*x_1

            sage: x[0] * vector([1,2]) < x[1] * vector([2,3])
            Traceback (most recent call last):
            ...
            ValueError: strict < is not allowed, use <= instead.

            sage: x[0] * vector([1,2]) > x[1] * vector([2,3])
            Traceback (most recent call last):
            ...
            ValueError: strict > is not allowed, use >= instead.

        TESTS::

            sage: lt = x[0] * vector([1,2])
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(10, lt, operator.le)
            Coercion on left operand via
                Conversion map:
                  From: Integer Ring
                  To:   Tensor product of Vector space of dimension 2 over Real 
                        Double Field and Linear functions over Real Double Field
            Arithmetic performed after coercions.
            Result lives in Tensor product of Vector space of dimension 2 over 
            Real Double Field and Linear functions over Real Double Field
            Tensor product of Vector space of dimension 2 over Real Double Field 
            and Linear functions over Real Double Field
        
            sage: operator.le(10, lt)
            (10.0, 10.0) <= (1.0, 2.0)*x_0
            sage: lt <= 1
            (1.0, 2.0)*x_0 <= (1.0, 1.0)
            sage: lt >= 1
            (1.0, 1.0) <= (1.0, 2.0)*x_0
            sage: 1 <= lt
            (1.0, 1.0) <= (1.0, 2.0)*x_0
            sage: 1 >= lt
            (1.0, 2.0)*x_0 <= (1.0, 1.0)
       """
        from sage.numerical.linear_tensor_constraints import \
            LinearTensorConstraintsParent
        LT = left.parent()
        LC = LinearTensorConstraintsParent(LT)
        left = LT(left)
        right = LT(right)
        if op == Py_LT:
            raise ValueError("strict < is not allowed, use <= instead.")
        elif op == Py_EQ:
            return LC(left, right, True)
        elif op == Py_GT:
            raise ValueError("strict > is not allowed, use >= instead.")
        elif op == Py_LE:
            return LC(left, right, False)
        elif op == Py_NE:
            raise ValueError("inequality != is not allowed, use one of <=, ==, >=.")
        elif op == Py_GE:
            return LC(right, left, False)
        else:
            assert(False)   # unreachable

    def __hash__(self):
        r"""
        Return a hash.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: lt0 = p[0] * vector([1,2])
            sage: lt0.__hash__()   # random output
            103987752
            sage: d = {}
            sage: d[lt0] = 3
        """
        # see _cmp_() if you want to change the hash function
        return hash_by_id(<void *> self)

    def __cmp__(left, right):
        """
        Implement comparison of two linear functions.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: f = p[0] * vector([1,2])
            sage: v0 = vector([0, 0])
            sage: v1 = vector([1, 1])
            sage: cmp(f, f)
            0
            sage: abs(cmp(f, f+v0))     # since we are comparing by id()
            1
            sage: abs(cmp(f, f+v1))
            1
            sage: len(set([f, f]))
            1
            sage: len(set([f, f+v0]))
            2
            sage: len(set([f, f+v1]))
            2
        """
        # Note: if you want to implement smarter comparison, you also
        # need to change __hash__(). The comparison function must
        # satisfy cmp(x,y)==0 => hash(x)==hash(y)
        if left is right:
            return 0
        if <size_t><void*>left < <size_t><void*>right:
            return -1
        else:
            return 1


