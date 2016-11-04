r"""
Arrays of non-negative integers with specified slice sums.

Given a sequence `x = (x_1,...,x_n)` of weak compositions with

.. MATH::

    x_i = (x_{i1}, \dotsc, x_{im_i}),

an integer array of shape `x` is an `n`-dimensional array

.. MATH::

    a_{r_1 r_2 \dotsb r_n}, 1\leq r_j \leq x_j

of non-negative integers such that the slice sums

.. MATH::

    \sum_{r_i = j} a_{r_1 r_2 \dotsb r_n} = x_{ij}

for `1\leq i\leq n` and `1\leq j\leq x_{ij}`.

When `n=2`, these coincide with integer matrices.

.. SEEALSO::

    :mod:`sage.combinat.integer_matrices`.

AUTHOR:

- Amritanshu Prasad (2016-11-04): initial implementation

**Relation to Kronecker coefficients**

Given a sequence of partitions `x_1,...,x_n` of a fixed integer `N`, their
Kronecker coefficient `g_{x_1...x_n}` is defined as the multiplicity of
the trivial representation in tensor product

.. MATH::

    V_{x_1}\otimes ... \otimes V_{x_n},

where `V_{x_i}` is the irreducible representation of the symmetric group
corresponding to `x_i`.

Using Young's rule one can show that the number of integer arrays with
slice sums `x_1,...,x_n` equals

.. MATH::

    \sum_{y_1,...,y_n} g_{y_1...y_n}K_{y_1 x_1}...K_{y_n x_n},

the sum being over all `n`-tuples of partitions of `N`, and where `K_{y_i x_i}`
is the number of semistandard Young tableaux of shape `y_i` and type `x_i`.

When `n=2`, the Kronecker coefficient `g_{x_1 x_2}` is `1` is `x_1=x_2` and `0`
otherwise. In this case, the Robinson-Schencted-Knuth correspondence is a
bijectivization of this combinatorial identity (see [Prasad15]_).

We test the identity in an example with `n=3`::

    sage: from itertools import product
    sage: def kronecker_coefficient(partns):
    ....:     S = SymmetricFunctions(QQ).schur()
    ....:     pr = reduce(lambda x, y: x.itensor(y), [S[la] for la in partns])
    ....:     return pr.coefficient(Partition([sum(partns[0])]))
    ....:
    sage: len(list(IntegerArrays([[2, 1], [2, 1], [2, 1]]))) == sum([kronecker_coefficient(partns)*
    ....: prod([symmetrica.kostka_number(partns[i], [2,1]) for i in range(3)], 1) for partns in 
    ....: product(Partitions(3), repeat=3)])
    True

REFERENCES:

.. [Prasad15] Prasad, A.  *Representation Theory: A Combinatorial Viewpoint*.
   Cambridge University Press, 2016.


"""
#*****************************************************************************
#       Copyright (C) 2016 Amritanshu Prasad <amri@imsc.res.in>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import product
from numpy import array
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.composition import Composition
from sage.numerical.mip import MixedIntegerLinearProgram

class IntegerArray(Element):
    r"""
    An array of non-negative integers.

    INPUT:

    - ``arr`` a numpy array.

    EXAMPLE::

        sage: from numpy import array
        sage: IntegerArray(array([[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]]))
        [[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]]
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, arr):
        r"""
        Create an integer array.

        EXAMPLES::

            sage: IntegerArray([[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]])
            [[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]]

        The parent class is the class of integer arrays with its slice sums.

            sage: IA = IntegerArray([[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]])
            sage: IA.parent().slice_sums()
            ([10, 4], [9, 5], [6, 4, 4])
        """
        arr = array(arr)
        d = len(arr.shape)
        def noi(n, i):
            L = range(n)
            L.pop(i)
            return tuple(L)
        slice_sums = [Composition(arr.sum(axis=noi(d, axis))) for axis in range(d)]
        P = IntegerArrays(slice_sums)
        return P(arr)

    def __init__(self, parent, arr):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: elt = IntegerArray([[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]])
            sage: TestSuite(elt).run()
        """
        self._arr = arr
        Element.__init__(self, parent)

    def __repr__(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: IntegerArray([[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]])
            [[[3, 2, 1], [2, 1, 1]], [[1, 1, 1], [0, 0, 1]]]
        """
        return str(self._arr.tolist())

    def __hash__(self):
        r"""
        TESTS::

            sage: hash(IntegerArray([]))
            11648069979105038
            sage: hash(IntegerArray([1]))
            8209412804330245758
            sage: from numpy import array
            sage: hash(IntegerArray(array(range(12)).reshape([2, 2, 3])))
            -784913444396038158
        """
        return hash(str(self))

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: IA1 = IntegerArray([])
            sage: IA2 = IntegerArray([1])
            sage: IA3 = IntegerArray([2])
            sage: IA4 = IntegerArray([1])
            sage: IA1 == IA2
            False
            sage: IA2 == IA3
            False
            sage: IA2 == IA4
            True
        """
        if not isinstance(other, IntegerArray):
            return False
        elif self.slice_sums() != other.slice_sums():
            return False
        else:
            b = self._arr == other._arr
            if isinstance(b, bool):
                return b
            else:
                return b.all()

    def __ne__(self, other):
        r"""
        TESTS::

            sage: IA1 = IntegerArray([])
            sage: IA2 = IntegerArray([1, 1])
            sage: IA3 = IntegerArray([2, 1])
            sage: IA4 = IntegerArray([1, 1])
            sage: IA1 != IA4
            True
            sage: IA1 != IA2
            True
            sage: IA2 != IA3
            True
            sage: IA2 != IA4
            False
        """
        if not isinstance(other, IntegerArray):
            return True
        elif self.slice_sums() != other.slice_sums():
            return True
        else:
            b = self._arr != other._arr
        if isinstance(b, bool):
            return b
        else:
            return b.any()

    def slice_sums(self):
        """
        Return the slice sums of ``self``.

        EXAMPLES::

            sage: IntegerArray([[[2], [1]], [[3], [1]], [[1], [3]]]).slice_sums()
            ([3, 4, 4], [6, 5], [11])
        """
        return self.parent().slice_sums()

    def array(self):
        """
        Return numpy array corresponding to ``self``.

        EXAMPLES::

            sage: A = IntegerArray([[[2], [1]], [[3], [1]], [[1], [3]]])
            sage: print A.array()
            [[[2]
              [1]]
            <BLANKLINE>
             [[3]
              [1]]
            <BLANKLINE>
             [[1]
              [3]]]
        """
        return self._arr

class IntegerArrays(UniqueRepresentation, Parent):
    r"""
    The class of integer arrays with slice sums given by ``slice_sums``.

    INPUT:

    - ``slice_sums`` -- list of vectors with non-negative integer entries.

    EXAMPLES:

        sage: IAC = IntegerArrays([[2, 1], [3], [1, 1, 1]])
        sage: for a in IAC:
        ....:     print a
        [[[0, 1, 1]], [[1, 0, 0]]]
        [[[1, 0, 1]], [[0, 1, 0]]]
        [[[1, 1, 0]], [[0, 0, 1]]]
    """
    @staticmethod
    def __classcall_private__(cls, slice_sums):
        r"""
        Create the class of integer arrays with slice sums given by
        ``slice_sums``.

        EXAMPLES::

            sage: IA1 = IntegerArrays([[2, 1], [3], [1, 1, 1]])
            sage: IA2 = IntegerArrays([[2, 1], [3], [1, 1, 1]])
            sage: IA1 is IA2
            True

        """
        slice_sums = tuple(map(Composition, slice_sums))
        return super(IntegerArrays, cls).__classcall__(cls, slice_sums)

    def __init__(self, slice_sums):
        r"""
        Initialize ``self``.

        TESTS::

            sage: IA = IntegerArrays([[2, 1], [3], [1, 1, 1]])
            sage: TestSuite(IA).run()
        """
        Parent.__init__(self, category = FiniteEnumeratedSets())
        self._slice_sums = slice_sums

    def _element_constructor_(self, arr):
        """
        Construct an element of ``self``.

        INPUT:

        - ``arr`` -- a numpy array of integers

        EXAMPLES::

            sage: IAC = IntegerArrays([[2, 1], [3], [1, 1, 1]])
            sage: elt = IAC([[[1, 0, 1]], [[0, 1, 0]]]); elt
            [[[1, 0, 1]], [[0, 1, 0]]]
            sage: elt.parent() is IAC
            True
        """
        arr = array(arr)
        return self.element_class(self, arr)

    Element = IntegerArray

    def __iter__(self):
        r"""
        Iterate over ``self``.

        """
        slice_sums = self._slice_sums
        p = MixedIntegerLinearProgram()
        a = p.new_variable(integer=True, nonnegative=True)
        d = len(slice_sums) # dimension of the array
        shape = map(len, slice_sums) # shape of the array
        entries = product(*map(range, shape)) # enumerator for entries of the array

        sums = [[0]*sh for sh in shape] # initializing all the sums of variable that will be constrained

        # for each entry of array, add the corresponding variable of p to sums

        for entry in entries:
            for sh in range(d):
                sums[sh][entry[sh]] += a[entry]

        # add constraints to p

        for sh in range(d):
            for slice in range(shape[sh]):
                p.add_constraint(sums[sh][slice] == slice_sums[sh][slice])

        # enumerate integer points in the polyhedron of p

        for pt in p.polyhedron().integral_points():
            yield self.element_class(self, array(pt).reshape(shape))

    def slice_sums(self):
        """
        Return slice sums of ``self``.
        """
        return self._slice_sums
