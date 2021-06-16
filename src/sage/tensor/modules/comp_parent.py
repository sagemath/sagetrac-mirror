r"""
Collection of components as indexed sets of ring elements w.r.t to certain
indices and symmetries
"""

#******************************************************************************
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

from sage.structure.parent import Parent
from .comp_element import Components_generic, ComponentsWithSym, ComponentsFullySym, ComponentsFullyAntiSym
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer

class CompParent(Parent, UniqueRepresentation):
    r"""

    """

    Element = Components_generic

    def __init__(self, nb_indices):
        r"""

        """
        Parent.__init__(self)
        self._nid = nb_indices

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: C = CompParent(2)
            sage: C._repr_()
            'Parent of 2-index components without symmetry'

        """
        prefix, suffix = self._repr_symmetry()
        return f"Parent of " + prefix + f"{self._nid}-index components" + suffix

    def _repr_symmetry(self):
        r"""
        Return a prefix and a suffix string describing the symmetry of ``self``.
        """
        return "", ""

    def _element_constructor_(self, *args, **kwargs):
        r"""
        Construct an indexed set of components w.r.t. ``frame`` over the ring
        ``ring``.
        """
        if isinstance(args[0], Components_generic):
            c = args[0]
            if not isinstance(c, Components_generic):
                raise TypeError("cannot coerce {} into an element of {}".format(c, self))
            else:
                return c
        if len(args) != 2:
            raise ValueError("{} missing {} required positional arguments".format(type(self), len(args)))
        ring = args[0]
        frame = args[1]
        start_index = kwargs.pop('start_index', 0)
        output_formatter = kwargs.pop('output_formatter', None)
        return self.element_class(self, ring, frame,
                                  start_index = start_index,
                                  output_formatter=output_formatter)

    def _check_indices(self, ind, ranges):
        r"""
        Check the validity of a list of indices and returns a tuple from it

        INPUT:

        - ``ind`` -- tuple of indices

        OUTPUT:

        - a tuple containing valid indices

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: cp = CompParent(2)
            sage: r = range(3)
            sage: ranges = (r, r)
            sage: cp._check_indices((0,1), ranges)
            (0, 1)
            sage: cp._check_indices([0,1], ranges)
            (0, 1)
            sage: cp._check_indices([2,1], ranges)
            (2, 1)
            sage: cp._check_indices([2,3], ranges)
            Traceback (most recent call last):
            ...
            IndexError: index out of range: 3 not in range(0, 3)
            sage: cp._check_indices([1], ranges)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 1 are provided
            sage: cp._check_indices([1,2,3], ranges)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 3 are provided

        """
        if len(ind) != self._nid:
            raise ValueError(("wrong number of indices: {} expected,"
                             " while {} are provided").format(self._nid, len(ind)))
        for i, range in zip(ind, ranges):
            if i not in range:
                raise IndexError(f"index out of range: {i} not in {range}")
        return tuple(ind)

    def index_generator(self, ranges):
        r"""
        Generator of indices.

        INPUT:

        - ``ranges`` -- a tuple of ranges for the indices

        OUTPUT:

        - an iterable index

        EXAMPLES:

        Indices on a 3-dimensional vector space::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: V = VectorSpace(QQ,3)
            sage: r = range(V.dimension())
            sage: cp = CompParent(1)
            sage: ranges = (r,)
            sage: list(cp.index_generator(ranges))
            [(0,), (1,), (2,)]
            sage: ranges = (range(1, 4),)
            sage: list(cp.index_generator(ranges))
            [(1,), (2,), (3,)]
            sage: cp = CompParent(2)
            sage: ranges = (r, r)
            sage: list(cp.index_generator(ranges))
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]

        """
        # this is really just itertools.product
        if len(ranges) != self._nid:
            raise TypeError('need a range for every index')
        ind = [r.start for r in ranges]
        ind_end = [r.start for r in ranges]
        ind_end[0] = ranges[0].stop
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(self._nid-1, -1, -1):
                if ind[pos] != ranges[pos].stop - 1:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = ranges[pos].stop # end point reached
                    else:
                        ind[pos] = ranges[pos].start
                        ret = 1

    def non_redundant_index_generator(self, ranges):
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

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: V = VectorSpace(QQ, 3)
            sage: r = range(V.dimension())
            sage: cp = CompParent(2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]
            sage: r_start_1 = range(1, V.dimension() + 1)
            sage: list(cp.non_redundant_index_generator((r_start_1, r_start_1)))
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1),
             (3, 2), (3, 3)]
        """
        yield from self.index_generator(ranges)

    @cached_method
    def symmetrize(self, *pos):
        r"""

        """
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the "
                                 "total number of indices")
        n_sym = len(pos) # number of indices involved in the symmetry
        if n_sym == self._nid:
            return CompParentFullySym(self._nid)
        return CompParentWithSym(self._nid, sym=pos)

    @cached_method
    def antisymmetrize(self, *pos):
        r"""

        """
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the "
                                 "total number of indices")
        n_sym = len(pos)  # number of indices involved in the antisymmetry
        if n_sym == self._nid:
            return CompParentFullyAntiSym(self._nid)
        return CompParentWithSym(self._nid, antisym=pos)

class CompParentWithSym(CompParent):
    r"""

    """

    @staticmethod
    def __classcall_private__(cls, nb_indices, sym=None, antisym=None):
        r"""
        Determine the correct class to return based upon the input. In
        particular, convert lists of symmetries to tuples and delegate to
        parent with correct symmetry.

        TESTS::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentWithSym
            sage: CompParentWithSym(3) is CompParent(3)
            True
            sage: CompParentWithSym(5, sym=(3, 4),) is CompParentWithSym(5, sym=((4, 3),))
            True
            sage: CompParentWithSym(5, sym=((1, 2), (3, 4))) is CompParentWithSym(5, sym=((4, 3), (2, 1)))
            True
        """
        if sym:
            if isinstance(sym[0], (int, Integer)):
                # a single symmetry is provided as a tuple or a range object;
                # it is converted to a 1-item list:
                sym = [tuple(sym)]
            sym_list = []
            for isym in sym:
                if len(isym) < 2:
                    raise IndexError("at least two index positions must be " +
                                     "provided to define a symmetry")
                for i in isym:
                    if i < 0 or i > nb_indices:
                        raise IndexError("invalid index position: " + str(i) +
                                         " not in [0," + str(nb_indices) + "]")
                sym_list.append(tuple(sorted(isym)))
            sym = tuple(sorted(sym_list))
        else:
            sym = ()

        if antisym:
            if isinstance(antisym[0], (int, Integer)):
                # a single antisymmetry is provided as a tuple or a range
                # object; it is converted to a 1-item list:
                antisym = [tuple(antisym)]
            antisym_list = []
            for isym in antisym:
                if len(isym) < 2:
                    raise IndexError("at least two index positions must be " +
                                     "provided to define an antisymmetry")
                for i in isym:
                    if i < 0 or i > nb_indices:
                        raise IndexError("invalid index position: " + str(i) +
                                         " not in [0," + str(self._nid-1) + "]")
                antisym_list.append(tuple(sorted(isym)))
            antisym = tuple(sorted(antisym_list))
        else:
            antisym = ()

        # Final consistency check:
        index_list = []
        for isym in sym:
            index_list += list(isym)
        for isym in antisym:
            index_list += list(isym)
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("incompatible lists of symmetries: the same " +
                             "index position appears more then once")

        if not sym and not antisym:
            return CompParent(nb_indices)

        if sym:
            if any(len(s) == nb_indices for s in sym):
                return CompParentFullySym(nb_indices)

        if antisym:
            if any(len(s) == nb_indices for s in antisym):
                return CompParentFullyAntiSym(nb_indices)

        return super(cls, CompParentWithSym).__classcall__(cls, nb_indices,
                                                           sym, antisym)

    Element = ComponentsWithSym

    def __init__(self, nb_indices, sym, antisym):
        r"""
        TESTS::


        """
        super().__init__(nb_indices)
        self._sym = sym
        self._antisym = antisym

    def _repr_symmetry(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentWithSym
            sage: cp = CompParent(2)
            sage: cp._repr_()
            'Parent of 2-index components without symmetry'

            sage: from sage.tensor.modules.comp_parent import CompParentWithSym
            sage: cp = CompParentWithSym(4, sym=(0,1))
            sage: cp._repr_()
            'Parent of 4-index components, with symmetry on the index positions (0, 1)'

            sage: cp = CompParentWithSym(4, sym=((0,1),), antisym=((2,3),))
            sage: cp._repr_()
            'Parent of 4-index components, with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)'
        """
        description = ""
        for isym in self._sym:
            description += ", with symmetry on the index positions " + \
                           str(tuple(isym))
        for isym in self._antisym:
            description += ", with antisymmetry on the index positions " + \
                           str(tuple(isym))
        return "", description

    @cached_method(key=lambda self, indices, ranges: (tuple(indices), tuple(ranges)))
    def _ordered_indices(self, indices, ranges):
        r"""
        Given a sequence of indices, return a sequence of indices with the indices
        at the positions of symmetries or antisymmetries being ordered,
        as well as some antisymmetry indicator.

        INPUT:

        - ``indices`` -- tuple of indices

        OUTPUT:

        - a pair ``(s,ind)`` where ``ind`` is a tuple that differs from the
          original list of indices by a reordering at the positions of
          symmetries and antisymmetries and

          * ``s = 0`` if the value corresponding to ``indices`` vanishes by
            antisymmetry (repeated indices); `ind` is then set to ``None``
          * ``s = 1`` if the value corresponding to ``indices`` is the same as
            that corresponding to ``ind``
          * ``s = -1`` if the value corresponding to ``indices`` is the
            opposite of that corresponding to ``ind``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParentWithSym
            sage: cp = CompParentWithSym(4, sym=((0,1),), antisym=((2,3),))
            sage: r = range(3)
            sage: ranges = (r, r, r, r)
            sage: cp._ordered_indices([0,1,1,2], ranges)
            (1, (0, 1, 1, 2))
            sage: cp._ordered_indices([1,0,1,2], ranges)
            (1, (0, 1, 1, 2))
            sage: cp._ordered_indices([0,1,2,1], ranges)
            (-1, (0, 1, 1, 2))
            sage: cp._ordered_indices([0,1,2,2], ranges)
            (0, None)

        """
        from sage.combinat.permutation import Permutation
        ind = list(self._check_indices(indices, ranges))
        for isym in self._sym:
            indsym = []
            for pos in isym:
                indsym.append(ind[pos])
            indsym_ordered = sorted(indsym)
            for k, pos in enumerate(isym):
                ind[pos] = indsym_ordered[k]
        sign = 1
        for isym in self._antisym:
            indsym = []
            for pos in isym:
                indsym.append(ind[pos])
            # Returns zero if some index appears twice:
            if len(indsym) != len(set(indsym)):
                return (0, None)
            # From here, all the indices in indsym are distinct and we need
            # to determine whether they form an even permutation of their
            # ordered series
            indsym_ordered = sorted(indsym)
            for k, pos in enumerate(isym):
                ind[pos] = indsym_ordered[k]
            if indsym_ordered != indsym:
                # Permutation linking indsym_ordered to indsym:
                #  (the +1 is required to fulfill the convention of Permutation)
                perm = [indsym.index(i) +1 for i in indsym_ordered]
                #c#     Permutation(perm).signature()
                sign *= Permutation(perm).signature()
        ind = tuple(ind)
        return (sign, ind)

    def non_redundant_index_generator(self, ranges):
        r"""
        Generator of indices, with only ordered indices in case of symmetries,
        so that only non-redundant indices are generated.

        OUTPUT:

        - an iterable index

        EXAMPLES:

        Indices on a 2-dimensional space::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentWithSym, \
            ....:  CompParentFullySym, CompParentFullyAntiSym
            sage: V = VectorSpace(QQ, 2)
            sage: r = range(V.dimension())

            sage: cp = CompParentFullySym(2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 0), (0, 1), (1, 1)]
            sage: r_start_1 = range(1, V.dimension() + 1)
            sage: list(cp.non_redundant_index_generator((r_start_1, r_start_1)))
            [(1, 1), (1, 2), (2, 2)]

            sage: cp = CompParentFullyAntiSym(2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 1)]

        Indices on a 3-dimensional space::

            sage: V = VectorSpace(QQ, 3)
            sage: r = range(V.dimension())

            sage: cp = CompParentFullySym(2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
            sage: r_start_1 = range(1, V.dimension() + 1)
            sage: list(cp.non_redundant_index_generator((r_start_1, r_start_1)))
            [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]

            sage: cp = CompParentFullyAntiSym(2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 1), (0, 2), (1, 2)]

            sage: cp = CompParentWithSym(3, sym=(1,2))  # symmetry on the last two indices
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2),
             (0, 2, 2), (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 1, 1),
             (1, 1, 2), (1, 2, 2), (2, 0, 0), (2, 0, 1), (2, 0, 2),
             (2, 1, 1), (2, 1, 2), (2, 2, 2)]

            sage: cp = CompParentWithSym(3, antisym=(1,2))  # antisymmetry on the last two indices
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 0, 1), (0, 0, 2), (0, 1, 2), (1, 0, 1), (1, 0, 2), (1, 1, 2),
             (2, 0, 1), (2, 0, 2), (2, 1, 2)]

            sage: cp = CompParentFullySym(3)
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2),
             (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]

            sage: cp = CompParentFullyAntiSym(3)
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 1, 2)]

        Indices on a 4-dimensional space::

            sage: V = VectorSpace(QQ, 4)
            sage: r = range(V.dimension())

            sage: cp = CompParent(1)
            sage: list(cp.non_redundant_index_generator((r,)))
            [(0,), (1,), (2,), (3,)]
            sage: cp = CompParentFullyAntiSym(2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: cp = CompParentFullyAntiSym(3)
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
            sage: cp = CompParentFullyAntiSym(4)
            sage: list(cp.non_redundant_index_generator((r, r, r, r)))
            [(0, 1, 2, 3)]
            sage: cp = CompParentFullyAntiSym(5)
            sage: list(cp.non_redundant_index_generator((r, r, r, r, r)))  # nothing since c is identically zero in this case (for 5 > 4)
            []

        """
        ind = [r.start for r in ranges]
        ind_end = [r.start for r in ranges]
        ind_end[0] = ranges[0].stop
        while ind != ind_end:
            ordered = True
            for isym in self._sym:
                for k in range(len(isym)-1):
                    if ind[isym[k+1]] < ind[isym[k]]:
                        ordered = False
                        break
            for isym in self._antisym:
                for k in range(len(isym)-1):
                    if ind[isym[k+1]] <= ind[isym[k]]:
                        ordered = False
                        break
            if ordered:
                yield tuple(ind)
            ret = 1
            for pos in range(self._nid -1, -1, -1):
                if ind[pos] != ranges[pos].stop - 1:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = ranges[pos].stop # end point reached
                    else:
                        ind[pos] = ranges[pos].start
                        ret = 1

    @cached_method
    def swap_adjacent_indices(self, pos1, pos2, pos3):
        r"""
        Swap two adjacent sets of indices.

        This method is essentially required to reorder the covariant and
        contravariant indices in the computation of a tensor product.

        The symmetries are preserved and the corresponding indices are adjusted
        consequently.

        INPUT:

        - ``pos1`` -- position of the first index of set 1 (with the convention
          position=0 for the first slot)
        - ``pos2`` -- position of the first index of set 2 = 1 + position of
          the last index of set 1 (since the two sets are adjacent)
        - ``pos3`` -- 1 + position of the last index of set 2

        OUTPUT:

        - Components with index set 1 permuted with index set 2.

        EXAMPLES:

        Swap of the index in position 0 with the pair of indices in position
        (1,2) in a set of components antisymmetric with respect to the indices
        in position (1,2)::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 3, antisym=(1,2))
            sage: c[0,0,1], c[0,0,2], c[0,1,2] = (1,2,3)
            sage: c[1,0,1], c[1,0,2], c[1,1,2] = (4,5,6)
            sage: c[2,0,1], c[2,0,2], c[2,1,2] = (7,8,9)
            sage: c[:]
            [[[0, 1, 2], [-1, 0, 3], [-2, -3, 0]],
             [[0, 4, 5], [-4, 0, 6], [-5, -6, 0]],
             [[0, 7, 8], [-7, 0, 9], [-8, -9, 0]]]
            sage: c1 = c.swap_adjacent_indices(0,1,3)
            sage: c.parent()._antisym   # c is antisymmetric with respect to the last pair of indices...
            ((1, 2),)
            sage: c1.parent()._antisym  #...while c1 is antisymmetric with respect to the first pair of indices
            ((0, 1),)
            sage: c[0,1,2]
            3
            sage: c1[1,2,0]
            3
            sage: c1[2,1,0]
            -3

        """
        # The symmetries:
        lpos = list(range(self._nid))
        new_lpos = lpos[:pos1] + lpos[pos2:pos3] + lpos[pos1:pos2] + lpos[pos3:]
        res_sym = []
        for s in self._sym:
            new_s = [new_lpos.index(pos) for pos in s]
            res_sym.append(tuple(sorted(new_s)))
        res_antisym = []
        for s in self._antisym:
            new_s = [new_lpos.index(pos) for pos in s]
            res_antisym.append(tuple(sorted(new_s)))
        return CompParentWithSym(self._nid, sym=res_sym,
                                 antisym=res_antisym)

    @cached_method
    def common_symmetries(self, other):
        r"""
        Return a collection of components with common symmetries.
        """
        common_sym = []
        for isym in self._sym:
            for osym in other._sym:
                com = tuple(set(isym).intersection(set(osym)))
                if len(com) > 1:
                    common_sym.append(com)
        common_antisym = []
        for isym in self._antisym:
            for osym in other._antisym:
                com = tuple(set(isym).intersection(set(osym)))
                if len(com) > 1:
                    common_antisym.append(com)
        if common_sym != [] or common_antisym != []:
            # convert to tuples
            result = CompParentWithSym(self._nid,
                                       sym=common_sym, antisym=common_antisym)
        else:
            # no common symmetry -> result is collection of generic components:
            result = CompParent(self._nid)

        return result

    @cached_method
    def tensor_product_sym(self, other):
        r"""
        Return the collection of components with the symmetries of tensor
        products of components in ``self`` with components in ``other``.
        """
        sym = list(self._sym)
        antisym = list(self._antisym)
        if isinstance(other, CompParentWithSym):
            if other._sym != []:
                for s in other._sym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    sym.append(ns)
            if other._antisym != []:
                for s in other._antisym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    antisym.append(ns)
        return CompParentWithSym(self._nid,
                                 sym=sym, antisym=antisym)

    @cached_method
    def contract_sym(self, pos1, pos2):
        r"""
        Return the collection of components built from contracting ``pos1``
        with ``pos2``.
        """
        pass

class CompParentFullySym(CompParentWithSym):
    r"""

    """

    Element = ComponentsFullySym

    def __init__(self, nb_indices):
        r"""

        """
        CompParentWithSym.__init__(self, nb_indices,
                                   sym=(tuple(range(nb_indices)),),
                                   antisym=())

    def _repr_symmetry(self):
        return "Fully symmetric ", ""

class CompParentFullyAntiSym(CompParentWithSym):
    r"""

    """

    Element = ComponentsFullyAntiSym

    def __init__(self, nb_indices):
        r"""

        """
        CompParentWithSym.__init__(self, nb_indices,
                                   sym=(),
                                   antisym=(tuple(range(nb_indices)),))

    def _repr_symmetry(self):
        return "Fully antisymmetric ", ""
