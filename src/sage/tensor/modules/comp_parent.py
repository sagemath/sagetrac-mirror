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
from .comp_element import Components_generic, CompWithSym, CompFullySym, CompFullyAntiSym
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer

class CompParent(Parent, UniqueRepresentation):
    r"""

    """

    Element = Components_generic

    @staticmethod
    def __classcall_private__(cls, dim, nb_indices, start_index=0):
        r"""
        Determine the correct class to return based upon the input.

        TODO: Remove this method as soon as :trac:`14167` is solved.

        TESTS::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: CompParent(ZZ, 4, 2) is CompParent(ZZ, 4, 2, start_index=0)
            True

        """
        return super(cls, CompParentWithSym).__classcall__(cls, dim,
                                            nb_indices, start_index=start_index)

    def __init__(self, dim, nb_indices, start_index=0):
        r"""

        """
        Parent.__init__(self)
        self._dim = dim
        self._nid = nb_indices
        self._sindex = start_index

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: C = CompParent(3, 2)
            sage: C._repr_()

        """
        description = "Collection of "
        description += str(self._nid)
        if self._nid == 1:
            description += "-index"
        else:
            description += "-indices"
        description += " components"
        return description

    def _element_constructor_(self, *args, **kwargs):
        r"""
        Construct an indexed set of components w.r.t. ``frame`` over the ring
        ``ring``.
        """
        if isinstance(args[0], Components_generic):
            c = args[0]
            if not isinstance(c, Components_generic):
                TypeError("cannot coerce {} into an element of {}".format(c, self))
            else:
                return c
        else:
            if len(args) != 2:
                ValueError("{} missing {} required positional arguments".format(type(self), len(args)))
            ring = args[0]
            frame = args[1]
            output_formatter = kwargs.pop('output_formatter', None)
        return self.element_class(self, ring, frame,
                                  output_formatter=output_formatter)

    @cached_method
    def _check_indices(self, ind):
        r"""
        Check the validity of a list of indices and returns a tuple from it

        INPUT:

        - ``ind`` -- tuple of indices

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
            IndexError: index out of range: 3 not in [0, 2]
            sage: c._check_indices(1)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 1 are provided
            sage: c._check_indices([1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 3 are provided

        """
        if len(ind) != self._nid:
            raise ValueError(("wrong number of indices: {} expected,"
                             " while {} are provided").format(self._nid, len(ind)))
        si = self._sindex
        imax = self._dim - 1 + si
        for k in range(self._nid):
            i = ind[k]
            if i < si or i > imax:
                raise IndexError("index out of range: " +
                                 "{} not in [{}, {}]".format(i, si, imax))
        return ind

    @cached_method
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
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(self._nid)]
        ind_end = [si for k in range(self._nid)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(self._nid-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    @cached_method
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
        for ind in self.index_generator():
            yield ind

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
            return CompParentFullySym(self._nid, start_index=self._sindex)
        return CompParentWithSym(self._nid, start_index=self._sindex, sym=pos)

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
            return CompParentFullyAntiSym(self._nid, start_index=self._sindex)
        return CompParentWithSym(self._nid, start_index=self._sindex,
                                 antisym=pos)

class CompParentWithSym(CompParent):
    r"""

    """

    @staticmethod
    def __classcall_private__(cls, dim, nb_indices, start_index=0, sym=None,
                              antisym=None):
        r"""
        Determine the correct class to return based upon the input. In
        particular, convert lists of symmetries to tuples and delegate to
        parent with correct symmetry.

        TESTS:



        """
        if sym is None and antisym is None:
            raise TypeError("for components without symmetries, please use "
                            "'CompParent' instead")
        if sym is not None:
            if any(len(s) == nb_indices for s in sym):
                return CompParentFullySym(dim, nb_indices,
                                          start_index=start_index)
        if antisym is not None:
            if any(len(s) == nb_indices for s in antisym):
                return CompParentFullyAntiSym(dim, nb_indices,
                                              start_index=start_index)
        sym = tuple(sorted(sym))
        antisym = tuple(sorted(antisym))
        return super(cls, CompParentWithSym).__classcall__(cls, dim, nb_indices,
                                                       start_index=start_index,
                                                       sym=sym, antisym=antisym)

    Element = CompWithSym

    def __init__(self, dim, nb_indices, start_index=0, sym=None,
                 antisym=None):
        r"""

        """
        CompParent.__init__(self, dim, nb_indices, start_index)
        self._sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):
                # a single symmetry is provided as a tuple or a range object;
                # it is converted to a 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) < 2:
                    raise IndexError("at least two index positions must be " +
                                     "provided to define a symmetry")
                for i in isym:
                    if i<0 or i>self._nid-1:
                        raise IndexError("invalid index position: " + str(i) +
                                         " not in [0," + str(self._nid-1) + "]")
                self._sym.append(tuple(isym))
        self._antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):
                # a single antisymmetry is provided as a tuple or a range
                # object; it is converted to a 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) < 2:
                    raise IndexError("at least two index positions must be " +
                                     "provided to define an antisymmetry")
                for i in isym:
                    if i<0 or i>self._nid-1:
                        raise IndexError("invalid index position: " + str(i) +
                                         " not in [0," + str(self._nid-1) + "]")
                self._antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self._sym:
            index_list += isym
        for isym in self._antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("incompatible lists of symmetries: the same " +
                             "index position appears more then once")

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: C = CompParent(3, 2)
            sage: C._repr_()

            sage: from sage.tensor.modules.comp_parent import CompParentSym
            sage: CompParentSym(3, 4, sym=(0,1))

            sage: CompWithSym(3, 4, sym=(0,1), antisym=(2,3))


        """
        description = "Collection of "
        description += str(self._nid)
        if self._nid == 1:
            description += "-index"
        else:
            description += "-indices"
        description += " components"
        for isym in self._sym:
            description += ", with symmetry on the index positions " + \
                           str(tuple(isym))
        for isym in self._antisym:
            description += ", with antisymmetry on the index positions " + \
                           str(tuple(isym))
        return description

    @cached_method
    def _ordered_indices(self, indices):
        r"""
        Given a set of indices, return a set of indices with the indices
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

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            sage: c._ordered_indices([0,1,1,2])
            (1, (0, 1, 1, 2))
            sage: c._ordered_indices([1,0,1,2])
            (1, (0, 1, 1, 2))
            sage: c._ordered_indices([0,1,2,1])
            (-1, (0, 1, 1, 2))
            sage: c._ordered_indices([0,1,2,2])
            (0, None)

        """
        from sage.combinat.permutation import Permutation
        ind = list(self._check_indices(indices))
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

    @cached_method
    def non_redundant_index_generator(self):
        r"""
        Generator of indices, with only ordered indices in case of symmetries,
        so that only non-redundant indices are generated.

        OUTPUT:

        - an iterable index

        EXAMPLES:

        Indices on a 2-dimensional space::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:  CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 2)
            sage: c = CompFullySym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 0), (0, 1), (1, 1)]
            sage: c = CompFullySym(QQ, V.basis(), 2, start_index=1)
            sage: list(c.non_redundant_index_generator())
            [(1, 1), (1, 2), (2, 2)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 1)]

        Indices on a 3-dimensional space::

            sage: V = VectorSpace(QQ, 3)
            sage: c = CompFullySym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
            sage: c = CompFullySym(QQ, V.basis(), 2, start_index=1)
            sage: list(c.non_redundant_index_generator())
            [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 1), (0, 2), (1, 2)]
            sage: c = CompWithSym(QQ, V.basis(), 3, sym=(1,2))  # symmetry on the last two indices
            sage: list(c.non_redundant_index_generator())
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2),
             (0, 2, 2), (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 1, 1),
             (1, 1, 2), (1, 2, 2), (2, 0, 0), (2, 0, 1), (2, 0, 2),
             (2, 1, 1), (2, 1, 2), (2, 2, 2)]
            sage: c = CompWithSym(QQ, V.basis(), 3, antisym=(1,2))  # antisymmetry on the last two indices
            sage: list(c.non_redundant_index_generator())
            [(0, 0, 1), (0, 0, 2), (0, 1, 2), (1, 0, 1), (1, 0, 2), (1, 1, 2),
             (2, 0, 1), (2, 0, 2), (2, 1, 2)]
            sage: c = CompFullySym(QQ, V.basis(), 3)
            sage: list(c.non_redundant_index_generator())
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2),
             (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: list(c.non_redundant_index_generator())
            [(0, 1, 2)]

        Indices on a 4-dimensional space::

            sage: V = VectorSpace(QQ, 4)
            sage: c = Components(QQ, V.basis(), 1)
            sage: list(c.non_redundant_index_generator())
            [(0,), (1,), (2,), (3,)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: list(c.non_redundant_index_generator())
            [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: list(c.non_redundant_index_generator())
            [(0, 1, 2, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 5)
            sage: list(c.non_redundant_index_generator())  # nothing since c is identically zero in this case (for 5 > 4)
            []

        """
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(self._nid)]
        ind_end = [si for k in range(self._nid)]
        ind_end[0] = imax+1
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
            for pos in range(self._nid-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
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

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: V = VectorSpace(QQ, 3)
            sage: c = CompWithSym(QQ, V.basis(), 3, antisym=(1,2))
            sage: c[0,0,1], c[0,0,2], c[0,1,2] = (1,2,3)
            sage: c[1,0,1], c[1,0,2], c[1,1,2] = (4,5,6)
            sage: c[2,0,1], c[2,0,2], c[2,1,2] = (7,8,9)
            sage: c[:]
            [[[0, 1, 2], [-1, 0, 3], [-2, -3, 0]],
             [[0, 4, 5], [-4, 0, 6], [-5, -6, 0]],
             [[0, 7, 8], [-7, 0, 9], [-8, -9, 0]]]
            sage: c1 = c.swap_adjacent_indices(0,1,3)
            sage: c._antisym   # c is antisymmetric with respect to the last pair of indices...
            [(1, 2)]
            sage: c1._antisym  #...while c1 is antisymmetric with respect to the first pair of indices
            [(0, 1)]
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
        return CompParentWithSym(self._dim, self._nid,
                                 start_index=self._sindex, sym=res_sym,
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
            result = CompParentWithSym(self._dim, self._nid,
                                       start_index=self._sindex,
                                       sym=common_sym, antisym=common_antisym)
        else:
            # no common symmetry -> result is collection of generic components:
            result = CompParent(self._dim, self._nid, start_index=self._sindex)

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
        return CompParentWithSym(self.self._dim, self._nid,
                                 start_index=self._sindex, sym=sym,
                                 antisym=antisym)

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

    Element = CompFullySym

    def __init__(self, dim, nb_indices, start_index=0):
        r"""

        """
        CompParentWithSym.__init__(self, dim, nb_indices,
                                   start_index=start_index,
                                   sym=range(nb_indices))

class CompParentFullyAntiSym(CompParentWithSym):
    r"""

    """

    Element = CompFullyAntiSym

    def __init__(self, dim, nb_indices, start_index=0):
        r"""

        """
        CompParentWithSym.__init__(self, dim, nb_indices,
                                   start_index=start_index,
                                   antisym=range(nb_indices))