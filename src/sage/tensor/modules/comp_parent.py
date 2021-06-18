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

import operator

from sage.structure.parent import Parent
from sage.structure.coerce_actions import ActedUponAction
from sage.categories.action import PrecomposedAction
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import get_coercion_model
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.modules.module import Module
from sage.categories.all import FreeModules
from .comp_element import Components_base
from .comp_element_dict import (Components_dict, ComponentsWithSym_dict,
                                ComponentsFullySym_dict, ComponentsFullyAntiSym_dict)

class CompParent(Module, UniqueRepresentation):
    r"""

    """

    Element = Components_dict

    def __init__(self, base_ring, nb_indices, category=None):
        r"""

        """
        if category is None:
            category = FreeModules(base_ring.category()).WithBasis()
        Module.__init__(self, base_ring, category=category)
        self._nid = nb_indices

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: C = CompParent(QQ, 2)
            sage: C._repr_()
            'Parent of 2-index components over Rational Field'

        """
        prefix, suffix = self._repr_symmetry()
        return f"Parent of " + prefix + f"{self._nid}-index components over {self.base_ring()}" + suffix

    def _repr_symmetry(self):
        r"""
        Return a prefix and a suffix string describing the symmetry of ``self``.
        """
        return "", ""

    def sym_antisym(self):
        return (), ()

    def _element_constructor_(self, x, **kwargs):
        r"""
        Construct an indexed set of components w.r.t. ``frame``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import (
            ....:     CompParent, CompParentWithSym, CompParentFullySym, CompParentFullyAntiSym)
            sage: cp123 = CompParentWithSym(QQ, 5, sym=(1, 2, 3)); cp123
            Parent of 5-index components over Rational Field, with symmetry on the index positions (1, 2, 3)
            sage: t = cp123(frame=range(6)); t
            5-index components w.r.t. range(0, 6), with symmetry on the index positions (1, 2, 3)
            sage: t[0, 4, 3, 5, 0] = 17
            sage: t.display('t')
            t_03450 = 17
            t_03540 = 17
            t_04350 = 17
            t_04530 = 17
            t_05340 = 17
            t_05430 = 17
            sage: cp = cp123.ambient(); cp
            Parent of 5-index components over Rational Field
            sage: lifted_t = cp(t); lifted_t
            5-index components w.r.t. range(0, 6)
            sage: lifted_t.display('lifted_t')
            lifted_t_03450 = 17
            lifted_t_03540 = 17
            lifted_t_04350 = 17
            lifted_t_04530 = 17
            lifted_t_05340 = 17
            lifted_t_05430 = 17
        """
        items = ()
        if isinstance(x, Components_base):
            if x.parent() is self:
                return +x  # make a copy
            # Coerce to our element class
            items = ((ind, x[[ind]])
                     for ind in self.non_redundant_index_generator(x._ranges)
                     if x[[ind]])
            frame = x._frame
            start_index = x._sindex
            output_formatter = x._output_formatter
        else:
            frame = kwargs.pop('frame', ()) # () is for coercion from 0
            start_index = kwargs.pop('start_index', 0)
            output_formatter = kwargs.pop('output_formatter', None)
            if x == 0:
                # coercion from 0
                pass
            else:
                raise TypeError('cannot construct Components from {x}')
        result = self.element_class(self, frame=frame,
                                    start_index=start_index,
                                    output_formatter=output_formatter)
        for ind, value in items:
            result[[ind]] = value
        return result

    def _check_indices(self, ind, ranges):
        r"""
        Check the validity of a list of indices and returns a tuple from it

        INPUT:

        - ``ind`` -- tuple of indices

        OUTPUT:

        - a tuple containing valid indices

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent
            sage: cp = CompParent(SR, 2)
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
            sage: cp = CompParent(QQ, 1)
            sage: ranges = (r,)
            sage: list(cp.index_generator(ranges))
            [(0,), (1,), (2,)]
            sage: ranges = (range(1, 4),)
            sage: list(cp.index_generator(ranges))
            [(1,), (2,), (3,)]
            sage: cp = CompParent(QQ, 2)
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
            sage: cp = CompParent(QQ, 2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]
            sage: r_start_1 = range(1, V.dimension() + 1)
            sage: list(cp.non_redundant_index_generator((r_start_1, r_start_1)))
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1),
             (3, 2), (3, 3)]
        """
        yield from self.index_generator(ranges)

    def _coerce_map_from_(self, other_parent):
        """
        Check whether a coercion from ``other_parent`` to ``self`` is possible.

        TESTS::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentWithSym
            sage: cp = CompParent(QQ, 2)
            sage: cp.coerce_map_from(CompParent(ZZ, 2))
            Coercion map:
              From: Parent of 2-index components over Integer Ring
              To:   Parent of 2-index components over Rational Field

            sage: cp.coerce_map_from(CompParentWithSym(ZZ, 2, sym=(1, 2)))
            Coercion map:
            From: Parent of Fully symmetric 2-index components over Integer Ring
            To:   Parent of 2-index components over Rational Field
        """
        # Version for parent with no symmetries.  TODO: Split out an abstract base class.
        if (isinstance(other_parent, CompParent)
            and self.base_ring().has_coerce_map_from(other_parent.base_ring())
            and self._nid == other_parent._nid):
            if self.common_symmetries(other_parent):
                return True
        return super()._coerce_map_from_(other_parent)

    @cached_method
    def symmetrize(self, *pos):
        r"""
        Return the parent of components in ``self`` that are symmetric in ``pos``.

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
            return CompParentFullySym(self.base_ring(), self._nid)
        return CompParentWithSym(self.base_ring(), self._nid, sym=pos)

    @cached_method
    def antisymmetrize(self, *pos):
        r"""
        Return the parent of components in ``self`` that are antisymmetric in ``pos``.

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
            return CompParentFullyAntiSym(self.base_ring(), self._nid)
        return CompParentWithSym(self.base_ring(), self._nid, antisym=pos)

    @cached_method
    def common_symmetries(self, other):
        if isinstance(other, CompParentWithSym):
            return other.common_symmetries(self)
        return self

    @cached_method
    def __add__(self, other):
        r"""
        Sum (= pushout) of two parents

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import (
            ....:     CompParent, CompParentWithSym, CompParentFullySym, CompParentFullyAntiSym)
            sage: from sage.categories.pushout import pushout
            sage: A1 = CompParentWithSym(QQ, 5, sym=(1, 2, 3))
            sage: A2 = CompParentWithSym(QQ, 5, sym=(2, 3, 4))
            sage: A1 + A2
            Parent of 5-index components over Rational Field, with symmetry on the index positions (2, 3)
            sage: pushout(A1, A2)
            Parent of 5-index components over Rational Field, with symmetry on the index positions (2, 3)
        """
        if not (isinstance(other, CompParent)
                and self._nid == other._nid
                and self.base_ring() == other.base_ring()):
            raise TypeError('can only add CompParent instances with the same number of indices and the same base ring')
        return self.common_symmetries(other)

    _pushout_ = __add__

    @cached_method
    def tensor(*parents, **kwargs):
        r"""
        Return the tensor product

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import (
            ....:     CompParent, CompParentWithSym, CompParentFullySym, CompParentFullyAntiSym)
            sage: A = CompParentWithSym(QQ, 3)
            sage: B = CompParentFullySym(ZZ, 2)
            sage: C = CompParentFullyAntiSym(QQ, 3)
            sage: A.tensor(B, C)
            Parent of 8-index components over Rational Field,
             with symmetry on the index positions (3, 4),
             with antisymmetry on the index positions (5, 6, 7)
        """
        tp_sym = []
        tp_antisym = []
        tp_nid = 0
        for p in parents:
            sym, antisym = p.sym_antisym()
            tp_sym.extend(tuple(i + tp_nid for i in clique)
                          for clique in sym)
            tp_antisym.extend(tuple(i + tp_nid for i in clique)
                              for clique in antisym)
            tp_nid += p._nid
        cm = get_coercion_model()
        tp_base_ring = cm.common_parent(*[p.base_ring() for p in parents])
        return CompParentWithSym(tp_base_ring, tp_nid,
                                 sym=tp_sym, antisym=tp_antisym)


class CompParentWithSym(CompParent):
    r"""

    """

    @staticmethod
    def __classcall_private__(cls, base_ring, nb_indices, sym=None, antisym=None):
        r"""
        Determine the correct class to return based upon the input. In
        particular, convert lists of symmetries to tuples and delegate to
        parent with correct symmetry.

        TESTS::

            sage: from sage.tensor.modules.comp_parent import (
            ....:     CompParent, CompParentWithSym, CompParentFullySym, CompParentFullyAntiSym)
            sage: CompParentWithSym(QQ, 3) is CompParent(QQ, 3)
            True
            sage: CompParentWithSym(SR, 5, sym=(3, 4),) is CompParentWithSym(SR, 5, sym=((4, 3),))
            True
            sage: CompParentWithSym(SR, 5, sym=((1, 2), (3, 4))) is CompParentWithSym(SR, 5, sym=((4, 3), (2, 1)))
            True
            sage: CompParentWithSym(RDF, 3, sym=((1, 3, 2))) is CompParentFullySym(RDF, 3)
            True
            sage: CompParentWithSym(CDF, 3, antisym=((2, 1, 3))) is CompParentFullyAntiSym(CDF, 3)
            True
            sage: CompParentWithSym(ZZ, 0) is CompParent(ZZ, 0)
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
            return CompParent(base_ring, nb_indices)

        if sym:
            if any(len(s) == nb_indices for s in sym):
                return CompParentFullySym(base_ring, nb_indices)

        if antisym:
            if any(len(s) == nb_indices for s in antisym):
                return CompParentFullyAntiSym(base_ring, nb_indices)

        return super(cls, CompParentWithSym).__classcall__(cls, base_ring, nb_indices,
                                                           sym, antisym)

    Element = ComponentsWithSym_dict

    def __init__(self, base_ring, nb_indices, sym, antisym):
        r"""
        TESTS::


        """
        category = FreeModules(base_ring.category()).WithBasis().Subobjects()
        super().__init__(base_ring, nb_indices, category=category)
        self._sym = sym
        self._antisym = antisym

    @cached_method
    def ambient(self):
        r"""
        Return the enclosing ``CompParent`` (without symmetries).

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParentWithSym
            sage: cp = CompParentWithSym(SR, 4, sym=(0,1)); cp
            Parent of 4-index components over Symbolic Ring, with symmetry on the index positions (0, 1)
            sage: cp.ambient()
            Parent of 4-index components over Symbolic Ring
        """
        return CompParent(self.base_ring(), self._nid)

    def sym_antisym(self):
        return self._sym, self._antisym

    def _repr_symmetry(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentWithSym
            sage: cp = CompParent(SR, 2)
            sage: cp._repr_()
            'Parent of 2-index components over Symbolic Ring'

            sage: from sage.tensor.modules.comp_parent import CompParentWithSym
            sage: cp = CompParentWithSym(SR, 4, sym=(0,1))
            sage: cp._repr_()
            'Parent of 4-index components over Symbolic Ring, with symmetry on the index positions (0, 1)'

            sage: cp = CompParentWithSym(SR, 4, sym=((0,1),), antisym=((2,3),))
            sage: cp._repr_()
            'Parent of 4-index components over Symbolic Ring, with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)'
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
            sage: cp = CompParentWithSym(SR, 4, sym=((0,1),), antisym=((2,3),))
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

            sage: cp = CompParentFullySym(SR, 2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 0), (0, 1), (1, 1)]
            sage: r_start_1 = range(1, V.dimension() + 1)
            sage: list(cp.non_redundant_index_generator((r_start_1, r_start_1)))
            [(1, 1), (1, 2), (2, 2)]

            sage: cp = CompParentFullyAntiSym(SR, 2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 1)]

        Indices on a 3-dimensional space::

            sage: V = VectorSpace(QQ, 3)
            sage: r = range(V.dimension())

            sage: cp = CompParentFullySym(SR, 2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
            sage: r_start_1 = range(1, V.dimension() + 1)
            sage: list(cp.non_redundant_index_generator((r_start_1, r_start_1)))
            [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]

            sage: cp = CompParentFullyAntiSym(SR, 2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 1), (0, 2), (1, 2)]

            sage: cp = CompParentWithSym(SR, 3, sym=(1,2))  # symmetry on the last two indices
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2),
             (0, 2, 2), (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 1, 1),
             (1, 1, 2), (1, 2, 2), (2, 0, 0), (2, 0, 1), (2, 0, 2),
             (2, 1, 1), (2, 1, 2), (2, 2, 2)]

            sage: cp = CompParentWithSym(SR, 3, antisym=(1,2))  # antisymmetry on the last two indices
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 0, 1), (0, 0, 2), (0, 1, 2), (1, 0, 1), (1, 0, 2), (1, 1, 2),
             (2, 0, 1), (2, 0, 2), (2, 1, 2)]

            sage: cp = CompParentFullySym(SR, 3)
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2),
             (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]

            sage: cp = CompParentFullyAntiSym(SR, 3)
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 1, 2)]

        Indices on a 4-dimensional space::

            sage: V = VectorSpace(QQ, 4)
            sage: r = range(V.dimension())

            sage: cp = CompParent(SR, 1)
            sage: list(cp.non_redundant_index_generator((r,)))
            [(0,), (1,), (2,), (3,)]
            sage: cp = CompParentFullyAntiSym(SR, 2)
            sage: list(cp.non_redundant_index_generator((r, r)))
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: cp = CompParentFullyAntiSym(SR, 3)
            sage: list(cp.non_redundant_index_generator((r, r, r)))
            [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
            sage: cp = CompParentFullyAntiSym(SR, 4)
            sage: list(cp.non_redundant_index_generator((r, r, r, r)))
            [(0, 1, 2, 3)]
            sage: cp = CompParentFullyAntiSym(SR, 5)
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
        return CompParentWithSym(self.base_ring(), self._nid,
                                 sym=res_sym, antisym=res_antisym)

    @cached_method
    def common_symmetries(self, other):
        r"""
        Return a collection of components with common symmetries.
        """
        if not isinstance(other, CompParentWithSym):
            return self
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
        return CompParentWithSym(self.base_ring(), self._nid,
                                 sym=common_sym, antisym=common_antisym)

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

    @staticmethod
    def __classcall_private__(cls, base_ring, nb_indices):
        r"""
        Determine the correct class to return based upon the input.

        TESTS::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentFullySym
            sage: CompParentFullySym(SR, 0) is CompParent(SR, 0)
            True
            sage: CompParentFullySym(SR, 1) is CompParent(SR, 1)
            True
            sage: CompParentFullySym(SR, 2) is CompParent(SR, 2)
            False
        """
        if nb_indices <= 1:
            return CompParent(base_ring, nb_indices)
        return super(cls, CompParentFullySym).__classcall__(cls, base_ring, nb_indices)

    Element = ComponentsFullySym_dict

    def __init__(self, base_ring, nb_indices):
        r"""

        """
        super().__init__(base_ring, nb_indices,
                         sym=(tuple(range(nb_indices)),),
                         antisym=())

    def _repr_symmetry(self):
        return "Fully symmetric ", ""

class CompParentFullyAntiSym(CompParentWithSym):
    r"""

    """

    @staticmethod
    def __classcall_private__(cls, base_ring, nb_indices):
        r"""
        Determine the correct class to return based upon the input.

        TESTS::

            sage: from sage.tensor.modules.comp_parent import CompParent, CompParentFullyAntiSym
            sage: CompParentFullyAntiSym(SR, 0) is CompParent(SR, 0)
            True
            sage: CompParentFullyAntiSym(SR, 2) is CompParent(SR, 2)
            False
        """
        if nb_indices == 0:
            return CompParent(base_ring, nb_indices)
        if nb_indices == 1:
            # TODO: A special-case parent for the trivial "fully antisymmetric 1-index tensors"
            # could be helpful here
            pass
        return super(cls, CompParentFullyAntiSym).__classcall__(cls, base_ring, nb_indices)

    Element = ComponentsFullyAntiSym_dict

    def __init__(self, base_ring, nb_indices):
        r"""

        """
        super().__init__(base_ring, nb_indices,
                         sym=(),
                         antisym=(tuple(range(nb_indices)),))

    def _repr_symmetry(self):
        return "Fully antisymmetric ", ""
