r"""
Implementation of Components using Dictionaries
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

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.parallel.decorate import parallel
from sage.parallel.parallelism import Parallelism
from operator import itemgetter
from .comp_element import Components_base

class Components_dict(Components_base):

    def __init__(self, parent, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import Components
            sage: Components(ZZ, [1,2,3], 2)
            2-index components w.r.t. [1, 2, 3]

        """
        super().__init__(parent, frame, start_index, output_formatter)
        self._comp = {} # the dictionary of components, with the index tuples
                        # as keys

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        EXAMPLES:

        Copy of a set of components with a single index::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = -2, 1, 5
            sage: b = a.copy() ; b
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: b[:]
            [-2, 1, 5]
            sage: b == a
            True
            sage: b is a  # b is a distinct object
            False

        """
        result = self._new_instance()
        for ind, val in self._comp.items():
            if isinstance(val, SageObject) and hasattr(val, 'copy'):
                result._comp[ind] = val.copy()
            else:
                result._comp[ind] = val
        return result

    def _del_zeros(self):
        r"""
        Deletes all the zeros in the dictionary :attr:`_comp`

        NB: The use case of this method must be rare because zeros are not
            stored in :attr:`_comp`.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._comp = {(0,1): 3, (0,2): 0, (1,2): -5, (2,2): 0}  # enforcing zero storage
            sage: c._del_zeros()
            sage: c._comp
            {(0, 1): 3, (1, 2): -5}

        """
        # The zeros are first searched; they are deleted in a second stage, to
        # avoid changing the dictionary while it is read
        zeros = []
        for ind, value in self._comp.items():
            if value == 0:
                zeros.append(ind)
        for ind in zeros:
            del self._comp[ind]

    def __getitem__(self, args):
        r"""
        Returns the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the
          components `T_{ij...}` (for a 2-index object, a matrix is returned)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c[1,2]    # unset components are zero
            0
            sage: c.__getitem__((1,2))
            0
            sage: c.__getitem__([1,2])
            0
            sage: c[1,2] = -4
            sage: c[1,2]
            -4
            sage: c.__getitem__((1,2))
            -4
            sage: c[:]
            [ 0  0  0]
            [ 0  0 -4]
            [ 0  0  0]
            sage: c.__getitem__(slice(None))
            [ 0  0  0]
            [ 0  0 -4]
            [ 0  0  0]

        """
        no_format = self._output_formatter is None
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            no_format = True
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            return self._get_list(indices, no_format, format_type)
        else:
            ind = self._check_indices(indices)
            if ind in self._comp:
                if no_format:
                    return self._comp[ind]
                elif format_type is None:
                    return self._output_formatter(self._comp[ind])
                else:
                    return self._output_formatter(self._comp[ind], format_type)
            else:  # if the value is not stored in self._comp, it is zero:
                if no_format:
                    return self.base_ring().zero()
                elif format_type is None:
                    return self._output_formatter(self.base_ring().zero())
                else:
                    return self._output_formatter(self.base_ring().zero(),
                                                 format_type)

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object); if ``[:]`` is provided, all the
          components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]`` (``slice(None)``)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c.__setitem__((0,1), -4)
            sage: c[:]
            [ 0 -4  0]
            [ 0  0  0]
            [ 0  0  0]
            sage: c[0,1] = -4
            sage: c[:]
            [ 0 -4  0]
            [ 0  0  0]
            [ 0  0  0]
            sage: c.__setitem__(slice(None), [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: c[:]
            [0 1 2]
            [3 4 5]
            [6 7 8]

        """
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            self._set_list(indices, format_type, value)
        else:
            ind = self._check_indices(indices)
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if zero_value:
                # if the component has been set previously, it is deleted,
                # otherwise nothing is done (zero components are not stored):
                if ind in self._comp:
                    del self._comp[ind]
            else:
                if format_type is None:
                    self._comp[ind] = self.base_ring()(value)
                else:
                    self._comp[ind] = self.base_ring()({format_type: value})
                    # NB: the writing
                    #   self._comp[ind] = self._ring(value, format_type)
                    # is not allowed when ring is an algebra and value some
                    # element of the algebra's base ring, cf. the discussion at
                    # http://trac.sagemath.org/ticket/16054

    def is_zero(self):
        r"""
        Return ``True`` if all the components are zero and ``False`` otherwise.

        EXAMPLES:

        A just-created set of components is initialized to zero::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: c = Components(QQ, V.basis(), 1)
            sage: c.is_zero()
            True
            sage: c[:]
            [0, 0, 0]
            sage: c[0] = 1 ; c[:]
            [1, 0, 0]
            sage: c.is_zero()
            False
            sage: c[0] = 0 ; c[:]
            [0, 0, 0]
            sage: c.is_zero()
            True

        It is equivalent to use the operator == to compare to zero::

            sage: c == 0
            True
            sage: c != 0
            False

        Comparing to a nonzero number is meaningless::

            sage: c == 1
            Traceback (most recent call last):
            ...
            TypeError: cannot compare a set of components to a number

        """
        if not self._comp:
            return True

        #!# What follows could be skipped since _comp should not contain
        # any zero value
        # In other words, the full method should be
        #   return self.comp == {}
        for val in self._comp.values():
            if not (val == 0):
                return False
        return True

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the components represented by ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 1)
            sage: c[:] = 5, 0, -4
            sage: a = c.__neg__() ; a
            1-index components w.r.t. [1, 2, 3]
            sage: a[:]
            [-5, 0, 4]
            sage: a == -c
            True

        """
        result = self._new_instance()
        for ind, val in self._comp.items():
             result._comp[ind] = - val
        return result

    def _add_(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT:

        - components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__add__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [5, 5, 3]
            sage: s == a+b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = a.__add__(b) ; s_par
            1-index components w.r.t. [1, 2, 3]
            sage: s_par[:]
            [5, 5, 3]
            sage: s_par == s
            True
            sage: b.__add__(a) == s  # test of commutativity of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components_base):
            raise TypeError("the second argument for the addition must be " +
                            "an instance of Components")
        if isinstance(other, ComponentsWithSym_dict):
            return other + self     # to deal properly with symmetries
        if not other._frame:
            return +self
        if not self._frame:
            return +other
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._nid != self._nid:
            raise ValueError("the two sets of components do not have the " +
                             "same number of indices")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        # Initialization of the result to self.copy(), so that there remains
        # only to add other:
        result = self.copy()
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = [ind for ind in other._comp]
            ind_step = max(1, int(len(ind_list)/nproc/2))
            local_list = lol(ind_list, ind_step)
            # list of input parameters
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_sum(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    partial.append([ind, a[[ind]]+b[[ind]]])
                return partial

            for ii, val in paral_sum(listParalInput):
                for jj in val:
                    result[[jj[0]]] = jj[1]

        else:
            # Sequential computation
            for ind, val in other._comp.items():
                result[[ind]] += val

        return result

    def __mul__(self, other):
        r"""
        Component tensor product.

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT:

        - the tensor product of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__mul__(b) ; s
            2-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [  4   5   6]
            [  0   0   0]
            [-12 -15 -18]
            sage: s == a*b
            True
            sage: t = b*a
            sage: aa = a*a ; aa
            Fully symmetric 2-index components w.r.t. [1, 2, 3]
            sage: aa[:]
            [ 1  0 -3]
            [ 0  0  0]
            [-3  0  9]

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__mul__(b) ; s_par
            2-index components w.r.t. [1, 2, 3]
            sage: s_par[:]
            [  4   5   6]
            [  0   0   0]
            [-12 -15 -18]
            sage: s_par == s
            True
            sage: b.__mul__(a) == t  # test of parallel comput.
            True
            sage: a.__mul__(a) == aa  # test of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Tensor product with a set of symmetric components::

            sage: s = b*aa; s
            3-index components w.r.t. [1, 2, 3], with symmetry on the index positions (1, 2)
            sage: s[:]
            [[[4, 0, -12], [0, 0, 0], [-12, 0, 36]],
             [[5, 0, -15], [0, 0, 0], [-15, 0, 45]],
             [[6, 0, -18], [0, 0, 0], [-18, 0, 54]]]
            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = b*aa; s_par
            3-index components w.r.t. [1, 2, 3], with symmetry on the index positions (1, 2)
            sage: s_par[:] == s[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if not isinstance(other, Components_base):
            raise TypeError("the second argument for the tensor product " +
                            "must be an instance of Components")
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        from .comp import Components
        if isinstance(other, ComponentsWithSym_dict):
            sym = []
            for s in other.parent()._sym:
                ns = tuple(s[i]+self._nid for i in range(len(s)))
                sym.append(ns)
            antisym = []
            for s in other.parent()._antisym:
                ns = tuple(s[i]+self._nid for i in range(len(s)))
                antisym.append(ns)
            result = Components(self.base_ring(), self._frame, self._nid + other._nid,
                                 self._sindex, self._output_formatter, sym,
                                 antisym)
        elif self._nid == 1 and other._nid == 1:
            if self is other:  # == would be dangerous here
                from .comp import CompFullySym
                # The result is symmetric:
                result = CompFullySym(self.base_ring(), self._frame, 2, self._sindex,
                                      self._output_formatter)
                # The loop below on self._comp.items() and
                # other._comp.items() cannot be used in the present case
                # (it would not deal correctly with redundant indices)
                # So we use a loop specific to the current case and return the
                # result:

                nproc = Parallelism().get('tensor')
                if nproc != 1:
                    # Parallel computation
                    lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
                    ind_list = [ind for ind in result.non_redundant_index_generator()]
                    ind_step = max(1, int(len(ind_list)/nproc))
                    local_list = lol(ind_list, ind_step)
                    # list of input parameters:
                    listParalInput = [(self, ind_part) for ind_part in local_list]

                    @parallel(p_iter='multiprocessing',ncpus=nproc)
                    def paral_mul(a, local_list_ind):
                        partial = []
                        for ind in local_list_ind:
                            partial.append([ind, a[[ind[0]]]*a[[ind[1]]]])
                        return partial

                    for ii,val in paral_mul(listParalInput):
                        for jj in val:
                            result[[jj[0]]] = jj[1]
                else:
                    # Sequential computation
                    for ind in result.non_redundant_index_generator():
                        result[[ind]] = self[[ind[0]]] * self[[ind[1]]]
                return result
            else:
                result = Components(self.base_ring(), self._frame, 2, self._sindex,
                                    self._output_formatter)
        else:
            result = Components(self.base_ring(), self._frame, self._nid + other._nid,
                                self._sindex, self._output_formatter)
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = [ind for ind  in self._comp]
            ind_step = max(1, int(len(ind_list)/nproc))
            local_list = lol(ind_list, ind_step)
            # list of input parameters:
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_mul(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    for ind_o, val_o in b._comp.items():
                        partial.append([ind + ind_o, a._comp[ind]*val_o])
                return partial

            for ii,val in paral_mul(listParalInput):
                for jj in val:
                    result._comp[jj[0]] = jj[1]
        else:
            # Sequential computation
            for ind_s, val_s in self._comp.items():
                for ind_o, val_o in other._comp.items():
                    result._comp[ind_s + ind_o] = val_s * val_o
        return result

    def _rmul_(self, other):
        r"""
        Reversed scalar multiplication (multiplication on the left by ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__rmul__(2) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [2, 0, -6]
            sage: s == 2*a
            True
            sage: a.__rmul__(0) == 0
            True

        """
        # Left multiplication by a "scalar":
        result = self._new_instance()
        if other == 0:
            return result   # because a just created Components is zero
        for ind, val in self._comp.items():
            result._comp[ind] = other * val
        return result

    _lmul_ = _rmul_

    def __truediv__(self, other):
        r"""
        Division (by a scalar).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(QQ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__truediv__(3) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [1/3, 0, -1]
            sage: s == a/3
            True
            sage: 3*s == a
            True

        """
        if isinstance(other, Components_base):
            raise NotImplementedError("division by an object of type " +
                                      "Components not implemented")
        result = self._new_instance()
        for ind, val in self._comp.items():
            result._comp[ind] = val / other
        return result

    def trace(self, pos1, pos2):
        r"""
        Index contraction.

        INPUT:

        - ``pos1`` -- position of the first index for the contraction (with the
          convention position=0 for the first slot)
        - ``pos2`` -- position of the second index for the contraction

        OUTPUT:

        - set of components resulting from the (pos1, pos2) contraction

        EXAMPLES:

        Self-contraction of a set of components with 2 indices::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: c.trace(0,1)
            15
            sage: c[0,0] + c[1,1] + c[2,2]  # check
            15

        Three self-contractions of a set of components with 3 indices::

            sage: v = Components(QQ, V.basis(), 1)
            sage: v[:] = (-1,2,3)
            sage: a = c*v ; a
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s = a.trace(0,1) ; s  # contraction on the first two indices
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [-15, 30, 45]
            sage: [sum(a[j,j,i] for j in range(3)) for i in range(3)]  # check
            [-15, 30, 45]
            sage: s = a.trace(0,2) ; s[:]  # contraction on the first and last indices
            [28, 32, 36]
            sage: [sum(a[j,i,j] for j in range(3)) for i in range(3)]  # check
            [28, 32, 36]
            sage: s = a.trace(1,2) ; s[:] # contraction on the last two indices
            [12, 24, 36]
            sage: [sum(a[i,j,j] for j in range(3)) for i in range(3)]  # check
            [12, 24, 36]

        """
        if self._nid < 2:
            raise ValueError("contraction can be performed only on " +
                             "components with at least 2 indices")
        if pos1 < 0 or pos1 > self._nid - 1:
            raise IndexError("pos1 out of range")
        if pos2 < 0 or pos2 > self._nid - 1:
            raise IndexError("pos2 out of range")
        if pos1 == pos2:
            raise IndexError("the two positions must differ for the " +
                             "contraction to be meaningful")
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 2:
            res = 0
            for i in range(si, nsi):
                res += self[[i,i]]
            return res
        else:
            # More than 2 indices
            from .comp import Components
            result = Components(self.base_ring(), self._frame, self._nid - 2,
                                self._sindex, self._output_formatter)
            if pos1 > pos2:
                pos1, pos2 = (pos2, pos1)
            for ind, val in self._comp.items():
                if ind[pos1] == ind[pos2]:
                    # there is a contribution to the contraction
                    ind_res = ind[:pos1] + ind[pos1+1:pos2] + ind[pos2+1:]
                    result[[ind_res]] += val
            return result

    def swap_adjacent_indices(self, pos1, pos2, pos3):
        r"""
        Swap two adjacent sets of indices.

        This method is essentially required to reorder the covariant and
        contravariant indices in the computation of a tensor product.

        INPUT:

        - ``pos1`` -- position of the first index of set 1 (with the convention
          ``position=0`` for the first slot)
        - ``pos2`` -- position of the first index of set 2 equals 1 plus the
          position of the last index of set 1 (since the two sets are adjacent)
        - ``pos3`` -- 1 plus position of the last index of set 2

        OUTPUT:

        - Components with index set 1 permuted with index set 2.

        EXAMPLES:

        Swap of the two indices of a 2-index set of components::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: c1 = c.swap_adjacent_indices(0,1,2)
            sage: c[:], c1[:]
            (
            [1 2 3]  [1 4 7]
            [4 5 6]  [2 5 8]
            [7 8 9], [3 6 9]
            )

        Swap of two pairs of indices on a 4-index set of components::

            sage: d = c*c1 ; d
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: d1 = d.swap_adjacent_indices(0,2,4)
            sage: d[0,1,1,2]
            16
            sage: d1[1,2,0,1]
            16
            sage: d1[0,1,1,2]
            24
            sage: d[1,2,0,1]
            24

        """
        result = self._new_instance()
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            result._comp[new_ind] = val
            # the above writing is more efficient than result[new_ind] = val
            # it does not work for the derived class ComponentsWithSym_dict, but for the
            # latter, the function ComponentsWithSym_dict.swap_adjacent_indices will be
            # called and not the present function.
        return result

#******************************************************************************

class ComponentsWithSym_dict(Components_dict):
    r"""
    Indexed set of ring elements forming some components with respect to a
    given "frame", with symmetries or antisymmetries regarding permutations
    of the indices.

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities,
    such as connection coefficients or structure coefficients.

    Subclasses of :class:`ComponentsWithSym_dict` are

    * :class:`ComponentsFullySym_dict` for fully symmetric components.
    * :class:`ComponentsFullyAntiSym_dict` for fully antisymmetric components.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an instance of ``ring`` and
      the second one, if any, some format specification.
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the indices: each symmetry is described by a tuple containing the
      positions of the involved indices, with the convention ``position=0``
      for the first slot; for instance:

        * ``sym = (0, 1)`` for a symmetry between the 1st and 2nd indices
        * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
          indices and a symmetry between the 2nd, 4th and 5th indices.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the indices, with the same convention as for ``sym``

    EXAMPLES:

    Symmetric components with 2 indices::

        sage: from sage.tensor.modules.comp import Components, CompWithSym
        sage: V = VectorSpace(QQ,3)
        sage: c = CompWithSym(QQ, V.basis(), 2, sym=(0,1))  # for demonstration only: it is preferable to use CompFullySym in this case
        sage: c[0,1] = 3
        sage: c[:]  # note that c[1,0] has been set automatically
        [0 3 0]
        [3 0 0]
        [0 0 0]

    Antisymmetric components with 2 indices::

        sage: c = CompWithSym(QQ, V.basis(), 2, antisym=(0,1))  # for demonstration only: it is preferable to use CompFullyAntiSym in this case
        sage: c[0,1] = 3
        sage: c[:]  # note that c[1,0] has been set automatically
        [ 0  3  0]
        [-3  0  0]
        [ 0  0  0]

    Internally, only non-redundant components are stored::

        sage: c._comp
        {(0, 1): 3}

    Components with 6 indices, symmetric among 3 indices (at position
    `(0, 1, 5)`) and antisymmetric among 2 indices (at position `(2, 4)`)::

        sage: c = CompWithSym(QQ, V.basis(), 6, sym=(0,1,5), antisym=(2,4))
        sage: c[0,1,2,0,1,2] = 3
        sage: c[1,0,2,0,1,2]  # symmetry between indices in position 0 and 1
        3
        sage: c[2,1,2,0,1,0]  # symmetry between indices in position 0 and 5
        3
        sage: c[0,2,2,0,1,1]  # symmetry between indices in position 1 and 5
        3
        sage: c[0,1,1,0,2,2]  # antisymmetry between indices in position 2 and 4
        -3

    Components with 4 indices, antisymmetric with respect to the first pair of
    indices as well as with the second pair of indices::

        sage: c = CompWithSym(QQ, V.basis(), 4, antisym=[(0,1),(2,3)])
        sage: c[0,1,0,1] = 3
        sage: c[1,0,0,1]  # antisymmetry on the first pair of indices
        -3
        sage: c[0,1,1,0]  # antisymmetry on the second pair of indices
        -3
        sage: c[1,0,1,0]  # consequence of the above
        3

    .. RUBRIC:: ARITHMETIC EXAMPLES

    Addition of a symmetric set of components with a non-symmetric one: the
    symmetry is lost::

        sage: V = VectorSpace(QQ, 3)
        sage: a = Components(QQ, V.basis(), 2)
        sage: a[:] = [[1,-2,3], [4,5,-6], [-7,8,9]]
        sage: b = CompWithSym(QQ, V.basis(), 2, sym=(0,1))  # for demonstration only: it is preferable to declare b = CompFullySym(QQ, V.basis(), 2)
        sage: b[0,0], b[0,1], b[0,2] = 1, 2, 3
        sage: b[1,1], b[1,2] = 5, 7
        sage: b[2,2] = 11
        sage: s = a + b ; s
        2-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:], s[:]
        (
        [ 1 -2  3]  [ 1  2  3]  [ 2  0  6]
        [ 4  5 -6]  [ 2  5  7]  [ 6 10  1]
        [-7  8  9], [ 3  7 11], [-4 15 20]
        )
        sage: a + b == b + a
        True

    Addition of two symmetric set of components: the symmetry is preserved::

        sage: c = CompWithSym(QQ, V.basis(), 2, sym=(0,1)) # for demonstration only: it is preferable to declare c = CompFullySym(QQ, V.basis(), 2)
        sage: c[0,0], c[0,1], c[0,2] = -4, 7, -8
        sage: c[1,1], c[1,2] = 2, -4
        sage: c[2,2] = 2
        sage: s = b + c ; s
        Fully symmetric 2-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: b[:], c[:], s[:]
        (
        [ 1  2  3]  [-4  7 -8]  [-3  9 -5]
        [ 2  5  7]  [ 7  2 -4]  [ 9  7  3]
        [ 3  7 11], [-8 -4  2], [-5  3 13]
        )
        sage: b + c == c + b
        True

    Check of the addition with counterparts not declared symmetric::

        sage: bn = Components(QQ, V.basis(), 2)
        sage: bn[:] = b[:]
        sage: bn == b
        True
        sage: cn = Components(QQ, V.basis(), 2)
        sage: cn[:] = c[:]
        sage: cn == c
        True
        sage: bn + cn == b + c
        True

    Addition of an antisymmetric set of components with a non-symmetric one:
    the antisymmetry is lost::

        sage: d = CompWithSym(QQ, V.basis(), 2, antisym=(0,1))  # for demonstration only: it is preferable to declare d = CompFullyAntiSym(QQ, V.basis(), 2)
        sage: d[0,1], d[0,2], d[1,2] = 4, -1, 3
        sage: s = a + d ; s
        2-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], d[:], s[:]
        (
        [ 1 -2  3]  [ 0  4 -1]  [ 1  2  2]
        [ 4  5 -6]  [-4  0  3]  [ 0  5 -3]
        [-7  8  9], [ 1 -3  0], [-6  5  9]
        )
        sage: d + a == a + d
        True

    Addition of two antisymmetric set of components: the antisymmetry is preserved::

        sage: e = CompWithSym(QQ, V.basis(), 2, antisym=(0,1))  # for demonstration only: it is preferable to declare e = CompFullyAntiSym(QQ, V.basis(), 2)
        sage: e[0,1], e[0,2], e[1,2] = 2, 3, -1
        sage: s = d + e ; s
        Fully antisymmetric 2-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: d[:], e[:], s[:]
        (
        [ 0  4 -1]  [ 0  2  3]  [ 0  6  2]
        [-4  0  3]  [-2  0 -1]  [-6  0  2]
        [ 1 -3  0], [-3  1  0], [-2 -2  0]
        )
        sage: e + d == d + e
        True

    """
    def __init__(self, parent, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: C = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            sage: TestSuite(C).run()

        """
        super().__init__(parent, frame, start_index=start_index,
                         output_formatter=output_formatter)

    def _ordered_indices(self, ind):
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
        parent = self.parent()
        return parent._ordered_indices(ind, self._ranges)

    def __getitem__(self, args):
        r"""
        Return the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the components
          `T_{ij...}` (for a 2-index object, a matrix is returned).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            sage: c.__getitem__((0,1,1,2)) # uninitialized components are zero
            0
            sage: c[0,1,1,2] = 5
            sage: c.__getitem__((0,1,1,2))
            5
            sage: c.__getitem__((1,0,1,2))
            5
            sage: c.__getitem__((0,1,2,1))
            -5
            sage: c[0,1,2,1]
            -5

        """
        no_format = self._output_formatter is None
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            no_format = True
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            return self._get_list(indices, no_format, format_type)
        else:
            sign, ind = self._ordered_indices(indices)
            if (sign == 0) or (ind not in self._comp): # the value is zero:
                if no_format:
                    return self.base_ring().zero()
                elif format_type is None:
                    return self._output_formatter(self.base_ring().zero())
                else:
                    return self._output_formatter(self.base_ring().zero(),
                                                 format_type)
            else: # non zero value
                if no_format:
                    if sign == 1:
                        return self._comp[ind]
                    else: # sign = -1
                        return -self._comp[ind]
                elif format_type is None:
                    if sign == 1:
                        return self._output_formatter(self._comp[ind])
                    else: # sign = -1
                        return self._output_formatter(-self._comp[ind])
                else:
                    if sign == 1:
                        return self._output_formatter(
                                                 self._comp[ind], format_type)
                    else: # sign = -1
                        return self._output_formatter(
                                                -self._comp[ind], format_type)

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) ; if ``[:]`` is provided, all the
          components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: c.__setitem__((1,2), 5)
            sage: c[:]
            [0 0 0]
            [0 0 5]
            [0 5 0]
            sage: c = CompWithSym(ZZ, [1,2,3], 2, antisym=(0,1))
            sage: c.__setitem__((1,2), 5)
            sage: c[:]
            [ 0  0  0]
            [ 0  0  5]
            [ 0 -5  0]
            sage: c.__setitem__((2,2), 5)
            Traceback (most recent call last):
            ...
            ValueError: by antisymmetry, the component cannot have a nonzero value for the indices (2, 2)

        """
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            self._set_list(indices, format_type, value)
        else:
            sign, ind = self._ordered_indices(indices)
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if sign == 0:
                if not zero_value:
                    raise ValueError("by antisymmetry, the component cannot " +
                                     "have a nonzero value for the indices " +
                                     str(indices))
                if ind in self._comp:
                    del self._comp[ind]  # zero values are not stored
            elif zero_value:
                if ind in self._comp:
                    del self._comp[ind]  # zero values are not stored
            else:
                if format_type is None:
                    if sign == 1:
                        self._comp[ind] = self.base_ring()(value)
                    else:   # sign = -1
                        self._comp[ind] = -self.base_ring()(value)
                else:
                    if sign == 1:
                        self._comp[ind] = self.base_ring()({format_type: value})
                    else:   # sign = -1
                        self._comp[ind] = -self.base_ring()({format_type: value})

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
        new_parent = self.parent().swap_adjacent_indices(pos1, pos2, pos3)
        result = self._new_instance(parent=new_parent)
        # The values:
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            result[new_ind] = val
        return result

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT:

        - components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: a = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: b[0,1], b[2,2] = 2, -3
            sage: s = a.__add__(b) ; s  # the symmetry is kept
            Fully symmetric 2-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [ 0  6  0]
            [ 6  0  5]
            [ 0  5 -3]
            sage: s == a + b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = a + b ; s_par
            Fully symmetric 2-index components w.r.t. [1, 2, 3]
            sage: s_par[:] == s[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Addition with different symmetries::

            sage: c = CompWithSym(ZZ, [1,2,3], 2, antisym=(0,1))
            sage: c[0,1], c[0,2] = 3, 7
            sage: s = a.__add__(c) ; s  # the symmetry is lost
            2-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [ 0  7  7]
            [ 1  0  5]
            [-7  5  0]

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = a + c ; s_par
            2-index components w.r.t. [1, 2, 3]
            sage: s_par[:] == s[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if isinstance(other, (int, Integer)) and other == 0:
            # TODO: This case can be removed once we define _add_ instead of __add__
            return +self
        if not isinstance(other, Components_base):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        if not other._frame:
            return +self
        if not self._frame:
            return +other
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._nid != self._nid:
            raise ValueError("the two sets of components do not have the " +
                             "same number of indices")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        from .comp import Components
        if isinstance(other, ComponentsWithSym_dict):
            # Are the symmetries of the same type ?
            self_sym = self.parent()._sym
            other_sym = other.parent()._sym
            diff_sym = set(self_sym).symmetric_difference(set(other_sym))
            self_antisym = self.parent()._antisym
            other_antisym = other.parent()._antisym
            diff_antisym = \
                set(self_antisym).symmetric_difference(set(other_antisym))
            if diff_sym == set() and diff_antisym == set():
                # The symmetries/antisymmetries are identical:
                result = self.copy()
                nproc = Parallelism().get('tensor')
                if nproc != 1:
                    # Parallel computation
                    lol = lambda lst, sz: [lst[i:i+sz] for i in
                                           range(0, len(lst), sz)]
                    ind_list = [ind for ind in other._comp]
                    ind_step = max(1, int(len(ind_list)/nproc/2))
                    local_list = lol(ind_list, ind_step)
                    # list of input parameters
                    listParalInput = [(self, other, ind_part) for ind_part in
                                      local_list]

                    @parallel(p_iter='multiprocessing', ncpus=nproc)
                    def paral_sum(a, b, local_list_ind):
                        partial = []
                        for ind in local_list_ind:
                            partial.append([ind, a[[ind]]+b[[ind]]])
                        return partial

                    for ii, val in paral_sum(listParalInput):
                        for jj in val:
                            result[[jj[0]]] = jj[1]

                else:
                    # Sequential computation
                    for ind, val in other._comp.items():
                        result[[ind]] += val
                return result
            else:
                # The symmetries/antisymmetries are different: only the
                # common ones are kept
                common_sym = []
                for isym in self.parent()._sym:
                    for osym in other.parent()._sym:
                        com = tuple(set(isym).intersection(set(osym)))
                        if len(com) > 1:
                            common_sym.append(com)
                common_antisym = []
                for isym in self.parent()._antisym:
                    for osym in other.parent()._antisym:
                        com = tuple(set(isym).intersection(set(osym)))
                        if len(com) > 1:
                            common_antisym.append(com)
                result = Components(self.base_ring(), self._frame, self._nid,
                                    self._sindex, self._output_formatter,
                                    common_sym, common_antisym)
        else:
            # other has no symmetry at all:
            result = Components(self.base_ring(), self._frame, self._nid,
                                self._sindex, self._output_formatter)
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = [ind for ind in result.non_redundant_index_generator()]
            ind_step = max(1, int(len(ind_list)/nproc/2))
            local_list = lol(ind_list, ind_step)
            # definition of the list of input parameters
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_sum(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    partial.append([ind, a[[ind]]+b[[ind]]])
                return partial

            for ii,val in paral_sum(listParalInput):
                for jj in val:
                    result[[jj[0]]] = jj[1]
        else:
            # Sequential computation
            for ind in result.non_redundant_index_generator():
                result[[ind]] = self[[ind]] + other[[ind]]
        return result

    def __mul__(self, other):
        r"""
        Component tensor product.

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT:

        - the tensor product of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: a = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: b[0,1], b[2,2] = 2, -3
            sage: s1 = a.__mul__(b) ; s1
            4-index components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: s1[1,0,0,1]
            8
            sage: s1[1,0,0,1] == a[1,0] * b[0,1]
            True
            sage: s1 == a*b
            True
            sage: c = CompWithSym(ZZ, [1,2,3], 2, antisym=(0,1))
            sage: c[0,1], c[0,2] = 3, 7
            sage: s2 = a.__mul__(c) ; s2
            4-index components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s2[1,0,2,0]
            -28
            sage: s2[1,0,2,0] == a[1,0] * c[2,0]
            True
            sage: s2 == a*c
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s1_par = a.__mul__(b) ; s1_par
            4-index components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: s1_par[1,0,0,1]
            8
            sage: s1_par[:] == s1[:]
            True
            sage: s2_par = a.__mul__(c) ; s2_par
            4-index components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s2_par[1,0,2,0]
            -28
            sage: s2_par[:] == s2[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if not isinstance(other, Components_base):
            raise TypeError("the second argument for the tensor product " +
                            "be an instance of Components")
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        sym = list(self.parent()._sym)
        antisym = list(self.parent()._antisym)
        if isinstance(other, ComponentsWithSym_dict):
            if other.parent()._sym != []:
                for s in other.parent()._sym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    sym.append(ns)
            if other.parent()._antisym != []:
                for s in other.parent()._antisym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    antisym.append(ns)
        from .comp import Components
        result = Components(self.base_ring(), self._frame, self._nid + other._nid,
                            self._sindex, self._output_formatter, sym, antisym)
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = [ind for ind in self._comp]
            ind_step = max(1, int(len(ind_list)/nproc))
            local_list = lol(ind_list, ind_step)
            # list of input parameters:
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_mul(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    for ind_o, val_o in b._comp.items():
                        partial.append([ind + ind_o, a._comp[ind]*val_o])
                return partial

            for ii,val in paral_mul(listParalInput):
                for jj in val:
                    result._comp[jj[0]] = jj[1]
        else:
            # Sequential computation
            for ind_s, val_s in self._comp.items():
                for ind_o, val_o in other._comp.items():
                    result._comp[ind_s + ind_o] = val_s * val_o
        return result

    def trace(self, pos1, pos2):
        r"""
        Index contraction, taking care of the symmetries.

        INPUT:

        - ``pos1`` -- position of the first index for the contraction (with
          the convention position=0 for the first slot)
        - ``pos2`` -- position of the second index for the contraction

        OUTPUT:

        - set of components resulting from the (pos1, pos2) contraction

        EXAMPLES:

        Self-contraction of symmetric 2-index components::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:   CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 3)
            sage: a = CompFullySym(QQ, V.basis(), 2)
            sage: a[:] = [[1,2,3],[2,4,5],[3,5,6]]
            sage: a.trace(0,1)
            11
            sage: a[0,0] + a[1,1] + a[2,2]
            11

        Self-contraction of antisymmetric 2-index components::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[1,2] = (3, -2, 1)
            sage: b.trace(0,1)  # must be zero by antisymmetry
            0

        Self-contraction of 3-index components with one symmetry::

            sage: v = Components(QQ, V.basis(), 1)
            sage: v[:] = (-2, 4, -8)
            sage: c = v*b ; c
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: s = c.trace(0,1) ; s
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [-28, 2, 8]
            sage: [sum(v[k]*b[k,i] for k in range(3)) for i in range(3)] # check
            [-28, 2, 8]
            sage: s = c.trace(1,2) ; s
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:] # is zero by antisymmetry
            [0, 0, 0]
            sage: c = b*v ; c
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: s = c.trace(0,1)
            sage: s[:]  # is zero by antisymmetry
            [0, 0, 0]
            sage: s = c.trace(1,2) ; s[:]
            [28, -2, -8]
            sage: [sum(b[i,k]*v[k] for k in range(3)) for i in range(3)]  # check
            [28, -2, -8]

        Self-contraction of 4-index components with two symmetries::

            sage: c = a*b ; c
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s = c.trace(0,1) ; s  # the symmetry on (0,1) is lost:
            Fully antisymmetric 2-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [  0  33 -22]
            [-33   0  11]
            [ 22 -11   0]
            sage: [[sum(c[k,k,i,j] for k in range(3)) for j in range(3)] for i in range(3)]  # check
            [[0, 33, -22], [-33, 0, 11], [22, -11, 0]]
            sage: s = c.trace(1,2) ; s  # both symmetries are lost by this contraction
            2-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [ 0  0  0]
            [-2  1  0]
            [-3  3 -1]
            sage: [[sum(c[i,k,k,j] for k in range(3)) for j in range(3)] for i in range(3)]  # check
            [[0, 0, 0], [-2, 1, 0], [-3, 3, -1]]

        """
        if self._nid < 2:
            raise TypeError("contraction can be performed only on " +
                            "components with at least 2 indices")
        if pos1 < 0 or pos1 > self._nid - 1:
            raise IndexError("pos1 out of range")
        if pos2 < 0 or pos2 > self._nid - 1:
            raise IndexError("pos2 out of range")
        if pos1 == pos2:
            raise IndexError("the two positions must differ for the " +
                             "contraction to take place")
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 2:
            res = 0
            for i in range(si, nsi):
                res += self[[i,i]]
            return res
        else:
            from .comp import Components
            # More than 2 indices
            if pos1 > pos2:
                pos1, pos2 = (pos2, pos1)
            # Determination of the remaining symmetries:
            sym_res = list(self.parent()._sym)
            for isym in self.parent()._sym:
                isym_res = list(isym)
                if pos1 in isym:
                    isym_res.remove(pos1)
                if pos2 in isym:
                    isym_res.remove(pos2)
                if len(isym_res) < 2:       # the symmetry is lost
                    sym_res.remove(isym)
                else:
                    sym_res[sym_res.index(isym)] = tuple(isym_res)
            antisym_res = list(self.parent()._antisym)
            for isym in self.parent()._antisym:
                isym_res = list(isym)
                if pos1 in isym:
                    isym_res.remove(pos1)
                if pos2 in isym:
                    isym_res.remove(pos2)
                if len(isym_res) < 2:       # the symmetry is lost
                    antisym_res.remove(isym)
                else:
                    antisym_res[antisym_res.index(isym)] = tuple(isym_res)
            # Shift of the index positions to take into account the
            # suppression of 2 indices:
            max_sym = 0
            for k in range(len(sym_res)):
                isym_res = []
                for pos in sym_res[k]:
                    if pos < pos1:
                        isym_res.append(pos)
                    elif pos < pos2:
                        isym_res.append(pos-1)
                    else:
                        isym_res.append(pos-2)
                max_sym = max(max_sym, len(isym_res))
                sym_res[k] = tuple(isym_res)
            max_antisym = 0
            for k in range(len(antisym_res)):
                isym_res = []
                for pos in antisym_res[k]:
                    if pos < pos1:
                        isym_res.append(pos)
                    elif pos < pos2:
                        isym_res.append(pos-1)
                    else:
                        isym_res.append(pos-2)
                max_antisym = max(max_antisym, len(isym_res))
                antisym_res[k] = tuple(isym_res)
            # Construction of the appropriate object in view of the
            # remaining symmetries:
            nid_res = self._nid - 2
            result = Components(self.base_ring(), self._frame, nid_res,
                                self._sindex, self._output_formatter,
                                sym=sym_res, antisym=antisym_res)
            # The contraction itself:
            for ind_res in result.non_redundant_index_generator():
                ind = list(ind_res)
                ind.insert(pos1, 0)
                ind.insert(pos2, 0)
                res = 0
                for i in range(si, nsi):
                    ind[pos1] = i
                    ind[pos2] = i
                    res += self[[ind]]
                result[[ind_res]] = res
            return result

    def symmetrize(self, *pos):
        r"""
        Symmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the
          symmetrization (with the convention ``position=0`` for the first
          slot); if none, the symmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`CompWithSym` describing the symmetrized
          components

        EXAMPLES:

        Symmetrization of 3-index components on a 3-dimensional space::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:   CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 3)
            sage: c[:] = [[[1,2,3], [4,5,6], [7,8,9]], [[10,11,12], [13,14,15], [16,17,18]], [[19,20,21], [22,23,24], [25,26,27]]]
            sage: cs = c.symmetrize(0,1) ; cs
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: s = cs.symmetrize() ; s
            Fully symmetric 3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: cs[:], s[:]
            ([[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
              [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
              [[13, 14, 15], [19, 20, 21], [25, 26, 27]]],
             [[[1, 16/3, 29/3], [16/3, 29/3, 14], [29/3, 14, 55/3]],
              [[16/3, 29/3, 14], [29/3, 14, 55/3], [14, 55/3, 68/3]],
              [[29/3, 14, 55/3], [14, 55/3, 68/3], [55/3, 68/3, 27]]])
            sage: s == c.symmetrize() # should be true
            True
            sage: s1 = cs.symmetrize(0,1) ; s1   # should return a copy of cs
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: s1 == cs    # check that s1 is a copy of cs
            True

        Let us now start with a symmetry on the last two indices::

            sage: cs1 = c.symmetrize(1,2) ; cs1
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: s2 = cs1.symmetrize() ; s2
            Fully symmetric 3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s2 == c.symmetrize()
            True

        Symmetrization alters pre-existing symmetries: let us symmetrize w.r.t.
        the index positions `(1, 2)` a set of components that is symmetric
        w.r.t. the index positions `(0, 1)`::

            sage: cs = c.symmetrize(0,1) ; cs
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: css = cs.symmetrize(1,2)
            sage: css # the symmetry (0,1) has been lost:
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: css[:]
            [[[1, 9/2, 8], [9/2, 8, 23/2], [8, 23/2, 15]],
             [[7, 21/2, 14], [21/2, 14, 35/2], [14, 35/2, 21]],
             [[13, 33/2, 20], [33/2, 20, 47/2], [20, 47/2, 27]]]
            sage: cs[:]
            [[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
             [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
             [[13, 14, 15], [19, 20, 21], [25, 26, 27]]]
            sage: css == c.symmetrize() # css differs from the full symmetrized version
            False
            sage: css.symmetrize() == c.symmetrize() # one has to symmetrize css over all indices to recover it
            True

        Another example of symmetry alteration: symmetrization over `(0, 1)` of
        a 4-index set of components that is symmetric w.r.t. `(1, 2, 3)`::

            sage: v = Components(QQ, V.basis(), 1)
            sage: v[:] = (-2,1,4)
            sage: a = v*s ; a
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2, 3)
            sage: a1 = a.symmetrize(0,1) ; a1 # the symmetry (1,2,3) has been reduced to (2,3):
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: a1.parent()._sym  # a1 has two distinct symmetries:
            ((0, 1), (2, 3))
            sage: a[0,1,2,0] == a[0,0,2,1]  # a is symmetric w.r.t. positions 1 and 3
            True
            sage: a1[0,1,2,0] == a1[0,0,2,1] # a1 is not
            False
            sage: a1[0,1,2,0] == a1[1,0,2,0] # but it is symmetric w.r.t. position 0 and 1
            True
            sage: a[0,1,2,0] == a[1,0,2,0] # while a is not
            False

        Partial symmetrization of 4-index components with an antisymmetry on
        the last two indices::

            sage: a = Components(QQ, V.basis(), 2)
            sage: a[:] = [[-1,2,3], [4,5,-6], [7,8,9]]
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[1,2] = (2, 4, 8)
            sage: c = a*b ; c
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (2, 3)
            sage: s = c.symmetrize(0,1) ; s  # symmetrization on the first two indices
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s[0,1,2,1] == (c[0,1,2,1] + c[1,0,2,1]) / 2 # check of the symmetrization
            True
            sage: s = c.symmetrize() ; s  # symmetrization over all the indices
            Fully symmetric 4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s == 0    # the full symmetrization results in zero due to the antisymmetry on the last two indices
            True
            sage: s = c.symmetrize(2,3) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (2, 3)
            sage: s == 0    # must be zero since the symmetrization has been performed on the antisymmetric indices
            True
            sage: s = c.symmetrize(0,2) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 2)
            sage: s != 0  # s is not zero, but the antisymmetry on (2,3) is lost because the position 2 is involved in the new symmetry
            True

        Partial symmetrization of 4-index components with an antisymmetry on
        the last three indices::

            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = (1, -2, 3)
            sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: b[0,1,2] = 4
            sage: c = a*b ; c
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2, 3)
            sage: s = c.symmetrize(0,1) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1),
               with antisymmetry on the index positions (2, 3)

        Note that the antisymmetry on `(1, 2, 3)` has been reduced to
        `(2, 3)` only::

            sage: s = c.symmetrize(1,2) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: s == 0 # because (1,2) are involved in the original antisymmetry
            True

        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the " \
                                 "total number of indices")
            pos = tuple(pos)
        pos_set = set(pos)
        # If the symmetry is already present, there is nothing to do:
        for isym in self.parent()._sym:
            if pos_set.issubset(set(isym)):
                return self.copy()

        #
        # Interference of the new symmetry with existing ones:
        #
        sym_res = [pos]  # starting the list of symmetries of the result
        for isym in self.parent()._sym:
            inter = pos_set.intersection(set(isym))
            # if len(inter) == len(isym), isym is included in the new symmetry
            # and therefore has not to be included in sym_res
            if len(inter) != len(isym):
                if len(inter) >= 1:
                    # some part of isym is lost
                    isym_set = set(isym)
                    for k in inter:
                        isym_set.remove(k)
                    if len(isym_set) > 1:
                        # some part of isym remains and must be included in sym_res:
                        isym_res = tuple(isym_set)
                        sym_res.append(isym_res)
                else:
                    # case len(inter)=0: no interference: the existing symmetry is
                    # added to the list of symmetries for the result:
                    sym_res.append(isym)
        #
        # Interference of the new symmetry with existing antisymmetries:
        #
        antisym_res = []  # starting the list of antisymmetries of the result
        zero_result = False
        for iasym in self.parent()._antisym:
            inter = pos_set.intersection(set(iasym))
            if len(inter) > 1:
                # If at least two of the symmetry indices are already involved
                # in the antisymmetry, the outcome is zero:
                zero_result = True
            elif len(inter) == 1:
                # some piece of antisymmetry is lost
                k = inter.pop()  # the symmetry index position involved in the
                                 # antisymmetry
                iasym_set = set(iasym)
                iasym_set.remove(k)
                if len(iasym_set) > 1:
                    iasym_res = tuple(iasym_set)
                    antisym_res.append(iasym_res)
                # if len(iasym_set) == 1, the antisymmetry is fully lost, it is
                # therefore not appended to antisym_res
            else:
                # case len(inter)=0: no interference: the antisymmetry is
                # added to the list of antisymmetries for the result:
                antisym_res.append(iasym)
        #
        # Creation of the result object
        #
        max_sym = 0
        for isym in sym_res:
            max_sym = max(max_sym, len(isym))
        from .comp import Components
        result = Components(self.base_ring(), self._frame, self._nid, self._sindex,
                            self._output_formatter, sym=sym_res,
                            antisym=antisym_res)
        if zero_result:
            return result   # since a just created instance is zero
        #
        # Symmetrization
        #
        n_sym = len(pos) # number of indices involved in the symmetry
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
        Antisymmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the antisymmetrization
          (with the convention ``position=0`` for the first slot); if none, the
          antisymmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`CompWithSym` describing the antisymmetrized
          components

        EXAMPLES:

        Antisymmetrization of 3-index components on a 3-dimensional space::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:  CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 3)
            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = (-2,1,3)
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[1,2] = (4,1,2)
            sage: c = a*b ; c   # tensor product of a by b
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            ([[[0, -8, -2], [8, 0, -4], [2, 4, 0]],
              [[0, 4, 1], [-4, 0, 2], [-1, -2, 0]],
              [[0, 12, 3], [-12, 0, 6], [-3, -6, 0]]],
             [[[0, 0, 0], [0, 0, 7/3], [0, -7/3, 0]],
              [[0, 0, -7/3], [0, 0, 0], [7/3, 0, 0]],
              [[0, 7/3, 0], [-7/3, 0, 0], [0, 0, 0]]])

        Check of the antisymmetrization::

            sage: all(s[i,j,k] == (c[i,j,k]-c[i,k,j]+c[j,k,i]-c[j,i,k]+c[k,i,j]-c[k,j,i])/6
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True

        Antisymmetrization over already antisymmetric indices does not change anything::

            sage: s1 = s.antisymmetrize(1,2) ; s1
            Fully antisymmetric 3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s1 == s
            True
            sage: c1 = c.antisymmetrize(1,2) ; c1
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: c1 == c
            True

        But in general, antisymmetrization may alter previous antisymmetries::

            sage: c2 = c.antisymmetrize(0,1) ; c2  # the antisymmetry (2,3) is lost:
            3-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: c2 == c
            False
            sage: c = s*a ; c
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1, 2)
            sage: s = c.antisymmetrize(1,3) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 2), with antisymmetry on the index positions (1, 3)
            sage: s.parent()._antisym  # the antisymmetry (0,1,2) has been reduced to (0,2), since 1 is involved in the new antisymmetry (1,3):
            ((0, 2), (1, 3))

        Partial antisymmetrization of 4-index components with a symmetry on
        the first two indices::

            sage: a = CompFullySym(QQ, V.basis(), 2)
            sage: a[:] = [[-2,1,3], [1,0,-5], [3,-5,4]]
            sage: b = Components(QQ, V.basis(), 2)
            sage: b[:] = [[1,2,3], [5,7,11], [13,17,19]]
            sage: c = a*b ; c
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: s = c.antisymmetrize(2,3) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1),
               with antisymmetry on the index positions (2, 3)

        Some check of the antisymmetrization::

            sage: all(s[2,2,i,j] == (c[2,2,i,j] - c[2,2,j,i])/2
            ....:     for i in range(3) for j in range(i,3))
            True

        The full antisymmetrization results in zero because of the symmetry on the
        first two indices::

            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s == 0
            True

        Similarly, the partial antisymmetrization on the first two indices results in zero::

            sage: s = c.antisymmetrize(0,1) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: s == 0
            True

        The partial antisymmetrization on the positions `(0, 2)` destroys
        the symmetry on `(0, 1)`::

            sage: s = c.antisymmetrize(0,2) ; s
            4-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 2)
            sage: s != 0
            True
            sage: s[0,1,2,1]
            27/2
            sage: s[1,0,2,1]  # the symmetry (0,1) is lost
            -2
            sage: s[2,1,0,1]  # the antisymmetry (0,2) holds
            -27/2

        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the " \
                                 "total number of indices")
            pos = tuple(pos)
        pos_set = set(pos)
        # If the antisymmetry is already present, there is nothing to do:
        for iasym in self.parent()._antisym:
            if pos_set.issubset(set(iasym)):
                return self.copy()
        #
        # Interference of the new antisymmetry with existing ones
        #
        antisym_res = [pos]  # starting the list of symmetries of the result
        for iasym in self.parent()._antisym:
            inter = pos_set.intersection(set(iasym))
            # if len(inter) == len(iasym), iasym is included in the new
            # antisymmetry and therefore has not to be included in antisym_res
            if len(inter) != len(iasym):
                if len(inter) >= 1:
                    # some part of iasym is lost
                    iasym_set = set(iasym)
                    for k in inter:
                        iasym_set.remove(k)
                    if len(iasym_set) > 1:
                        # some part of iasym remains and must be included in
                        # antisym_res:
                        iasym_res = tuple(iasym_set)
                        antisym_res.append(iasym_res)
                else:
                    # case len(inter)=0: no interference: the existing
                    # antisymmetry is added to the list of antisymmetries for
                    # the result:
                    antisym_res.append(iasym)
        #
        # Interference of the new antisymmetry with existing symmetries
        #
        sym_res = []  # starting the list of symmetries of the result
        zero_result = False
        for isym in self.parent()._sym:
            inter = pos_set.intersection(set(isym))
            if len(inter) > 1:
                # If at least two of the antisymmetry indices are already
                # involved in the symmetry, the outcome is zero:
                zero_result = True
            elif len(inter) == 1:
                # some piece of the symmetry is lost
                k = inter.pop()  # the antisymmetry index position involved in
                                 # the symmetry
                isym_set = set(isym)
                isym_set.remove(k)
                if len(isym_set) > 1:
                    isym_res = tuple(isym_set)
                    sym_res.append(isym_res)
                # if len(isym_set) == 1, the symmetry is fully lost, it is
                # therefore not appended to sym_res
            else:
                # case len(inter)=0: no interference: the symmetry is
                # added to the list of symmetries for the result:
                sym_res.append(isym)
        #
        # Creation of the result object
        #
        max_sym = 0
        for isym in antisym_res:
            max_sym = max(max_sym, len(isym))
        from .comp import Components
        result = Components(self.base_ring(), self._frame, self._nid, self._sindex,
                            self._output_formatter, sym=sym_res,
                            antisym=antisym_res)
        if zero_result:
            return result   # since a just created instance is zero
        #
        # Antisymmetrization
        #
        n_sym = len(pos) # number of indices involved in the antisymmetry
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

#******************************************************************************

class ComponentsFullySym_dict(ComponentsWithSym_dict):
    r"""
    Indexed set of ring elements forming some components with respect to a
    given "frame" that are fully symmetric with respect to any permutation
    of the indices.

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an instance of ``ring`` and
      the second one, if any, some format specification.

    EXAMPLES:

    Symmetric components with 2 indices on a 3-dimensional space::

        sage: from sage.tensor.modules.comp import CompFullySym, CompWithSym
        sage: V = VectorSpace(QQ, 3)
        sage: c = CompFullySym(QQ, V.basis(), 2)
        sage: c[0,0], c[0,1], c[1,2] = 1, -2, 3
        sage: c[:] # note that c[1,0] and c[2,1] have been updated automatically (by symmetry)
        [ 1 -2  0]
        [-2  0  3]
        [ 0  3  0]

    Internally, only non-redundant and non-zero components are stored::

        sage: c._comp  # random output order of the component dictionary
        {(0, 0): 1, (0, 1): -2, (1, 2): 3}

    Same thing, but with the starting index set to 1::

        sage: c1 = CompFullySym(QQ, V.basis(), 2, start_index=1)
        sage: c1[1,1], c1[1,2], c1[2,3] = 1, -2, 3
        sage: c1[:]
        [ 1 -2  0]
        [-2  0  3]
        [ 0  3  0]

    The values stored in ``c`` and ``c1`` are equal::

        sage: c1[:] == c[:]
        True

    but not ``c`` and ``c1``, since their starting indices differ::

        sage: c1 == c
        False

    Fully symmetric components with 3 indices on a 3-dimensional space::

        sage: a = CompFullySym(QQ, V.basis(), 3)
        sage: a[0,1,2] = 3
        sage: a[:]
        [[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
         [[0, 0, 3], [0, 0, 0], [3, 0, 0]],
         [[0, 3, 0], [3, 0, 0], [0, 0, 0]]]
        sage: a[0,1,0] = 4
        sage: a[:]
        [[[0, 4, 0], [4, 0, 3], [0, 3, 0]],
         [[4, 0, 3], [0, 0, 0], [3, 0, 0]],
         [[0, 3, 0], [3, 0, 0], [0, 0, 0]]]

    The full symmetry is preserved by the arithmetics::

        sage: b = CompFullySym(QQ, V.basis(), 3)
        sage: b[0,0,0], b[0,1,0], b[1,0,2], b[1,2,2] = -2, 3, 1, -5
        sage: s = a + 2*b ; s
        Fully symmetric 3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:], s[:]
        ([[[0, 4, 0], [4, 0, 3], [0, 3, 0]],
          [[4, 0, 3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [3, 0, 0], [0, 0, 0]]],
         [[[-2, 3, 0], [3, 0, 1], [0, 1, 0]],
          [[3, 0, 1], [0, 0, 0], [1, 0, -5]],
          [[0, 1, 0], [1, 0, -5], [0, -5, 0]]],
         [[[-4, 10, 0], [10, 0, 5], [0, 5, 0]],
          [[10, 0, 5], [0, 0, 0], [5, 0, -10]],
          [[0, 5, 0], [5, 0, -10], [0, -10, 0]]])

    It is lost if the added object is not fully symmetric::

        sage: b1 = CompWithSym(QQ, V.basis(), 3, sym=(0,1))  # b1 has only symmetry on index positions (0,1)
        sage: b1[0,0,0], b1[0,1,0], b1[1,0,2], b1[1,2,2] = -2, 3, 1, -5
        sage: s = a + 2*b1 ; s  # the result has the same symmetry as b1:
        3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with symmetry on the index positions (0, 1)
        sage: a[:], b1[:], s[:]
        ([[[0, 4, 0], [4, 0, 3], [0, 3, 0]],
          [[4, 0, 3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [3, 0, 0], [0, 0, 0]]],
         [[[-2, 0, 0], [3, 0, 1], [0, 0, 0]],
          [[3, 0, 1], [0, 0, 0], [0, 0, -5]],
          [[0, 0, 0], [0, 0, -5], [0, 0, 0]]],
         [[[-4, 4, 0], [10, 0, 5], [0, 3, 0]],
          [[10, 0, 5], [0, 0, 0], [3, 0, -10]],
          [[0, 3, 0], [3, 0, -10], [0, 0, 0]]])
        sage: s = 2*b1 + a ; s
        3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with symmetry on the index positions (0, 1)
        sage: 2*b1 + a == a + 2*b1
        True

    """
    def __init__(self, parent, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: C = CompFullySym(ZZ, (1,2,3), 2)
            sage: TestSuite(C).run()

        """
        super().__init__(parent, frame, start_index=start_index,
                         output_formatter=output_formatter)

    def __getitem__(self, args):
        r"""
        Return the component corresponding to the given indices of ``self``.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the
          components `T_{ij...}` (for a 2-index object, a matrix is returned)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c[0,1] = 4
            sage: c.__getitem__((0,1))
            4
            sage: c.__getitem__((1,0))
            4
            sage: c.__getitem__(slice(None))
            [0 4 0]
            [4 0 0]
            [0 0 0]

        """
        no_format = self._output_formatter is None
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            no_format = True
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]

        if isinstance(indices, slice):
            return self._get_list(indices, no_format, format_type)

        ind = self._ordered_indices(indices)[1]  # [0]=sign is not used
        if ind in self._comp: # non zero value
            if no_format:
                return self._comp[ind]
            elif format_type is None:
                return self._output_formatter(self._comp[ind])
            else:
                return self._output_formatter(self._comp[ind], format_type)

        # the value is zero
        if no_format:
            return self.base_ring().zero()
        elif format_type is None:
            return self._output_formatter(self.base_ring().zero())
        else:
            return self._output_formatter(self.base_ring().zero(),
                                             format_type)

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:

        - ``indices`` -- list of indices (possibly a single integer if
          self is a 1-index object) ; if [:] is provided, all the components
          are set.
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c.__setitem__((0,1), 4)
            sage: c[:]
            [0 4 0]
            [4 0 0]
            [0 0 0]
            sage: c.__setitem__((2,1), 5)
            sage: c[:]
            [0 4 0]
            [4 0 5]
            [0 5 0]
            sage: c.__setitem__(slice(None), [[1, 2, 3], [2, 4, 5], [3, 5, 6]])
            sage: c[:]
            [1 2 3]
            [2 4 5]
            [3 5 6]

        """
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            self._set_list(indices, format_type, value)
        else:
            ind = self._ordered_indices(indices)[1]  # [0]=sign is not used
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if zero_value:
                # if the component has been set previously, it is deleted,
                # otherwise nothing is done (zero components are not stored):
                if ind in self._comp:
                    del self._comp[ind]  # zero values are not stored
            else:
                if format_type is None:
                    self._comp[ind] = self.base_ring()(value)
                else:
                    self._comp[ind] = self.base_ring()({format_type: value})

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT:

        - components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: a = CompFullySym(ZZ, (1,2,3), 2)
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompFullySym(ZZ, (1,2,3), 2)
            sage: b[0,1], b[2,2] = 2, -3
            sage: s = a.__add__(b) ; s  # the symmetry is kept
            Fully symmetric 2-index components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  6  0]
            [ 6  0  5]
            [ 0  5 -3]
            sage: s == a + b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__add__(b) ; s_par
            Fully symmetric 2-index components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  6  0]
            [ 6  0  5]
            [ 0  5 -3]
            sage: s_par == s
            True
            sage: s_par == b.__add__(a)  # test of commutativity of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Addition with a set of components having different symmetries::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: c = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: c[0,1], c[0,2] = 3, 7
            sage: s = a.__add__(c) ; s  # the symmetry is lost
            2-index components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  7  7]
            [ 1  0  5]
            [-7  5  0]
            sage: s == a + c
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__add__(c) ; s_par
            2-index components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  7  7]
            [ 1  0  5]
            [-7  5  0]
            sage: s_par[:] == s[:]
            True
            sage: s_par == c.__add__(a)  # test of commutativity of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components_base):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        if not other._frame:
            return +self
        if not self._frame:
            return +other
        if isinstance(other, ComponentsFullySym_dict):
            if other._frame != self._frame:
                raise ValueError("the two sets of components are not defined " +
                                 "on the same frame")
            if other._nid != self._nid:
                raise ValueError("the two sets of components do not have the " +
                                 "same number of indices")
            if other._sindex != self._sindex:
                raise ValueError("the two sets of components do not have the " +
                                 "same starting index")
            # Initialization of the result to self.copy(), so that there
            # remains only to add other:
            result = self.copy()
            nproc = Parallelism().get('tensor')
            if nproc != 1:
                # parallel sum
                lol = lambda lst, sz: [lst[i:i+sz] for i
                                       in range(0, len(lst), sz)]
                ind_list = [ind for ind in other._comp]
                ind_step = max(1, int(len(ind_list)/nproc/2))
                local_list = lol(ind_list, ind_step)
                # definition of the list of input parameters
                listParalInput = [(self, other, ind_part) for ind_part
                                  in local_list]

                @parallel(p_iter='multiprocessing', ncpus=nproc)
                def paral_sum(a, b, local_list_ind):
                    partial = []
                    for ind in local_list_ind:
                        partial.append([ind, a[[ind]]+b[[ind]]])
                    return partial

                for ii,val in paral_sum(listParalInput):
                    for jj in val:
                        result[[jj[0]]] = jj[1]

            else:
                # sequential sum
                for ind, val in other._comp.items():
                    result[[ind]] += val
            return result
        else:
            return super().__add__(other)

#******************************************************************************

class ComponentsFullyAntiSym_dict(ComponentsWithSym_dict):
    r"""
    Indexed set of ring elements forming some components with respect to a
    given "frame" that are fully antisymmetric with respect to any permutation
    of the indices.

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an instance of ``ring`` and
      the second one, if any, some format specification.

    EXAMPLES:

    Antisymmetric components with 2 indices on a 3-dimensional space::

        sage: from sage.tensor.modules.comp import CompWithSym, CompFullyAntiSym
        sage: V = VectorSpace(QQ, 3)
        sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
        sage: c[0,1], c[0,2], c[1,2] = 3, 1/2, -1
        sage: c[:]  # note that all components have been set according to the antisymmetry
        [   0    3  1/2]
        [  -3    0   -1]
        [-1/2    1    0]

    Internally, only non-redundant and non-zero components are stored::

        sage: c._comp # random output order of the component dictionary
        {(0, 1): 3, (0, 2): 1/2, (1, 2): -1}

    Same thing, but with the starting index set to 1::

        sage: c1 = CompFullyAntiSym(QQ, V.basis(), 2, start_index=1)
        sage: c1[1,2], c1[1,3], c1[2,3] = 3, 1/2, -1
        sage: c1[:]
        [   0    3  1/2]
        [  -3    0   -1]
        [-1/2    1    0]

    The values stored in ``c`` and ``c1`` are equal::

        sage: c1[:] == c[:]
        True

    but not ``c`` and ``c1``, since their starting indices differ::

        sage: c1 == c
        False

    Fully antisymmetric components with 3 indices on a 3-dimensional space::

        sage: a = CompFullyAntiSym(QQ, V.basis(), 3)
        sage: a[0,1,2] = 3  # the only independent component in dimension 3
        sage: a[:]
        [[[0, 0, 0], [0, 0, 3], [0, -3, 0]],
         [[0, 0, -3], [0, 0, 0], [3, 0, 0]],
         [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]]

    Setting a nonzero value incompatible with the antisymmetry results in an
    error::

        sage: a[0,1,0] = 4
        Traceback (most recent call last):
        ...
        ValueError: by antisymmetry, the component cannot have a nonzero value for the indices (0, 1, 0)
        sage: a[0,1,0] = 0   # OK
        sage: a[2,0,1] = 3   # OK

    The full antisymmetry is preserved by the arithmetics::

        sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
        sage: b[0,1,2] = -4
        sage: s = a + 2*b ; s
        Fully antisymmetric 3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:], s[:]
        ([[[0, 0, 0], [0, 0, 3], [0, -3, 0]],
          [[0, 0, -3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -4], [0, 4, 0]],
          [[0, 0, 4], [0, 0, 0], [-4, 0, 0]],
          [[0, -4, 0], [4, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -5], [0, 5, 0]],
          [[0, 0, 5], [0, 0, 0], [-5, 0, 0]],
          [[0, -5, 0], [5, 0, 0], [0, 0, 0]]])

    It is lost if the added object is not fully antisymmetric::

        sage: b1 = CompWithSym(QQ, V.basis(), 3, antisym=(0,1))  # b1 has only antisymmetry on index positions (0,1)
        sage: b1[0,1,2] = -4
        sage: s = a + 2*b1 ; s  # the result has the same symmetry as b1:
        3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with antisymmetry on the index positions (0, 1)
        sage: a[:], b1[:], s[:]
        ([[[0, 0, 0], [0, 0, 3], [0, -3, 0]],
          [[0, 0, -3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -4], [0, 0, 0]],
          [[0, 0, 4], [0, 0, 0], [0, 0, 0]],
          [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -5], [0, -3, 0]],
          [[0, 0, 5], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]])
        sage: s = 2*b1 + a ; s
        3-index components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with antisymmetry on the index positions (0, 1)
        sage: 2*b1 + a == a + 2*b1
        True

    """
    def __init__(self, parent, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: C = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: TestSuite(C).run()

        """
        super().__init__(parent, frame, start_index=start_index,
                         output_formatter=output_formatter)

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT:

        - components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: a = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: b[0,1], b[0,2] = 2, -3
            sage: s = a.__add__(b) ; s  # the antisymmetry is kept
            Fully antisymmetric 2-index components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  6 -3]
            [-6  0  5]
            [ 3 -5  0]
            sage: s == a + b
            True
            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c[0,1], c[0,2] = 3, 7
            sage: s = a.__add__(c) ; s  # the antisymmetry is lost
            2-index components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  7  7]
            [-1  0  5]
            [ 7 -5  0]
            sage: s == a + c
            True

        Parallel computation::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: Parallelism().set('tensor', nproc=2)
            sage: a = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: b[0,1], b[0,2] = 2, -3
            sage: s_par = a.__add__(b) ; s_par  # the antisymmetry is kept
            Fully antisymmetric 2-index components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  6 -3]
            [-6  0  5]
            [ 3 -5  0]
            sage: s_par == a + b
            True
            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c[0,1], c[0,2] = 3, 7
            sage: s_par = a.__add__(c) ; s_par  # the antisymmetry is lost
            2-index components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  7  7]
            [-1  0  5]
            [ 7 -5  0]
            sage: s_par == a + c
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components_base):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        if not other._frame:
            return +self
        if not self._frame:
            return +other
        if isinstance(other, ComponentsFullyAntiSym_dict):
            if other._frame != self._frame:
                raise ValueError("the two sets of components are not defined " +
                                 "on the same frame")
            if other._nid != self._nid:
                raise ValueError("the two sets of components do not have the " +
                                 "same number of indices")
            if other._sindex != self._sindex:
                raise ValueError("the two sets of components do not have the " +
                                 "same starting index")
            # Initialization of the result to self.copy(), so that there remains
            # only to add other:
            result = self.copy()
            nproc = Parallelism().get('tensor')
            if nproc != 1:
                # Parallel computation
                lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
                ind_list = [ind for ind in other._comp]
                ind_step = max(1, int(len(ind_list)/nproc/2))
                local_list = lol(ind_list, ind_step)
                # list of input parameters
                listParalInput = [(self, other, ind_part) for ind_part in local_list]

                @parallel(p_iter='multiprocessing', ncpus=nproc)
                def paral_sum(a, b, local_list_ind):
                    partial = []
                    for ind in local_list_ind:
                        partial.append([ind, a[[ind]]+b[[ind]]])
                    return partial

                for ii, val in paral_sum(listParalInput):
                    for jj in val:
                        result[[jj[0]]] = jj[1]

            else:
                # Sequential computation
                for ind, val in other._comp.items():
                    result[[ind]] += val

            return result
        else:
            return super().__add__(other)

    def interior_product(self, other):
        r"""
        Interior product with another set of fully antisymmetric components.

        The interior product amounts to a contraction over all the `p` indices
        of ``self`` with the first `p` indices of ``other``, assuming that
        the number `q` of indices of ``other`` obeys `q\geq p`.

        .. NOTE::

            ``self.interior_product(other)`` yields the same result as
            ``self.contract(0,..., p-1, other, 0,..., p-1)``
            (cf. :meth:`~sage.tensor.modules.comp.Components.contract`), but
            ``interior_product`` is more efficient, the antisymmetry of ``self``
            being not used to reduce the computation in
            :meth:`~sage.tensor.modules.comp.Components.contract`.

        INPUT:

        - ``other`` -- fully antisymmetric components defined on the same frame
          as ``self`` and with a number of indices at least equal to that of
          ``self``

        OUTPUT:

        - base ring element (case `p=q`) or set of components (case `p<q`)
          resulting from the contraction over all the `p` indices of ``self``
          with the first `p` indices of ``other``

        EXAMPLES:

        Interior product of a set of components ``a`` with ``p`` indices with a
        set of components ``b`` with ``q`` indices on a 4-dimensional vector
        space.

        Case ``p=2`` and ``q=2``::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: V = VectorSpace(QQ, 4)
            sage: a = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: a[0,1], a[0,2], a[0,3] = -2, 4, 3
            sage: a[1,2], a[1,3], a[2,3] = 5, -3, 1
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[0,3] = 3, -4, 2
            sage: b[1,2], b[1,3], b[2,3] = 2, 5, 1
            sage: c = a.interior_product(b)
            sage: c
            -40
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case ``p=2`` and ``q=3``::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: b[0,1,2], b[0,1,3], b[0,2,3], b[1,2,3] = 3, -4, 2, 5
            sage: c = a.interior_product(b)
            sage: c[:]
            [58, 10, 6, 82]
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case ``p=2`` and ``q=4``::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: b[0,1,2,3] = 5
            sage: c = a.interior_product(b)
            sage: c[:]
            [  0  10  30  50]
            [-10   0  30 -40]
            [-30 -30   0 -20]
            [-50  40  20   0]
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case ``p=3`` and ``q=3``::

            sage: a = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: a[0,1,2], a[0,1,3], a[0,2,3], a[1,2,3] = 2, -1, 3, 5
            sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: b[0,1,2], b[0,1,3], b[0,2,3], b[1,2,3] = -2, 1, 4, 2
            sage: c = a.interior_product(b)
            sage: c
            102
            sage: c == a.contract(0, 1, 2, b, 0, 1, 2)
            True

        Case ``p=3`` and ``q=4``::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: b[0,1,2,3] = 5
            sage: c = a.interior_product(b)
            sage: c[:]
            [-150, 90, 30, 60]
            sage: c == a.contract(0, 1, 2, b, 0, 1, 2)
            True

        Case ``p=4`` and ``q=4``::

            sage: a = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: a[0,1,2,3] = 3
            sage: c = a.interior_product(b)
            sage: c
            360
            sage: c == a.contract(0, 1, 2, 3, b, 0, 1, 2, 3)
            True

        """
        from sage.functions.other import factorial
        # Sanity checks:
        if not isinstance(other, ComponentsFullyAntiSym_dict):
            raise TypeError("{} is not a fully antisymmetric ".format(other) +
                            "set of components")
        if other._frame != self._frame:
            raise ValueError("The {} are not defined on the ".format(other) +
                             "same frame as the {}".format(self))
        if other._nid < self._nid:
            raise ValueError("The {} have less indices than ".format(other) +
                             "the {}".format(self))
        # Number of indices of the result:
        res_nid = other._nid - self._nid
        #
        # Case of a scalar result
        #
        if res_nid == 0:
            res = 0
            for ind in self.non_redundant_index_generator():
                res += self[[ind]] * other[[ind]]
            return factorial(self._nid)*res
        #
        # Case of component result
        #
        from .comp import Components, CompFullyAntiSym
        if res_nid == 1:
            res = Components(self.base_ring(), self._frame, res_nid,
                             start_index=self._sindex,
                             output_formatter=self._output_formatter)
        else:
            res = CompFullyAntiSym(self.base_ring(), self._frame, res_nid,
                                   start_index=self._sindex,
                                   output_formatter=self._output_formatter)
        factorial_s = factorial(self._nid)
        for ind in res.non_redundant_index_generator():
            sm = 0
            for ind_s, cmp_s in self._comp.items():
                ind_o = ind_s + ind
                sm += cmp_s * other[[ind_o]]
            res[[ind]] = factorial_s*sm
        return res

#******************************************************************************

class ComponentsKroneckerDelta_dict(ComponentsFullySym_dict):
    r"""
    Kronecker delta `\delta_{ij}`.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method ``__getitem__``); ``output_formatter`` must
      take 1 or 2 arguments: the first argument must be an instance of
      ``ring`` and the second one, if any, some format specification

    EXAMPLES:

    The Kronecker delta on a 3-dimensional space::

        sage: from sage.tensor.modules.comp import KroneckerDelta
        sage: V = VectorSpace(QQ,3)
        sage: d = KroneckerDelta(QQ, V.basis()) ; d
        Kronecker delta of size 3x3
        sage: d[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]

    One can read, but not set, the components of a Kronecker delta::

        sage: d[1,1]
        1
        sage: d[1,1] = 2
        Traceback (most recent call last):
        ...
        TypeError: the components of a Kronecker delta cannot be changed

    Examples of use with output formatters::

        sage: d = KroneckerDelta(QQ, V.basis(), output_formatter=Rational.numerical_approx)
        sage: d[:]  # default format (53 bits of precision)
        [ 1.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000  1.00000000000000]
        sage: d[:,10] # format = 10 bits of precision
        [ 1.0 0.00 0.00]
        [0.00  1.0 0.00]
        [0.00 0.00  1.0]
        sage: d = KroneckerDelta(QQ, V.basis(), output_formatter=str)
        sage: d[:]
        [['1', '0', '0'], ['0', '1', '0'], ['0', '0', '1']]

    """
    def __init__(self, parent, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import KroneckerDelta
            sage: d = KroneckerDelta(ZZ, (1,2,3))
            sage: TestSuite(d).run(skip='_test_category')

        """
        super().__init__(parent, frame, start_index=start_index,
                         output_formatter=output_formatter)
        for i in range(self._sindex, self._dim + self._sindex):
            self._comp[(i,i)] = self.base_ring()(1)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import KroneckerDelta
            sage: KroneckerDelta(ZZ, (1,2,3))
            Kronecker delta of size 3x3

        """
        n = str(self._dim)
        return "Kronecker delta of size " + n + "x" + n

    def __setitem__(self, args, value):
        r"""
        Should not be used (the components of a Kronecker delta are constant)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import KroneckerDelta
            sage: d = KroneckerDelta(ZZ, (1,2,3))
            sage: d.__setitem__((0,0), 1)
            Traceback (most recent call last):
            ...
            TypeError: the components of a Kronecker delta cannot be changed

        """
        raise TypeError("the components of a Kronecker delta cannot be changed")
