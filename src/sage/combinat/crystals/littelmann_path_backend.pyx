r"""
Littelmann Paths Backend

This provides a Cython implementation of Littelmann paths using
lists with entries as dense vectors in the weight space given in
terms of the corresponding basis. This is designed to be as fast
and lightweight as possible.

AUTHORS:

- Travis Scrimshaw (2016): Initial implementation
"""

#****************************************************************************
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper cimport ElementWrapper
from sage.functions.other import floor
cimport cython

cdef class LittelmannPath(object):
    """
    TESTS::

        sage: C = crystals.LSPaths(['E',6],[1,0,0,0,0,0])
        sage: c = C.an_element()
        sage: TestSuite(c).run()
    """
    def __init__(self, list value):
        """
        Initialize ``self``.
        """
        self.value = value

    def __nonzero__(self):
        return bool(self.value)

    def __repr__(self):
        return repr(self.value)

    def __cmp__(self, other):
        return cmp(self.value, other.value)

    def __reduce__(self):
        return LittelmannPath, (self.value,)

    cpdef copy(self):
        """
        Make a (deep) copy of ``self``.
        """
        return LittelmannPath([list(v) for v in self.value])

    cpdef list endpoint(self):
        r"""
        Computes the endpoint of the path.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b.endpoint()
            Lambda[1] + Lambda[2]
            sage: b.f_string([1,2,2,1])
            (-Lambda[1] - Lambda[2],)
            sage: b.f_string([1,2,2,1]).endpoint()
            -Lambda[1] - Lambda[2]
            sage: b.f_string([1,2])
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            sage: b.f_string([1,2]).endpoint()
            0
            sage: b = C([])
            sage: b.endpoint()
            0
        """
        if not self.value:
            return None
        ret = [0]*len(self.value[0])
        for v in self.value:
            add_inplace_lists(ret, v)
        return ret

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef compress(self):
        r"""
        Merge consecutive positively parallel steps present in the path
        by mutating ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: Lambda = C.R.weight_space().fundamental_weights(); Lambda
            Finite family {1: Lambda[1], 2: Lambda[2]}
            sage: c = C(tuple([1/2*Lambda[1]+1/2*Lambda[2], 1/2*Lambda[1]+1/2*Lambda[2]]))
            sage: c.compress()
            (Lambda[1] + Lambda[2],)
        """
        if not self.value:
            return self
        cdef list curr = self.value[0]
        cdef int pos = 0
        cdef list v
        while pos < len(self.value) - 1:
            v = self.value[pos+1]
            if positively_parallel_weights(curr, v):
                add_inplace_lists(curr, v)
                self.value.pop(pos+1)
            else:
                curr = v
                pos += 1

    cdef split_step(self, int which_step, r):
        r"""
        Split indicated step into two parallel steps of relative
        lengths `r` and `1-r` by mutating ``self``.

        INPUT:

        - ``which_step`` -- a position in the tuple ``self``
        - ``r`` -- a rational number between 0 and 1

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b.split_step(0, 1/3)
            (1/3*Lambda[1] + 1/3*Lambda[2], 2/3*Lambda[1] + 2/3*Lambda[2])
        """
        cdef list v = self.value[which_step]
        cdef list s1,s2
        s1 = [r*a for a in v]
        r = r.parent().one() - r
        s2 = [r*a for a in v]
        self.value[which_step] = s1
        self.value.insert(which_step+1, s2)
        #self.value = self.value[:which_step] + [s1, s2] + self.value[which_step+1:]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef reflect_step(self, int which_step, int i, list root):
        r"""
        Apply the `i`-th simple reflection to the indicated step in ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b.reflect_step(0,1)
            (-Lambda[1] + 2*Lambda[2],)
            sage: b.reflect_step(0,2)
            (2*Lambda[1] - Lambda[2],)
        """
        c = -self.value[which_step][i]
        inplace_axpy(c, root, self.value[which_step])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef list _string_data(self, int i):
        r"""
        Computes the `i`-string data of ``self``.

        TESTS::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b._string_data(1)
            ()
            sage: b._string_data(2)
            ()
            sage: b.f(1)._string_data(1)
            ((0, -1, -1),)
            sage: b.f(1).f(2)._string_data(2)
            ((0, -1, -1),)
        """
        if not self.value:
            return []
        # Get the wet step data
        cdef list minima_pos = []
        cdef int ix
        ps = 0
        psmin = 0
        for ix in range(len(self.value)):
            # Compute the i-height of the step
            step = self.value[ix][i]
            ps += step
            if ps < psmin:
                minima_pos.append((ix, ps, step))
                psmin = ps
        return minima_pos

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef LittelmannPath dualize(self):
        r"""
        Return the dualized path of ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: for c in C:
            ....:   print c, c.dualize()
            (Lambda[1] + Lambda[2],) (-Lambda[1] - Lambda[2],)
            (-Lambda[1] + 2*Lambda[2],) (Lambda[1] - 2*Lambda[2],)
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]) (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            (Lambda[1] - 2*Lambda[2],) (-Lambda[1] + 2*Lambda[2],)
            (-Lambda[1] - Lambda[2],) (Lambda[1] + Lambda[2],)
            (2*Lambda[1] - Lambda[2],) (-2*Lambda[1] + Lambda[2],)
            (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2]) (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2])
            (-2*Lambda[1] + Lambda[2],) (2*Lambda[1] - Lambda[2],)
        """
        if not self.value:
            return self
        self.value.reverse()
        cdef int n = len(self.value[0])
        cdef int i, j
        cdef list row
        for i in range(len(self.value)):
            row = self.value[i]
            for j in range(n):
                row[j] = -row[j]
        return self

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef int epsilon(self, int i):
        r"""
        Returns the distance to the beginning of the `i`-string.

        This method overrides the generic implementation in the category of crystals
        since this computation is more efficient.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: [c.epsilon(1) for c in C]
            [0, 1, 0, 0, 1, 0, 1, 2]
            sage: [c.epsilon(2) for c in C]
            [0, 0, 1, 2, 1, 1, 0, 0]
        """
        if not self.value:
            return 0
        cdef int ix
        ps = 0
        psmin = 0
        for ix in range(len(self.value)):
            # Compute the i-height of the step
            ps += self.value[ix][i]
            if ps < psmin:
                psmin = ps
        return floor(-psmin)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef int phi(self, int i):
        r"""
        Returns the distance to the end of the `i`-string.

        This method overrides the generic implementation in the category of crystals
        since this computation is more efficient.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: [c.phi(1) for c in C]
            [1, 0, 0, 1, 0, 2, 1, 0]
            sage: [c.phi(2) for c in C]
            [1, 2, 1, 0, 0, 0, 0, 1]
        """
        if not self.value:
            return 0
        cdef int ix
        ps = 0
        psmin = 0
        for ix in range(len(self.value)-1,-1,-1):
            # Compute the i-height of the step
            ps -= self.value[ix][i]
            if ps < psmin:
                psmin = ps
        return floor(-psmin)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef e(self, int i, list root, int power=1, to_string_end=False):
        r"""
        Return the `i`-th crystal raising operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set of the underlying root system
        - ``root`` -- the ``i``-th simple root as a list
        - ``power`` -- positive integer; specifies the power of the raising operator
          to be applied (default: 1)
        - ``to_string_end`` -- boolean; if set to True, returns the dominant end of the
          `i`-string of ``self``. (default: False)
        - ``length_only`` -- boolean; if set to True, returns the distance to the dominant
          end of the `i`-string of ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C[2]; c
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            sage: c.e(1)
            sage: c.e(2)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.e(2,to_string_end=True)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.e(1,to_string_end=True)
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            sage: c.e(1,length_only=True)
            0
        """
        cdef list data = self._string_data(i)
        # compute the minimum i-height M on the path
        cdef int p
        if not data:
            M = 0
        else:
            M = data[len(data)-1][1]
        cdef int max_raisings = floor(-M)

        # set the power of e_i to apply
        if to_string_end:
            p = max_raisings
        else:
            p = power
        if p > max_raisings:
            return None

        cdef int ix = len(data) - 1
        cdef int j
        while ix >= 0 and data[ix][1] < M + p:
        # get the index of the current step to be processed
            j = data[ix][0]
            # find the i-height where the current step might need to be split
            if ix == 0:
                prev_ht = M + p
            else:
                prev_ht = min(data[ix-1][1], M+p)
            # if necessary split the step, then reflect the wet part
            if data[ix][1] - data[ix][2] > prev_ht:
                self.split_step(j, 1 - (prev_ht - data[ix][1]) / -data[ix][2])
                self.reflect_step(j + 1, i, root)
            else:
                self.reflect_step(j, i, root)
            ix = ix - 1
        self.compress()
        return self

    cpdef f(self, int i, list root, int power=1, to_string_end=False):
        r"""
        Return the `i`-th crystal lowering operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set of the underlying root system
        - ``power`` -- positive integer; specifies the power of the lowering operator
          to be applied (default: 1)
        - ``to_string_end`` -- boolean; if set to True, returns the anti-dominant end of the
          `i`-string of ``self``. (default: False)

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C.module_generators[0]
            sage: c.f(1)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.f(1,power=2)
            sage: c.f(2)
            (2*Lambda[1] - Lambda[2],)
            sage: c.f(2,to_string_end=True)
            (2*Lambda[1] - Lambda[2],)
            sage: c.f(2,length_only=True)
            1

            sage: C = crystals.LSPaths(['A',2,1],[-1,-1,2])
            sage: c = C.module_generators[0]
            sage: c.f(2,power=2)
            (Lambda[0] + Lambda[1] - 2*Lambda[2],)
        """
        dual_path = self.dualize().e(i, root, power, to_string_end)
        if dual_path is None:
            return None
        return dual_path.dualize()

    cpdef LittelmannPath s(self, int i, list root):
        r"""
        Compute the reflection of ``self`` along the `i`-string.

        This method is more efficient than the generic implementation since
        it uses powers of `e` and `f` in the Littelmann model directly.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C.module_generators[0]
            sage: c.s(1)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.s(2)
            (2*Lambda[1] - Lambda[2],)

            sage: C = crystals.LSPaths(['A',2,1],[-1,0,1])
            sage: c = C.module_generators[0]; c
            (-Lambda[0] + Lambda[2],)
            sage: c.s(2)
            (Lambda[1] - Lambda[2],)
            sage: c.s(1)
            (-Lambda[0] + Lambda[2],)
            sage: c.f(2).s(1)
            (Lambda[0] - Lambda[1],)
        """
        cdef int diff = self.phi(i) - self.epsilon(i)
        if diff >= 0:
            return self.f(i, root, power=diff)
        else:
            return self.e(i, root, power=-diff)


#####################################################################
## B(\infty)


cdef class InfinityLittelmannPath(LittelmannPath):
    def __init__(self, list value, list rho):
        LittelmannPath.__init__(self, value)
        self.rho = rho

    def __reduce__(self):
        return InfinityLittelmannPath, (self.value, self.rho)

    cpdef copy(self):
        """
        Make a (deep) copy of ``self``.
        """
        return InfinityLittelmannPath([list(v) for v in self.value], self.rho)

    def e(self, int i, list root, int power=1, to_string_end=False):
        r"""
        Return the `i`-th crystal raising operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set
        - ``power`` -- (default: 1) positive integer; specifies the
          power of the lowering operator to be applied
        - ``length_only`` -- (default: ``False``) boolean; if ``True``,
          then return the distance to the anti-dominant end of the
          `i`-string of ``self``

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['B',3,1])
            sage: mg = B.module_generator()
            sage: mg.e(0)
            sage: mg.e(1)
            sage: mg.e(2)
            sage: x = mg.f_string([1,0,2,1,0,2,1,1,0])
            sage: all(x.f(i).e(i) == x for i in B.index_set())
            True
            sage: all(x.e(i).f(i) == x for i in B.index_set() if x.epsilon(i) > 0)
            True

        TESTS:

        Check that this works in affine types::

            sage: B = crystals.infinity.LSPaths(['A',3,1])
            sage: mg = B.highest_weight_vector()
            sage: x = mg.f_string([0,1,2,3])
            sage: x.e_string([3,2,1,0]) == mg
            True

        We check that :meth:`epsilon` works::

            sage: B = crystals.infinity.LSPaths(['D',4])
            sage: mg = B.highest_weight_vector()
            sage: x = mg.f_string([1,3,4,2,4,3,2,1,4])
            sage: [x.epsilon(i) for i in B.index_set()]
            [1, 1, 0, 1]

        Check that :trac:`21671` is fixed::

            sage: B = crystals.infinity.LSPaths(['G',2])
            sage: len(B.subcrystal(max_depth=7))
            116
        """
        if LittelmannPath.e(self, i, root, power=power) is None:
            return None

        cdef int n = len(root)
        cdef list endpoint = self.endpoint()
        cdef int j
        cdef list zero = [0]*len(self.rho)
        # Do not count \delta
        if self.rho[-1] == 0:
            n -= 1

        if not positively_parallel_weights(self.value[-1], self.rho):
            self.value.append(list(self.rho))
            add_inplace_lists(endpoint, self.rho)

        while any(endpoint[j] < 1 for j in xrange(n)):
            add_inplace_lists(self.value[-1], self.rho)
            add_inplace_lists(endpoint, self.rho)
        mrho = [-val for val in self.rho]  # Negate rho
        while all(endpoint[j] > 1 for j in xrange(n)) and self.value[-1] != zero:
            add_inplace_lists(self.value[-1], mrho)
            add_inplace_lists(endpoint, mrho)
        while self.value[-1] == zero:
            self.value.pop()
        return self

    def f(self, int i, list root, int power=1, to_string_end=False):
        r"""
        Return the `i`-th crystal lowering operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set
        - ``power`` -- (default: 1) positive integer; specifies the
          power of the lowering operator to be applied
        - ``length_only`` -- (default: ``False``) boolean; if ``True``,
          then return the distance to the anti-dominant end of the
          `i`-string of ``self``

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['D',3,2])
            sage: mg = B.highest_weight_vector()
            sage: mg.f(1)
            (3*Lambda[0] - Lambda[1] + 3*Lambda[2],
             2*Lambda[0] + 2*Lambda[1] + 2*Lambda[2])
            sage: mg.f(2)
            (Lambda[0] + 2*Lambda[1] - Lambda[2],
             2*Lambda[0] + 2*Lambda[1] + 2*Lambda[2])
            sage: mg.f(0)
            (-Lambda[0] + 2*Lambda[1] + Lambda[2] - delta,
             2*Lambda[0] + 2*Lambda[1] + 2*Lambda[2])
        """
        self.dualize()
        if LittelmannPath.e(self, i, root, power) is None:
            return None
        self.dualize()

        cdef int n = len(root)
        cdef list endpoint = self.endpoint()
        cdef int j
        cdef list zero = [0]*len(self.rho)
        # Do not count \delta
        if self.rho[-1] == 0:
            n -= 1

        if not positively_parallel_weights(self.value[-1], self.rho):
            self.value.append(list(self.rho))
            add_inplace_lists(endpoint, self.rho)

        while any(endpoint[j] < 1 for j in xrange(n)):
            add_inplace_lists(self.value[-1], self.rho)
            add_inplace_lists(endpoint, self.rho)
        mrho = [-val for val in self.rho]  # Negate rho
        while all(endpoint[j] > 1 for j in xrange(n)) and self.value[-1] != zero:
            add_inplace_lists(self.value[-1], mrho)
            add_inplace_lists(endpoint, mrho)
        while self.value[-1] == zero:
            self.value.pop()
        return self


#####################################################################
## Helper functions

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline add_inplace_lists(list l, list a):
    for i, c in enumerate(a):
        l[i] += c

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline inplace_axpy(a, list x, list y):
    for i, c in enumerate(x):
        y[i] += a * c

@cython.boundscheck(False)
@cython.wraparound(False)
cdef bint positively_parallel_weights(list v, list w):
    """
    Check whether the vectors ``v`` and ``w`` are positive scalar
    multiples of each other.

    EXAMPLES::

        sage: from sage.combinat.crystals.littelmann_path import positively_parallel_weights
        sage: La = RootSystem(['A',5,2]).weight_space(extended=True).fundamental_weights()
        sage: rho = sum(La)
        sage: positively_parallel_weights(rho, 4*rho)
        True
        sage: positively_parallel_weights(4*rho, rho)
        True
        sage: positively_parallel_weights(rho, -rho)
        False
        sage: positively_parallel_weights(rho, La[1] + La[2])
        False
    """
    cdef int i
    for i,vi in enumerate(v):
        if vi == 0:
            continue
        return vi*w[i] > 0 and [vi*a for a in w] == [w[i]*a for a in v]
    return False

#####################################################################
## Element classes

cdef class LSPathElement(ElementWrapper):
    """
    TESTS::

        sage: C = crystals.LSPaths(['E',6],[1,0,0,0,0,0])
        sage: c = C.an_element()
        sage: TestSuite(c).run()
    """
    def __hash__(self):
        return hash(tuple([tuple(v) for v in self.value.value]))

    def __iter__(self):
        WLR = self._parent.weight_lattice_realization()
        I = WLR.basis().keys()
        for v in self.value.value:
            yield WLR._from_dict({I[i]: c for i,c in enumerate(v)})

    def endpoint(self):
        r"""
        Computes the endpoint of the path.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b.endpoint()
            Lambda[1] + Lambda[2]
            sage: b.f_string([1,2,2,1])
            (-Lambda[1] - Lambda[2],)
            sage: b.f_string([1,2,2,1]).endpoint()
            -Lambda[1] - Lambda[2]
            sage: b.f_string([1,2])
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            sage: b.f_string([1,2]).endpoint()
            0
            sage: b = C([])
            sage: b.endpoint()
            0
        """
        WLR = self._parent.weight_lattice_realization()
        if not self.value.endpoint():
            return WLR.zero()
        I = WLR.basis().keys()
        return WLR._from_dict({I[i]: c for i,c in enumerate(self.value.endpoint())})

    def compress(self):
        r"""
        Merges consecutive positively parallel steps present in the path.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: Lambda = C.R.weight_space().fundamental_weights(); Lambda
            Finite family {1: Lambda[1], 2: Lambda[2]}
            sage: c = C(tuple([1/2*Lambda[1]+1/2*Lambda[2], 1/2*Lambda[1]+1/2*Lambda[2]]))
            sage: c.compress()
            (Lambda[1] + Lambda[2],)
        """
        # TODO: Deprecate
        return self

    def split_step(self, which_step, r):
        r"""
        Splits indicated step into two parallel steps of relative lengths `r` and `1-r`.

        INPUT:

        - ``which_step`` -- a position in the tuple ``self``
        - ``r`` -- a rational number between 0 and 1

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b.split_step(0, 1/3)
            (1/3*Lambda[1] + 1/3*Lambda[2], 2/3*Lambda[1] + 2/3*Lambda[2])
        """
        # TODO: Deprecate
        assert 0 <= which_step and which_step <= len(self.value)
        v = self.value.value[which_step]
        data = self.value.value[:which_step] + [r*v,(1-r)*v] + self.value.value[which_step+1:]
        cdef LittelmannPath L = LittelmannPath(data)
        return self._parent.element_class(self._parent, L)

    def reflect_step(self, which_step, i):
        r"""
        Apply the `i`-th simple reflection to the indicated step in ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: b = C.module_generators[0]
            sage: b.reflect_step(0, 1)
            (-Lambda[1] + 2*Lambda[2],)
            sage: b.reflect_step(0, 2)
            (2*Lambda[1] - Lambda[2],)
        """
        assert 0 <= which_step and which_step <= len(self.value)
        i = self._parent._inverse_index_map[i]
        root = self._parent._simple_root_as_list(i)
        cdef LittelmannPath data = self.value.copy()
        data.reflect_step(which_step, i, root)
        return self._parent.element_class(self._parent, data)

    def epsilon(self, i):
        r"""
        Returns the distance to the beginning of the `i`-string.

        This method overrides the generic implementation in the category of crystals
        since this computation is more efficient.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: [c.epsilon(1) for c in C]
            [0, 1, 0, 0, 1, 0, 1, 2]
            sage: [c.epsilon(2) for c in C]
            [0, 0, 1, 2, 1, 1, 0, 0]
        """
        i = self._parent._inverse_index_map[i]
        return self.value.epsilon(i)

    def phi(self, i):
        r"""
        Returns the distance to the end of the `i`-string.

        This method overrides the generic implementation in the category of crystals
        since this computation is more efficient.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: [c.phi(1) for c in C]
            [1, 0, 0, 1, 0, 2, 1, 0]
            sage: [c.phi(2) for c in C]
            [1, 2, 1, 0, 0, 0, 0, 1]
        """
        i = self._parent._inverse_index_map[i]
        return self.value.phi(i)

    def e(self, i, power=1, to_string_end=False, length_only=False):
        r"""
        Returns the `i`-th crystal raising operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set of the underlying root system
        - ``power`` -- positive integer; specifies the power of the raising operator
          to be applied (default: 1)
        - ``to_string_end`` -- boolean; if set to True, returns the dominant end of the
          `i`-string of ``self``. (default: False)
        - ``length_only`` -- boolean; if set to True, returns the distance to the dominant
          end of the `i`-string of ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C[2]; c
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            sage: c.e(1)
            sage: c.e(2)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.e(2,to_string_end=True)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.e(1,to_string_end=True)
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            sage: c.e(1,length_only=True)
            0
        """
        if length_only:
            return self.epsilon(i)
        root = self._parent._simple_root_as_list(i)
        i = self._parent._inverse_index_map[i]
        ret = self.value.copy().e(i, root, power, to_string_end)
        if ret is None:
            return None
        return self._parent.element_class(self._parent, ret)

    def dualize(self):
        r"""
        Returns dualized path.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: for c in C:
            ....:     print("{} {}".format(c, c.dualize()))
            (Lambda[1] + Lambda[2],) (-Lambda[1] - Lambda[2],)
            (-Lambda[1] + 2*Lambda[2],) (Lambda[1] - 2*Lambda[2],)
            (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]) (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
            (Lambda[1] - 2*Lambda[2],) (-Lambda[1] + 2*Lambda[2],)
            (-Lambda[1] - Lambda[2],) (Lambda[1] + Lambda[2],)
            (2*Lambda[1] - Lambda[2],) (-2*Lambda[1] + Lambda[2],)
            (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2]) (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2])
            (-2*Lambda[1] + Lambda[2],) (2*Lambda[1] - Lambda[2],)
        """
        if not self.value:
            return self
        cdef LittelmannPath data = self.value.copy()
        data.dualize()
        return self._parent.element_class(self._parent, data)

    def f(self, i, power=1, to_string_end=False, length_only=False):
        r"""
        Returns the `i`-th crystal lowering operator on ``self``.

        INPUT:

        - ``i`` -- element of the index set of the underlying root system
        - ``power`` -- positive integer; specifies the power of the lowering operator
          to be applied (default: 1)
        - ``to_string_end`` -- boolean; if set to True, returns the anti-dominant end of the
          `i`-string of ``self``. (default: False)
        - ``length_only`` -- boolean; if set to True, returns the distance to the anti-dominant
          end of the `i`-string of ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C.module_generators[0]
            sage: c.f(1)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.f(1,power=2)
            sage: c.f(2)
            (2*Lambda[1] - Lambda[2],)
            sage: c.f(2,to_string_end=True)
            (2*Lambda[1] - Lambda[2],)
            sage: c.f(2,length_only=True)
            1

            sage: C = crystals.LSPaths(['A',2,1],[-1,-1,2])
            sage: c = C.module_generators[0]
            sage: c.f(2,power=2)
            (Lambda[0] + Lambda[1] - 2*Lambda[2],)
        """
        if length_only:
            return self.phi(i)
        root = self._parent._simple_root_as_list(i)
        i = self._parent._inverse_index_map[i]
        ret = self.value.copy().f(i, root, power, to_string_end)
        if ret is None:
            return None
        return self._parent.element_class(self._parent, ret)

    def s(self, i):
        r"""
        Computes the reflection of ``self`` along the `i`-string.

        This method is more efficient than the generic implementation since it uses
        powers of `e` and `f` in the Littelmann model directly.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C.module_generators[0]
            sage: c.s(1)
            (-Lambda[1] + 2*Lambda[2],)
            sage: c.s(2)
            (2*Lambda[1] - Lambda[2],)

            sage: C = crystals.LSPaths(['A',2,1],[-1,0,1])
            sage: c = C.module_generators[0]; c
            (-Lambda[0] + Lambda[2],)
            sage: c.s(2)
            (Lambda[1] - Lambda[2],)
            sage: c.s(1)
            (-Lambda[0] + Lambda[2],)
            sage: c.f(2).s(1)
            (Lambda[0] - Lambda[1],)
        """
        root = self._parent._simple_root_as_list(i)
        i = self._parent._inverse_index_map[i]
        return self._parent.element_class(self._parent, self.value.copy().s(i, root))

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: B = crystals.LSPaths(['A',1,1],[1,0])
            sage: b = B.highest_weight_vector()
            sage: b.f(0).weight()
            -Lambda[0] + 2*Lambda[1] - delta
        """
        return self.endpoint()

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        WLR = self._parent.weight_lattice_realization()
        I = WLR.basis().keys()
        return repr(tuple([WLR._from_dict({I[i]: c for i,c in enumerate(v)})
                           for v in self.value.value]))

    def _latex_(self):
        r"""
        Latex method for ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: c = C.module_generators[0]
            sage: c._latex_()
            \left[\Lambda_{1} + \Lambda_{2}\right]
        """
        WLR = self._parent.weight_lattice_realization()
        I = WLR.basis().keys()
        from sage.misc.latex import latex
        return latex([WLR._from_dict({I[i]: c for i,c in enumerate(v)})
                      for v in self.value.value])

