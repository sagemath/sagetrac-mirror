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
            sage: b.split_step(0,1/3)
            (1/3*Lambda[1] + 1/3*Lambda[2], 2/3*Lambda[1] + 2/3*Lambda[2])
        """
        cdef list v = self.value[which_step]
        cdef list s1,s2
        s1 = [r*a for a in v]
        r = 1 - r
        s2 = [r*a for a in v]
        self.value = self.value[:which_step] + [s1, s2] + self.value[which_step+1:]

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
    cpdef LittelmannPath dualize(self, inplace=True):
        r"""
        Return the dualized path of ``self``.

        EXAMPLES::

            sage: C = crystals.LSPaths(['A',2],[1,1])
            sage: for c in C:
            ...     print c, c.dualize()
            ...
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
        if not inplace:
            return LittelmannPath([[-c for c in v] for v in reversed(self.value)])
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


class InfinityCrystalOfLSPathsElement(LittelmannPath):
    def e(self, int i, list root, int power=1):
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
        """
        ret = LittelmannPath.e(self, i, root, power=power)
        if ret is None:
            return None
        WLR = self.parent().weight_lattice_realization()
        value = list(ret.value)
        endpoint = sum(p for p in value)
        rho = WLR.rho()
        h = WLR.simple_coroots()
        I = self.parent().index_set()

        if not positively_parallel_weights(value[-1], rho):
            value.append(rho)
            endpoint += rho

        while any(endpoint.scalar(alc) < 1 for alc in h):
            value[-1] += rho
            endpoint += rho
        while all(endpoint.scalar(alc) > 1 for alc in h):
            value[-1] -= rho
            endpoint -= rho
        while value[-1] == WLR.zero():
            value.pop()
        ret.value = tuple(value)
        return ret

    def f(self, int i, list root, int power=1):
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
        dual_path = self.dualize()
        dual_path = LittelmannPath.e(dual_path, i, root, power)
        if dual_path is None:
            return None
        ret = dual_path.dualize()
        WLR = self.parent().weight_lattice_realization()
        value = list(ret.value)
        endpoint = sum(p for p in value)
        rho = WLR.rho()
        h = WLR.simple_coroots()

        if not positively_parallel_weights(value[-1], rho):
            value.append(rho)
            endpoint += rho

        while any(endpoint.scalar(alc) < 1 for alc in h):
            value[-1] += rho
            endpoint += rho
        while all(endpoint.scalar(alc) > 1 for alc in h):
            value[-1] -= rho
            endpoint -= rho
        while value[-1] == WLR.zero():
            value.pop()
        ret.value = tuple(value)
        return ret

    def weight(self):
        """
        Return the weight of ``self``.

        .. TODO::

            This is a generic algorithm. We should find a better
            description and implement it.

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['E',6])
            sage: mg = B.highest_weight_vector()
            sage: f_seq = [1,4,2,6,4,2,3,1,5,5]
            sage: x = mg.f_string(f_seq)
            sage: x.weight()
            -3*Lambda[1] - 2*Lambda[2] + 2*Lambda[3] + Lambda[4] - Lambda[5]

            sage: al = B.cartan_type().root_system().weight_space().simple_roots()
            sage: x.weight() == -sum(al[i] for i in f_seq)
            True
        """
        WLR = self.parent().weight_lattice_realization()
        alpha = WLR.simple_roots()
        return -WLR.sum(alpha[i] for i in self.to_highest_weight()[1])

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        Let `\pi \in \mathcal{B}(\infty)`. Define

        .. MATH::

            \varphi_i(\pi) := \varepsilon_i(\pi) + \langle h_i,
            \mathrm{wt}(\pi) \rangle,

        where `h_i` is the `i`-th simple coroot and `\mathrm{wt}(\pi)`
        is the :meth:`weight` of `\pi`.

        INPUT:

        - ``i`` -- element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.LSPaths(['D',4])
            sage: mg = B.highest_weight_vector()
            sage: x = mg.f_string([1,3,4,2,4,3,2,1,4])
            sage: [x.phi(i) for i in B.index_set()]
            [-1, 4, -2, -3]
        """
        WLR = self.parent().weight_lattice_realization()
        h = WLR.simple_coroots()
        return self.epsilon(i) + WLR(self.weight()).scalar(h[i])


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
cdef positively_parallel_weights(list v, list w):
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

