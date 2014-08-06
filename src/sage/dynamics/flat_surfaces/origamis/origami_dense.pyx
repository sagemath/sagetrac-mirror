"""
Dense origamis

An *origami* is a pair `(r,u)` of permutations up to conjugacy. The computation
of canonic representative of an origami is performed by an independant C
program in normal_form.c

The two permutations r and u are contiguous in memory and for each
initialization we perform only one call to sage_malloc.

A *pillowcase cover* is a quadruple `(g_0, g_1, g_2, g_3)` of permutations such
that the product `g_0 g_1 g_2 g_3` is the identity. It is sometimes called a
4-constellations. The computation of canonic representative is performed by an
independant C program in noraml_form.c
"""

include "../../../ext/stdsage.pxi"
include "../../../ext/interrupt.pxi"
from cpython.list cimport *
from cpython.tuple cimport *

#cimport list
#cimport tuple

from cpython cimport bool

from sage.rings.integer cimport Integer, smallInteger
from sage.rings.integer import GCD_list

from libc.string cimport memset, memcpy
from libc.limits cimport UINT_MAX

cdef extern from "normal_form.h":
    int origami_normal_form(int *x, int *y, int *ren, unsigned int n)
    int pillowcase_cover_normal_form(int *g, int *ren, unsigned int n)
    inline int origami_diff(int *o1, int *o2, unsigned int n)
    inline int pillowcase_cover_diff(int *g1, int *g2, unsigned int n)

cdef extern from "lyapunov_exponents.h":
    struct origami_data
    origami_data * new_origami_data(size_t degree, size_t nb_vectors, int *pa, int *pb)
    void free_origami_data(origami_data *o)
    void lyapunov_exponents(origami_data *o, size_t nb_iterations, double *theta)

    struct origami_with_involution_data
    origami_with_involution_data * new_origami_with_involution_data(size_t degree, size_t nb_vectors_p, size_t nb_vectors_m, int *pa, int *pb, int *s)
    void free_origami_with_involution_data(origami_with_involution_data *o)
    void lyapunov_exponents_with_involution(origami_with_involution_data *o, size_t NB_ITERATIONS, double * ttheta)


cdef inline tuple array_to_tuple(int * x, unsigned int n):
# the following raise a SIGSEGV error
#    cdef int i
#    cdef tuple res = PyTuple_New(<Py_ssize_t> n)
#
#    for i from 0 <= i < n:
#        PyTuple_SetItem(res, i, PyInt_FromLong(x[i]))
#
#    return res
    cdef list res = [None]*n
    cdef int i

    for i from 0 <= i < n:
        res[i] = x[i]

    return tuple(res)

cdef inline tuple array_to_tuple_i(int * x, unsigned int n):
#    cdef int i
#    cdef tuple res = PyTuple_New(<Py_ssize_t> n)
#
#    for i from 0 <= i < n:
#        PyTuple_SetItem(res, x[i], PyInt_FromLong(i))
#
#    return res
    cdef list res = [None]*n
    cdef int i

    for i from 0 <= i < n:
        res[x[i]] = i

    return tuple(res)

cdef tuple projectivize_edges(l):
    r"""
    Action of `PGL(2,\ZZ)` or `PSL(2,\ZZ)` knowing the one of `GL(2,\ZZ)` or
    `SL(2,\ZZ)`.

    INPUT:

    - ``l`` - a list of dictionnaries
    """
    cdef int i
    cdef Origami_dense_pyx o1, oo1, o2, oo2, ooo
    cdef dict ll
    cdef set waiting
    cdef int *xx
    cdef int *yy

    o1 = l[0].iterkeys().next() # pick a random element in l[0] !!!!
    oo1 = o1.inverse()
    oo1._set_standard_form()

    # (trivial) case 1: -Id preserves o (and hence preserves pointwise the
    # orbit)
    # (rk: this case corresponds to orientation cover orbits)
    if o1 == oo1:
        return l

    # case 2: -Id does not preserve the orbit
    elif oo1 not in l[0]:
        ll = [{} for _ in xrange(len(l))]
        for o in l[0]:
            o1 = o
            o2 = o.inverse()
            o2._set_standard_form()
            if o2 < o1:
                o = o2
            
            for i in xrange(len(l)):
                oo1 = l[i][o1]
                oo2 = oo1.inverse()
                oo2._set_standard_form()
                if oo1 < oo2:
                    ll[i][o] = oo1
                else:
                    ll[i][o] == oo2

    # case 3: -Id preserve the orbit
    ll = [{} for _ in xrange(len(l))]
    waiting = set(l[0])
    while waiting:
        o1 = waiting.pop()
        o2 = o1.inverse()
        o2._set_standard_form()
        waiting.remove(o2)
        if o2 < o1:
            ooo = o1
            o1 = o2
            o2 = ooo

        for i in xrange(len(l)):
            oo1 = l[i][o1]
            oo2 = l[i][o2]
            if oo1 < oo2:
                ll[i][o1] = oo1
            else:
                ll[i][o1] = oo2

    return ll

def sl_orbit_from_gl_orbit(o,L,I):
    r"""
    Compute the sl2z orbit of the origami ``o`` knowing the action of gl2z.

    TODO: this has nothing to do with origamis... but rather to the action of
    SL(2,Z)/GL(2,Z)/PSL(2,Z).

    INPUT:

    - ``o`` - an origami

    - ``L`` - the action of the matrix l

    - ``I`` - the action of the matrix i

    EXAMPLES:

    On the following example, the SL(2,Z) action has two orbits whereas the
    GL(2,Z) action as only one::

        sage: l_edges = {0:1,1:9,2:8,3:0,4:13,5:10,6:5,7:4,8:11,9:3,10:6,11:2,12:7,13:12}
        sage: i_edges = {0:10,1:5,2:12,3:4,4:3,5:1,6:11,7:9,8:13,9:7,10:0,11:6,12:2,13:8}
        sage: from sage.dynamics.flat_surfaces.origamis.origami_dense import sl_orbit_from_gl_orbit
        sage: s0 = sl_orbit_from_gl_orbit(0, l_edges, i_edges); s0
        ({0: 1, 1: 9, 2: 8, 3: 0, 8: 11, 9: 3, 11: 2},
         {0: 11, 1: 0, 2: 9, 3: 8, 8: 2, 9: 3, 11: 1},
         {0: 2, 1: 8, 2: 0, 3: 9, 8: 1, 9: 3, 11: 11})
        sage: s0 == sl_orbit_from_gl_orbit(1, l_edges, i_edges)
        True
        sage: s4 = sl_orbit_from_gl_orbit(4, l_edges, i_edges)
        sage: s0 == s4
        False
    """
    l = {}
    r = {}
    s = {}
    ri = {}

    waiting = set([o])

    while True:
        # compute L orbit and put new guys in waiting
        ooo = o
        oo = L[o]
        while oo != o:
            l[ooo] = oo
            waiting.add(oo)
            ooo = oo
            oo = L[oo]
        l[ooo] = oo

        # compute R images while we do not find somebody new
        while waiting:
            oo = waiting.pop()
            o = I[L[I[oo]]]
            r[oo] = o
            ri[o] = oo
            if not (o in l):
                waiting.add(o)
                break
        else:
            break

    for o in l: # s = l ~r l
        s[o] = l[ri[l[o]]]

    return l,r,s

cdef class Origami_dense_pyx:
    r"""
    Origami or square tiled surface.

    An origami is a flat surface which is a covering of a one punctered torus.
    It can be described either by a couple of permutations up to conugacy (in
    the symmetric group) or by a subgroup of finite index of the free group on
    two generators.

    EXAMPLES::

        sage: Origami([2,1,3], [3,2,1])
        (1,2)(3)
        (1,3)(2)
    """
    #
    # Intialization and copy
    #

    def __cinit__(self):
        r"""
        TESTS::

            sage: o = Origami('(1,2)','(1,3)')
            sage: loads(dumps(o)) == o
            True
        """
        self._n = 0
        self._r = NULL
        self._u = NULL

    def __init__(self, r, u):
        r"""
        TESTS::

            sage: o = Origami([2,1,3],[1,3,2])
            sage: o == loads(dumps(o))
            True
        """
        cdef int i

        assert len(r) == len(u)

        self._n = len(r)
        self._r = <int *> sage_malloc(2*self._n*sizeof(int))
        self._u = self._r + self._n

        if self._r == NULL:
            raise MemoryError, "not able to allocate"
       
        for i from 0 <= i < self._n:
            self._r[i] = r[i]
            self._u[i] = u[i]

        self._l_edges = {}
        self._i_edges = {}

    def __dealloc__(self):
        if self._r != NULL: sage_free(self._r)

    cdef Origami_dense_pyx _new_c(self, int * rr_and_uu):
        r"""
        Return an origami with given permutations.

        Beware that we assume that the created origami is in the same orbit
        under the action of GL(2,Z)
        """
        cdef Origami_dense_pyx other = PY_NEW_SAME_TYPE(self)

        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__

        other._n = self._n
        other._r = rr_and_uu
        other._u = rr_and_uu + self._n
        
        other._l_edges = self._l_edges
        other._i_edges = self._i_edges

        return other

    def __copy__(self):
        r"""
        Return a copy of the origami

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: oo = copy(o)
            sage: o == oo
            True
            sage: o is oo
            False
        """
        cdef int * r = <int *> sage_malloc(2 * self._n * sizeof(int))
        cdef int * u = r + self._n

        memcpy(r,self._r, 2 * self._n * sizeof(int))

        return self._new_c(r)

    #
    # Comparisons
    #

    def __richcmp__(self,other,i):
        r"""
        Comparison

         0: <
         1: <=
         2: ==
         3: !=
         4: >
         5: >=
        
        EXAMPLES:

        First compare the number of squares::

            sage: o1 = Origami('(1,2)','(2,1)')
            sage: o2 = Origami('(1,2,3)','(1,3)')
            sage: o1 == o2
            False
            sage: o1 != o2
            True
            sage: (o1 < o2) and (o1 <= o2)
            True
            sage: (o1 > o2) or (o1 >= o2)
            False

        Then compare the permutations::

            sage: o1 = Origami('(1,2)','(1,2,3)')
            sage: o2 = Origami('(1,2)','(1,3,2)')
            sage: o1 == o2
            False
            sage: o1 != o2
            True
            sage: (o1 < o2) and (o1 <= o2)
            True
            sage: (o1 > o2) or (o1 >= o2)
            False
        """
        cdef Origami_dense_pyx s
        cdef Origami_dense_pyx o

        if not isinstance(other, Origami_dense_pyx):
            return NotImplemented

        s = <Origami_dense_pyx> self
        o = <Origami_dense_pyx> other

        # compare the number of squares
        if s._n != o._n:
            if i < 2: return s._n < o._n
            if i > 3: return s._n > o._n
            return i == 3

        # find the first index where self and other differ and make the
        # difference.
        test = origami_diff(s._r, o._r, s._n)

        if test == 0: # equality
            return i != 0 and i != 3 and i != 4
        else: # different
            if i < 2: return test < 0
            if i > 3: return test > 0
            return i == 3

    def __hash__(self):
        r"""
        Hash value for self

        TESTS::

            sage: h = []
            sage: from itertools import permutations
            sage: for p in permutations(range(5)):
            ....:     for q in permutations(range(5)):
            ....:         h.append(hash(Origami(p,q,as_tuple=True,check=False)))
            sage: len(h) == len(set(h))
            True
        """
        cdef int i, h=0, br=12, bu=37

        for i from 0 <= i < self._n:
            h += self._r[i]*br + self._u[i]*bu
            br *= 503
            bu *= 251

        return h

    #
    # Python access to attributes
    #

    def nb_squares(self):
        r"""
        Return the number of squares.

        EXAMPLES::

            sage: o = Origami('(1,2)', '(1,3)')
            sage: o.nb_squares()
            3
        """
        return self._n

    def r_tuple(self):
        r"""
        Return the right permutation of the origami as a tuple on {0,...,n-1}

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: o.r_tuple()
            (1, 0, 2)
        """
        return array_to_tuple(self._r, self._n)

    def r_inv_tuple(self):
        r"""
        Return the inverse of the right permutation as a tuple on {0,...,n-1}

        EXAMPLES::

            sage: o = Origami('(1,2,3)','(1,2)')
            sage: o.r_inv_tuple()
            (2, 0, 1)
        """
        return array_to_tuple_i(self._r, self._n)

    def u_tuple(self):
        r"""
        Return the up permutation of the origami as a tuple on {0,...,n-1}

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: o.u_tuple()
            (2, 1, 0)
        """
        return array_to_tuple(self._u, self._n)

    def u_inv_tuple(self):
        r"""
        Return the inverse of the up permutation as a tuple on {0,...,n-1}

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,2,3)')
            sage: o.u_inv_tuple()
            (2, 0, 1)
        """
        return array_to_tuple_i(self._u, self._n)

    def widths_and_heights(self):
        r"""
        Return the list of widths and heigths of cylinder.

        EXAMPLES::

            sage: Origami('(1,2)','(1,3)').widths_and_heights()
            [(1, 1), (2, 1)]
            sage: Origami('(1,2)(3,4)','(1,3,5)(2,4)').widths_and_heights()
            [(1, 1), (2, 2)]
            sage: Origami('(1,2)','(1,3,4)').widths_and_heights()
            [(1, 2), (2, 1)]
            sage: Origami('(1,2)(3,4)','(1,3,5,6)(2,4)').widths_and_heights()
            [(1, 2), (2, 2)]
        """
        cdef int * r = self._r
        cdef int * u = self._u
        cdef int * seen = <int *> sage_malloc(self._n * sizeof(int))
        cdef int w,h,i
        cdef list wh = []

        # compute the set of squares that are on a top of a cylinder
        # for each top we pick one square and record the width
        memset(seen, 0, self._n*sizeof(int))
        for i in range(self._n):
            if r[u[i]] != u[r[i]]:
                i = u[i]
                if not seen[i]:
                    w = 0
                    while seen[i] == 0:
                        seen[i] = 1
                        i = r[i]
                        w += 1
                    wh.append((w,i))

        for j,(w,i) in enumerate(wh):
            h = 1
            i = u[i]
            while seen[i] == 0:
                i = u[i]
                h += 1
            wh[j] = (smallInteger(w),smallInteger(h))

        sage_free(seen)
        return wh

    def period_generators(self):
        r"""
        Return a list of periods that generate the lattice of periods.

        EXAMPLES::

            sage: o = Origami('(1,3,6)(2,5,7)(4)', '(1,2,4,3,5,6,7)')
            sage: sorted(o.period_generators())
            [(-1, 2), (0, 1), (0, 2), (1, 0), (1, 0), (2, 0)]

            sage: r = '(1,17,6,18,10,8)(11,2,12,5,7,9)(15,16)(3,13)(14,4)'
            sage: u = '(1,11,15,3,14,10,7,6,12)(17,2,16,13,4,8,9,18,5)'
            sage: o = Origami(r,u)
            sage: o.stratum()
            H_2(2)
            sage: sorted(o.period_generators())
            [(-2, 2), (0, 2), (0, 3), (2, 0), (2, 0), (4, 0)]
        """
        # TODO: it is stupid as we do twice the job... in standard form we
        # already compute the singularities
        cdef int * r = self._r
        cdef int * u = self._u
        cdef int * memory = <int *> sage_malloc(2*self._n * sizeof(int))
        cdef int * br_sg = memory
        cdef int * i_to_tr = memory + self._n
         # array: i-> distance to the nearest square that has a singularity in
         # its top right corner on the left (or -1 if none)
        cdef int i,k
        cdef list periods = []

        # compute the set squares which have a singularity in either their top
        # right or bottom right corners
        memset(br_sg, 0, self._n*sizeof(int))
        for i in range(self._n):
            if r[u[i]] != u[r[i]]:
                br_sg[u[i]] = 1
#                print "tr/br at %d/%d"%(i,u[i])

        # now compute the horizontal saddles and build i_to_tr
#        print "compute horiz"
        memset(i_to_tr, -1, self._n*sizeof(int))
        for i in range(self._n):
            if br_sg[u[i]]:
                i_to_tr[i] = 0
                j = r[i]
                k = 1
                while br_sg[u[j]] == 0:
                    i_to_tr[j] = k
                    j = r[j]
                    k += 1
#                print "new (%d,%d)"%(k,0)
                periods.append((smallInteger(k),smallInteger(0)))

#        print "i_to_tr"
#        print "".join("%3d"%i_to_tr[i] for i in range(self._n))

        # now compute the period of vertical transversals in each cylinder
#        print "compute vert"
        for i in range(self._n):
            if br_sg[i] == 1:
                j = i
                k = 1
                while i_to_tr[j] == -1:
                    j = u[j]
                    k += 1
#                print "new (%d,%d)"%(-i_to_tr[j],k)
                periods.append((smallInteger(-i_to_tr[j]), k))

        sage_free(memory)
        return periods

    def lattice_of_periods(self):
        r"""
        Returns (a,t,u) where ((a,0),(t,u)) is a standard basis
        for the lattice of periods of self

        The lattice of periods of an origami is the sublattice of ZZ^2
        generated by the holonomy vectors of its saddle connections.
        Any sublattice of ZZ^2 has a standard basis consisting of
        a horizontal vector (a,0) and a nonhorizontal vector (t,u),
        where a, t, u are integers satisfying 0 <= t < a and 0 < t.

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: o.lattice_of_periods()
            (1, 0, 1)
            sage: r = '(1,2,3,4,5,6,7,8,9)(10,11,12)(13,14,15,16,17,18,19,20,21)'
            sage: u = '(1,14,21,19,8,3,10,5,12,4,11,6)(2,15,16,17,18,7)(9,13,20)'
            sage: oy = Origami(r,u)
            sage: oy.lattice_of_periods()
            (3, 2, 1)
        """
        cdef Integer ZZ_0 = smallInteger(0)
        cdef Integer ZZ_1 = smallInteger(1)
        cdef list periods = self.period_generators()
        cdef Integer a = GCD_list([x for x,y in periods if y == ZZ_0])
        cdef list w = sorted(set((h, t%a) for (t,h) in periods))
        cdef int i

        while len(w) > 1:
            while w[0][0] == ZZ_0:
                a = a.gcd(w[0][1])
                del w[0]
            if a == ZZ_1:
                if w[0] == (1,0):
                    return (ZZ_1,ZZ_0,ZZ_1)
                else:
                    return (ZZ_1,ZZ_0,GCD_list([x for x,y in w]))
            for i in range(1,len(w)):
                w[i] = ((w[i][0]%w[0][0], (w[i][1]-(w[i][0]//w[0][0])*w[0][1])%a))
            w = sorted(set(w))
        return (a,w[0][1],w[0][0])

    def is_reduced(self):
        r"""
        Test of reducibility

        An origami is reduced, if it is not a ramified cover of a bigger torus
        with only one ramification point. In other terms, it is equivalent to
        say that the period of the origami generates `\ZZ^2`.

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: o.is_reduced()
            True
            sage: o = Origami('(1,2,3,4)(5,6)','(1,5)(2,6)')
            sage: o.is_reduced()
            False
            sage: o = Origami('(1,2)(3,4)','(1,3,5,6)(2,4)')
            sage: o.is_reduced()
            False
            sage: o = Origami('(1,2,3,4)(5,6)','(1,5)(2,6)')
            sage: o.is_reduced()
            False
            sage: o = Origami('(1,2,3,4)(5,6)','(1,5)(2,6)')
            sage: o.is_reduced()
            False
        """
        return self.lattice_of_periods() == (smallInteger(1),smallInteger(0),smallInteger(1))

    #
    # standard form (canonic labels)
    #

    cpdef _set_standard_form(self,return_map=False):
        r"""
        Renumerote the origami in its standard form and the map associated to
        the renumerotation if ``return_map`` is set to True.

        INPUT:

        - ``return_map`` - boolean (default: False)

        EXAMPLES::

            sage: from sage.dynamics.flat_surfaces.origamis.origami_dense import Origami_dense_pyx
            sage: o = Origami_dense_pyx((1,0,2),(0,2,1))
            sage: o.r_tuple()
            (1, 0, 2)
            sage: o.u_tuple()
            (0, 2, 1)
            sage: o._set_standard_form()
            sage: o.r_tuple()
            (0, 2, 1)
            sage: o.u_tuple()
            (1, 0, 2)
        """
        cdef int *ren = <int *> sage_malloc(self._n * sizeof(int))
        cdef object m = None

        if self._n != 1:
            origami_normal_form(self._r,self._u,ren,self._n)

        if return_map:
            m = array_to_tuple(ren,self._n)

        sage_free(ren)
        return m

    # TODO: compute at the same time the lattice generated by the holonomies
    def to_standard_form(self, return_map=False):
        r"""
        Return an isomorphic origami in standard form.

        INPUT:

        - ``return_map`` - boolean (default: False) - if True return the
          associated mapping

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: oo,m = o.to_standard_form(return_map=True)
            sage: ~m * o.r() * m == oo.r()
            True
            sage: ~m * o.u() * m == oo.u()
            True


            sage: o = Origami('(1,2,3)(4,5)(6)(7,8,9)','(1,6,3,8,2,7,4,9)')
            sage: oo,m = o.to_standard_form(return_map=True)
            sage: ~m * o.r() * m == oo.r()
            True
            sage: ~m * o.u() * m == oo.u()
            True
        """
        oo = self.__copy__()

        if return_map:
            from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
            m = oo._set_standard_form(return_map=True)
            return oo, PermutationGroupElement([i+1 for i in m])
        else:
            oo._set_standard_form(return_map=False)
            return oo

    def relabel(self, return_map=False, inplace=False):
        r"""
        Relabel self

        INPUT:

        - ``return_map`` -- return the labelization

        - ``inplace`` -- modify self, the default is False. It might be
          dangerous to set it True as an origami is hashable and
          the hash is modified by this operation.

        EXAMPLES::

            sage: o = Origami('(1,2,3,4)(5,6)','(2,5)(3,6)')
            sage: o2,p = o.relabel(return_map=True)
            sage: old_r = o.r_tuple()
            sage: old_u = o.u_tuple()
            sage: new_r = o2.r_tuple()
            sage: new_u = o2.u_tuple()
            sage: all(p[old_r[i]] == new_r[p[i]] for i in xrange(6))
            True
            sage: all(p[old_u[i]] == new_u[p[i]] for i in xrange(6))
            True
        """
        if inplace:
            return self._set_standard_form(return_map=return_map)

        else:
            o = self.__copy__()
            m = o._set_standard_form(return_map=return_map)
            if return_map:
                return o,m
            return o



    #
    # GL(2,Z), SL(2,Z), PGL(2,Z) and PSL(2,Z) actions
    #

    cpdef inverse(self):
        r"""
        Return the origami `(r^{-1}, u^{-1})` which corresponds to the action of
        `-Id` on the origami.

        EXAMPLES::

            sage: o = Origami('(1,2,3)','(1,2)')
            sage: o
            (1,2,3)
            (1,2)(3)
            sage: o.inverse()
            (1,3,2)
            (1,2)(3)

        We check that it commutes with the inverse operation on cylinder
        diagrams::

            sage: o = Origami('(1,2,3)(5,6)', '(1,4)(3,5,7)')
            sage: os = o.inverse(); os
            (1,3,2)(4)(5,6)(7)
            (1,4)(2)(3,7,5)(6)

            sage: c = o.cylinder_diagram()
            sage: cs1 = o.cylinder_diagram().inverse()
            sage: cs1.relabel()
            sage: cs2 = os.cylinder_diagram()
            sage: cs2.relabel()
            sage: cs1 == cs2
            True
        """
        cdef int i
        cdef int *rr = <int *>sage_malloc(2*self._n*sizeof(int))
        cdef int *uu = rr + self._n
    
        for i from 0 <= i < self._n:
            rr[self._r[i]] = i
            uu[self._u[i]] = i

        return self._new_c(rr)

    def vertical_symmetry(self):
        r"""
        Return the origami `(r^{-1}, u)`.

        EXAMPLES::

            sage: o = Origami('(1,2,3)(5,6)', '(1,4)(3,5,7)')
            sage: ov = o.vertical_symmetry(); ov
            (1,3,2)(4)(5,6)(7)
            (1,4)(2)(3,5,7)(6)

        We check that it commutes with vertical symmetry of the cylinder
        diagram::

            sage: cv1 = o.cylinder_diagram().vertical_symmetry()
            sage: cv1.relabel()
            sage: cv2 = ov.cylinder_diagram()
            sage: cv2.relabel()
            sage: cv1 == cv2
            True
        """
        cdef int i
        cdef int * rr = <int *> sage_malloc(2*self._n*sizeof(int))
        cdef int * uu = rr + self._n

        for i in range(self._n):
            rr[self._r[i]] = i
            uu[i] = self._u[i]

        return self._new_c(rr)

    def horizontal_symmetry(self):
        r"""
        Return the origami `(r, u^{-1})`.

        EXAMPLES::

            sage: o = Origami('(1,2,3)(5,6)', '(1,4)(3,5,7)')
            sage: oh = o.horizontal_symmetry(); oh
            (1,2,3)(4)(5,6)(7)
            (1,4)(2)(3,7,5)(6)

        We check that it commutes with vertical symmetry of the cylinder
        diagram::

            sage: ch1 = o.cylinder_diagram().horizontal_symmetry()
            sage: ch1.relabel()
            sage: ch2 = oh.cylinder_diagram()
            sage: ch2.relabel()
            sage: ch1 == ch2
            True
        """
        cdef int i
        cdef int * rr = <int *> sage_malloc(2*self._n*sizeof(int))
        cdef int * uu = rr + self._n

        for i in range(self._n):
            rr[i] = self._r[i]
            uu[self._u[i]] = i

        return self._new_c(rr)

    cpdef mirror(self):
        r"""
        Return the origami (u,r)

        EXAMPLES::

            sage: o = Origami('(1,2,3,4)','(1,5)'); o
            (1,2,3,4)(5)
            (1,5)(2)(3)(4)
            sage: o.mirror()
            (1,5)(2)(3)(4)
            (1,2,3,4)(5)
        """
        cdef int i
        cdef int *rr = <int *>sage_malloc(2*self._n*sizeof(int))
        cdef int *uu = rr + self._n
    
        for i in range(self._n):
            rr[i] = self._u[i]
            uu[i] = self._r[i]

        return self._new_c(rr)

    cpdef horizontal_twist(self,width=1):
        r"""
        Return the origami `(r, ur^{-k})` which is obtained by the action of an
        horizontal twist of width ``k`` on this origami

        INPUT:

        - ``width`` - integer (default: 1) - the width of the twist

        EXAMPLES::

            sage: o = Origami('(1,2,3,4,5,6)','(1,7)')
            sage: o
            (1,2,3,4,5,6)(7)
            (1,7)(2)(3)(4)(5)(6)
            sage: o.horizontal_twist()
            (1,2,3,4,5,6)(7)
            (1,6,5,4,3,2,7)
            sage: o.horizontal_twist(-1)
            (1,2,3,4,5,6)(7)
            (1,2,3,4,5,6,7)
            sage: o.horizontal_twist(6) == o
            True
        """
        cdef int i,ii
        cdef int *rr = <int *>sage_malloc(2*self._n*sizeof(int))
        cdef int *uu = rr + self._n

        for i from 0 <= i < self._n:
            rr[i] = self._r[i]
            if width <= 0:
                ii = i
                for j from 0 <= j < -width:
                    ii = self._r[ii]
                uu[i] = self._u[ii]
            else:
                ii = self._r[i]
                for j from 1 <= j < width:
                    ii = self._r[ii]
                uu[ii] = self._u[i]

        return self._new_c(rr)

    cpdef vertical_twist(self,width=1):
        r"""
        Return the origami `(ru^{-k},u)` which is obtained by the action of a
        vertical twist of width `k` on this origami

        INPUT:

        - ``width`` - integer (default: 1) - the width of the twist

        EXAMPLES::

            sage: o = Origami('(1,2,3,4)','(4,5,6,7)')
            sage: o.vertical_twist()
            (1,2,3,4,7,6,5)
            (1)(2)(3)(4,5,6,7)
            sage: o.vertical_twist(-1)
            (1,2,3,4,5,6,7)
            (1)(2)(3)(4,5,6,7)
            sage: o.vertical_twist(4) == o
            True
        """
        cdef int i
        cdef int *rr = <int *>sage_malloc(2*self._n*sizeof(int))
        cdef int *uu = rr + self._n
    
        for i from 0 <= i < self._n:
            uu[i] = self._u[i]
            if width <= 0:
                ii = i
                for j from 0 <= j < -width:
                    ii = self._u[ii]
                rr[i] = self._r[ii]
            else:
                ii = self._u[i]
                for j from 1 <= j < width:
                    ii = self._u[ii]
                rr[ii] = self._r[i]

        return self._new_c(rr)

    cdef _compute_gl2z_edges(self):
        r"""
        Compute the action of `GL(2,\ZZ)`

        The generators of `GL(2,ZZ)` considered are

        I =
        0 1
        1 0

        L =
        1 1
        0 1

        EXAMPLES::

            sage: o = Origami((1,0,2),(0,2,1),as_tuple=True)
            sage: l,i = o.gl2z_edges() #indirect doctest
            sage: len(l)
            3
            sage: len(i)
            3
        """
        cdef dict l_edges = self._l_edges
        cdef dict i_edges = self._i_edges
        cdef bool VERBOSE=False
        cdef int i, n=self._n
        cdef size_t N = 2*self._n*sizeof(int)
        cdef set waiting = set([])
        cdef int * renum = <int *> sage_malloc(self._n * sizeof(int))
        cdef int *r = <int *> sage_malloc(N)
        cdef int *u = r+n
        cdef int * rr = <int *> sage_malloc(N)
        cdef int * uu = rr+n
        cdef int *rrr
        cdef int *uuu
        cdef Origami_dense_pyx o,oo,ooo

        # set o=self
        # and put o in normal form
        for i from 0 <= i < n:
            r[i] = self._r[i]
            u[i] = self._u[i]
        origami_normal_form(r,u,renum,n)
        o = self._new_c(r)
        waiting.add(o)

        # at each step o,r,u is set to the current origami
        # and rr,uu is memory allocated
        while True:
            if VERBOSE:
                print "loop from %s %s" %(str(o.r_tuple()),str(o.u_tuple()))
                print "len(l) = %d len(i) = %d" %(len(l_edges),len(i_edges))
            if o in l_edges:
                raise ValueError, "%s seen before" %str(o)

            # we compute backword the l-cusp of o
            for i from 0 <= i < n:
                rr[i] = r[i]
                uu[i] = r[u[i]]
            origami_normal_form(rr, uu, renum, n)

            ooo = o
            while origami_diff(r,rr,n):
                oo = o._new_c(rr)
                if VERBOSE:
                    print "new element in cups %s %s" %(str(oo.r_tuple()),str(oo.u_tuple()))
                waiting.add(oo)
                l_edges[oo] = ooo
                ooo = oo
                rrr = rr
                uuu = uu
                rr = <int *>sage_malloc(N)
                uu = rr+n
                for i from 0 <= i < n:
                    rr[i] = rrr[i]
                    uu[i] = rrr[uuu[i]]
                origami_normal_form(rr,uu,renum,n)

            l_edges[o] = ooo

            if VERBOSE:
                print "cups computed"
                print "len(l) = %d len(i) = %d" %(len(l_edges),len(i_edges))
                print "%d origami wait" %len(waiting)

            # then we create i-edges until we find a new guy
            # we set r,u to be available
            r = rr; u = uu
            while waiting:
                oo = waiting.pop()
                if VERBOSE:
                    print "try i-link from %s %s" %(str(oo.r_tuple()),str(oo.u_tuple()))
                rr = oo._r; uu = oo._u
                for i from 0 <= i < n:
                    r[i] = uu[i]
                    u[i] = rr[i]
                origami_normal_form(r,u,renum,n)
                if origami_diff(r,rr,n): # not symmetric under r <-> u
                    o = self._new_c(r)
                    if VERBOSE:
                        print "find new guy %s %s" %(str(o.r_tuple()),str(o.u_tuple()))
                    rr = <int *>sage_malloc(N)
                    uu = rr+n
                    i_edges[o] = oo
                    i_edges[oo] = o
                    if o in waiting: # we find a fake new guy
                        if VERBOSE:
                            print "he was there before"
                        waiting.remove(o)
                        r = rr; u = uu
                    else: # we find a real new guy
                        break
                else: # symmetric under r <-> u
                    if VERBOSE:
                        print "symmetric one"
                    i_edges[oo] = oo
                    rr = NULL
                    uu = NULL
            else:
                break

        sage_free(renum)
        sage_free(rr)

    def gl2z_edges(self):
        r"""
        Return a couple of dictionnaries ``(l_edges, i_edges)`` associated to
        the action of `GL(2,\ZZ)`

        The generators of `GL(2,ZZ)` considered are

        BEWARE: Do not modify the output dictionnary!

        .. MATH::

            L=\begin{pmatrix}1&1\\0&1\end{pmatrix}
            I=\begin{pmatrix}0&1\\1&0\end{pmatrix}
            \quad \text{and} \quad

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: l,i = o.gl2z_edges()
            sage: for oo in l: print "(%s,%s) -> (%s,%s)" %(oo.r(),oo.u(),l[oo].r(),l[oo].u())
            ((1,2,3),(2,3)) -> ((1,2,3),(2,3))
            ((2,3),(1,2)) -> ((2,3),(1,2,3))
            ((2,3),(1,2,3)) -> ((2,3),(1,2))

        TESTS::

            sage: o = Origami('(1,2,3,4)','(1,5)')
            sage: l,i = o.gl2z_edges()
            sage: ll,ii = o.gl2z_edges()
            sage: l is ll
            True
            sage: i is ii
            True
        """
        if not self._l_edges: #an empty dictionnary means that we do not have
                              #computed yet
            self._compute_gl2z_edges()

        return self._l_edges, self._i_edges

    def sl2z_edges(self):
        r"""
        Action of the matrices l,r,s

        L =
        1 1
        0 1

        R =
        1 0
        1 1

        S = L~RL =
        0 -1
        1 0

        EXAMPLES::

            sage: o = Origami('(1,2)','(2,3)')
            sage: l,r,s = o.sl2z_edges()
            sage: len(l)
            3
            sage: len(r)
            3
            sage: len(s)
            3
        """
        o = self.to_standard_form()
        L,I = self.gl2z_edges()
        return sl_orbit_from_gl_orbit(o,L,I)

    def pgl2z_edges(self):
        r"""
        Action of `PGL(2,\ZZ)`

        Projective action of the matrices

        L=
        1 1
        0 1

        and

        I=
        0 1
        1 0

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: l,i = o.pgl2z_edges()
            sage: len(l)
            3
            sage: len(i)
            3
        """
        return projectivize_edges(self.gl2z_edges())

    def psl2z_edges(self):
        r"""
        Return the action of `PSL(2,\ZZ)`

        The generators are

        L=
        1 1
        0 1

        R=
        0 1
        1 1

        and

        S=
        0 -1
        1 0

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: l,r,s = o.psl2z_edges()
            sage: len(l)
            3
            sage: len(r)
            3
            sage: len(s)
            3
        """
        return projectivize_edges(self.sl2z_edges())

    #
    # Lyapunov exponents
    #

    def lyapunov_exponents_approx(self,
            nb_iterations=0X5000, nb_experiments=4, only_mean=True,
            seed=None,
            nb_vectors=None, nb_vectors_p=None, nb_vectors_m=None,
            involution=None):
        r"""
        Approximation of Lyapunov exponents of the Kontsevich-Zorich cocycle.

        An origami defines a Teichmuller curve in the moduli space of
        translation surfaces. The Kontsevich-Zorich cocycle above this
        Teichmuller curve (for the Haar measure) has the following form

        .. MATH::

            1 = \lambda_1, \lambda_2, ... \lambda_g, -\lambda_g, ...,
            -\lambda_2, -\lambda_1 = -1

        This function return the approximations of `\lambda_2, \lambda_3, ...,
        \lambda_g)` as a list.

        INPUT:

        - ``nb_iterations`` - integer (default: 2**17) - the number of
          iterations performed in the algorithm

        - ``nb_experiments`` - integer (default: 4) - the number of experiments
          to perform

        - ``only_mean`` - boolean (default: ``True``) - if ``True``, returns the
          list of mean exponents, otherwise returns a list of lists.

        - ``nb_vectors`` - integer (default: genus-1) - the number of vectors to
          consider

        - ``involution`` - permutation or boolean - if ``True`` or an inan involution for the
          origami with
          derivative either 1.

        - ``nb_vectors_p``, ``nb_vectors_m`` - if involution is not None, then
          it will be interpreted as the number of + and - vectors to consider.

        EXAMPLES::

            sage: o = Origami('(1,2)(3,4)(5,6)','(2,3)(4,5)')
            sage: lexp = o.lyapunov_exponents_approx()
            sage: lexp # random
            [0.665250480769510, 0.332194948308155]
            sage: 0.6 < lexp[0] < 0.7 and 0.3 < lexp[1] < 0.4
            True

            sage: o = Origami('(1,2)(3,4)(5,6)(7,8)(9,10)','(2,3)(4,5)(6,7)(8,9)')
            sage: s = SymmetricGroup(10)('(1,10)(2,9)(3,8)(4,7)(5,6)')
            sage: o.lyapunov_exponents_approx(involution=s)  # random
            ([0.600372348286643, 0.199902392953331],
             [0.800242827363281, 0.399695139521823])
        """
        #TODO: use the seed to init.

        if involution is None:
            if nb_vectors is None:
                nb_vectors = self.genus()-1
            nb_vectors = int(nb_vectors)
            assert nb_vectors >= 0, "nb_vectors should be >= 0"

            return self.pyx_lyapunov_exponents_approx(nb_iterations,
                    nb_experiments, nb_vectors, only_mean)

        else:
            if involution is True:
                A = self.automorphism_group()
                assert A.order() == 2, "The automorphism group is not of order 2"
                involution = A.list()[1]
            else:
                from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
                if not isinstance(involution, PermutationGroupElement):
                    involution = PermutationGroupElement(involution)

            assert involution * self.r() * involution == self.r(), "srs is different from s"
            assert involution * self.u() * involution == self.u(), "sus is different from u"

            if nb_vectors_p is None or nb_vectors_m is None:
                from sage.groups.perm_gps.permgroup import PermutationGroup
                oo = self.quotient(PermutationGroup([involution]))
                if nb_vectors_p is None:
                    nb_vectors_p = oo.genus()-1
                if nb_vectors_m is None:
                    nb_vectors_m = self.genus()-oo.genus()

            nb_vectors_p = int(nb_vectors_p)
            assert nb_vectors_p >= 0
            nb_vectors_m = int(nb_vectors_m)
            assert nb_vectors_m >= 0

            involution = [x-1 for x in involution.tuple()]
            involution.extend(range(len(involution),self.nb_squares()))

            return self.pyx_lyapunov_exponents_approx_with_involution(
                    involution,
                    nb_iterations, nb_experiments,
                    nb_vectors_p, nb_vectors_m,
                    only_mean)

    cdef pyx_lyapunov_exponents_approx_with_involution(self, involution,
            nb_iterations, nb_experiments,
            nb_vectors_p, nb_vectors_m,
            only_mean):
        import sys
        from time import time
        cdef origami_with_involution_data * o
        cdef double * theta
        n_p = nb_vectors_p
        n_m = nb_vectors_m

        if n_p == 1: n_p=2
        if n_m == 1: n_m=2
        n = n_p + n_m

        res = [[] for _ in xrange(n)]
        theta = <double *> sage_malloc((n+1)*sizeof(double))
        s = <int *> sage_malloc((self.nb_squares()) * sizeof(int))

        for i in xrange(self.nb_squares()):
            s[i] = involution[i]

        o = new_origami_with_involution_data(
                self.nb_squares(), # degree
                n_p,               # nb_vectors_p
                n_m,               # nb_vectors_m
                self._r,           # r
                self._u,           # u
                s)                 # involution

        from sage.rings.real_mpfr import RealField
        R = RealField()
        for _ in xrange(nb_experiments):
            lyapunov_exponents_with_involution(o, nb_iterations, theta)
            for i in xrange(n):
                res[i].append(R(theta[i+1] / (2*theta[0])))

        free_origami_with_involution_data(o)
        sage_free(s)
        sage_free(theta)

        if only_mean:
            rres = []
            for i in xrange(n):
                rres.append(sum(res[i]) / nb_experiments)

            return rres[:nb_vectors_p], rres[n_p:n_p+nb_vectors_m]

        return res[:nb_vectors_p], res[n_p:n_p+nb_vectors_m]

    cdef pyx_lyapunov_exponents_approx(self,
            nb_iterations, nb_experiments,
            nb_vectors, only_mean):
        import sys
        cdef origami_data * o
        cdef double * theta
        n = max(2,nb_vectors)

        res = [[] for _ in xrange(n)]
        theta = <double *> sage_malloc((n+1)*sizeof(double))

        o = new_origami_data(
                self.nb_squares(),
                n,
                self._r,
                self._u)

        from sage.rings.real_mpfr import RealField
        R = RealField()
        for _ in xrange(nb_experiments):
            lyapunov_exponents(o, nb_iterations,theta)
            for i in xrange(n):
                res[i].append(R(theta[i+1] / (2*theta[0])))

        free_origami_data(o)
        sage_free(theta)

        if only_mean:
            rres = []
            for i in xrange(n):
                rres.append(sum(res[i]) / nb_experiments)

            return rres[:nb_vectors]

        return res[:nb_vectors]


cdef gl2z_orbits(origamis,int n, int limit):
    r"""
    Compute the action of `GL(2,\ZZ)` on the set ``origamis`` whose elements
    should be in normal form and have same number of squares ``n``. An optional
    argument ``limit`` can be used in order to stop the computation if too much
    origami are found (and then ``None`` is returned). If limit is a non
    positive number then the function does not take care of it.

    The generators of `GL(2,ZZ)` considered are

    I =
    0 1
    1 0

    L =
    1 1
    0 1
    """
    cdef list orbits = []
    cdef dict l_edges
    cdef dict i_edges
    cdef int i
    cdef size_t N = 2*n*sizeof(int)
    cdef set waiting = set([])
    cdef int * renum = <int *> sage_malloc(n * sizeof(int))
    cdef int *r = <int *> sage_malloc(N)
    cdef int *u = r+n
    cdef int *rr=NULL
    cdef int *uu=NULL
    cdef int *rrr=NULL
    cdef int *uuu=NULL
    cdef Origami_dense_pyx o,oo,ooo
    cdef bool VERBOSE=False

    while origamis:
        o = origamis.pop()
        if o._l_edges: # TODO: if data is here... use it instead of clear!
            o._l_edges.clear()
            o._i_edges.clear()
        rr = <int *> sage_malloc(N)
        uu = rr+n
        if VERBOSE: 
            print "pop origami\n r=%s\n u=%s"%(str(o.r()),str(o.u()))
            print "check:"
            print "  l = %s at %d"%(o._l_edges, id(o._l_edges))
            print "  i = %s at %d"%(o._i_edges, id(o._i_edges))
        l_edges = o._l_edges
        i_edges = o._i_edges

        for i from 0 <= i < n:
            r[i] = o._r[i]
            u[i] = o._u[i]
        waiting.add(o)

        # at each step o,r,u is set to the current origami
        # rr,uu is memory allocated
        # r,u   pointed by o
        while True:
            if VERBOSE: 
                print "start new cusp..."
            if o in l_edges:
                raise ValueError("%s seen before" %str(o))
    
            # we compute backward the l-orbit of o
            for i from 0 <= i < n:
                rr[i] = r[i]
                uu[i] = r[u[i]]
            origami_normal_form(rr, uu, renum, n)
    
            ooo = o
            while origami_diff(r,rr,n):
                oo = o._new_c(rr)
                if VERBOSE: 
                    print " new elt in cusp"
                    print " r = %s"%oo.r()
                    print " u = %s"%oo.u()
                    print " go"
                if oo in origamis:
                    if VERBOSE:
                        print " remove origami in the set"
                    origamis.remove(oo)

                waiting.add(oo)
                l_edges[oo] = ooo
                ooo = oo
                rrr = rr
                uuu = uu
                rr = <int *>sage_malloc(N)
                uu = rr+n
                for i from 0 <= i < n:
                    rr[i] = rrr[i]
                    uu[i] = rrr[uuu[i]]
                origami_normal_form(rr,uu,renum,n)
    
            # at this point rr is memory allocated

            l_edges[o] = ooo
            if limit > 0 and len(l_edges) > limit: # test the size
                if VERBOSE: 
                    print "oversize"
                l_edges.clear()
                i_edges.clear()
                waiting.clear()
                r = rr; u = uu
                break
   
            # then we create i-edges until we find a new guy
            # we set r,u to be available
            if VERBOSE: print "end of cusp, apply symmetry"
            r = rr; u = uu
            while waiting:
                if VERBOSE:
                    print " new try..."
                oo = waiting.pop()
                rr = oo._r; uu = oo._u
                for i from 0 <= i < n:
                    r[i] = uu[i]
                    u[i] = rr[i]
                origami_normal_form(r,u,renum,n)
                if origami_diff(r,rr,n): # not symmetric under r <-> u
                    o = o._new_c(r)
                    rr = <int *>sage_malloc(N)
                    uu = rr+n
                    i_edges[o] = oo
                    i_edges[oo] = o
                    if o in waiting: # we find a fake new guy
                        if VERBOSE: print " ...was already there"
                        waiting.remove(o)
                        r = rr; u = uu
                    else: # we find a real new guy
                        if VERBOSE:
                            print "go elsewhere"
                        if o in origamis:
                            origamis.remove(o)
                        break
                else: # symmetric under r <-> u
                    if VERBOSE:
                        print " ...symmetric guy"
                    i_edges[oo] = oo
                    rr = NULL
                    uu = NULL
            else:
                break

        if l_edges: # append if we do not quit because of oversize        
            orbits.append((l_edges,i_edges))
            if VERBOSE: 
                    print "new orbit of size %d"%(len(l_edges))
                    print "check: waiting=",waiting

    sage_free(renum)
    sage_free(rr)

    return orbits

cpdef sl2z_orbits(origamis,int n,int limit):
    r"""
    Action of the matrices l,r,s

    L =
    1 1
    0 1

    R =
    1 0
    1 1

    S = L~RL =
    0 -1
    1 0
    """
    slorbits = []
    glorbits = gl2z_orbits(origamis,n,limit)

    for L,I in glorbits:
        o = L.iterkeys().next()
        l,r,s = sl_orbit_from_gl_orbit(o,L,I)
        slorbits.append((l,r,s))
        if len(l) != len(L):
            l,r,s = sl_orbit_from_gl_orbit(o.mirror().relabel(),L,I)
            slorbits.append((l,r,s))

    return slorbits


cdef class PillowcaseCover_dense_pyx:
    def __cinit__(self):
        r"""
        TESTS::

            sage: o = Origami('(1,2)','(1,3)')
            sage: loads(dumps(o)) == o
            True
        """
        self._n = 0
        self._g = NULL
        self._l_edges = None
        self._i_edges = None

    def __init__(self, tuple g0, tuple g1, tuple g2, tuple g3):
        r"""
        TESTS::

            sage: o = Origami([2,1,3],[1,3,2])
            sage: o == loads(dumps(o))
            True
        """
        cdef int i

        self._n = len(g0)
        self._g = <int *> sage_malloc(4*self._n*sizeof(int))

        if self._g == NULL:
            raise MemoryError, "not able to allocate"
       
        for i from 0 <= i < self._n:
            self._g[i]     = g0[i]
            self._g[i+1*self._n] = g1[i]
            self._g[i+2*self._n] = g2[i]
            self._g[i+3*self._n] = g3[i]

        self._l_edges = {}
        self._i_edges = {}

    def __dealloc__(self):
        if self._g != NULL: sage_free(self._g)

    cdef PillowcaseCover_dense_pyx _new_c(self, int * g):
        r"""
        Return an origami with given permutations.

        Beware that we assume that the created origami is in the same orbit
        under the action of PGL(2,Z)
        """
        cdef PillowcaseCover_dense_pyx other = PY_NEW_SAME_TYPE(self)

        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__

        other._n = self._n
        other._g = g
        
        other._l_edges = self._l_edges
        other._i_edges = self._i_edges

        return other
    
    def __copy__(self):
        r"""
        Return a copy of the origami

        EXAMPLES::

            sage: o = Origami('(1,2)','(1,3)')
            sage: oo = copy(o)
            sage: o == oo
            True
            sage: o is oo
            False
        """
        cdef int * g = <int *> sage_malloc(4 * self._n * sizeof(int))

        memcpy(g,self._g, 4 * self._n * sizeof(int))

        return self._new_c(g)

    def g_tuple(self, int i):
        r"""
        Return a tuple.

        EXAMPLES::

            sage: p = PillowcaseCover([2,1,3,4],[3,2,1,4],[4,2,3,1])
            sage: p.g_tuple(0)
            (1, 0, 2, 3)
            sage: p.g_tuple(1)
            (2, 1, 0, 3)
            sage: p.g_tuple(2)
            (3, 1, 2, 0)
            sage: p.g_tuple(3)
            (3, 0, 1, 2)
        """
        if i < 0 or i > 3:
            raise ValueError
        return array_to_tuple(self._g+i*self._n,self._n)

    def degree(self):
        r"""
        The degree of the covering.
        
        EXAMPLES::
        
            sage: p = PillowcaseCover([2,1,3,4],[3,2,1,4],[4,2,3,1])
            sage: p.degree()
            4
        """
        return self._n
