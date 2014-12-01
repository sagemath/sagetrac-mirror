"""
A Cython implementation of elements of path algebras.

AUTHORS:

- Simon King (2014-12-04)

"""

#*****************************************************************************
#     Copyright (C) 2014 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "algebra_elements.pxi"
from sage.misc.cachefunc import cached_method
from sage.misc.misc import repr_lincomb

## TODO
# Move representation from pxi to here and repr/latex_monomial to the parent
### (can only be done after making PathAlgebraElement the default element class)

cdef class PathAlgebraElement(RingElement):
    def __cinit__(self, *args, **kwds):
        self.data = NULL
    def __dealloc__(self):
        homog_poly_free(self.data)
    def __init__(self, S, data, order="negdegrevlex"):
        self._hash = -1
        if 'order' in data:
            order = data.pop('order')
        if order=="negdegrevlex":
            self.cmp_terms = negdegrevlex
        elif order=="degrevlex":
            self.cmp_terms = degrevlex
        elif order=="revlex":
            self.cmp_terms = revlex
        else:
            raise ValueError("Unknown term order '{}'".format(order))
        cdef QuiverPath tmp = None
        RingElement.__init__(self, S)
        cdef dict homog = {}
        cdef list L
        for tmp, c in data.iteritems():
            homog.setdefault((tmp.initial_vertex(),tmp.terminal_vertex()),[]).append((tmp,c))
        cdef path_homog_poly_t *HP
        for (s,e),L in sorted(homog.iteritems(), reverse=True):
            HP = homog_poly_init_list(s,e,L,self.cmp_terms,0,-1)
            HP.nxt = self.data
            self.data = HP

    def __reduce__(self):
        return path_algebra_element_unpickle, (self._parent, homog_poly_pickle(self.data))

    cdef list _sorted_items_for_printing(self):
        """
        Return list of pairs ``(M,c)``, where ``c`` is a coefficient and ``M``
        will be passed to ``self.parent()._repr_monomial`` resp. to
        ``self.parent()._latex_monomial``, providing the indices of the
        algebra generators occurring in the monomial.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: X         # indirect doctest
            5*e_0 + 3*e_2 + a + 2*b + 3*c
            sage: latex(X)  # indirect doctest
            5e_0 + 3e_2 + a + 2b + 3c

        """
        cdef path_homog_poly_t *H = self.data
        cdef list L, L_total
        cdef size_t i
        cdef path_term_t * T
        L_total = []
        cdef list vertices = self._parent.quiver().vertices()
        cdef unsigned int offset = len(vertices)
        while H != NULL:
            L = []  # data for a single component (given by start- and endpoints)
            T = H.poly.lead
            while T!=NULL:
                if T.mon.path.length:
                    L.append(([offset+biseq_getitem(T.mon.path,i) for i in range(T.mon.path.length)],
                              <object>(T.coef)))
                else:
                    L.append(([vertices.index(H.start)], <object>(T.coef)))
                T = T.nxt
            if len(L) != H.poly.nterms:
                print "Term count of polynomial is wrong, got",len(L), "expected", H.poly.nterms
            L_total.extend(L)
            H = H.nxt
        return L_total

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: X         # indirect doctest
            5*e_0 + 3*e_2 + a + 2*b + 3*c

        """
        return repr_lincomb(self._sorted_items_for_printing(), strip_one=True,
                            scalar_mult=self.parent()._print_options['scalar_mult'],
                            repr_monomial = self._parent._repr_monomial
                            )

    def _latex_(self):
        """
        Latex string representation.

        EXAMPLES::

            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('a+2*b+3*c+5*e_0+3*e_2', A.gens_dict())
            sage: latex(X)  # indirect doctest
            5e_0 + 3e_2 + a + 2b + 3c

        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult       = self.parent()._print_options['scalar_mult'],
                            latex_scalar_mult = self.parent()._print_options['latex_scalar_mult'],
                            latex_scalar_mult = None,
                            repr_monomial = self._parent._latex_monomial,
                            is_latex=True, strip_one = True)

    # Basic properties

    def __nonzero__(self):
        return self.data != NULL

    def __len__(self):
        cdef size_t l = 0
        cdef path_homog_poly_t *H = self.data
        while H != NULL:
            l += H.poly.nterms
            H = H.nxt
        return l

    cpdef ssize_t degree(self) except -2:
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef ssize_t deg = -1
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                if deg == -1:
                    deg = T.mon.path.length
                elif deg != T.mon.path.length:
                    raise ValueError("Element is not homogeneous.")
                T = T.nxt
            H = H.nxt
        return deg

    def is_homogeneous(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef ssize_t deg = -1
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                if deg == -1:
                    deg = T.mon.path.length
                elif deg != T.mon.path.length:
                    return False
                T = T.nxt
            H = H.nxt
        return True

    cpdef dict monomial_coefficients(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        cdef dict D = {}
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                D[tmp] = <object>T.coef
                T = T.nxt
            H = H.nxt
        return D

    cpdef list coefficients(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef list L = []
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                L.append(<object>T.coef)
                T = T.nxt
            H = H.nxt
        return L

    cpdef list monomials(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_homog_poly_t *out
        cdef path_term_t *T
        cdef object one = self.base_ring().one()
        cdef list L = []
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                out = homog_poly_create(H.start, H.end)
                out.poly.lead = term_create_keep_mon(one, mon_copy(T.mon))
                out.poly.lead.nxt = NULL
                out.poly.nterms = 1
                L.append(self._new_(out))
                T = T.nxt
            H = H.nxt
        return L

    cpdef list terms(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_homog_poly_t *out
        cdef path_term_t *T
        cdef object one = self.base_ring().one()
        cdef list L = []
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                out = homog_poly_create(H.start, H.end)
                out.poly.lead = term_copy(T)
                out.poly.lead.nxt = NULL
                out.poly.nterms = 1
                L.append(self._new_(out))
                T = T.nxt
            H = H.nxt
        return L

    cpdef list support(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        cdef list L = []
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                L.append(tmp)
                T = T.nxt
            H = H.nxt
        return L

    def support_of_term(self):
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        if self.data != NULL and self.data.nxt == NULL:
            if self.data.poly.lead != NULL:
                tmp = sample._new_(self.data.start, self.data.end)
                biseq_init_copy(tmp._path, self.data.poly.lead.mon.path)
                return tmp
        raise ValueError("{} is not a single term".format(self))

    def __iter__(self):
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                yield (tmp, <object>T.coef)
                T = T.nxt
            H = H.nxt

    cdef inline PathAlgebraElement _new_(self, path_homog_poly_t *h):
        cdef PathAlgebraElement out = PY_NEW(type(self))
        out._parent = self._parent
        out.cmp_terms = self.cmp_terms
        out.data = h
        out._hash = -1
        return out

    def __copy__(self):
        return self._new_(homog_poly_copy(self.data))

    def __getitem__(self, k):
        cdef path_homog_poly_t *H
        cdef path_term_t *T
        cdef path_mon_t *kM
        cdef PathAlgebraElement out
        cdef QuiverPath K
        if isinstance(k, tuple):
            H = homog_poly_get_predecessor_of_component(self.data,k[0],k[1])
            if H == NULL:
                if self.data.start == k[0] and self.data.end == k[1]:
                    out = self._new_(homog_poly_create(self.data.start, self.data.end))
                    out.data.nxt = NULL
                    poly_icopy(out.data.poly, self.data.poly)
                else:
                    return self._new_(NULL)
            else:
                if H.nxt == NULL or H.nxt.start != k[0] or H.nxt.end != k[1]:
                    return self._new_(NULL)
                out = self._new_(homog_poly_create(H.nxt.start, H.nxt.end))
                out.data.nxt = NULL
                poly_icopy(out.data.poly, H.nxt.poly)
            return out
        elif isinstance(k, QuiverPath):
            if self.data == NULL:
                return self.base_ring().zero()
            K = k
            H = homog_poly_get_predecessor_of_component(self.data, K._start, K._end)
            if H == NULL:
                if self.data.start != K._start or self.data.end != K._end:
                    return self.base_ring().zero()
                H = self.data
            else:
                H = H.nxt
            if H == NULL:
                return self.base_ring().zero()
            # Now, H points to the component that belongs to K
            kM = mon_create_keep(K._path, 0, -1)
            T = H.poly.lead
            while T != NULL:
                if self.cmp_terms(T.mon, kM) == 0:
                    return <object>T.coef
                T = T.nxt
        return self.base_ring().zero()

    cpdef object coefficient(self, QuiverPath P):
        if self.data == NULL:
            return self.base_ring().zero()
        H = homog_poly_get_predecessor_of_component(self.data, P._start, P._end)
        if H == NULL:
            if self.data.start != P._start or self.data.end != P._end:
                return self.base_ring().zero()
            H = self.data
        else:
            H = H.nxt
        if H == NULL:
            return self.base_ring().zero()
        # Now, H points to the component that belongs to K
        pM = mon_create_keep(P._path, 0, -1)
        T = H.poly.lead
        while T != NULL:
            if self.cmp_terms(T.mon, pM) == 0:
                return <object>T.coef
            T = T.nxt
        return self.base_ring().zero()        

    def filtration_by_vertices(self):
        cdef path_homog_poly_t * H = self.data
        cdef PathAlgebraElement out
        cdef list C = []
        while H != NULL:
            out = self._new_(homog_poly_create(H.start, H.end))
            out.data.nxt = NULL
            poly_icopy(out.data.poly, H.poly)
            C.append((out, H.start, H.end))
            H = H.nxt
        return C

    ####
    ## Arithmetics
    # Hash and Comparison
    def __hash__(self):
        if self._hash==-1:
            self._hash = hash(frozenset(self.monomial_coefficients().items()))
        return self._hash

    def __cmp__(left, right):
        """
        TESTS::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(GF(3))
            sage: x = PathAlgebraElement(A, sage_eval('b*e*b*e+4*b*e+e_0', A.gens_dict()))
            sage: y = PathAlgebraElement(A, sage_eval('a*c+d*c*b*f', A.gens_dict()))

        BUG::

            sage: x*(y+x) == x*y+x*x
            False

        """
        return (<Element>left)._cmp(right)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        cdef PathAlgebraElement other = right
        cdef PathAlgebraElement self = left
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2 = other.data
        cdef int c
        while H1 != NULL and H2 != NULL:
            c = cmp(H1.start, H2.start)
            if c != 0:
                return c
            c = cmp(H1.end, H2.end)
            if c != 0:
                return c
            c = poly_cmp(H1.poly, H2.poly, self.cmp_terms)
            if c != 0:
                return c
            H1 = H1.nxt
            H2 = H2.nxt
        if H1 == NULL:
            if H2 == NULL:
                return 0
            return -1
        return 1

    # negation
    cpdef ModuleElement _neg_(self):
        return self._new_(homog_poly_neg(self.data))

    # addition
    cpdef ModuleElement _add_(self, ModuleElement other):
        """
        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(GF(3))
            sage: x = PathAlgebraElement(A, sage_eval('b*e*b*e+4*b*e+e_0', A.gens_dict()))
            sage: y = PathAlgebraElement(A, sage_eval('a*c', A.gens_dict()))
            sage: x+y

        """
        cdef PathAlgebraElement right = other
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2 = right.data
        cdef path_poly_t *P
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *tmp
        cdef int c
        while True:
            if H1 == NULL:
                if out == NULL:
                    if H2 == NULL:
                        return self._new_(NULL)                    
                    return self._new_(homog_poly_copy(H2))
                else:
                    if H2 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H2)
                    return self._new_(out)
            elif H2 == NULL:
                if out == NULL:
                    if H1 == NULL:
                        return self._new_(NULL)                    
                    return self._new_(homog_poly_copy(H1))
                else:
                    if H1 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H1)
                    return self._new_(out)
            else:
                c = cmp((H1.start, H1.end), (H2.start, H2.end))
                if c == 1:
                    if out == NULL:
                        out = homog_poly_create(H2.start, H2.end)
                        poly_icopy(out.poly, H2.poly)
                        tmp = out
                    else:
                        tmp.nxt = homog_poly_create(H2.start, H2.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H2.poly)
                    H2 = H2.nxt
                elif c == -1:
                    if out == NULL:
                        out = homog_poly_create(H1.start, H1.end)
                        poly_icopy(out.poly, H1.poly)
                        tmp = out
                    else:
                        tmp.nxt = homog_poly_create(H1.start, H1.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H1.poly)
                    H1 = H1.nxt
                else:
                    # start- and endpoints match
                    P = poly_add(H1.poly, H2.poly, self.cmp_terms)
                    if P.lead != NULL:
                        if out == NULL:
                            out = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = out
                        else:
                            tmp.nxt = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = tmp.nxt
                    else:
                        poly_free(P)
                    H1 = H1.nxt
                    H2 = H2.nxt

    cpdef ModuleElement _sub_(self, ModuleElement other):
        """
        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(GF(3))
            sage: x = PathAlgebraElement(A, sage_eval('b*e*b*e+4*b*e+e_0', A.gens_dict()))
            sage: y = PathAlgebraElement(A, sage_eval('a*c', A.gens_dict()))
            sage: x+y

        """
        cdef PathAlgebraElement right = other
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2 = right.data
        cdef path_poly_t *P
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *tmp
        cdef int c
        while True:
            if H1 == NULL:
                if out == NULL:
                    if H2 == NULL:
                        return self._new_(NULL)                    
                    return self._new_(homog_poly_copy(H2))
                else:
                    if H2 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H2)
                    return self._new_(out)
            elif H2 == NULL:
                if out == NULL:
                    if H1 == NULL:
                        return self._new_(NULL)                    
                    return self._new_(homog_poly_copy(H1))
                else:
                    if H1 != NULL:
                        # If out is not NULL then tmp isn't either
                        tmp.nxt = homog_poly_copy(H1)
                    return self._new_(out)
            else:
                c = cmp((H1.start, H1.end), (H2.start, H2.end))
                if c == 1:
                    if out == NULL:
                        out = homog_poly_create(H2.start, H2.end)
                        poly_icopy(out.poly, H2.poly)
                        tmp = out
                    else:
                        tmp.nxt = homog_poly_create(H2.start, H2.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H2.poly)
                    H2 = H2.nxt
                elif c == -1:
                    if out == NULL:
                        out = homog_poly_create(H1.start, H1.end)
                        poly_icopy(out.poly, H1.poly)
                        tmp = out
                    else:
                        tmp.nxt = homog_poly_create(H1.start, H1.end)
                        tmp = tmp.nxt
                        poly_icopy(tmp.poly, H1.poly)
                    H1 = H1.nxt
                else:
                    # start- and endpoints match
                    P = poly_sub(H1.poly, H2.poly, self.cmp_terms)
                    if P.lead != NULL:
                        if out == NULL:
                            out = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = out
                        else:
                            tmp.nxt = homog_poly_init_poly(H1.start, H1.end, P)
                            tmp = tmp.nxt
                    else:
                        poly_free(P)
                    H1 = H1.nxt
                    H2 = H2.nxt


#    # TODO: REMOVE
#    def iadd_mul(self, coef, QuiverPath L, PathAlgebraElement P, QuiverPath R):
#        """
#        from sage.quivers.algebra_elements import PathAlgebraElement
#        A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
#        DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().inject_variables()
#        x = PathAlgebraElement(A, sage_eval('a*c+7*c*b+c*b*e*b*e*b+3*c*b*e*b', A.gens_dict()))
#        y = PathAlgebraElement(A, sage_eval('a*c+d*c*b*f', A.gens_dict()))
#        y.iadd_mul(2, d, x, f)
#
#        sage: x = PathAlgebraElement(A, sage_eval('b*e*b*e+4*b*e+e_0', A.gens_dict()))
#        sage: copy(x).iadd_mul(1,b*e,x,e_0)  # correct, but memory leak
#        e_0 + 5*b*e + 5*b*e*b*e + b*e*b*e*b*e
#
#
#        """
#        coef = self._parent.base_ring()(coef)
#        cdef path_homog_poly_t *prev = NULL
#        cdef path_homog_poly_t *H2 = P.data
#        if not coef:
#            return
#        while H2 != NULL:
#            if H2.start == L._end and H2.end == R._start:
#                break
#            H2 = H2.nxt
#        if H2 == NULL:
#            return self
#        if self.data == NULL:
#            self.data = homog_poly_create(L._start, R._end)
#            poly_iadd_lrmul(self.data.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
#            return self
#        prev = homog_poly_get_predecessor_of_component(self.data, L._start, R._end)
#        cdef path_homog_poly_t *H1
#        if prev == NULL:
#            H1  = self.data
#            if H1!=NULL and H1.start == L._start and H1.end == R._end:
#                # We need to override H1
#                poly_iadd_lrmul(H1.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
#                # Do we need to delete H1?
#                if H1.poly.lead == NULL:
#                    self.data = H1.nxt
#                    sage_free(H1.poly)
#                    sage_free(H1)
#                return self
#            self.data = homog_poly_create(L._start, R._end)
#            poly_iadd_lrmul(self.data.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
#            if self.data.poly.lead == NULL:
#                sage_free(self.data.poly)
#                sage_free(self.data)
#                self.data = H1
#            else:
#                self.data.nxt = H1
#        else:
#            H1 = prev.nxt
#            if H1.start == L._start and H1.end == R._end:
#                # We need to override H1
#                poly_iadd_lrmul(H1.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
#                # Do we need to delete H1?
#                if H1.poly.lead == NULL:
#                    prev.nxt = H1.nxt
#                    sage_free(H1.poly)
#                    sage_free(H1)
#                return self
#            prev.nxt = homog_poly_create(L._start, R._end)
#            poly_iadd_lrmul(prev.nxt.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
#            if prev.nxt.poly.lead == NULL:
#                sage_free(prev.nxt.poly)
#                sage_free(prev.nxt)
#                prev.nxt = H1
#            else:
#                prev.nxt.nxt = H1
#        return self

## (scalar) multiplication

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: X = sage_eval('3*a+3*b+3*c+3*e_0+3*e_2', A.gens_dict())
            sage: x = PathAlgebraElement(A, X)
            sage: x*2

        ::

            sage: X = sage_eval('3*a+4*b+3*c+4*e_0+3*e_2', A.gens_dict())
            sage: x = PathAlgebraElement(A, X)
            sage: z = (x^3)[0,2]
            sage: z*5   # crash

        """
        cdef path_homog_poly_t * out = homog_poly_scale(self.data, right)
        if out.poly.nterms == 0:
            homog_poly_free(out)
            return self._new_(NULL)
        return self._new_(out)

    cpdef ModuleElement _rmul_(self, RingElement right):
        cdef path_homog_poly_t * out = homog_poly_scale(self.data, right)
        if out.poly.nterms == 0:
            homog_poly_free(out)
            return self._new_(NULL)
        return self._new_(out)

    def __div__(self, x):
        """
        Division by coefficients

        """
        if isinstance(self, PathAlgebraElement):
            return (<PathAlgebraElement>self)._new_(homog_poly_scale((<PathAlgebraElement>self).data, ~(self.base_ring()( x ))))
        return (<PathAlgebraElement>x)._new_(homog_poly_scale((<PathAlgebraElement>x).data, ~(x.base_ring()( self ))))

    cpdef RingElement _mul_(self, RingElement  other):
        """
        Former leaks or bugs...

        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = PathAlgebraElement(A, b*e*b*e+4*b*e+e_0)
            sage: y = PathAlgebraElement(A, a*c+5*f*e)
            sage: x*y
            sage: y*x
            sage: y*y
            sage: x*x

        ::

            sage: x = PathAlgebraElement(A, sage_eval('b*e*b*e+4*b*e+e_0', A.gens_dict()))
            sage: y = PathAlgebraElement(A, sage_eval('a*c+d*c*b*f', A.gens_dict()))
            sage: x*(y+x) == x*y+x*x
            True

        ::

            sage: X = sage_eval('a+2*b+3*c+e_0+e_2', A.gens_dict())
            sage: x = PathAlgebraElement(A, X)
            sage: x^3 == PathAlgebraElement(A, X^3)
            True

        """
        cdef PathAlgebraElement right = other
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2
        cdef path_term_t *T2
        cdef path_poly_t *P
        cdef path_homog_poly_t *out_orig = NULL
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *nxt
        cdef int c
        while H1 != NULL:
            H2 = right.data
            while H2 != NULL:
                if H2.start == H1.end:
                    out = homog_poly_get_predecessor_of_component(out_orig, H1.start, H2.end)
                    if out == NULL:
                        if out_orig == NULL:
                            out_orig = homog_poly_create(H1.start, H2.end)
                        else:
                            if out_orig.start != H1.start or out_orig.end != H2.end:
                                nxt = out_orig
                                out_orig = homog_poly_create(H1.start, H2.end)
                                out_orig.nxt = nxt
                    else:
                        if out.nxt==NULL or out.nxt.start != H1.start or out.nxt.end != H2.end:
                            nxt = out.nxt
                            out.nxt = homog_poly_create(H1.start, H2.end)
                            out.nxt.nxt = nxt
                    T2 = H2.poly.lead
                    # now, either out==NULL, and we need to put the product
                    # into out_orig; or out!=NULL, and we need to put the
                    # product into out.nxt
                    if out == NULL:
                        while T2 != NULL:
                            #assert poly_is_sane(out_orig.poly)
                            #print "out==0, T2!=0", H1.start, H1.end, H2.start,H2.end,
                            poly_iadd_lmul(out_orig.poly, <object>T2.coef, H1.poly, T2.mon.path,
                                           self.cmp_terms, 0, 0, -1)
                            T2 = T2.nxt
                    else:
                        while T2 != NULL:
                            #assert poly_is_sane(out.nxt.poly)
                            #print "out!=0, T2!=0", H1.start, H1.end, H2.start, H2.end, 
                            poly_iadd_lmul(out.nxt.poly, <object>T2.coef, H1.poly, T2.mon.path,
                                           self.cmp_terms, 0, 0, -1)
                            T2 = T2.nxt                                
                H2 = H2.nxt
            H1 = H1.nxt
        while out_orig != NULL and out_orig.poly.lead == NULL:
            tmp = out_orig.nxt
            sage_free(out_orig.poly)
            sage_free(out_orig)
            out_orig = tmp
        if out_orig == NULL:
            return self._new_(NULL)
        tmp = out_orig
        while tmp.nxt != NULL:
            if tmp.nxt.poly.lead == NULL:
                nxt = tmp.nxt.nxt
                sage_free(tmp.nxt.poly)
                sage_free(tmp.nxt)
                tmp.nxt = nxt
            else:
                tmp = tmp.nxt
        return self._new_(out_orig)

cpdef PathAlgebraElement path_algebra_element_unpickle(P, list data):
    cdef PathAlgebraElement out = PY_NEW(P.element_class)
    out._parent = P
    order = P.order_string()
    if order=="negdegrevlex":
        out.cmp_terms = negdegrevlex
    elif order=="degrevlex":
        out.cmp_terms = degrevlex
    elif order=="revlex":
        out.cmp_terms = revlex
    else:
        raise ValueError("Unknown term order '{}'".format(order))
    out.data = homog_poly_unpickle(data)
    out._hash = -1
    return out
