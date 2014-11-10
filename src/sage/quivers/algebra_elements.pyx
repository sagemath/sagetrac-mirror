include "algebra_elements.pxi"

cdef class PathAlgebraElement(RingElement):
    def __cinit__(self, *args, **kwds):
        self.data = NULL
    def __dealloc__(self):
        homog_poly_free(self.data)
    def __init__(self, S, data, order="negdegrevlex"):
        if order=="negdegrevlex":
            self.cmp_terms = negdegrevlex
        elif order=="degrevlex":
            self.cmp_terms = degrevlex
        elif order=="revlex":
            self.cmp_terms = revlex
        else:
            raise ValueError("Unknown term order '{}'".format(order))
        self._monomial_coefficients = None
        cdef QuiverPath tmp = None
        RingElement.__init__(self, S)
        cdef dict homog = {}
        cdef list L
        for tmp, c in data.monomial_coefficients().iteritems():
            homog.setdefault((tmp.initial_vertex(),tmp.terminal_vertex()),[]).append((tmp,c))
        cdef path_homog_poly_t *HP
        for (s,e),L in sorted(homog.iteritems()):
            HP = homog_poly_init_list(s,e,L,self.cmp_terms,0,-1)
            HP.nxt = self.data
            self.data = HP
    # String representation
    def _repr_monomial(self, m):
        # m is [list, pos, mid], where the list gives the nb of arrows, pos
        # gives the component in the module, and mid gives the length of the
        # left factor in a two-sided module.
        cdef list arrows = self._parent._quiver.edge_labels()
        cdef list L
        if m[0]:
            L = [arrows[n] for n in (m[0][:m[2]] if m[2]!=-1 else m[0])]
            if m[2]!=-1:
                L.append('I_{}'.format(m[1]))
                L.extend([arrows[n] for n in m[0][m[2]:]])
        elif m[2]!=-1:
            L = ['e_{}I_{}'.format(m[3],m[2])]
        else:
            L = ['e_{}'.format(m[3])]
        return '*'.join(L)

    def __repr__(self):
        return repr_homog_poly(self.data, self._repr_monomial)

    def monomial_coefficients(self):
        if self._monomial_coefficients is not None:
            return self._monomial_coefficients
        cdef path_homog_poly_t *H = self.data
        cdef path_term_t *T
        cdef QuiverPath sample = self._parent.semigroup().gen(0)
        cdef QuiverPath tmp
        self._monomial_coefficients = {}
        while H!=NULL:
            T = H.poly.lead
            while T!=NULL:
                tmp = sample._new_(H.start, H.end)
                biseq_init_copy(tmp._path, T.mon.path)
                self._monomial_coefficients[tmp] = <object>T.coef
                T = T.nxt
            H = H.nxt
        return self._monomial_coefficients

    def __iter__(self):
        if self._monomial_coefficients is not None:
            return self._monomial_coefficients.iteritems()
        return self.monomial_coefficients().iteritems()

    cdef inline PathAlgebraElement _new_(self, path_homog_poly_t *h):
        cdef PathAlgebraElement out = PY_NEW(type(self))
        out._parent = self._parent
        out.cmp_terms = self.cmp_terms
        out.data = h
        return out

    def __copy__(self):
        return self._new_(homog_poly_copy(self.data))

    cpdef ModuleElement _add_(self, ModuleElement other):
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

    def iadd_mul(self, coef, QuiverPath L, PathAlgebraElement P, QuiverPath R):
        """
        from sage.quivers.algebra_elements import PathAlgebraElement
        A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
        DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().inject_variables()
        x = PathAlgebraElement(A, sage_eval('a*c+7*c*b+c*b*e*b*e*b+3*c*b*e*b', A.gens_dict()))
        y = PathAlgebraElement(A, sage_eval('a*c+d*c*b*f', A.gens_dict()))
        y.iadd_mul(2, d, x, f)

        """
        coef = self._parent.base_ring()(coef)
        cdef path_homog_poly_t *prev = NULL
        cdef path_homog_poly_t *H2 = P.data
        if not coef:
            return
        while H2 != NULL:
            if H2.start == L._end and H2.end == R._start:
                break
            H2 = H2.nxt
        if H2 == NULL:
            return self
        if self.data == NULL:
            # self is completely null
            self.data = homog_poly_create(L._start, R._end)
            poly_iadd_lrmul(self.data.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
            return self
        prev = homog_poly_get_predecessor_of_component(self.data, L._start, R._end)
        cdef path_homog_poly_t *H1
        if prev == NULL:
            H1  = self.data
            if H1.start == L._start and H1.end == R._end:
                # We need to override H1
                poly_iadd_lrmul(H1.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
                # Do we need to delete H1?
                if H1.poly.lead == NULL:
                    self.data = H1.nxt
                    sage_free(H1.poly)
                    sage_free(H1)
                return self
            self.data = homog_poly_create(L._start, R._end)
            poly_iadd_lrmul(self.data.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
            if self.data.poly.lead == NULL:
                sage_free(self.data.poly)
                sage_free(self.data)
                self.data = H1
            else:
                self.data.nxt = H1
        else:
            H1 = prev.nxt
            if H1.start == L._start and H1.end == R._end:
                # We need to override H1
                poly_iadd_lrmul(H1.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
                # Do we need to delete H1?
                if H1.poly.lead == NULL:
                    prev.nxt = H1.nxt
                    sage_free(H1.poly)
                    sage_free(H1)
                return self
            prev.nxt = homog_poly_create(L._start, R._end)
            poly_iadd_lrmul(prev.nxt.poly, coef, L._path, H2.poly, R._path, self.cmp_terms, 0, 0, -1)
            if prev.nxt.poly.lead == NULL:
                sage_free(prev.nxt.poly)
                sage_free(prev.nxt)
                prev.nxt = H1
            else:
                prev.nxt.nxt = H1
        return self

    cpdef RingElement _mul_(self, RingElement  other):
        cdef PathAlgebraElement right = other
        cdef path_homog_poly_t *H1 = self.data
        cdef path_homog_poly_t *H2
        cdef path_poly_t *P
        cdef path_homog_poly_t *out = NULL
        cdef path_homog_poly_t *prev = NULL
        cdef int c
        while H1 != NULL:
            H2 = right.data
            
            H1 = H1.nxt
