include "algebra_elements.pxi"

## TODO
# remove inplace changes
# Move representation from pxi to here and repr/latex_monomial to the parent
# neg
# sub
# getitem/coefficient
# support?
# monomials
# terms
# coefficients
# _acted_upon_/_lmul_,_rmul_
# __div__


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
        for (s,e),L in sorted(homog.iteritems(), reverse=True):
            HP = homog_poly_init_list(s,e,L,self.cmp_terms,0,-1)
            HP.nxt = self.data
            self.data = HP

    def __nonzero__(self):
        return self.data != NULL

    def __len__(self):
        cdef size_t l = 0
        cdef path_homog_poly_t *H = self.data
        while H != NULL:
            l += H.poly.nterms
            H = H.nxt
        return l

    # String representation
    # TODO: Move the following two methods to the parent.
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

    def _latex_monomial(self, m):
        # m is [list, pos, mid], where the list gives the nb of arrows, pos
        # gives the component in the module, and mid gives the length of the
        # left factor in a two-sided module.
        cdef list arrows = self._parent._quiver.edge_labels()
        cdef list L
        if m[0]:
            L = [arrows[n] for n in (m[0][:m[2]] if m[2]!=-1 else m[0])]
            if m[2]!=-1:
                L.append('I_{%s}'%(m[1]))
                L.extend([arrows[n] for n in m[0][m[2]:]])
        elif m[2]!=-1:
            L = ['e_{%s}\\cdot I_{%s}'%(m[3],m[2])]
        else:
            L = ['e_{%s}'%(m[3])]
        return '\\cdot '.join(L)

    cdef list _sorted_items_for_printing(self):
        """
        Return list of pairs (M,c), where c is a coefficient and M is data
        that can be understood by self.parent()._repr_monomial or
        self.parent()._latex_monomial.
        """
        cdef dict = self.parent().print_options()
        cdef path_homog_poly_t *H = self.data
        cdef list L, L_total
        cdef size_t i
        cdef path_term_t * T
        L_total = []
        while H != NULL:
            L = []
            T = H.poly.lead
            while T!=NULL:
                L.append((([biseq_getitem(T.mon.path,i) for i in range(T.mon.path.length)],
                       T.mon.pos, T.mon.mid, H.start),<object>(T.coef)))
                T = T.nxt
            if len(L) != H.poly.nterms:
                print "Term count of polynomial is wrong, got",len(L), "expected", H.poly.nterms
            L_total.extend(L)
            H = H.nxt
        return L_total

    def _repr_(self):
        return repr_lincomb(self._sorted_items_for_printing(), strip_one=True,
                            repr_monomial = self._repr_monomial,
                            scalar_mult       = '*',
                            latex_scalar_mult = None,
                            )

    def _latex_(self):
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult       = '*',
                            latex_scalar_mult = None,
                            repr_monomial = self._latex_monomial,
                            is_latex=True, strip_one = True)

    cpdef dict monomial_coefficients(self):
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

    def __getitem__(self, k):
        cdef path_homog_poly_t * H
        cdef PathAlgebraElement out
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

    def components(self):
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
    # Comparison

    def __hash__(self):
        cdef long out = 0
        cdef path_homog_poly_t *H = self.data
        while H != NULL:
            out = (out<<5) | (out>>(sizeof(long)-5))
            out += (H.start)<<7
            out += (H.end)<<11
            out += poly_hash(H.poly)
            H = H.nxt
        if out == -1:
            return -2
        return out

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

## multiplication

    cpdef RingElement _mul_(self, RingElement  other):
        """
        EXAMPLES::

            sage: from sage.quivers.algebra_elements import PathAlgebraElement
            sage: A = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup().algebra(ZZ.quo(15))
            sage: A.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: x = PathAlgebraElement(A, b*e*b*e+4*b*e+e_0)
            sage: y = PathAlgebraElement(A, a*c+5*f*e)
            sage: x*y   # no leak
            sage: y*x   # no leak
            sage: y*y   # no leak
            sage: x*x   # leak

        BUG::

            sage: x = PathAlgebraElement(A, sage_eval('b*e*b*e+4*b*e+e_0', A.gens_dict()))
            sage: y = PathAlgebraElement(A, sage_eval('a*c+d*c*b*f', A.gens_dict()))
            sage: x*(y+x) == x*y+x*x
            False

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
                            poly_iadd_lmul(out_orig.poly, <object>T2.coef, H1.poly, T2.mon.path,
                                           self.cmp_terms, 0, 0, -1)
                            T2 = T2.nxt
                    else:
                        while T2 != NULL:
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
