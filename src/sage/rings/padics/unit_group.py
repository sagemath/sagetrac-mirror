r"""
Unit groups of p-adic fields

EXAMPLES::

   sage: x = polygen(QQ)
   sage: f = x^6 + 3*x^3 + 3
   sage: K.<a> = Qp(3).ext(f)
   sage: U = pAdicUnitGroup(K)
   sage: 
"""

from sage.groups.abelian_gps.values import AbelianGroupWithValues_class, AbelianGroupWithValuesElement
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.arith.srange import xsrange
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.misc.misc_c import prod
from sage.databases.cremona import cremona_letter_code
from copy import copy
from collections import Counter

class pAdicUnit(AbelianGroupWithValuesElement):
    """
    """
    def _div_(left, right):
        """
        Divide ``left`` by ``right``

        TESTS::


        """
        G = left.parent() # = right.parent() by coercion model
        value = left.value() / right.value()
        prec = value.precision_relative()
        exponents = [ (x-y)%order if order !=0 else x-y
                      for x, y, order in
                      zip(left._exponents, right._exponents, G.gens_orders(prec)) ]
        return G.element_class(G, exponents, value)

    def _mul_(left, right):
        """
        Multiply ``left`` and ``right``
        """
        G = left.parent()
        value = left.value() * right.value()
        prec = value.precision_relative()
        exponents = [ (x+y)%order if order != 0 else x+y
                      for x, y, order in
                      zip(left._exponents, right._exponents, G.gens_orders(prec)) ]
        return G.element_class(G, exponents, value)

    def __pow__(self, n):
        """
        Exponentiation.
        """
        G = self.parent()
        value = self.value()**n
        prec = value.precision_relative()
        exponents = [ (n*x)%order if order != 0 else n*x
                      for x, order in
                      zip(self._exponents, G.gens_orders(prec)) ]
        return G.element_class(G, exponents, value)

    def _repr_(self):
        return repr(self.value())

class pAdicUnitGroup(AbelianGroupWithValues_class, UniqueRepresentation):
    Element = pAdicUnit
    @staticmethod
    def __classcall__(cls, K, gens=None, additive_basis=None, e0_basis=None, prec=None, names=None, check=True):
        field = K.is_field()
        f = K.absolute_f()
        k = K.residue_field()
        V = k.vector_space()
        def check_basis(basis):
            if len(basis) != f:
                raise ValueError("Wrong number of elements in basis")
            basis = tuple(K(c) for c in basis)
            if matrix([V(k(c)) for c in basis]).rank() != f:
                raise ValueError("Does not reduce to a basis for the residue field")
            return basis
        p = K.prime()
        q = k.cardinality()
        e = K.absolute_e()
        n = K.absolute_degree()
        a = K(k.gen()).lift_to_precision()
        if prec is None:
            prec = K.precision_cap()
        eps = -K(p).unit_part()
        qthzeta, s, _ = K._primitive_qth_root_of_unity(infinity)
        mup_case = (s >= 1) # and prec >= e/(p-1))
        if mup_case:
            v, u = e.val_unit(p)
            e0 = u // (p-1)
        if check and gens is not None:
            gens = tuple(K(gen) for gen in gens)
            expected_len = n
            if q != 2:
                expected_len += 1
                if gens[0].multiplicative_order() != q-1:
                    raise ValueError("First entry must be a multiplicative generator modulo the maximal ideal")
            if mup_case:
                expected_len += 1
                pzeta_index, ordinate = (0, "First") if q == 2 else (1, "Second")
                if gens[pzeta_index].multiplicative_order() != p**s:
                    raise ValueError(ordinate + " entry must be a p^%s root of unity"%s)
            if field:
                expected_len += 1
                if gens[-1].valuation() != 1:
                    raise ValueError("Last entry must be a uniformizer")
            if len(gens) != expected_len:
                raise ValueError("Provided generators have length %s but should have length %s"%(len(gens), expected_len))
        if additive_basis is None:
            additive_basis = tuple(a**i for i in range(f))
        elif check:
            additive_basis = check_basis(additive_basis)
        if e0_basis is None:
            if mup_case:
                if s == v + 1:
                    e0_basis = [(qthzeta - 1) >> e0]
                else:
                    e0_basis = [K(k(eps).nth_root(e//e0)).lift_to_precision(prec)]
                vec = V(k(e0_basis[0]))
                j = 0
                while vec[j] == 0:
                    j += 1
                e0_basis += [a**i for i in range(f) if i != j]
                e0_basis = tuple(e0_basis)
        elif check and mup_case:
            e0_basis = check_basis(e0_basis)
            if mup_case and (e0_basis[0]**(e//e0) - eps).valuation() == 0:
                raise ValueError("First entry in basis survives in level e/(p-1)")
        elif not mup_case:
            e0_basis = None
        return super(pAdicUnitGroup, cls).__classcall__(cls, K, gens, additive_basis, e0_basis, qthzeta, s, prec, names)
    def __init__(self, K, gens, additive_basis, e0_basis, qthzeta, s, prec, names):
        self._prec_cap = prec
        self._field = K.is_field()
        self.K = K.fraction_field()
        self.OK = K.integer_ring()
        self.k = k = K.residue_field()
        self.V = V = k.vector_space()
        self.p = p = K.prime()
        self.q = q = k.cardinality()
        self.n = n = K.absolute_degree()
        self.e = e = K.absolute_e()
        self.f = f = K.absolute_f()
        self.pi = pi = K.uniformizer()
        v, u = e.val_unit(p)
        self.mu0 = mu0 = v + 1
        self.e0 = e0 = u // (p-1) # only used in the mup_case, when u divisible by p-1
        self.eps = eps = -K(p).unit_part()
        self.cutoff = cutoff = p*e/(p-1)
        self.cutoff_floor = cutoff_floor = cutoff.floor()
        self._resgen = K(k.gen()).lift_to_precision(prec)
        self._mup_case = mup_case = (s >= 1) # and prec >= e/(p-1))
        self._additive_basis = additive_basis
        self._e0_basis = e0_basis

        self._pbasis = []
        orders = [ZZ(0)] * n
        for level in self.fundamental_levels():
            if mup_case and level == e0:
                self._pbasis.extend([1 + w * pi**level for w in e0_basis])
            else:
                self._pbasis.extend([1 + w * pi**level for w in additive_basis])
        if mup_case:
            if gens is None and s != mu0 and prec < cutoff + (mu0 - s)*e:
                # This may never be reached, since in the low precision case
                # an error may already have been raised by K._primitive_qth_root_of_unity
                # (or even the Krasner check, once that's enabled).
                # If prec is too small, we will run into problems in the
                # Smith normal form when the smallest order element won't be the root of unity
                raise NotImplementedError("Need more precision")
            self._wstar = wstar = self._find_wstar()
            self._eta_star = eta_star = 1 + wstar * pi**cutoff_floor
            self._pbasis.append(eta_star)
            orders = [p**s] + orders
        if q != 2:
            orders = [q-1] + orders
        if self._field:
            orders.append(ZZ(0))

        if gens is None:
            self._reszeta = reszeta = k.multiplicative_generator()
            reszeta = K.teichmuller(reszeta, prec)
            if mup_case:
                gens, self._to_gens = self._pbasis_to_gens(reszeta, qthzeta, s)
            else:
                gens = [reszeta] + self._pbasis
                self._to_gens = None
            if self._field:
                gens += (pi,)
        else:
            if q == 2:
                self._reszeta = k(1)
            else:
                self._reszeta = k(gens[0])
            if self._field:
                self.pi = gens[-1]
            self._to_gens = self._given_to_gens(gens, s)
        self._gens_data = self._compute_gens_data(gens)
        if names is None:
            names = self._name_the_gens(gens, s)

        AbelianGroupWithValues_class.__init__(self, tuple(orders), tuple(names), gens, K)

    def _find_wstar(self):
        # In this case, the p-th powering map is not surjective mapping from
        # 1 + a*pi^v + O(pi^{v+1}) to
        # 1 + (a^p - eps*a)*pi^{pv} + O(pi^{pv+1})
        # when v = e/(p-1).
        # eta_star generates the cokernel, of order p.
        p, eps, K, k, f, prec = self.p, self.eps, self.K, self.k, self.f, self._prec_cap
        R = k['x']
        a = k.gen()
        x = R.gen()
        eps_bar = k(eps)
        for i in range(f):
            poly = x**p - eps_bar*x - a**i
            try:
                poly.any_root()
            except ValueError:
                return K(a**i).lift_to_precision(prec)
        else: # this shouldn't happen
            raise RuntimeError("pth power map at cutoff level surjective")

    def _pbasis_to_gens(self, reszeta, qthzeta, s):
        p, e, n, e0, mu0 = self.p, self.e, self.n, self.e0, self.mu0
        prec = self._prec_cap
        S = self._pbasis
        i0 = self._index(e0)
        b = S[i0]**(p**mu0)
        if b == 1:
            # We have full pth-power roots of unity already, so the change of basis is just a swap
            gens = copy(S)
            gens.insert(0,gens.pop(i0))
            change_of_basis = identity_matrix(n+1)
            change_of_basis.swap_rows(0, i0)
            return gens, change_of_basis
        # We don't have full pth-power roots of unity, so
        # we have to modify the basis using Smith normal form

        rel = self._dlog(b)
        ind = 0
        M = []
        for level in self.fundamental_levels():
            order = self._compute_order(level, prec)
            for i in range(self.f):
                if level == e0 and i == 0:
                    row = rel
                    row[ind] = -p**mu0
                else:
                    row = [ZZ(0)] * len(rel)
                    row[ind] = order
                M.append(row)
                ind += 1
        D, change_of_basis, _ = matrix(M).transpose().smith_form()
        reverse_change = (~change_of_basis).change_ring(ZZ)
        gens = [prod(c**d for c,d in zip(S, reverse_change.column(i))) for i in range(self.n+1)]
        # We reorder so that the generators are sorted by level
        # (aside from having the first be the root of unity)
        top_row = change_of_basis.row(0)
        L = zip(gens[1:], change_of_basis.rows()[1:])
        L.sort(key=lambda u: (u[0]-1).valuation())
        gens, rows = zip(*L)
        change_of_basis = matrix((top_row,) + rows)
        # We now need to change the first generator to the actual root of unity
        # (it could differ in high order terms)
        level_coords = self._dlog(qthzeta)
        coords = change_of_basis * vector(level_coords)
        if coords[0] % p**s != 1:
            raise NotImplementedError("Need more precision")
        for j in range(1,len(coords)):
            change_of_basis.add_multiple_of_row(j,0,-coords[j])
        # Should probably reduce change of basis modulo appropriate powers of p

        gens = (reszeta, qthzeta) + gens
        return gens, change_of_basis

    def _given_to_gens(self, gens):
        # Strip things other than principal units
        if self.q != 2:
            gens = gens[1:]
        if self._field:
            gens = gens[:-1]
        reverse_change = []
        for gen in gens:
            reverse_change.append(self._dlog(gen))
        try:
            # TODO: May only be invertible modulo a power of p...
            change_of_basis = (~matrix(reverse_change)).change_ring(ZZ)
        except TypeError:
            raise ValueError("Generators provided do not span the torsion-free quotient")
        # Should probably reduce change of basis modulo appropriate powers of p
        return change_of_basis

    def _std_name(self, level, i):
        if self.f == 1:
            return 'u%s'%level
        else:
            return 'u%s_%s'%(level, i)

    def _name_the_gens(self, gens, s):
        names = []
        unif = gens[-1] # only used when self._field is True
        if self.q != 2:
            names.append('zeta_u')
            gens = gens[1:]
        if self._mup_case:
            names.append('zeta_p%s'%('' if s == 1 else s))
            name_counter = Counter()
            gens = gens[1:]
            if self._field:
                gens = gens[:-1]
            for gen in gens:
                # level, u = (gen-1).val_unit() #FIXME
                level, u = (gen-1).valuation(), (gen-1).unit_part()
                i = self.k(u).polynomial().degree()
                if u == self._resgen**i:
                    names.append(self._std_name(level, i))
                else:
                    nlevel = name_counter[level]
                    names.append('u%s_%s'%(level, cremona_letter_code(nlevel)))
                    name_counter[level] += 1
        else:
            for level in self.fundamental_levels():
                names.extend([self._std_name(level, i) for i in range(f)])
        if self._field:
            if unif == K.uniformizer():
                names.append(K._uniformizer_print())
            else:
                names.append('pi')
        return names

    def _compute_gens_data(self, gens):
        p, e0, mu0, cutoff_floor = self.p, self.e0, self.mu0, self.cutoff_floor
        gens_data = []
        if self._to_gens is None:
            # Using standard pbasis
            for level in self.fundamental_levels():
                gens_data.append((level, (cutoff_floor // level).exact_log(p)))
        else:
            if self.q != 2:
                gens = gens[1:]
            if self._field:
                gens = gens[:-1]
            for gen in gens:
                level = (gen-1).valuation()
                max_s1 = (cutoff_floor // level).exact_log(p)
                extra_pi = ZZ(0)
                if self._mup_case:
                    flevel, s1, s2, s3 = self._fundamental_level(level)
                    if flevel == e0 and level < cutoff_floor:
                        cutoff_pow = gen**(p**max_s1)
                        extra_pi = (cutoff_pow-1).valuation() - cutoff_floor
                gens_data.append((level, max_s1, extra_pi))
        return gens_data

    def _dlog_simple(self, b):
        """
        INPUT:

        ``b`` -- a principal unit in the ring or field

        OUTPUT:

        A list of exponents expressing b as a product of powers of the ``_pbasis``
        """
        k, p = self.k, self.p
        S = self._pbasis
        decomp = [ZZ(0)] * self.n
        if self._mup_case:
            decomp.append(ZZ(0))
        while b != 1:
            z = b - 1
            level = z.valuation()
            reducer, flevel, s1, s2, s3, star_case, tmp_basis = self._additive_reducer(level)
            ppow = p**(s1+s2+s3)
            m = self._index(flevel)
            vec = reducer * self.V(k(z.unit_part()))
            for i, c in enumerate(vec):
                if star_case and i == 0:
                    cp = ZZ(c) * p**s3
                    decomp[-1] += cp
                    b_new = b / S[-1]**cp
                else:
                    cp = ZZ(c) * ppow
                    decomp[m+i] += cp
                    b_new = b / S[m+i]**cp
            if (b_new - 1).valuation() <= level:
                raise RuntimeError("Level did not increase")
            b = b_new
        return decomp

    def _dlog(self, b, algorithm=None):
        return self._dlog_simple(b)

    def discrete_log(self, b, algorithm=None):
        K = self.K
        b = K(b)
        prec = b.precision_relative()
        gens = self._values
        v = b.valuation()
        if v != 0:
            if not self._field:
                raise ValueError("Not a unit")
            b /= gens[-1]**v
        if self.q != 2:
            u = self.k(b)
            if u != 1:
                b /= K.teichmuller(u, prec)
            ulog = u.log(self._reszeta)
        coords = self._dlog(b, algorithm)
        if self._to_gens is not None:
            coords = list(self._to_gens * vector(coords))
        orders = self.gens_orders(prec)
        if self.q != 2:
            coords = [ulog] + coords
        for i in range(len(coords)):
            coords[i] = coords[i] % orders[i]
        if self._field:
            coords.append(v)
        return coords

    def _element_constructor_(self, x):
        if isinstance(x, (list, tuple)):
            if len(x) != len(self._values):
                raise ValueError("Length of input must match the number of generators (%s)"%(len(self._values)))
            coords = x
            x = prod(u**c for u,c in zip(coords, self._values))
        else:
            coords = None
        if isinstance(x, pAdicUnit):
            L = x.value().parent()
            if L is self.K or L is self.OK:
                coords = x._exponents
            x = x.value()
        x = self.K(x) if self._field else self.OK(x)
        if coords is None:
            coords = self.discrete_log(x)
        return self.element_class(self, coords, x)

    def fundamental_levels(self):
        p = self.p
        for i in xsrange(1, self.cutoff):
            if not p.divides(i):
                yield ZZ(i)

    def _fundamental_level(self, level):
        """
        Returns the corresponding fundamental level, together with the numbers of p-th powers used.
        """
        p, e, cutoff, cutoff_floor = self.p, self.e, self.cutoff, self.cutoff_floor
        if level > cutoff:
            s3 = 1 + (level - cutoff_floor - 1) // e
            level -= e * s3
        else:
            s3 = 0
        if level == cutoff:
            s2 = 1
            level -= e
        else:
            s2 = 0
        s1, level = level.val_unit(p)
        return level, s1, s2, s3

    def _additive_reducer(self, level):
        """
        Returns various data used in computing discrete logs

        INPUT:

        - ``level`` -- a positive integer

        OUTPUT:

        - ``reducer`` -- a change of basis matrix between the standard basis
           of ``(1+P^v)/(1+P^(v+1))`` and the appropriate power of the
           fundamental level.
        - ``flevel`` -- the fundamental level corresponding to ``level``
        - ``s1`` -- the number of `p`th powers needed to go from the fundamental level to the cutoff.
        - ``s2`` -- `0` or 
        """
        k, V = self.k, self.V
        p, eps = self.p, k(self.eps)
        flevel, s1, s2, s3 = self._fundamental_level(level)
        if self._mup_case and flevel == self.e0:
            basis = self._e0_basis
        else:
            basis = self._additive_basis
        basis = [k(w)**(p**s1) for w in basis]
        if s2:
            basis = [w**p - eps*w for w in basis]
        basis = [(-eps)**s3 * w for w in basis]
        # If reduction passes through the cutoff, we have to use eta_star rather than eta_e0,1
        star_case = (self._mup_case and flevel == self.e0 and level >= self.cutoff)
        if star_case:
            basis[0] = (-eps)**s3 * k(self._wstar)
        reducer = matrix([V(w) for w in basis]).transpose().inverse()
        return reducer, flevel, s1, s2, s3, star_case, basis

    def _index(self, level, i=0):
        """
        Returns the index in the list of generators given a fundamental level and a shift (corresponding to fixing one of the `f` generators at that level).

        INPUT:

        - ``level`` -- a fundamental level
        - ``i`` -- an integer between `0` and `f-1`

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^6 + 3*x^3 + 3
            sage: K.<w> = Qp(3).ext(f)
            sage: U = pAdicUnitGroup(K)
            sage: [U._index(level) for level in [1,2,4,5,7,8]]
            [0, 1, 2, 3, 4, 5]
        """
        return ((level-1) - (level-1) // self.p) * self.f + i

    def _compute_order(self, level, prec, max_s1=None, extra_pi=0):
        """
        Compute the order of a generator `u` of specified level at a given precision.

        INPUT:

        - ``level`` -- The valuation of `u-1`
        - ``prec`` -- The precision at which a power of `u` is considered `1`.
        - ``max_s1`` -- ``log_p(floor(cutoff/level))``, computed if not provided.
        - ``extra_pi`` -- Extra valuation obtained when pth powering past the cutoff.
          Can only be nonzero for fundamental level `e_0`.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^6 + 3*x^3 + 3
            sage: K.<w> = Qp(3).ext(f)
            sage: U = pAdicUnitGroup(K)
            sage: U._compute_order(2, 12)
            9
            sage: (1 + w^2)^9
            1 + 2*w^12 + ...
        """
        p, e, cutoff = self.p, self.e, self.cutoff_floor
        if level >= prec:
            return ZZ(1)
        elif prec <= cutoff:
            return p**(1 + ((prec-1) // level).exact_log(p))
        if max_s1 is None:
            max_s1 = (cutoff // level).exact_log(p)
        if prec <= cutoff + extra_pi:
            return p**max_s1
        else:
            return p**(1 + max_s1 + (prec - 1 - extra_pi - level * p**max_s1)//e)

    def gens_orders(self, prec=None):
        """
        Returns the orders of the generators.

        INPUT:

        - ``prec`` -- Non-negative integer or ``None``.  If given,
          will return the order of the generators modulo the ``prec``
          power of the uniformizer.

        OUTPUT:

        A tuple of integers, each giving the order of the corresponding
        entry of the :meth:`gens` method.  As usual, zero indicates
        infinite order.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^6 + 3*x^3 + 3
            sage: K.<w> = Qp(3).ext(f)
            sage: U = pAdicUnitGroup(K)

        If no precision given, returns the mathematical orders of the elements::

            sage: U.gens_orders()
            (2, 3, 0, 0, 0, 0, 0, 0, 0)

        When a precision is given, the order of each generator
        up to that precision is returned::

            sage: U.gens_orders(9)
            (2, 3, 9, 9, 9, 3, 3, 3, 0)
            sage: U.gen(2)^9
            1 + 2*w^13 + ...
        """
        if prec is None:
            # We're modeling a group that has infinite order elements, but elements have
            # finite order at finite precision.  When no precision is specified, we return
            # the mathematical orders
            return self._gens_orders
        # Now the user has requested the orders of the generators in (1+P)/(1+P^k)
        p, e = self.p, self.e
        resorder = ZZ(1) if prec == 0 else self.q-1
        if prec == 0:
            resorder = ZZ(1)
        else:
            resorder = self.q - 1
        orders = []
        if self._to_gens is None:
            for level, max_s1 in self._gens_data:
                order = self._compute_order(level, prec, max_s1)
                orders.extend([order] * self.f)
        else:
            for level, max_s1, extra_pi in self._gens_data:
                orders.append(self._compute_order(level, prec, max_s1, extra_pi))
        if self.q != 2:
            orders = [resorder] + orders
        if self._field:
            orders = orders + [ZZ(0)]
        return tuple(orders)

    @property
    def _pbasis_gens(self):
        """
        Stores the generators in the basis of 1+P as elements of this group.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^6 + 3*x^3 + 3
            sage: K.<w> = Qp(3).ext(f)
            sage: U = pAdicUnitGroup(K)
            sage: B = U._pbasis_gens; B
            [1 + 2*w + O(w^120),
             1 + w^2 + O(w^120),
             1 + w^4 + O(w^120),
             1 + w^5 + O(w^120),
             1 + w^7 + O(w^120),
             1 + w^8 + O(w^120),
             1 + w^9 + O(w^120)]
            sage: all(b.parent() is U for b in B)
            True
        """
        return [self(u) for u in self._pbasis]

    def gens_level(self, level):
        """
        Returns a list of generators for the group of units ``u`` with
        ``(u-1).valuation() >= level``.

        INPUT:

        - ``level`` -- a non-negative integer

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^6 + 3*x^3 + 3
            sage: K.<w> = Qp(3).ext(f)
            sage: U = pAdicUnitGroup(K)
            sage: U.gens_level(5)
            [1 + w^6 + 2*w^8 + 2*w^10 + w^11 + w^13 + O(w^120),
             1 + 2*w^10 + w^12 + w^13 + 2*w^14 + w^17 + O(w^120),
             1 + w^5 + O(w^120),
             1 + w^7 + O(w^120),
             1 + w^8 + O(w^120),
             1 + w^9 + O(w^120)]
        """
        if level < 0:
            raise ValueError("level must be non-negative")
        elif level == 0:
            gens = self._pbasis_gens
            if self.q != 2:
                gens = [self(self._values[0])] + gens
        else:
            p = self.p
            exponents = [self._compute_order(flevel, level) for flevel in self.fundamental_levels() for _ in range(self.f)]
            if self._mup_case:
                if level > self.cutoff / p:
                    exponents[self._index(self.e0)] = None
                exponents.append(self._compute_order(self.cutoff_floor, level))
            gens = [u**a for u,a in zip(self._pbasis_gens, exponents) if a is not None]
        return gens

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: f = x^6 + 3*x^3 + 3
            sage: K.<w> = Qp(3).ext(f)
            sage: U = pAdicUnitGroup(K)
            sage: U # indirect doctest
            Unit group with structure Z/2 x Z/3 x Z x Z3^6 of 3-adic Eisenstein Extension Field in w defined by x^6 + 3*x^3 + 3
        """
        structure = ""
        if self.q != 2:
            structure += ' x Z/%d' % (self.q-1)
            ps = self._gens_orders[1]
        else:
            ps = self._gens_orders[0]
        if ps > 1:
            structure += ' x Z/%d' % ps
        if self._field:
            structure += ' x Z'
            ring = self.K
        else:
            ring = self.OK
        structure += ' x Z%d^%d' % (self.p, self.n)
        structure = structure[3:]
        return "Unit group with structure %s of %s" % (structure, ring)

