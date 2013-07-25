r"""

Contains classes and reduction algorithms for Hilbert modular groups.

AUTHORS:

- Fredrik Stromberg (2013): copied from stoll_cremona.py
"""

#*****************************************************************************
#       Copyright (C) 2013 Florian Bouyer <f.j.s.c.bouyer@gmail.com>,
#                          Marco Streng <marco.streng@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject

def _HilbertFundDomData(K, embs=None, group='sl'):
    """
    Create an object containing the data for a fundamental domain
    for the action of SL_2(O_K) on (CC\setminus RR)^g for a totally real
    field K of degree g.

    Warning: assumes deg(K) <= 2 and cl(K) = 1.

    Input:

     - ``K`` -- A totally real number field.

     - ``embs`` (default=None) -- None, or the complete set
       of embeddings of K into a real field that should be used,
       or the precision to use for creating such a set of embeddings.

     - ``group`` (default='sl') -- a string: 'sl' for SL_2(O_K),
       'gl' for GL_2(O_K), 'plus' for 
       the subgroup GL_2(O_K)^+ of elements of totally positive
       determinant in GL_2(O_K).
       
    EXAMPLES::
        
        sage: from sage.schemes.hyperelliptic_curves.stoll_cremona import _HilbertFundDomData
        sage: K = QuadraticField(41,'a')
        sage: K(K.unit_group().gen(1)).norm()
        -1
        sage: z = [0.01*CC.gen()+23.34,0.02+0.03*CC.gen()]
        sage: h = _HilbertFundDomData(K, group='sl')
        sage: r = h.reduce(z); r
        (
        [    9/2*a - 11/2 -125/2*a + 801/2]                                     
        
        [   -3/2*a - 17/2         2*a - 13], [1.03241154280086 +
        15.5033759282630*I, -1.33856972357549 + 0.0496716993976175*I]
        )
        sage: prod([e.imag() for e in r[1]]), max([log(abs(e.imag())) for e in r[1]])
        (0.770079028756939, 2.74105780203325)
        sage: h = _HilbertFundDomData(K, group='gl')
        sage: r = h.reduce(z); r
        (
        [15/2*a + 99/2  5/2*a - 37/2]                                           
        
        [-3/2*a - 17/2      2*a - 13], [1.83799715597054 + 0.242181137053011*I,
        -2.13249730501399 - 3.17976469236058*I]
        )
        sage: prod([e.imag() for e in r[1]]), max([log(abs(e.imag())) for e in r[1]])
        (-0.770079028756901, 1.15680719794121)
    """
    if K.degree() == 1:
        return _HilbertFundDomData_classical(embs)
    if K.degree() == 2:
        return _HilbertFundDomData_quadratic(K, embs, group)
    return _HilbertFundDomData_generic(K, embs, group)



class _HilbertFundDomData_generic(SageObject):
    """
    An object containing the data for a fundamental domain
    for the action of SL_2(O_K) on HH^g for a totally real
    field K of degree g.
    """

    _K = None
    _embs = None
    _RR = None

    def __init__(self, K, embs, group):
        """
        Create the object.
        """
        self._K = K
        if embs == None:
            embs = K.embeddings(RR)
            self._RR = RR
        elif embs in ZZ:
            self._RR = RealField(embs)
            embs = K.embeddings(self._RR)
        else:
            self._RR = embs[0].codomain()
        self._embs = embs

    def to_additive(self, x):
        """
        Move x in RR^g to an additive fundamental domain.

        Returns an pair (a, x') with a in O_K and x' = x-a in the fundamental domain.
        """
        raise NotImplementedError, "Sorry reduction to a Hilbert fundamental domain is not implemented for g>2"

    def step(self, z, bound_on_c=None):
        """
        Do one reduction round.

        Returns a pair (M, z') with z' = M(z) nearer the fundamental domain if possible.
        """
        M1 = self.up(z, bound_on_c=bound_on_c)
        z = self.act(M1, z)
        M2 = self.to_multiplicative([e.imag() for e in z])[0]
        z = self.act(M2, z)
        t = self.to_additive([e.real() for e in z])[0]
        M3 = Matrix([[1, -t], [0, 1]])
        z = self.act(M3, z)
        return M3*M2*M1, z

    def reduce(self, z, max_rounds = 0, bound_on_c=None):
        """
        Do reduction rounds until reduced.

        Returns a pair (M, z') with z' = M(z) in the fundamental domain.

        max_rounds = 0: go on until finishes, risking an infinite loop from bugs or rounding errors
        max_rounds > 0: do at most that number of rounds
        """
        M = Matrix([[1,0],[0,1]])
        rounds = 0
        while true:
            rounds = rounds + 1
            M1, z = self.step(z, bound_on_c=bound_on_c)
            M = M1 * M
            if M1 == 1 or rounds == max_rounds:
                return M, z

    def act(self, M, z):
        """
        Return M(z)
        """
        a = M[0,0]
        b = M[0,1]
        c = M[1,0]
        d = M[1,1]
        embs = self._embs
        return [((embs[k](a)*z[k] + embs[k](b))/(embs[k](c)*z[k] + embs[k](d))) for k in range(len(embs))]

    def _up_cd_lll(self, z):
        """
        As _up_cd, but instead of making prod(Im(z)) maximal
        and log(Im(z)) in a fundamental domain, try to make min(Im(z)) large
        using LLL.
        """
        # The new imaginary part Im(z) / |cz+d|^2
        # In particular, max(|cz+d|^2/Im(z)) must be as small as possible,
        # and less than 1 / mim(Im(z)).
        # Let v = ((cz+d)/sqrt(Im(z)))_i.
        # Then we want |v|_infty to be as small as possible and less than
        # 1 / min(sqrt(Im(z))).
        # Note, g*|v|_2 >= g*|v|_infty^2 >= |v|_2, so it is necessary that
        # |v|_2 < g/(min(Im(z))) and sufficient that |v|_2 < 1/(min(Im(z)))
        g = self._K.degree()
        embs = self._embs
        m = 2**(z[0].parent().precision())
        O = self._K.maximal_order()
        vectors = ( [ [(embs[i](b)) for i in range(g)] +
        [0 for i in range(g)] for b in O.basis()] +
        [ [(embs[i](b)*z[i].real()) for i in range(g)] +
        [ (embs[i](b)*z[i].imag()) for i in range(g)]
        for b in O.basis()])
        id = identity_matrix(2*g)
        vectors = [id[j].list() +
                   [(m*vectors[j][i]/(z[i%g].imag().abs().sqrt())).round()
                      for i in range(len(vectors[j]))] for j in range(2*g)]
#        c = []
#        for i in range(len(vectors)):
#            for j in range(i, len(vectors)):
#                c.append(sum([vectors[i][k]*vectors[j][k] for k in range(2*g)]))
#        Q = QuadraticForm(ZZ, 2*g, c)
#        l = Q.short_vector_list_up_to_length(g/min([z.imag().abs()]))
#        for x in l[1:]:
        M = Matrix(vectors).LLL()
        d = sum([M[0][i]*O.basis()[i] for i in range(g)])
        c = sum([M[0][g+i]*O.basis()[i] for i in range(g)])
        B = c * O + d * O
        if B == 1*O:
            return c, d
        if B.is_principal():
            g = B.gens_reduced()[0]
            return c/g, d/g
        # TODO: if B is not principal, find another candidate pair c, d

    def up(self, z, bound_on_c=None):
        """
        Move z in HH^g, i.e., increase the norm of its imaginary part,
        if possible.

        Returns either None or a matrix M in SL_2(O_K) that such that
        Mz has a larger norm of its imaginary part.

        Assumes K has class number one.

        EXAMPLES::

            sage: from sage.schemes.hyperelliptic_curves.stoll_cremona import _HilbertFundDomData
            sage: h = _HilbertFundDomData(QuadraticField(15))
            sage: h.up([123+CC.0*0.05, 456+CC.0*0.0234])
            [          0          -1]
            [          1 -43*a - 290]
            sage: h.act(_, [123+CC.0*0.05, 456+CC.0*0.0234])
            [2.14072851429831 + 0.231823024326897*I, 1.85425166669223 + 0.0806070738044104*I]
            h.reduce([123+CC.0*0.05, 456+CC.0*0.0234])
            (
            [      2*a + 3 -708*a - 2158]
            [           -2    86*a + 579], [-0.0404815976687969 + 3.15206362011738*I, -0.618912376653362 + 2.90579848705465*I]
            )
        """
        cd = self._up_cd(z, bound_on_c=bound_on_c)
        if cd == None:
            return Matrix([[1,0],[0,1]])
        c, d = cd
        a = d.inverse_mod(c)
        b = (a*d-1)/c
        return Matrix([[a,b],[c,d]])

    def _up_cd(self, z):
        """
        Move z up in HH^g, i.e., increase the norm of its imaginary part,
        if possible. Actually, not implemented in this base class: just calls
        _up_cd_lll instead.
        """
        return self._up_cd_lll(z)


class _HilbertFundDomData_quadratic(_HilbertFundDomData_generic):
    """
    As _HilbertFundDomData, but only for quadratic fields.
    """

    _eps = None
    _epsmat = None
    _epslog = None
    _RR = None

    def __init__(self, K, embs, group):
        """
        Create the object.
        """
        _HilbertFundDomData_generic.__init__(self, K, embs, group)
        eps = K(K.unit_group().gens()[1])
        if group == 'sl' or (group == 'plus' and eps.norm() == -1):
            self._epsmat = Matrix([[eps, 0], [0, eps**-1]])
            eps = eps**2
        elif group in ['plus', 'gl']:
            if self._embs[0](eps) < 0:
                eps = -eps
            self._epsmat = Matrix([[eps, 0], [0, 1]])
        else:
            raise ValueError, "Unknown group in _HilbertFundDomData: %s" %group
        self._eps = eps
        self._epslog = self._RR(self._embs[1](eps).abs().log())

    def to_additive(self, x):
        """
        Move x in RR^g to an additive fundamental domain.

        Returns an pair (a, x') with a in O_K and x' = x-a in the fundamental domain.

        Examples::

            sage: from sage.schemes.hyperelliptic_curves.stoll_cremona import _HilbertFundDomData
            sage: h = _HilbertFundDomData(QuadraticField(5))
            sage: h.to_additive([3.2,7.9])
            (a + 6, [-0.563932022500210, -0.336067977499789])
        """
        embs = self._embs
        K = self._K
        om = K.maximal_order().gens()[1]
        y = ((x[1]-x[0]) / (embs[1](om) - embs[0](om))).round()
        a = y * om
        x = [x[0] - y*embs[0](om), x[1] - y*embs[1](om)]
        y = ((x[1]+x[0]) / 2).round()
        a = a + y
        x = [x[0] - y, x[1] - y]
        return a, x

    def to_multiplicative(self, y):
        """
        Move y in RR_+^g to a multiplicative fundamental domain.

        Returns a pair (M, y') with M a diagonal matrix and M*y = y' and y' in a fundamental domain.

        Examples::

            sage: from sage.schemes.hyperelliptic_curves.stoll_cremona import _HilbertFundDomData
            sage: h = _HilbertFundDomData(QuadraticField(15))
            sage: h.to_multiplicative([1,45678])
            (
            [-a + 4      0]
            [     0  a + 4], [61.9838667696635, 736.933695500892]
            )
            sage: h.act(_[0], [CC.0, CC.0*45678])
            [61.9838667696594*I, 736.933695500891*I]
            sage: h = _HilbertFundDomData(QuadraticField(15), sl=False)
            sage: h.to_multiplicative([1,45678])
            (
            [-63*a + 244           0]
            [          0           1], [487.997950815860, 93.6028520695257]
            )
            sage: h.act(_[0], [CC.0, CC.0*45678])
            [487.997950811067*I, 93.6028520686063*I]
        """
        embs = self._embs
        K = self._K
        eps = self._eps
        l = self._epslog
        e = ((y[1].abs().log() - y[0].abs().log()) / (2*l)).round()
        a = eps**e
        y = [y[0]/embs[0](a), y[1]/embs[1](a)]
        return self._epsmat**-e, y


    def _up_cd(self, z, bound_on_c=None):
        """
        Move z in HH^g, i.e., increase the norm of its imaginary part,
        if possible.

        Returns either None or a matrix M in SL_2(O_K) that such that
        Mz has a larger norm of its imaginary part.

        Assumes K has class number one.

        bound_on_c None means try a small one first, double each time until a proven one is reached
        """
        embs = self._embs
        cd = self._up_cd_lll(z)
        K = self._K
        O = K.maximal_order()
        if not cd is None:
            (c, d) = cd
            cz_plus_d = [embs[i](c)*z[i] + embs[i](d) for i in range(2)]
            if abs(prod(cz_plus_d)) < 1:
                assert c*O+d*O == 1*O
                return c, d
        # Assuming class number 1 and the above did not work, we are pretty
        # high up in the fundamental domain, so the search ranges below are not
        # that large.
        # TODO: two different notions of reduced are tested, so something
        # may be missed here. In that case, we may be high up after all.

        om = O.gens()[1]

        # M = (a, b; c, d) we need norm(cz+d)<1
        # and we want such c and d only up to units.
        # We have 1 > norm(|cz+d|) >= |norm(c)|*norm(im(z)),
        # so |norm(c)| < norm(im(z))^-1.
        # We can find all such c up to units by enumerating ideals
        # and finding generators.
        # Given c, we want a coprime d with norm(|cx+d|)^2 minimal.
        x = [z[0].real(), z[1].real()]
        y = [z[0].imag(), z[1].imag()]
        norm_im_z = self._RR(abs(prod(y)))
        bound = ZZ((norm_im_z**-1).floor())
#        if bound_on_c is None:
#            bound_to_try = 1
#            while bound_to_try < 2*bound:
#                m = self._up_cd(z, bound_on_c = bound_to_try)
#                if m != None:
#                    return m
#                bound_to_try = bound_to_try * 2
#            return None
        if (not bound_on_c is None) and bound > bound_on_c:
            bound = bound_on_c
        ids = K.ideals_of_bdd_norm(bound)
        for n in ids:
            for C in ids[n]:
                if C.is_principal():
                    eps = K(K.unit_group().gens()[1])
                    c = C.gens_reduced()[0]
                    c_times_y = [embs[0](c)*y[0], embs[1](c)*y[1]]
                    e = ((self._RR(c_times_y[0]/c_times_y[1]).abs().log()) / 
                         (2*(self._RR(embs[0](eps)).abs().log())))
#                    print c, e
                    c = c / eps ** (e.round())
#                    print c
                        
                    # Next, we try tofind d with |norm(cz+d)|<1.
                    # We start by moving to an additive fundamental domain,
                    # and do a and thorough search if that does not work
                    # (which it usually does not).
                    c_times_x = [embs[0](c)*x[0], embs[1](c)*x[1]]
                    d0, xp = self.to_additive(c_times_x)
                    d0 = -d0
                    # now xp = c*x + d_0 is in the fundamental domain
                    cz_plus_d = [embs[0](c)*z[0] + embs[0](d0), embs[1](c)*z[1] + embs[1](d0)]
                    # Actually, the following check already uses class number != 1, because
                    # otherwise the check is supposed to be < Norm(c*O+d0*O) or something.

                    if abs(prod(cz_plus_d)) < 1:
                        if c*O+d0*O == 1*O:
                           return c, d0
                        # In the class number > 1 case, we may want to do something
                        # if the test c*O+d0*O == 1*O fails.
                                
                    # We have 1 > norm(|cz+d|^2) =

                    #            |norm(c)|^2*norm(im(z))^2
                    #          + |c_1|^2*im(z_1)^2 *|c_0x_0+d_0|^2
                    #          + |c_0|^2*im(z_0)^2 *|c_1x_1+d_1|^2
                    #          +  norm(|cx+d|)^2,
                    # define x_0' = c_0x_0+d_0
                    # and    x_1' = c_1x_1+d_1, then they have the following
                    # absolute value bounds:
                    #
                    bound_on_x0p = ((1 - n**2 * norm_im_z**2)
                                     / (embs[1](c)**2 * y[1]**2)).sqrt()
                    bound_on_x1p = ((1 - n**2 * norm_im_z**2)
                                     / (embs[0](c)**2 * y[0]**2)).sqrt()
                    d0_interval = [-embs[0](c) * x[0] - bound_on_x0p,
                                   -embs[0](c) * x[0] + bound_on_x0p]
                    d1_interval = [-embs[1](c) * x[1] - bound_on_x1p,
                                   -embs[1](c) * x[1] + bound_on_x1p]
                    # Now write d = e + f*om = (g + f*delta)/2,
                    # where g = 2*e + f*trace(om),
                    # and om = (trace(om) + delta) / 2.
                    delta = 2*om - om.trace()
                    assert embs[0](delta) < 0
                    assert embs[1](delta) > 0
                    delta_num = embs[1](delta)
                    # Then g = d0 + d1 is in the following interval:
                    g_interval = [d0_interval[i] + d1_interval[i]
                                  for i in [0,1]]
                    g_start = sum([-embs[0](c) * x[i] for i in [0,1]]).round()
                    g_bound = (bound_on_x0p + bound_on_x1p + 3).floor()
                    t = om.trace()
                    for i in xrange(g_bound):
                        # TODO: this range can be way too big, something needs
                        #       to be done about that.
                        if t == 0 or not ((g_start + i) % 2):
                            for s in ([-1,1] if i else [1]):
                                g = g_start + s * i
                                # now f = (-g + 2*d)/delta, and we know
                                # bounds on each embedding:
                                f0_lower = (-g + 2*d0_interval[1]) / -delta_num
                                f1_lower = (-g + 2*d1_interval[0]) /  delta_num
                                f0_upper = (-g + 2*d0_interval[0]) / -delta_num
                                f1_upper = (-g + 2*d1_interval[1]) /  delta_num
                                f_lower = max([f0_lower, f1_lower]).ceil()
                                f_upper = min([f0_upper, f1_upper]).floor()
                                for f in range(f_lower, f_upper+1):
                                    if not ((g - f*t) % 2):
                                        d = (g+f*delta)/2
                                        assert d in O
                                        cz_plus_d = [embs[0](c)*z[0] +
                                                     embs[0](d),
                                                     embs[1](c)*z[1] +
                                                     embs[1](d)]
                                        if abs(prod(cz_plus_d)) < 1:
                                            if c*O+d*O == 1*O:
                                                return c, d
                                            # In the class number > 1 case,
                                            # we may want to do something
                                            # if the test c*O+d0*O == 1*O
                                            # fails. For now, the code just
                                            # skips c,d in that case.
                                        # 
                    # # if we write d = e + f*om, then d0+d1 = 2*e + f*tr(om)
                    # # is in the following interval:
                    # 
                    # d0-d1 = f*(om_0 - om_1) is in the following interval
                    # interval = [d0_interval[0]-d1_interval[1], d0_interval[1]-d1_interval[0]]
                    # f_interval = [interval[i] / (embs[0](om)-embs[1](om)) for i in [0,1]]
                    # f_interval.sort()
                    # for f in range(int(self._RR(f_interval[0]).ceil()), int(self._RR(f_interval[1]).floor()+1)):
                        # e_lower = self._RR(max([d0_interval[0]-f*embs[0](om), d1_interval[0]-f*embs[1](om)])).ceil()
                        # e_upper = self._RR(min([d0_interval[1]-f*embs[0](om), d1_interval[1]-f*embs[1](om)])).floor()
                        # for e in range(int(e_lower), int(e_upper+1)):
                            # d = e + f*om
                            # cz_plus_d = [embs[0](c)*z[0] + embs[0](d), embs[1](c)*z[1] + embs[1](d)]
                            # if abs(prod(cz_plus_d)) < 1:
                                # if c*O+d*O == 1*O:
                                    # return c, d
                                # # In the class number > 1 case, we may want to do something
                                # # if the test c*O+d0*O == 1*O fails.


class _HilbertFundDomData_classical(_HilbertFundDomData_generic):
    """
    As _HilbertFundDomData, but only for the classical upper half
    plane and the action of SL_2(ZZ).
    """

    _eps = None

    def __init__(self, embs=None):
        """
        Create the object.
        """
        self._K = QQ

    def to_additive(self, x):
        """
        Move x in RR^1 to an additive fundamental domain.
        """
        a = x[0].round()
        x = [x[0]-a]
        return a, x

    def to_multiplicative(self, y):
        """
        Returns (identitymatrix(2), y).
        """
        return Matrix([[1,0],[0,1]]), y


    def up(self, z):
        """
        Move z up in HH, i.e., increase its imaginary part,
        if possible by the transformation z -> -1/z

        Returns either 1 or a matrix M in SL_2(ZZ) such that
        Mz has a larger imaginary part.

        Assumes the real part of x is between -1/2 and 1/2
        """
        if z[0].abs() < 1:
            return Matrix(QQ, [[0,-1],[1,0]])
        return Matrix(QQ, [[1,0],[0,1]])

    def step(self, z):
        """
        Do one reduction round.

        Returns a pair (M, z') with z' = M(z) nearer the fundamental domain if possible.
        """
        M1 = self.up(z)
        z = self.act(M1, z)
        t = self.to_additive([e.real() for e in z])[0]
        M3 = Matrix([[1, -t], [0, 1]])
        z = self.act(M3, z)
        return M3*M1, z


    def act(self, M, z):
        """
        Return M(z)
        """
        a = M[0,0]
        b = M[0,1]
        c = M[1,0]
        d = M[1,1]
        z = z[0]
        return [(a*z+b)/(c*z+d)]


