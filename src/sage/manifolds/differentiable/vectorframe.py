r"""
Vector frames

The class :class:`VectorFrame` implements vector frames on differentiable
manifolds over `\RR`.
By *vector frame*, it is meant a field `e` on some open domain `U` of a
manifold `S` endowed with a differentiable mapping `\Phi: U\rightarrow V` to a
parallelizable domain `V` of a manifold `M` such that for each `p\in U`,
`e(p)` is a vector basis of the tangent space `T_{\Phi(p)}M`.

The standard case of a vector frame *on* `U` corresponds to `S=M`, `U=V`
and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
immersion and `\Phi` being a curve in `V` (`U` is then an open interval
of `\RR`).

A derived class of :class:`VectorFrame` is :class:`CoordFrame`; it regards the
vector frames associated with a chart, i.e. the so-called *coordinate bases*.

The vector frame duals, i.e. the coframes, are implemented via the class
:class:`CoFrame`. The derived class :class:`CoordCoFrame` is devoted to
coframes deriving from a chart.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

REFERENCES:

- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)

EXAMPLES:

    Setting a vector frame on a 3-dimensional manifold::

        sage: M = DiffManifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e') ; e
        vector frame (M, (e_0,e_1,e_2))
        sage: latex(e)
        \left(M, \left(e_0,e_1,e_2\right)\right)

    The first frame defined on a manifold is its default frame; in the present
    case it is the coordinate frame defined when introducing the chart c_xyz::

        sage: M.default_frame()
        coordinate frame (M, (d/dx,d/dy,d/dz))

    The default frame can be changed via the method
    :meth:`~sage.manifolds.differentiable.manifold.DiffManifold.set_default_frame`::

        sage: M.set_default_frame(e)
        sage: M.default_frame()
        vector frame (M, (e_0,e_1,e_2))

    The elements of a vector frame are vector fields on the manifold::

        sage: e._vec
        (vector field 'e_0' on the 3-dimensional manifold 'M', vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M')

    Each element can be accessed by its index::

        sage: e[0]
        vector field 'e_0' on the 3-dimensional manifold 'M'

    The index range depends on the starting index defined on the manifold::

        sage: M = DiffManifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e')
        sage: e._vec
        (vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M', vector field 'e_3' on the 3-dimensional manifold 'M')
        sage: e[1], e[2], e[3]
        (vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M', vector field 'e_3' on the 3-dimensional manifold 'M')

    Let us check that the vector fields e(i) are the frame vectors from
    their components w.r.t. to the frame e::

        sage: e[1].comp(e)[:]
        [1, 0, 0]
        sage: e[2].comp(e)[:]
        [0, 1, 0]
        sage: e[3].comp(e)[:]
        [0, 0, 1]

    Defining a vector frame on a manifold automatically creates the dual
    coframe, which bares the same name (here e)::

        sage: M.coframes()
        [coordinate coframe (M, (dx,dy,dz)), coframe (M, (e^1,e^2,e^3))]
        sage: f = M.coframes()[1] ; f
        coframe (M, (e^1,e^2,e^3))

    Each element of the coframe is a 1-form::

        sage: f[1], f[2], f[3]
        (1-form 'e^1' on the 3-dimensional manifold 'M',
        1-form 'e^2' on the 3-dimensional manifold 'M',
        1-form 'e^3' on the 3-dimensional manifold 'M')
        sage: latex(f[1]), latex(f[2]), latex(f[3])
        (e^1, e^2, e^3)

    Let us check that the coframe (e^i) is indeed the dual of the vector
    frame (e_i)::

        sage: f[1](e[1]) # the 1-form e^1 applied to the vector field e_1
        scalar field 'e^1(e_1)' on the 3-dimensional manifold 'M'
        sage: f[1](e[1]).expr() # the explicit expression of e^1(e_1)
        1
        sage: f[1](e[1]).expr(), f[1](e[2]).expr(), f[1](e[3]).expr()
        (1, 0, 0)
        sage: f[2](e[1]).expr(), f[2](e[2]).expr(), f[2](e[3]).expr()
        (0, 1, 0)
        sage: f[3](e[1]).expr(), f[3](e[2]).expr(), f[3](e[3]).expr()
        (0, 0, 1)

    The coordinate frame associated to spherical coordinates of the
    sphere `S^2`::

        sage: M = DiffManifold(2, 'S^2', start_index=1)
        sage: c_spher.<th,ph> = M.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: b = M.default_frame() ; b
        coordinate frame (S^2, (d/dth,d/dph))
        sage: b[1]
        vector field 'd/dth' on the 2-dimensional manifold 'S^2'
        sage: b[2]
        vector field 'd/dph' on the 2-dimensional manifold 'S^2'

    The orthonormal frame constructed from the coordinate frame::

        sage: change_frame = M.automorphism_field()
        sage: change_frame[:] = [[1,0], [0, 1/sin(th)]]
        sage: e = b.new_frame(change_frame, 'e') ; e
        vector frame (S^2, (e_1,e_2))
        sage: e[1][:]
        [1, 0]
        sage: e[2][:]
        [0, 1/sin(th)]

    The change-of-frame matrices::

        sage: M.change_of_frame(c_spher.frame(), e)
        field of tangent-space automorphisms on the 2-dimensional manifold 'S^2'
        sage: M.change_of_frame(c_spher.frame(), e)[:]
        [        1         0]
        [        0 1/sin(th)]
        sage: M.change_of_frame(e, c_spher.frame())
        field of tangent-space automorphisms on the 2-dimensional manifold 'S^2'
        sage: M.change_of_frame(e, c_spher.frame())[:]
        [      1       0]
        [      0 sin(th)]

"""

#******************************************************************************
#       Copyright (C) 2013, 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_basis import FreeModuleBasis, \
                                                              FreeModuleCoBasis
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule

class VectorFrame(FreeModuleBasis):
    r"""
    Vector frame on a differentiable manifold over `\RR`.

    By *vector frame*, it is meant a field `e` on some open domain `U` of a
    manifold `S` endowed with a differentiable mapping `\Phi: U\rightarrow V`
    to a parallelizable domain `V` of a manifold `M` such that for each
    `p\in U`, `e(p)` is a vector basis of the tangent space `T_{\Phi(p)}M`.

    The standard case of a vector frame *on* `U` corresponds to `S=M`, `U=V`
    and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `V` (`U` is then an open interval
    of `\RR`).

    For each instanciation of a vector frame, a coframe is automatically
    created, as an instance of the class :class:`CoFrame`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U\subset S` with values on `\Phi(U)\subset V \subset M`
    - ``symbol`` -- a letter (of a few letters) to denote a
      generic vector of the frame; can be set to None if the parameter
      ``from_frame`` is filled.
    - ``latex_symbol`` -- (default: ``None``) symbol to denote a generic vector of
      the frame; if ``None``, the value of ``symbol`` is used.
    - ``from_frame`` -- (default: ``None``) vector frame `\tilde e` on the codomain
      `V` of the destination map `\Phi`; the frame `e` = ``self`` is then
      constructed so that `\forall p \in U, e(p) = \tilde e(\Phi(p))`


    EXAMPLES:

    Setting a vector frame on a 3-dimensional manifold::

        sage: DiffManifold._clear_cache_() # for doctests only
        sage: M = DiffManifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e') ; e
        vector frame (M, (e_0,e_1,e_2))
        sage: latex(e)
        \left(M, \left(e_0,e_1,e_2\right)\right)

    The LaTeX symbol can be specified::

        sage: e = M.vector_frame('E', r"\epsilon")
        sage: latex(e)
        \left(M, \left(\epsilon_0,\epsilon_1,\epsilon_2\right)\right)

    Example with a non-trivial mapping `\Phi`: vector frame along a curve::

        sage: R.<t> = RealLine()
        sage: U = R.open_interval(-1, 1)
        sage: Phi = U.diff_map(M, [cos(t), sin(t), t], name='Phi', latex_name=r'\Phi')
        sage: Phi
        Curve 'Phi' in the 3-dimensional manifold 'M'
        sage: f = U.vector_frame('f', dest_map=Phi) ; f
        vector frame ((-1, 1), (f_0,f_1,f_2)) with values on the 3-dimensional
         manifold 'M'
        sage: f.domain()
        Real interval (-1, 1)
        sage: p = U(0, name='p') ; p
        point 'p' on field R of real numbers
        sage: f.at(p)
        Basis (f_0,f_1,f_2) on the tangent space at point 'Phi(p)' on
         3-dimensional manifold 'M'

    """
    def __init__(self, vector_field_module, symbol, latex_symbol=None,
                 from_frame=None):
        from sage.manifolds.differentiable.manifold import DiffManifold
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._dest_map = vector_field_module._dest_map
        self._from_frame = from_frame
        self._manifold = self._domain._manifold
        if symbol is None:
            if from_frame is None:
                raise TypeError("Some frame symbol must be provided.")
            symbol = 'X'  # provisory symbol
        FreeModuleBasis.__init__(self, vector_field_module,
                                 symbol, latex_symbol=latex_symbol)
        # Redefinition of the name and the LaTeX name:
        if from_frame is None:
            self._name = "(" + self._domain._name + ", " + self._name + ")"
            self._latex_name = r"\left(" + self._domain._latex_name + ", " + \
                          self._latex_name + r"\right)"
        else:
            if not from_frame._domain.is_subset(self._dest_map._codomain):
                raise ValueError("The domain of the frame 'from_frame' is " +
                                 "not included in the codomain of the " +
                                 "destination map.")
            n = self._fmodule.rank()
            for i in range(n):
                self._vec[i]._name = from_frame._vec[i]._name
                self._vec[i]._latex_name = from_frame._vec[i]._latex_name
            self._name = "(" + self._domain._name + ", (" + \
                        ",".join([self._vec[i]._name for i in range(n)]) + "))"
            self._latex_name = r"\left(" + self._domain._latex_name + \
                        r" ,\left(" + \
                        ",".join([self._vec[i]._latex_name for i in range(n)])+ \
                        r"\right)\right)"
            self._symbol = from_frame._symbol
            self._latex_symbol = from_frame._latex_symbol
            # Names of the dual coframe:
            self_dual = self.dual_basis()
            from_dual = from_frame.dual_basis()
            for i in range(n):
                self_dual._form[i]._name = from_dual._form[i]._name
                self_dual._form[i]._latex_name = from_dual._form[i]._latex_name
            self_dual._name = "(" + self._domain._name + ", (" + \
                  ",".join([self_dual._form[i]._name for i in range(n)]) + "))"
            self_dual._latex_name = r"\left(" + self._domain._latex_name + \
                r" ,\left(" + \
                ",".join([self_dual._form[i]._latex_name for i in range(n)])+ \
                r"\right)\right)"
        # The frame is added to the domain's set of frames, as well as to all
        # the superdomains' sets of frames; moreover the first defined frame
        # is considered as the default one
        dest_map = self._dest_map
        for sd in self._domain._supersets:
            for other in sd._frames:
                if repr(self) == repr(other):
                    raise ValueError("the " + str(self) + " already exists " +
                                     "on the " + str(sd))
            sd._frames.append(self)
            sd._top_frames.append(self)
            if sd._def_frame is None:
                sd._def_frame = self
            if isinstance(sd, DiffManifold):
                # Initialization of the zero elements of tensor field modules:
                if dest_map in sd._vector_field_modules:
                    xsd = sd._vector_field_modules[dest_map]
                    if not isinstance(xsd, FiniteRankFreeModule):
                        for t in xsd._tensor_modules.itervalues():
                            t(0).add_comp(self)
                            # (since new components are initialized to zero)
        if dest_map is self._domain._identity_map:
            # The frame is added to the list of the domain's covering frames:
            self._domain._set_covering_frame(self)
        #
        # Dual coframe
        self._coframe = self.dual_basis()  # self._coframe = a shortcut for
                                           # self._dual_basis
        #
        # Derived quantities:
        self._structure_coef = None
        # Initialization of the set of frames that are restrictions of the
        # current frame to subdomains of the frame domain:
        self._subframes = set([self])
        # Initialization of the set of frames which the current frame is a
        # restriction of:
        self._superframes = set([self])
        #
        self._restrictions = {} # dict. of the restrictions of self to
                               # subdomains of self._domain, with the
                               # subdomains as keys
        # NB: set(self._restrictions.itervalues()) is identical to self._subframes


    ###### Methods that must be redefined by derived classes of FreeModuleBasis ######

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "vector frame " + self._name
        if self._dest_map is not self._domain._identity_map:
            description += " with values on the " + str(self._dest_map._codomain)
        return description


    def _init_dual_basis(self):
        r"""
        Construct the basis dual to ``self``.

        OUTPUT:

        - instance of :class:`CoFrame` representing the dual of
          ``self``

        """
        return CoFrame(self, self._symbol, latex_symbol=self._latex_symbol)

    def _new_instance(self, symbol, latex_symbol=None):
        r"""
        Construct a new vector frame on the same vector field module
        as ``self``.

        INPUT:

        - ``symbol`` -- (string) a letter (of a few letters) to denote a
          generic element of the vector frame
        - ``latex_symbol`` -- (string; default: ``None``) symbol to denote a
          generic element of the vector frame; if ``None``, the value of
          ``symbol`` is used.

        OUTPUT:

        - instance of :class:`VectorFrame`

        """
        return VectorFrame(self._fmodule, symbol, latex_symbol=latex_symbol)

    ###### End of methods redefined by derived classes ######

    def domain(self):
        r"""
        Return the domain on which ``self`` is defined.
        """
        return self._domain

    def coframe(self):
        r"""
        Return the dual coframe.
        """
        return self._coframe

    def new_frame(self, change_of_frame, symbol, latex_symbol=None):
        r"""
        Define a new vector frame from the current one.

        The new vector frame is defined on the same domain as ``self`` from
        a field of automorphisms.

        INPUT:

        - ``change_of_frame`` -- instance of
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
          describing the automorphism `P` that relates the current frame
          `(e_i)` (described by ``self``) to the new frame `(n_i)` according
          to `n_i = P(e_i)`
        - ``symbol`` -- a letter (of a few letters) to denote a generic vector
          of the frame
        - ``latex_symbol`` -- (default: ``None``) symbol to denote a generic vector
          of the frame; if ``None``, the value of ``symbol`` is used.

        OUTPUT:

        - the new frame `(n_i)`, as an instance of :class:`VectorFrame`

        EXAMPLES:

        Frame resulting from a pi/3-rotation in the Euclidean plane::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2,'R^2')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e') ; M.set_default_frame(e)
            sage: M._frame_changes
            {}
            sage: rot = M.automorphism_field()
            sage: rot[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: n = e.new_frame(rot, 'n')
            sage: n[0][:]
            [1/2*sqrt(3), 1/2]
            sage: n[1][:]
            [-1/2, 1/2*sqrt(3)]
            sage: a =  M.change_of_frame(e,n)
            sage: a[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a == rot
            True
            sage: a is rot
            False
            sage: a._components # random (dictionary output)
            {vector frame (R^2, (e_0,e_1)): 2-indices components w.r.t. vector frame (R^2, (e_0,e_1)),
            vector frame (R^2, (n_0,n_1)): 2-indices components w.r.t. vector frame (R^2, (n_0,n_1))}
            sage: a.comp(n)[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a1 = M.change_of_frame(n,e)
            sage: a1[:]
            [1/2*sqrt(3)         1/2]
            [       -1/2 1/2*sqrt(3)]
            sage: a1 == rot.inverse()
            True
            sage: a1 is rot.inverse()
            False
            sage: e[0].comp(n)[:]
            [1/2*sqrt(3), -1/2]
            sage: e[1].comp(n)[:]
            [1/2, 1/2*sqrt(3)]

        """
        the_new_frame = self.new_basis(change_of_frame, symbol,
                                       latex_symbol=latex_symbol)
        for sdom in self._domain._supersets:
            sdom._frame_changes[(self, the_new_frame)] = \
                            self._fmodule._basis_changes[(self, the_new_frame)]
            sdom._frame_changes[(the_new_frame, self)] = \
                            self._fmodule._basis_changes[(the_new_frame, self)]
        return the_new_frame

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some subdomain of ``self._domain``.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- subdomain `V` of the current frame domain `U`

        OUTPUT:

        - the restriction of ``self`` to `V`, as an instance of
          :class:`VectorFrame`.

        EXAMPLE:

        Restriction of a frame defined on `\RR^2` to the unit disk::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: a = M.automorphism_field()
            sage: a[:] = [[1-y^2,0], [1+x^2, 2]]
            sage: e = c_cart.frame().new_frame(a, 'e') ; e
            vector frame (R^2, (e_0,e_1))
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1})
            sage: e_U = e.restrict(U) ; e_U
            vector frame (U, (e_0,e_1))

        The vectors of the restriction have the same symbols as those of the
        original frame::

            sage: e_U[0].display()
            e_0 = (-y^2 + 1) d/dx + (x^2 + 1) d/dy
            sage: e_U[1].display()
            e_1 = 2 d/dy

        They are actually the restrictions of the original frame vectors::

            sage: e_U[0] is e[0].restrict(U)
            True
            sage: e_U[1] is e[1].restrict(U)
            True

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("The provided domain is not a subdomain of " +
                                 "the current frame's domain.")
            sdest_map = self._dest_map.restrict(subdomain)
            res = VectorFrame(subdomain.vector_field_module(sdest_map,
                                                            force_free=True),
                              self._symbol, latex_symbol=self._latex_symbol)
            for dom in subdomain._supersets:
                if dom is not subdomain:
                    dom._top_frames.remove(res)  # since it was added by
                                                 # VectorFrame constructor
            n = self._fmodule.rank()
            new_vectors = list()
            for i in range(n):
                vrest = self._vec[i].restrict(subdomain)
                for j in self._fmodule.irange():
                    vrest.add_comp(res)[j] = 0
                vrest.add_comp(res)[i] = 1
                new_vectors.append(vrest)
            res._vec = tuple(new_vectors)
            # Update of superframes and subframes:
            res._superframes.update(self._superframes)
            for sframe in self._superframes:
                sframe._subframes.add(res)
                sframe._restrictions[subdomain] = res # includes sframe = self
        return self._restrictions[subdomain]

    def structure_coef(self):
        r"""
        Evaluate the structure coefficients associated to the vector frame.

        `n` being the manifold's dimension, the structure coefficients of the
        vector frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}`
        defined by

        .. MATH::

            [e_i, e_j] = C^k_{\ \, ij} e_k

        OUPUT:

        - the structure coefficients `C^k_{\ \, ij}`, as an instance of
          :class:`~sage.tensor.modules.comp.CompWithSym`
          with 3 indices ordered as `(k,i,j)`.

        EXAMPLE:

        Structure coefficients of the orthonormal frame associated to
        spherical coordinates in the Euclidean space `\RR^3`::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(3, 'R^3', '\RR^3', start_index=1)
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: ch_frame = M.automorphism_field()
            sage: ch_frame[1,1], ch_frame[2,2], ch_frame[3,3] = 1, 1/r, 1/(r*sin(th))
            sage: M.frames()
            [coordinate frame (R^3, (d/dr,d/dth,d/dph))]
            sage: e = c_spher.frame().new_frame(ch_frame, 'e')
            sage: e[1][:]  # components of e_1 in the manifold's default frame (d/dr, d/dth, d/dth)
            [1, 0, 0]
            sage: e[2][:]
            [0, 1/r, 0]
            sage: e[3][:]
            [0, 0, 1/(r*sin(th))]
            sage: c = e.structure_coef() ; c
            3-indices components w.r.t. vector frame (R^3, (e_1,e_2,e_3)), with antisymmetry on the index positions (1, 2)
            sage: c[:]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, -1/r, 0], [1/r, 0, 0], [0, 0, 0]],
             [[0, 0, -1/r], [0, 0, -cos(th)/(r*sin(th))], [1/r, cos(th)/(r*sin(th)), 0]]]
            sage: c[2,1,2]  # C^2_{12}
            -1/r
            sage: c[3,1,3]  # C^3_{13}
            -1/r
            sage: c[3,2,3]  # C^3_{23}
            -cos(th)/(r*sin(th))

        """
        from sage.tensor.modules.comp import CompWithSym
        if self._structure_coef is None:
            fmodule = self._fmodule
            self._structure_coef = CompWithSym(self._fmodule._ring, self, 3,
                                                   start_index=fmodule._sindex,
                                    output_formatter=fmodule._output_formatter,
                                                                 antisym=(1,2))
            si = fmodule._sindex
            nsi = si + fmodule.rank()
            for k in range(si,nsi):
                ce_k = self._coframe._form[k-si]
                for i in range(si, nsi):
                    e_i = self._vec[i-si]
                    for j in range(i+1, nsi):
                        e_j = self._vec[j-si]
                        self._structure_coef[[k,i,j]] = ce_k(e_j.lie_der(e_i))
        return self._structure_coef

    def at(self, point):
        r"""
        Return the value of the frame at a given point, i.e. a basis of the
        tangent vector space.

        INPUT:

        - ``point`` -- (instance of
          :class:`~sage.manifolds.differentiable.point.ManifoldPoint`) point `p` in
          the domain `U` of ``self`` (denoted `e` hereafter)

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the basis `e(p)` of the tangent vector space
          `T_{\Phi(p)} M`, where `\Phi: U \rightarrow V\subset M` is
          the differentiable mapping associated with `e` (possibly
          `\Phi = \mathrm{Id}_U`)

        EXAMPLES:

        Basis of a tangent space to a 2-dimensional manifold::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((-1,2), name='p')
            sage: e = X.frame() ; e
            coordinate frame (M, (d/dx,d/dy))
            sage: ep = e.at(p) ; ep
            Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: type(ep)
            <class 'sage.tensor.modules.free_module_basis.FreeModuleBasis'>
            sage: ep[0]
            tangent vector d/dx at point 'p' on 2-dimensional manifold 'M'
            sage: ep[1]
            tangent vector d/dy at point 'p' on 2-dimensional manifold 'M'

        Note that the symbols used to denote the vectors are same as those
        for the vector fields of the frame. At this stage, ep is the unique
        basis on the tangent space at p::

            sage: Tp = p.tangent_space()
            sage: Tp.bases()
            [Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M']

        Let us consider a vector frame that is a not a coordinate one::

            sage: aut = M.automorphism_field()
            sage: aut[:] = [[1+y^2, 0], [0, 2]]
            sage: f = e.new_frame(aut, 'f') ; f
            vector frame (M, (f_0,f_1))
            sage: fp = f.at(p) ; fp
            Basis (f_0,f_1) on the tangent space at point 'p' on 2-dimensional manifold 'M'

        There are now two bases on the tangent space::

            sage: Tp.bases()
            [Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M',
             Basis (f_0,f_1) on the tangent space at point 'p' on 2-dimensional manifold 'M']

        Moreover, the changes of bases in the tangent space have been computed
        from the known relation between the frames e and f (automorphism field
        aut defined above)::

            sage: Tp.change_of_basis(ep, fp)
            Automorphism of the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: Tp.change_of_basis(ep, fp).display()
            5 d/dx*dx + 2 d/dy*dy
            sage: Tp.change_of_basis(fp, ep)
            Automorphism of the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: Tp.change_of_basis(fp, ep).display()
            1/5 d/dx*dx + 1/2 d/dy*dy

        The dual bases::

            sage: e.coframe()
            coordinate coframe (M, (dx,dy))
            sage: ep.dual_basis()
            Dual basis (dx,dy) on the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: ep.dual_basis() is e.coframe().at(p)
            True
            sage: f.coframe()
            coframe (M, (f^0,f^1))
            sage: fp.dual_basis()
            Dual basis (f^0,f^1) on the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: fp.dual_basis() is f.coframe().at(p)
            True

        """
        # Case of a non-trivial destination map
        if self._from_frame is not None:
            if self._dest_map.is_identity():  #!# probably not necessary
                raise ValueError("the destination map should not be the " +
                                 "identity")
            ambient_point = self._dest_map(point)
            return self._from_frame.at(ambient_point)
        # Determination of the tangent space:
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))
        if self._dest_map.is_identity():
            ambient_point = point
        else:
            ambient_point = self._dest_map(point)
        ts = ambient_point.tangent_space()
        # If the basis has already been constructed, it is simply returned:
        ts_frame_bases = ts._frame_bases
        if self in ts_frame_bases:
            return ts_frame_bases[self]
        for frame in ts_frame_bases:
            if self in frame._subframes or self in frame._superframes:
                return ts_frame_bases[frame]
        # If this point is reached, the basis has to be constructed from
        # scratch:
        basis = ts.basis(symbol=self._symbol, latex_symbol=self._latex_symbol)
        # Names of basis vectors set to those of the frame vector fields:
        n = ts.dim()
        for i in range(n):
            basis._vec[i]._name = self._vec[i]._name
            basis._vec[i]._latex_name = self._vec[i]._latex_name
        basis._name = "(" + \
                ",".join([basis._vec[i]._name for i in range(n)]) + ")"
        basis._latex_name = r"\left(" + \
             ",".join([basis._vec[i]._latex_name for i in range(n)])+ \
             r"\right)"
        basis._symbol = basis._name
        basis._latex_symbol = basis._latex_name
        # Names of cobasis linear forms set to those of the coframe
        # 1-forms:
        coframe = self.coframe()
        cobasis = basis.dual_basis()
        for i in range(n):
            cobasis._form[i]._name = coframe._form[i]._name
            cobasis._form[i]._latex_name = coframe._form[i]._latex_name
        cobasis._name = "(" + \
             ",".join([cobasis._form[i]._name for i in range(n)]) + ")"
        cobasis._latex_name = r"\left(" + \
          ",".join([cobasis._form[i]._latex_name for i in range(n)])+ \
          r"\right)"
        ts_frame_bases[self] = basis
        # Update of the change of bases in the tangent space:
        for frame_pair, automorph in self._domain._frame_changes.iteritems():
            frame1 = frame_pair[0] ; frame2 = frame_pair[1]
            if frame1 is self:
                fr2 = None
                for frame in ts_frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr2 is not None:
                    basis1 = basis
                    basis2 = ts_frame_bases[fr2]
                    auto = ts.automorphism()
                    for frame, comp in automorph._components.iteritems():
                        bas = None
                        if frame is frame1:
                            bas = basis1
                        if frame is frame2:
                            bas = basis2
                        if bas is not None:
                            cauto = auto.add_comp(bas)
                            for ind, val in comp._comp.iteritems():
                                cauto._comp[ind] = val(point)
                    ts._basis_changes[(basis1, basis2)] = auto
            if frame2 is self:
                fr1 = None
                for frame in ts_frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                if fr1 is not None:
                    basis1 = ts_frame_bases[fr1]
                    basis2 = basis
                    auto = ts.automorphism()
                    for frame, comp in automorph._components.iteritems():
                        bas = None
                        if frame is frame1:
                            bas = basis1
                        if frame is frame2:
                            bas = basis2
                        if bas is not None:
                            cauto = auto.add_comp(bas)
                            for ind, val in comp._comp.iteritems():
                                cauto._comp[ind] = val(point)
                    ts._basis_changes[(basis1, basis2)] = auto
        return basis

    def along(self, mapping):
        r"""
        Return the vector frame deduced from ``self`` via a mapping, the
        codomain of which is included in the domain of ``self``.

        If ``self`` is the vector frame `e` on the open set `V` and if
        `\Phi: U \rightarrow V` is a differentiable mapping from an open
        set `U` of some manifold `M` to `V`, the returned object is
        a vector frame `\tilde e` along `U` with values on `V` such that

        .. MATH::

           \forall p \in U,\  \tilde e(p) = e(\Phi(p)).

        INPUT:

        - ``mapping`` -- differentiable mapping `\Phi: U \rightarrow V`

        OUTPUT:

        - vector frame `\tilde e` along `U` defined above.

        EXAMPLE:

        Vector frame along a curve::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: Phi = R.diff_map(M, {X: (cos(t), t)}, name='Phi',
            ....:                  latex_name=r'\Phi') ; Phi
            Curve 'Phi' in the 2-dimensional manifold 'M'
            sage: e = X.frame() ; e
            coordinate frame (M, (d/dx,d/dy))
            sage: te = e.along(Phi) ; te
            vector frame (R, (d/dx,d/dy)) with values on the 2-dimensional manifold 'M'


        Check of the formula `\tilde e(p) = e(\Phi(p))`::

            sage: p = R(pi) ; p
            point on field R of real numbers
            sage: te[0].at(p) == e[0].at(Phi(p))
            True
            sage: te[1].at(p) == e[1].at(Phi(p))
            True

        The result is cached::

            sage: te is e.along(Phi)
            True

        """
        dom = self._domain
        if mapping.codomain().is_subset(dom):
            rmapping = mapping
        else:
            rmapping = None
            for doms, rest in mapping._restrictions.iteritems():
                if doms[1].is_subset(dom):
                    rmapping = rest
                    break
            else:
                raise ValueError("the codomain of {} is not ".format(mapping)
                             + "included in the domain of {}".format(self))
        vmodule = rmapping.domain().vector_field_module(dest_map=rmapping)
        return vmodule.basis(from_frame=self)


#******************************************************************************

class CoordFrame(VectorFrame):
    r"""
    Coordinate frame on a differentiable manifold over `\RR`.

    By *coordinate frame*, it is meant a vector frame on a manifold `M` that
    is associated to a coordinate system (chart) on `M`.

    INPUT:

    - ``chart`` -- the chart defining the coordinates

    EXAMPLES:

    The coordinate frame associated to spherical coordinates of the
    sphere `S^2`::

        sage: DiffManifold._clear_cache_() # for doctests only
        sage: M = DiffManifold(2, 'S^2', start_index=1)
        sage: M.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
        chart (S^2, (th, ph))
        sage: b = M.default_frame()
        sage: b
        coordinate frame (S^2, (d/dth,d/dph))
        sage: b[1]
        vector field 'd/dth' on the 2-dimensional manifold 'S^2'
        sage: b[2]
        vector field 'd/dph' on the 2-dimensional manifold 'S^2'
        sage: latex(b)
        \left(S^2 ,\left(\frac{\partial}{\partial {\theta} },\frac{\partial}{\partial {\phi} }\right)\right)

    """
    def __init__(self, chart):
        from sage.misc.latex import latex
        from sage.manifolds.differentiable.chart import DiffChart
        if not isinstance(chart, DiffChart):
            raise TypeError("the first argument must be a chart")
        self._chart = chart
        VectorFrame.__init__(self,
                             chart._domain.vector_field_module(force_free=True),
                             symbol='X')
        # In the above:
        # - force_free=True ensures that a free module is constructed in case
        #   it is the first call to the vector field module on chart._domain
        # - 'X' is a provisory symbol
        n = self._manifold._dim
        for i in range(n):
            self._vec[i]._name = "d/d" + str(self._chart._xx[i])
            self._vec[i]._latex_name = r"\frac{\partial}{\partial" + \
                                     latex(self._chart._xx[i]) + r"}"
        self._name = "(" + self._domain._name + ", (" + \
                    ",".join([self._vec[i]._name for i in range(n)]) + "))"
        self._latex_name = r"\left(" + self._domain._latex_name + r" ,\left(" + \
                       ",".join([self._vec[i]._latex_name for i in range(n)])+ \
                       r"\right)\right)"
        self._symbol = self._name
        self._latex_symbol = self._latex_name


    ###### Methods that must be redefined by derived classes of FreeModuleBasis ######

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "coordinate frame " + self._name

    def _init_dual_basis(self):
        r"""
        Construct the basis dual to ``self``.

        OUTPUT:

        - instance of :class:`CoordCoFrame` representing the dual of
          ``self``

        """
        return CoordCoFrame(self)

    ###### End of methods redefined by derived classes ######

    def chart(self):
        r"""
        Return the chart defining this coordinate frame.
        """
        return self._chart

    def structure_coef(self):
        r"""
        Returns the structure coefficients associated to the coordinate frame.

        `n` being the manifold's dimension, the structure coefficients of the
        frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}`
        defined by

        .. MATH::

            [e_i, e_j] = C^k_{\ \, ij} e_k

        In the present case, where `(e_i)` is a coordinate frame,
        `C^k_{\ \, ij}=0`.

        OUPUT:

        - the structure coefficients `C^k_{\ \, ij}`, as a vanishing instance
          of :class:`~sage.tensor.modules.comp.CompWithSym` with 3 indices
          ordered as `(k,i,j)`

        EXAMPLE:

        Structure coefficients of the coordinate frame associated to
        spherical coordinates in the Euclidean space `\RR^3`::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_spher = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: b = M.default_frame() ; b
            coordinate frame (R^3, (d/dr,d/dth,d/dph))
            sage: c = b.structure_coef() ; c
            3-indices components w.r.t. coordinate frame (R^3, (d/dr,d/dth,d/dph)), with antisymmetry on the index positions (1, 2)
            sage: c == 0
            True

        """
        from sage.tensor.modules.comp import CompWithSym
        if self._structure_coef is None:
            self._structure_coef = CompWithSym(self._fmodule._ring, self, 3,
                                               start_index=self._fmodule._sindex,
                                output_formatter=self._fmodule._output_formatter,
                                                                 antisym=(1,2))
            # A just created CompWithSym is zero
        return self._structure_coef


#******************************************************************************

class CoFrame(FreeModuleCoBasis):
    r"""
    Coframe on a differentiable manifold over `\RR`.

    By *coframe*, it is meant a field `f` on some open domain `U` of a
    manifold `S` endowed with a differentiable mapping `\Phi: U\rightarrow V`
    to a parallelizable domain `V` of a manifold `M` such that for each
    `p\in U`, `f(p)` is a basis of the vector space dual to the tangent space
    `T_{\Phi(p)}M`.

    The standard case of a coframe *on* `U` corresponds to `S=M`, `U=V`
    and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `V` (`U` is then an open interval
    of `\RR`).

    INPUT:

    - ``frame`` -- the vector frame dual to the coframe
    - ``symbol`` -- a letter (of a few letters) to denote a generic 1-form in
      the coframe
    - ``latex_symbol`` -- (default: ``None``) symbol to denote a generic 1-form in
      the coframe; if ``None``, the value of ``symbol`` is used.

    EXAMPLES:

    Coframe on a 3-dimensional manifold::

        sage: DiffManifold._clear_cache_() # for doctests only
        sage: M = DiffManifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: v = M.vector_frame('v')
        sage: from sage.manifolds.differentiable.vectorframe import CoFrame
        sage: e = CoFrame(v, 'e') ; e
        coframe (M, (e^1,e^2,e^3))

    Instead of importing CoFrame in the global namespace, the coframe can be
    obtained by means of the method
    :meth:`~sage.tensor.modules.free_module_basis.FreeModuleBasis.dual_basis`,
    but the symbol is then the same as the frame::

        sage: a = v.dual_basis() ; a
        coframe (M, (v^1,v^2,v^3))
        sage: a[1] == e[1]
        True
        sage: a[1] is e[1]
        False
        sage: e[1].display(v)
        e^1 = v^1

    The 1-forms composing the coframe are obtained via the () operator::

        sage: e[1], e[2], e[3]
        (1-form 'e^1' on the 3-dimensional manifold 'M',
         1-form 'e^2' on the 3-dimensional manifold 'M',
         1-form 'e^3' on the 3-dimensional manifold 'M')

    Checking that e is the dual of v::

        sage: e[1](v[1]).expr(), e[1](v[2]).expr(), e[1](v[3]).expr()
        (1, 0, 0)
        sage: e[2](v[1]).expr(), e[2](v[2]).expr(), e[2](v[3]).expr()
        (0, 1, 0)
        sage: e[3](v[1]).expr(), e[3](v[2]).expr(), e[3](v[3]).expr()
        (0, 0, 1)

    """
    def __init__(self, frame, symbol, latex_symbol=None):
        self._domain = frame._domain
        self._manifold = self._domain._manifold
        FreeModuleCoBasis.__init__(self, frame, symbol,
                                   latex_symbol=latex_symbol)
        # Redefinition of the name and the LaTeX name:
        self._name = "(" + self._domain._name + ", " + self._name + ")"
        self._latex_name = r"\left(" + self._domain._latex_name + ", " + \
                          self._latex_name + r"\right)"
        # The coframe is added to the domain's set of coframes, as well as to
        # all the superdomains' sets of coframes
        for sd in self._domain._supersets:
            for other in sd._coframes:
                if repr(self) == repr(other):
                    raise ValueError("The " + str(self) + " already exist on" +
                                     " the " + str(sd))
            sd._coframes.append(self)


    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "coframe " + self._name

    def at(self, point):
        r"""
        Return the value of the coframe at a given point on the manifold, i.e.
        a basis of the dual of the tangent vector space.

        INPUT:

        - ``point`` -- (instance of
          :class:`~sage.manifolds.differentiable.point.ManifoldPoint`) point `p` in the
          domain of ``self`` (denoted `f` hereafter)

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleCoBasis`
          representing the basis `f(p)` of the vector space `T_p^* M`,
          dual of the tangent space to `M` at `p`
          (`M` being the manifold on which ``self`` is defined)

        EXAMPLES:

        Cobasis of a tangent space on a 2-dimensional manifold::

            sage: from sage.manifolds.differentiable.tangentspace import TangentSpace # for doctests only
            sage: TangentSpace._clear_cache_() ; DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((-1,2), name='p')
            sage: f = X.coframe() ; f
            coordinate coframe (M, (dx,dy))
            sage: fp = f.at(p) ; fp
            Dual basis (dx,dy) on the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: type(fp)
            <class 'sage.tensor.modules.free_module_basis.FreeModuleCoBasis'>
            sage: fp[0]
            Linear form dx on the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: fp[1]
            Linear form dy on the tangent space at point 'p' on 2-dimensional manifold 'M'
            sage: fp is X.frame().at(p).dual_basis()
            True

        """
        return self._basis.at(point).dual_basis()

#******************************************************************************

class CoordCoFrame(CoFrame):
    r"""
    Coordinate coframe on a differentiable manifold over `\RR`.

    By *coordinate coframe*, it is meant the n-tuple of the differentials of
    the coordinates of some chart on the manifold.

    INPUT:

    - ``coord_frame`` -- coordinate frame dual to the coordinate coframe

    EXAMPLES:

    Coordinate coframe on a 3-dimensional manifold::

        sage: DiffManifold._clear_cache_() # for doctests only
        sage: M = DiffManifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: M.frames()
        [coordinate frame (M, (d/dx,d/dy,d/dz))]
        sage: M.coframes()
        [coordinate coframe (M, (dx,dy,dz))]
        sage: dX = M.coframes()[0] ; dX
        coordinate coframe (M, (dx,dy,dz))

    The 1-forms composing the coframe are obtained via the () operator::

        sage: dX[1]
        1-form 'dx' on the 3-dimensional manifold 'M'
        sage: dX[2]
        1-form 'dy' on the 3-dimensional manifold 'M'
        sage: dX[3]
        1-form 'dz' on the 3-dimensional manifold 'M'
        sage: dX[1][:]
        [1, 0, 0]
        sage: dX[2][:]
        [0, 1, 0]
        sage: dX[3][:]
        [0, 0, 1]

    The coframe is the dual of the coordinate frame::

        sage: e = c_xyz.frame() ; e
        coordinate frame (M, (d/dx,d/dy,d/dz))
        sage: dX[1](e[1]).expr(), dX[1](e[2]).expr(), dX[1](e[3]).expr()
        (1, 0, 0)
        sage: dX[2](e[1]).expr(), dX[2](e[2]).expr(), dX[2](e[3]).expr()
        (0, 1, 0)
        sage: dX[3](e[1]).expr(), dX[3](e[2]).expr(), dX[3](e[3]).expr()
        (0, 0, 1)

    Each 1-form of a coordinate coframe is closed::

        sage: dX[1].exterior_der()
        2-form 'ddx' on the 3-dimensional manifold 'M'
        sage: dX[1].exterior_der() == 0
        True

    """
    def __init__(self, coord_frame):
        from sage.misc.latex import latex
        if not isinstance(coord_frame, CoordFrame):
            raise TypeError("The first argument must be a coordinate frame.")
        CoFrame.__init__(self, coord_frame, 'X') # 'X' = provisory symbol
        self._chart = coord_frame._chart
        n = self._manifold._dim
        for i in range(n):
            self._form[i]._name = "d" + str(self._chart._xx[i])
            self._form[i]._latex_name = r"\mathrm{d}" + latex(self._chart._xx[i])
        self._name = "(" + self._domain._name + ", (" + \
                    ",".join([self._form[i]._name for i in range(n)]) + "))"
        self._latex_name = r"\left(" + self._domain._latex_name + \
                          r" ,\left(" + \
                ",".join([self._form[i]._latex_name for i in range(n)])+ \
                          r"\right)\right)"

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "coordinate coframe " + self._name
