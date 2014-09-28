r"""
Affine connections

The class :class:`AffConnection` implements affine connections on 
differentiable manifolds over `\RR`. 

A subclass of :class:`AffConnection` is :class:`LeviCivitaConnection` for
Levi-Civita connections associated to pseudo-Riemannian metrics. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014) : initial version

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

from sage.structure.sage_object import SageObject
from domain import Domain

class AffConnection(SageObject):
    r"""
    Affine connection on a differentiable manifold.

    Given a differentiable manifold `M` and denoting by `\mathcal{X}(M)` 
    the `C^\infty(M)`-module of vector fields on `M`, an *affine connection* 
    on `M` is an operator

    .. MATH::

        \begin{array}{cccc}
        \nabla: & \mathcal{X}(M)\times \mathcal{X}(M) & \longrightarrow & 
                 \mathcal{X}(M) \\
                & (u,v) & \longmapsto & \nabla_u v
        \end{array}
                
    that
    
    - is a bilinear when considering `\mathcal{X}(M)` as a
      vector space over `\RR`
    - is `C^\infty(M)`-linear w.r.t. the first argument:
      `\forall f\in C^\infty(M),\ \nabla_{fu} v = f\nabla_u v` 
    - obeys Leibniz rule w.r.t. the second argument:
      `\forall f\in C^\infty(M),\ \nabla_u (f v) = \mathrm{d}f(u)\, v + f  \nabla_u v`

    The affine connection `\nabla` gives birth to the *covariant derivative
    operator* acting on tensor fields, denoted by the same symbol:
    
    .. MATH::

        \begin{array}{cccc}
        \nabla: &  T^{(k,l)}(M) & \longrightarrow & T^{(k,l+1)}(M)\\
                & t & \longmapsto & \nabla t
        \end{array}
        
    where `T^{(k,l)}(M)` stands for the `C^\infty(M)`-module of tensor fields
    of type `(k,l)` on `M` (with the convention `T^{(0,0)}(M):=C^\infty(M)`).
    For a vector field `v`,  the covariant derivative `\nabla v` is a 
    type-(1,1) tensor field such that 
    
    .. MATH::
    
        \forall u \in\mathcal{X}(M), \   \nabla_u v = \nabla v(., u) 
        
    More generally for any tensor field `t\in T^{(k,l)}(M)`, we have

    .. MATH::
    
        \forall u \in\mathcal{X}(M), \   \nabla_u t = \nabla t(\ldots, u) 
        
    
    .. NOTE::
    
        The above convention means that, in terms of index notation, 
        the "derivation index" in `\nabla t` is the *last* one:
        
        .. MATH::
        
            \nabla_c t^{a_1\ldots a_k}_{\quad\quad b_1\ldots b_l} =
                (\nabla t)^{a_1\ldots a_k}_{\quad\quad b_1\ldots b_l c}


    INPUT:
    
    - ``domain`` -- the manifold domain on which the connection is defined
      (must be an instance of class 
      :class:`~sage.geometry.manifolds.domain.Domain`)
    - ``name`` -- name given to the affine connection
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the affine 
      connection

    EXAMPLES:
    
    Affine connection on a 3-dimensional manifold::
    
        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: nab = M.aff_connection('nabla', r'\nabla') ; nab
        affine connection 'nabla' on the 3-dimensional manifold 'M'
        
    A just-created connection has no connection coefficients::
    
        sage: nab._coefficients
        {}

    The connection coefficients relative to the manifold's default frame 
    [here `(\partial/\partial x, \partial/\partial y, \partial/\partial z)`],
    are created by providing the relevant indices inside square brackets::
    
        sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz 
        sage: nab._coefficients
        {coordinate frame (M, (d/dx,d/dy,d/dz)): 3-indices components w.r.t. coordinate frame (M, (d/dx,d/dy,d/dz))}
        
    Unset components are initialized to zero::
    
        sage: nab[:] # list of coefficients relative to the manifold's default vector frame
        [[[0, x^2, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, y*z], [0, 0, 0]]]
      
    The treatment of connection coefficients in vector frames different 
    from the manifold's default one is similar to that of tensor components; 
    see therefore the class 
    :class:`~sage.geometry.manifolds.tensorfield.TensorField` for the 
    documentation. 
    
    Action on a scalar field::
    
        sage: f = M.scalar_field(x^2 - y^2, name='f')
        sage: Df = nab(f) ; Df
        1-form 'df' on the 3-dimensional manifold 'M'
        sage: Df[:]
        [2*x, -2*y, 0]
        sage: Df == f.differential()  # a necessary condition for any affine connection
        True

    A generic affine connection has some torsion::
    
        sage: DDf = nab(Df) ; DDf
        tensor field 'nabla df' of type (0,2) on the 3-dimensional manifold 'M'
        sage: DDf.antisymmetrize()[:] # nabla does not commute on scalar fields:
        [   0 -x^3    0]
        [ x^3    0    0]
        [   0    0    0]
        
    Let us check the standard formula
    
    .. MATH::
            
        \nabla_j \nabla_i \, f - \nabla_i \nabla_j \, f = T^k_{\ \, ij} \nabla_k \, f , 
        
    where the `T^k_{\ \, ij}`'s are the components of the connection's 
    torsion tensor::
            
        sage: 2*DDf.antisymmetrize() == nab.torsion().contract(0,Df)
        True

    The connection acting on a vector field::

        sage: v = M.vector_field('v')
        sage: v[:] = (y*z, x*z, x*y)
        sage: Dv = nab(v) ; Dv
        field of endomorphisms 'nabla v' on the 3-dimensional manifold 'M'
        sage: Dv[:]
        [            0 (x^2*y + 1)*z             y]
        [            z             0             x]
        [            y             x       x*y*z^2]
    
    Another example: connection on a non-parallelizable 2-dimensional manifold::
    
        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U') ; V = M.open_domain('V') 
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: c_xyW = c_xy.restrict(W) ; c_uvW = c_uv.restrict(W)
        sage: eUW = c_xyW.frame() ; eVW = c_uvW.frame()
        sage: nab = M.aff_connection('nabla', r'\nabla')
        
    The connection is first defined on the domain U by means of its 
    coefficients w.r.t. the frame eU (the manifold's default frame)::
    
        sage: nab[0,0,0], nab[1,0,1] = x, x*y

    The coefficients w.r.t the frame eV are deduced by continuation of the 
    coefficients w.r.t. the frame eVW on the domain `W=U\cap V`::
    
        sage: for i in M.irange():
        ....:     for j in M.irange():
        ....:         for k in M.irange():
        ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()
        ....:             

    At this stage, the connection is fully defined on all the manifold::
    
        sage: nab.coef(eU)[:]
        [[[x, 0], [0, 0]], [[0, x*y], [0, 0]]]
        sage: nab.coef(eV)[:]
        [[[1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v, -1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v],
          [1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v, -1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v]],
         [[-1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v, 1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v],
          [-1/16*u^2 + 1/16*v^2 + 1/8*u + 1/8*v, 1/16*u^2 - 1/16*v^2 + 1/8*u + 1/8*v]]]
    
    We may let it act on a vector field defined globally on `M`::
    
        sage: a = M.vector_field('a')
        sage: a[eU,:] = [-y,x]
        sage: a[eV,0] = a[eVW,0,c_uvW].expr()
        sage: a[eV,1] = a[eVW,1,c_uvW].expr()
        sage: a.view(eU)
        a = -y d/dx + x d/dy
        sage: a.view(eV)
        a = v d/du - u d/dv
        sage: da = nab(a) ; da
        field of endomorphisms 'nabla a' on the 2-dimensional manifold 'M'
        sage: da.view(eU)
        nabla a = -x*y d/dx*dx - d/dx*dy + d/dy*dx - x*y^2 d/dy*dy
        sage: da.view(eV)
        nabla a = (-1/16*u^3 + 1/16*u^2*v + 1/16*(u + 2)*v^2 - 1/16*v^3 - 1/8*u^2) d/du*du + (1/16*u^3 - 1/16*u^2*v - 1/16*(u - 2)*v^2 + 1/16*v^3 - 1/8*u^2 + 1) d/du*dv + (1/16*u^3 - 1/16*u^2*v - 1/16*(u - 2)*v^2 + 1/16*v^3 - 1/8*u^2 - 1) d/dv*du + (-1/16*u^3 + 1/16*u^2*v + 1/16*(u + 2)*v^2 - 1/16*v^3 - 1/8*u^2) d/dv*dv

    A few tests::
    
        sage: nab(a.restrict(V)) == da.restrict(V)
        True
        sage: nab.restrict(V)(a) == da.restrict(V)
        True
        sage: nab.restrict(V)(a.restrict(U)) == da.restrict(W)
        True
        sage: nab.restrict(U)(a.restrict(V)) == da.restrict(W)
        True

    """
    def __init__(self, domain, name, latex_name=None):
        if not isinstance(domain, Domain):
            raise TypeError("The first argument must be a domain.")
        self._manifold = domain._manifold
        self._domain = domain
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._coefficients = {}  # dict. of connection coefficients, with the
                                 # vector frames as keys
        # Initialization of derived quantities:
        AffConnection._init_derived(self) 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "affine connection"
        if self._name is not None:
            description += " '%s'" % self._name
        description += " on the " + str(self._domain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        self._restrictions = {} # dict. of restrictions of ``self`` some 
                                # subdomains, with the subdomains as keys
        self._torsion = None
        self._riemann = None
        self._ricci = None
        self._connection_forms = {}
        self._torsion_forms = {}
        self._curvature_forms = {}

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._restrictions.clear()
        self._torsion = None
        self._riemann = None
        self._ricci = None
        self._connection_forms.clear()
        self._torsion_forms.clear()
        self._curvature_forms.clear()

    def domain(self):
        r"""
        Return the domain on which the affine connection is defined.
        
        OUTPUT:
        
        - instance of class :class:`~sage.geometry.manifolds.domain.OpenDomain` 
          representing the manifold's open subset on which ``self`` is defined. 
        
        EXAMPLES::
        
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab.domain()
            3-dimensional manifold 'M'
            sage: U = M.open_domain('U', coord_def={c_xyz: x>0})
            sage: nabU = U.aff_connection('D')
            sage: nabU.domain()
            open domain 'U' on the 3-dimensional manifold 'M'

        """
        return self._domain


    def _new_coef(self, frame): 
        r"""
        Create the connection coefficients w.r.t. the given frame. 
        
        This method, to be called by :meth:`coef`, must be redefined by derived 
        classes to adapt the output to the relevant subclass of 
        :class:`~sage.tensor.modules.comp.Components`.
        
        """
        from sage.tensor.modules.comp import Components
        from scalarfield import ScalarField
        return Components(frame._domain.scalar_field_algebra(), frame, 3, 
                          start_index=self._manifold._sindex,
                          output_formatter=ScalarField.function_chart)
        
    def coef(self, frame=None):
        r"""
        Return the connection coefficients relative to the given frame.
        
        `n` being the manifold's dimension, the connection coefficients 
        relative to the vector frame `(e_i)` are the `n^3` scalar fields 
        `\Gamma^k_{\ \, ij}` defined by 
        
        .. MATH::
            
            \nabla_{e_j} e_i = \Gamma^k_{\ \, ij} e_k
        
        
        If the connection coefficients are not known already, they are computed
        from the above formula. 
        
        INPUT:
        
        - ``frame`` -- (default: None) vector frame relative to which the 
          connection coefficients are required; if none is provided, the 
          domain's default frame is assumed
 
        OUTPUT: 
        
        - connection coefficients relative to the frame ``frame``, as an 
          instance of the class :class:`~sage.tensor.modules.comp.Components` 
          with 3 indices ordered as `(k,i,j)`
        
        EXAMPLES:      
 
        Connection coefficient of an affine connection on a 3-dimensional
        manifold::
    
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: nab.coef()
            3-indices components w.r.t. coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: print type(nab.coef())
            <class 'sage.tensor.modules.comp.Components'>
            sage: M.default_frame()
            coordinate frame (M, (d/dx,d/dy,d/dz))
            sage: nab.coef() is nab.coef(c_xyz.frame())
            True
            sage: nab.coef()[:]  # full list of coefficients:
            [[[0, x^2, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, y*z], [0, 0, 0]]]

        """
        if frame is None: 
            frame = self._domain._def_frame
        if frame not in self._coefficients:
            # the coefficients must be computed
            #
            # Check whether frame is a subframe of a frame in which the 
            # coefficients are already known:
            for oframe in self._coefficients:
                if frame in oframe._subframes:
                    self._coefficients[frame] = self._new_coef(frame)
                    comp_store = self._coefficients[frame]._comp
                    ocomp_store = self._coefficients[oframe]._comp
                    for ind, value in ocomp_store.iteritems():
                        comp_store[ind] = value.restrict(frame._domain)
                    break
            else:
                # If not, the coefficients must be computed from scratch:
                manif = self._manifold
                ev = frame  # the vector frame
                ef = ev._coframe # the dual frame
                gam = self._new_coef(ev)
                for k in manif.irange():
                    for i in manif.irange():
                        for j in manif.irange():
                            gam[[k,i,j]] = self(ev[i])(ef[k],ev[j])
                self._coefficients[frame] = gam
        return self._coefficients[frame]
        

    def set_coef(self, frame=None):
        r"""
        Return the connection coefficients in a given frame for assignment.
        
        See method :meth:`coef` for details about the definition of the 
        connection coefficents. 

        The connection coefficients with respect to other frames are deleted, 
        in order to avoid any inconsistency. To keep them, use the method 
        :meth:`add_coef` instead.

        INPUT:
        
        - ``frame`` -- (default: None) vector frame in which the connection 
          coefficients are defined; if none is provided, the domain's default 
          frame is assumed.

        OUTPUT: 
        
        - connection coefficients in the given frame, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components`; if such 
          connection coefficients did not exist previously, they are created. 
          See method :meth:`coef` for the storage convention of the connection 
          coefficents. 
        

        """
        if frame is None: frame = self._domain._def_frame
        if frame not in self._coefficients:
            if frame not in self._domain._frames:
                raise ValueError("The vector frame " + frame +
                                 " has not been defined on the " + 
                                 str(self._domain))
            self._coefficients[frame] = self._new_coef(frame)
        self._del_derived() # deletes the derived quantities
        self.del_other_coef(frame)
        return self._coefficients[frame]

    def add_coef(self, frame=None):
        r"""
        Return the connection coefficients in a given frame for assignment, 
        keeping the coefficients in other frames. 
        
        See method :meth:`coef` for details about the definition of the 
        connection coefficents. 

        To delete the connection coefficients in other frames, use the method 
        :meth:`set_coef` instead.
        
        INPUT:
        
        - ``frame`` -- (default: None) vector frame in which the connection 
          coefficients are defined; if none is provided, the domain's default 
          frame is assumed.
          
        .. WARNING::
        
            If the connection has already coefficients in other frames, it 
            is the user's responsability to make sure that the coefficients
            to be added are consistent with them. 
         
        OUTPUT: 
        
        - connection coefficients in the given frame, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components`; if such 
          connection coefficients did not exist previously, they are created. 
          See method :meth:`coef` for the storage convention of the connection 
          coefficents. 
        

        """
        if frame is None: frame = self._domain._def_frame
        if frame not in self._coefficients:
            if frame not in self._domain._frames:
                raise ValueError("The vector frame " + frame +
                                 " has not been defined on the " + 
                                 str(self._domain))
            self._coefficients[frame] = self._new_coef(frame)
        self._del_derived() # deletes the derived quantities
        return self._coefficients[frame]


    def del_other_coef(self, frame=None):
        r"""
        Delete all the coefficients but those corresponding to ``frame``.
        
        """
        if frame is None: frame = self._domain._def_frame
        if frame not in self._coefficients:
            raise ValueError("The coefficients w.r.t. the vector frame " + 
                             frame + " have not been defined.")
        to_be_deleted = []
        for other_frame in self._coefficients:
            if other_frame != frame:
                to_be_deleted.append(other_frame)
        for other_frame in to_be_deleted:
            del self._coefficients[other_frame]

    def __getitem__(self, indices):
        r"""
        Return the connection coefficient w.r.t. the domain default frame 
        corresponding to the given indices. 

        INPUT:
        
        - ``indices`` -- list of indices
    
        """
        return self.coef()[indices]
        
    def __setitem__(self, indices, value):
        r"""
        Set the connection coefficient w.r.t. the domain default frame 
        corresponding to the given indices.

        INPUT:
        
        - ``indices`` -- list of indices
    
        """        
        self.set_coef()[indices] = value

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some subdomain.
        
        If such restriction has not been defined yet, it is constructed here.

        INPUT:
        
        - ``subdomain`` -- open subset `U` of ``self._domain`` (must be an 
          instance of :class:`~sage.geometry.manifolds.domain.OpenDomain`)
          
        OUTPUT:
        
        - instance of :class:`AffConnection` representing the restriction.

        EXAMPLE:
        
        Restriction of a connection on a 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[2,1,1] = x^2, x+y
            sage: nab[:]
            [[[0, x^2], [0, 0]], [[x + y, 0], [0, 0]]]
            sage: U = M.open_domain('U', coord_def={c_xy: x>0})
            sage: nabU = nab.restrict(U) ; nabU
            affine connection 'nabla' on the open domain 'U' on the 2-dimensional manifold 'M'
            sage: nabU.domain()
            open domain 'U' on the 2-dimensional manifold 'M'
            sage: nabU[:]
            [[[0, x^2], [0, 0]], [[x + y, 0], [0, 0]]]
            
        The result is cached::
        
            sage: nab.restrict(U) is nabU
            True

        ...until the connection is modified::
        
            sage: nab[1,2,2] = -y
            sage: nab.restrict(U) is nabU
            False
            sage: nab.restrict(U)[:]
            [[[0, x^2], [0, -y]], [[x + y, 0], [0, 0]]]

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subdomain(self._domain):
                raise ValueError("The provided domain is not a subdomain of " + 
                                 "the current connection's domain.")
            resu = AffConnection(subdomain, name=self._name, 
                                 latex_name=self._latex_name)
            for frame in self._coefficients:
                for sframe in subdomain._top_frames: 
                    if sframe in frame._subframes:
                        comp_store = self._coefficients[frame]._comp
                        scoef = resu._new_coef(sframe)
                        scomp_store = scoef._comp
                        # the coefficients of the restriction are evaluated 
                        # index by index:
                        for ind, value in comp_store.iteritems():
                            scomp_store[ind] = value.restrict(sframe._domain)
                        resu._coefficients[sframe] = scoef
            if self._torsion is not None:
                resu._torsion = self._torsion.restrict(subdomain)
            if self._riemann is not None:
                resu._riemann = self._riemann.restrict(subdomain)
            if self._ricci is not None:
                resu._ricci = self._ricci.restrict(subdomain)
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]
        
    def _common_frame(self, other):
        r"""
        Find a common vector frame for the coefficients of ``self`` and
        the components of  ``other``. 
        
        In case of multiple common frames, the domain's default frame is 
        privileged. 
        
        INPUT:
        
        - ``other`` -- a tensor field on parallelizable domain, as an 
          instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorFieldParal`
        
        OUPUT:
        
        - common frame; if no common frame is found, None is returned. 
        
        """
        # The domain of search is restricted to other._domain: 
        dom = other._domain
        # 1/ Does each object have components on the domain's default frame ? 
        def_frame = dom._def_frame
        if def_frame in self._coefficients and def_frame in other._components:
            return def_frame
        # 2/ Search for a common frame among the existing components, i.e. 
        #    without performing any component transformation. 
        #    -------------------------------------------------------------
        for frame in self._coefficients:
            if frame in other._components:
                return frame
        # 3/ Search for a common frame among the subframes of self's frames: 
        #    --------------------------------------------------------------
        for frame in self._coefficients:
            for oframe in other._components:
                if oframe in frame._subframes:
                    self.coef(oframe) # update the coefficients of self in oframe
                    return oframe
        #
        # 4/ Search for a common frame via one component transformation 
        #    ----------------------------------------------------------
        # If this point is reached, it is necessary to perform at least 
        # one component transformation to get a common frame
        for frame in self._coefficients:
            for oframe in other._components:
                if (oframe, frame) in dom._frame_changes:
                    other.comp(frame, from_basis=oframe)
                    return frame
        # 5/ Search for a common frame via one component transformation to
        #    a subframe of self's frames: 
        #    -------------------------------------------------------------
        for frame in self._coefficients:
            for oframe in other._components:
                for sframe in frame._subframes:
                    if (oframe, sframe) in dom._frame_changes:
                        self.coef(sframe)
                        other.comp(sframe, from_basis=oframe)
                        return sframe
        #
        # If this point is reached, no common frame could be found, even at 
        # the price of a component transformation:
        return None

    def __call__(self, tensor):
        r"""
        Action of the connection on a tensor field.
        
        INPUT:
        
        - ``tensor`` -- a tensor field `T`, of type `(k,\ell)`
        
        OUTPUT:
        
        - tensor field `\nabla T`. 
          
        """
        from tensorfield import TensorFieldParal
        from utilities import format_unop_txt, format_unop_latex
        dom_resu = self._domain.intersection(tensor._domain)
        tensor_r = tensor.restrict(dom_resu)
        if tensor_r._tensor_type == (0,0):  # scalar field case
            return tensor_r.differential()   
        if isinstance(tensor_r, TensorFieldParal):
            return self._derive_paral(tensor_r)
        resu_rst = []
        for dom, rst in tensor_r._restrictions.iteritems():
            # the computation is performed only if dom is not a subdomain
            # of another restriction:
            for odom in tensor_r._restrictions:
                if dom in odom._subdomains and dom is not odom:
                    break
            else:
                # dom is a not a subdomain and the computation is performed:
                resu_rst.append(self.__call__(rst))
        tensor_type_resu = (tensor_r._tensor_type[0], tensor_r._tensor_type[1]+1)
        name_resu = format_unop_txt(self._name + ' ', tensor_r._name)
        latex_name_resu=format_unop_latex(self._latex_name + ' ', 
                                                            tensor_r._latex_name) 
        vmodule = dom_resu.vector_field_module()
        resu = vmodule.tensor(tensor_type_resu, name=name_resu,
                              latex_name=latex_name_resu,
                              sym=resu_rst[0]._sym, 
                              antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def _derive_paral(self, tensor):
        r"""
        Action of the connection on a tensor field on a parallelizable domain. 
        
        INPUT:
        
        - ``tensor`` -- a tensor field `T`, of type `(k,\ell)`
        
        OUTPUT:
        
        - tensor field `\nabla T`. 
          
        """
        from scalarfield import ScalarField
        from sage.tensor.modules.comp import Components, CompWithSym
        from utilities import format_unop_txt, format_unop_latex
        manif = self._manifold
        tdom = tensor._domain
        frame = self._common_frame(tensor)
        if frame is None:
            raise ValueError("No common frame found for the computation.")
        # Component computation in the common frame:
        tc = tensor._components[frame]
        gam = self._coefficients[frame]
        if tensor._sym == [] and tensor._antisym == []:
            resc = Components(tdom.scalar_field_algebra(), frame,
                              tensor._tensor_rank+1, 
                              start_index=self._manifold._sindex,
                              output_formatter=ScalarField.function_chart)
        else:
            resc = CompWithSym(tdom.scalar_field_algebra(), frame,
                              tensor._tensor_rank+1, 
                              start_index=self._manifold._sindex,
                              output_formatter=ScalarField.function_chart,
                              sym=tensor._sym, antisym=tensor._antisym)
        n_con = tensor._tensor_type[0]
        n_cov = tensor._tensor_type[1]
        for ind in resc.non_redundant_index_generator():
            p = ind[-1]  # derivation index
            ind0 = ind[:-1]
            rsum = frame[p](tc[[ind0]])
            # loop on contravariant indices:
            for k in range(n_con): 
                for i in manif.irange():
                    indk = list(ind0)
                    indk[k] = i  
                    rsum += gam[[ind0[k], i, p]] * tc[[indk]]
            # loop on covariant indices:
            for k in range(n_con, tensor._tensor_rank): 
                for i in manif.irange():
                    indk = list(ind0)
                    indk[k] = i  
                    rsum -= gam[[i, ind0[k], p]] * tc[[indk]]
            resc[[ind]] = rsum
        # Resulting tensor field
        return tdom.vector_field_module().tensor_from_comp((n_con, n_cov+1),
                        resc, 
                        name=format_unop_txt(self._name + ' ', tensor._name),
                        latex_name=format_unop_latex(self._latex_name + ' ', 
                                                        tensor._latex_name) )

    def torsion(self):
        r""" 
        Return the connection's torsion tensor.
        
        The torsion tensor is the tensor field `T` of type (1,2) defined by

        .. MATH::
            
            T(\omega, u, v) = \left\langle \omega, \nabla_u v - \nabla_v u
                - [u, v] \right\rangle
        
        for any 1-form  `\omega`  and any vector fields `u` and `v`. 
                  
        OUTPUT:
        
        - the torsion tensor `T`, as an instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
        
        EXAMPLES:
        
        Torsion of an affine connection on a 3-dimensional manifold::
    
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: t = nab.torsion() ; t
            tensor field of type (1,2) on the 3-dimensional manifold 'M'
            sage: t.symmetries()
            no symmetry;  antisymmetry: (1, 2)
            sage: t[:]
            [[[0, -x^2, 0], [x^2, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, -y*z], [0, y*z, 0]]]

        The torsion expresses the lack of commutativity of two successive 
        derivatives of a scalar field::
        
            sage: f = M.scalar_field(x*z^2 + y^2 - z^2, name='f')
            sage: DDf = nab(nab(f)) ; DDf
            tensor field 'nabla df' of type (0,2) on the 3-dimensional manifold 'M'
            sage: DDf.antisymmetrize()[:]  # two successive derivatives do not commute:
            [             0   -1/2*x^2*z^2              0]
            [   1/2*x^2*z^2              0 -(x - 1)*y*z^2]
            [             0  (x - 1)*y*z^2              0]
            sage: 2*DDf.antisymmetrize() == nab.torsion().contract(0,nab(f))
            True

        The above identity is the standard formula
         
        .. MATH::
            
            \nabla_j \nabla_i \, f - \nabla_i \nabla_j \, f = T^k_{\ \, ij} \nabla_k \, f , 
        
        where the `T^k_{\ \, ij}`'s are the components of the torsion tensor. 

        The result is cached::
        
            sage: nab.torsion() is t
            True

        ...as long as the connection remains unchanged::
        
            sage: nab[2,1,3] = 1+x    # changing the connection
            sage: nab.torsion() is t  # a new computation of the torsion has been made
            False
            sage: (nab.torsion() - t).view()
            (-x - 1) d/dy*dx*dz + (x + 1) d/dy*dz*dx

        Another example: torsion of some connection on a non-parallelizable 
        2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U') ; V = M.open_domain('V') 
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: c_xyW = c_xy.restrict(W) ; c_uvW = c_uv.restrict(W)
            sage: eUW = c_xyW.frame() ; eVW = c_uvW.frame()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[0,0,0], nab[0,1,0], nab[1,0,1] = x, x-y, x*y
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         for k in M.irange():
            ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()
            ....: 
            sage: t = nab.torsion() ; t
            tensor field of type (1,2) on the 2-dimensional manifold 'M'
            sage: t.parent()
            module T^(1,2)(M) of type-(1,2) tensors fields on the 2-dimensional manifold 'M'          
            sage: t[eU,:]
            [[[0, x - y], [-x + y, 0]], [[0, -x*y], [x*y, 0]]]
            sage: t[eV,:]
            [[[0, 1/8*u^2 - 1/8*v^2 - 1/2*v], [-1/8*u^2 + 1/8*v^2 + 1/2*v, 0]],
             [[0, -1/8*u^2 + 1/8*v^2 - 1/2*v], [1/8*u^2 - 1/8*v^2 + 1/2*v, 0]]]

        Check of the torsion formula::
        
            sage: f = M.scalar_field({c_xy: (x+y)^2, c_uv: u^2}, name='f')
            sage: DDf = nab(nab(f)) ; DDf
            tensor field 'nabla df' of type (0,2) on the 2-dimensional manifold 'M'
            sage: DDf.antisymmetrize().view(eU)
            (-x^2*y - (x + 1)*y^2 + x^2) dx/\dy
            sage: DDf.antisymmetrize().view(eV)
            (1/8*u^3 - 1/8*u*v^2 - 1/2*u*v) du/\dv
            sage: 2*DDf.antisymmetrize() == nab(f).contract(nab.torsion())
            True


        """
        if self._torsion is None:
            manif = self._manifold
            resu = self._domain.tensor_field(1, 2, antisym=(1,2))
            for frame, gam in self._coefficients.iteritems():
                sc = frame.structure_coef()
                res = resu.add_comp(frame)
                for k in manif.irange():
                    for i in manif.irange():
                         for j in manif.irange(start=i+1):
                             res[[k,i,j]] = gam[[k,j,i]] - gam[[k,i,j]] - \
                                            sc[[k,i,j]]
            self._torsion = resu
        return self._torsion 

    def riemann(self):
        r""" 
        Return the connection's Riemann curvature tensor.

        The Riemann curvature tensor is the tensor field `R` of type (1,3) 
        defined by

        .. MATH::
            
            R(\omega, u, v, w) = \left\langle \omega, \nabla_u \nabla_v w
                - \nabla_v \nabla_u w - \nabla_{[u, v]} w \right\rangle
        
        for any 1-form  `\omega`  and any vector fields `u`, `v` and `w`. 

        OUTPUT:
        
        - the Riemann curvature tensor `R`, as an instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
        
        EXAMPLES:
        
        Curvature of an affine connection on a 3-dimensional manifold::
            
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla') ; nab
            affine connection 'nabla' on the 3-dimensional manifold 'M'
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz 
            sage: r = nab.riemann() ; r
            tensor field of type (1,3) on the 3-dimensional manifold 'M'
            sage: r.parent()
            free module T^(1,3)(M) of type-(1,3) tensors fields on the 3-dimensional manifold 'M'
            sage: r.symmetries()
            no symmetry;  antisymmetry: (2, 3)
            sage: r[:]
            [[[[0, 2*x, 0], [-2*x, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, z], [0, -z, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]]]

        Another example: Riemann curvature tensor of some connection on a 
        non-parallelizable 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U') ; V = M.open_domain('V') 
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: c_xyW = c_xy.restrict(W) ; c_uvW = c_uv.restrict(W)
            sage: eUW = c_xyW.frame() ; eVW = c_uvW.frame()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[0,0,0], nab[0,1,0], nab[1,0,1] = x, x-y, x*y
            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         for k in M.irange():
            ....:             nab.add_coef(eV)[i,j,k] = nab.coef(eVW)[i,j,k,c_uvW].expr()
            ....: 
            sage: r = nab.riemann() ; r
            tensor field of type (1,3) on the 2-dimensional manifold 'M'
            sage: r.parent()
            module T^(1,3)(M) of type-(1,3) tensors fields on the 2-dimensional manifold 'M'
            sage: r.view(eU)
            (x^2*y - x*y^2) d/dx*dx*dx*dy + (-x^2*y + x*y^2) d/dx*dx*dy*dx + d/dx*dy*dx*dy - d/dx*dy*dy*dx - (x^2 - 1)*y d/dy*dx*dx*dy + (x^2 - 1)*y d/dy*dx*dy*dx + (-x^2*y + x*y^2) d/dy*dy*dx*dy + (x^2*y - x*y^2) d/dy*dy*dy*dx
            sage: r.view(eV)
            (1/32*u^3 - 1/32*u*v^2 - 1/32*v^3 + 1/32*(u^2 + 4)*v - 1/8*u - 1/4) d/du*du*du*dv + (-1/32*u^3 + 1/32*u*v^2 + 1/32*v^3 - 1/32*(u^2 + 4)*v + 1/8*u + 1/4) d/du*du*dv*du + (1/32*u^3 - 1/32*u*v^2 + 3/32*v^3 - 1/32*(3*u^2 - 4)*v - 1/8*u + 1/4) d/du*dv*du*dv + (-1/32*u^3 + 1/32*u*v^2 - 3/32*v^3 + 1/32*(3*u^2 - 4)*v + 1/8*u - 1/4) d/du*dv*dv*du + (-1/32*u^3 + 1/32*u*v^2 + 5/32*v^3 - 1/32*(5*u^2 + 4)*v + 1/8*u - 1/4) d/dv*du*du*dv + (1/32*u^3 - 1/32*u*v^2 - 5/32*v^3 + 1/32*(5*u^2 + 4)*v - 1/8*u + 1/4) d/dv*du*dv*du + (-1/32*u^3 + 1/32*u*v^2 + 1/32*v^3 - 1/32*(u^2 + 4)*v + 1/8*u + 1/4) d/dv*dv*du*dv + (1/32*u^3 - 1/32*u*v^2 - 1/32*v^3 + 1/32*(u^2 + 4)*v - 1/8*u - 1/4) d/dv*dv*dv*du

        """
        if self._riemann is None:
            manif = self._manifold
            resu = self._domain.tensor_field(1, 3, antisym=(2,3))
            for frame, gam in self._coefficients.iteritems():
                # The computation is performed only on the top frames:
                for oframe in self._coefficients:
                    if frame in oframe._subframes and frame is not oframe:
                        break
                else:
                    # frame in not a subframe and the computation is performed:
                    sc = frame.structure_coef()
                    gam_gam = gam.contract(1, gam, 0)
                    gam_sc = gam.contract(2, sc, 0)
                    res = resu.add_comp(frame)
                    for i in manif.irange():
                        for j in manif.irange():
                            for k in manif.irange():
                                # antisymmetry of the Riemann tensor taken into 
                                # account by l>k: 
                                for l in manif.irange(start=k+1):
                                    res[i,j,k,l] = frame[k](gam[[i,j,l]]) - \
                                                   frame[l](gam[[i,j,k]]) + \
                                                   gam_gam[[i,k,j,l]] -  \
                                                   gam_gam[[i,l,j,k]] -  \
                                                   gam_sc[[i,j,k,l]]
            self._riemann = resu
        return self._riemann 
        

    def ricci(self):
        r""" 
        Return the connection's Ricci tensor.
        
        The Ricci tensor is the tensor field `Ric` of type (0,2) 
        defined from the Riemann curvature tensor `R` by 

        .. MATH::
            
            Ric(u, v) = R(e^i, u, e_i, v)
        
        for any vector fields `u` and `v`, `(e_i)` being any vector frame and
        `(e^i)` the dual coframe. 
        
        OUTPUT:
        
        - the Ricci  tensor `Ric`, as an instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
        
        EXAMPLES:
        
        Ricci tensor of an affine connection on a 3-dimensional manifold::
            
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla') ; nab
            affine connection 'nabla' on the 3-dimensional manifold 'M'
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz 
            sage: r = nab.ricci() ; r
            tensor field of type (0,2) on the 3-dimensional manifold 'M'
            sage: r[:]
            [  0 2*x   0]
            [  0  -z   0]
            [  0   0   0]
       
        """
        if self._ricci is None:
            self._ricci = self.riemann().trace(0,2)
        return self._ricci 
        
    def connection_form(self, i, j, frame=None):
        r"""
        Return the connection 1-form corresponding to the given index and
        vector frame.
        
        The connection 1-forms with respect to the frame `(e_i)` are the 
        `n^2` 1-forms `\omega^i_{\ \, j}` defined by 
        
        .. MATH::
        
            \nabla_v e_j = \langle \omega^i_{\ \, j}, v \rangle
                \, e_i
                
        for any vector `v`. 
        
        The components of `\omega^i_{\ \, j}` in the coframe `(e^i)` dual to 
        `(e_i)` are nothing but the connection coefficients `\Gamma^i_{\ \, jk}`
        relative to the frame `(e_i)`:

        .. MATH::
        
            \omega^i_{\ \, j} = \Gamma^i_{\ \, jk} e^k
        
        
        INPUT:
        
        - ``i``, ``j`` -- indices identifying the 1-form `\omega^i_{\ \, j}`
        - ``frame`` -- (default: None) vector frame relative to which the 
          connection 1-forms are defined; if none is provided, the default 
          frame of the connection's domain is assumed. 
          
        OUTPUT:
        
        - the 1-form `\omega^i_{\ \, j}`, as an instance of 
          :class:`~sage.geometry.manifolds.diffform.OneForm`
        
        EXAMPLES:
        
        Connection 1-forms on a 3-dimensional manifold::
        
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.connection_form(1,1)  # connection 1-form (i,j)=(1,1) w.r.t. M's default frame
            1-form 'nabla connection 1-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.connection_form(1,1)[:]
            [x*y*z, x^2, -y*z]
            
        Connection 1-forms w.r.t. a non-holonomic frame::
        
            sage: ch_basis = M.automorphism_field()
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = M.default_frame().new_frame(ch_basis, 'e')
            sage: e[1][:], e[2][:], e[3][:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: nab.connection_form(1,1,e)
            1-form 'nabla connection 1-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.connection_form(1,1,e).comp(e)[:]
            [x*y^2*z, (x^2*y + 1)*z/y, -x*y*z]
            
        Check of the formula `\omega^i_{\ \, j} = \Gamma^i_{\ \, jk} e^k`::
        
            sage: #... on the manifold's default frame (d/dx, d/dy, d:dz)
            sage: dx = M.default_frame().coframe() ; dx
            coordinate coframe (M, (dx,dy,dz))
            sage: check = []
            sage: for i in M.irange():
            ...       for j in M.irange():
            ...           check.append( nab.connection_form(i,j) == sum( nab[[i,j,k]]*dx[k] for k in M.irange() ) )
            ...
            sage: check
            [True, True, True, True, True, True, True, True, True]
            sage: #... on the frame e
            sage: ef = e.coframe() ; ef
            coframe (M, (e^1,e^2,e^3))
            sage: check = []
            sage: for i in M.irange():
            ...       for j in M.irange():
            ...           s = nab.connection_form(i,j,e).comp(c_xyz.frame(), from_basis=e) 
            ...           check.append( nab.connection_form(i,j,e) == sum( nab.coef(e)[[i,j,k]]*ef[k] for k in M.irange() ) )
            ...
            sage: check
            [True, True, True, True, True, True, True, True, True]

        Check of the formula `\nabla_v e_j = \langle \omega^i_{\ \, j}, v \rangle e_i`::
        
            sage: v = M.vector_field()
            sage: v[:] = (x*y, z^2-3*x, z+2*y)
            sage: b = M.default_frame()
            sage: for j in M.irange():  # check on M's default frame
            ...       nab(b[j]).contract(v) == sum( nab.connection_form(i,j)(v)*b[i] for i in M.irange())
            True
            True
            True
            sage: for j in M.irange():  # check on frame e
            ...       nab(e[j]).contract(v) == sum( nab.connection_form(i,j,e)(v)*e[i] for i in M.irange())
            True
            True
            True

        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._connection_forms:
            forms = {}
            frame_dom = frame._domain
            for i1 in self._manifold.irange():
                for j1 in self._manifold.irange():
                    name = self._name + " connection 1-form (" + str(i1) + \
                           "," + str(j1) + ")"
                    latex_name = r"\omega^" + str(i1) + r"_{\ \, " + str(j1) + \
                                 "}"
                    omega = frame_dom.one_form(name=name, 
                                               latex_name=latex_name)
                    comega = omega.set_comp(frame)
                    for k in self._manifold.irange():
                        comega[k] = self.coef(frame)[[i1,j1,k]]
                    forms[(i1,j1)] = omega
            self._connection_forms[frame] = forms
        return  self._connection_forms[frame][(i,j)] 
                    
    def torsion_form(self, i, frame=None):
        r"""
        Return the torsion 2-form corresponding to the given index and
        vector frame.
        
        The torsion 2-forms with respect to the frame `(e_i)` are the 
        `n` 2-forms `\theta^i` defined by 
        
        .. MATH::
        
            \theta^i(u,v) = T(e^i, u, v)
            
        where `T` is the connection's torsion tensor (cf. :meth:`torsion`),
        `(e^i)` is the coframe dual to `(e_i)` and `(u,v)` is a generic pair of 
        vectors.
        
        INPUT:
        
        - ``i`` -- index identifying the 2-form `\theta^i`
        - ``frame`` -- (default: None) vector frame relative to which the 
          torsion 2-forms are defined; if none is provided, the default frame
          of the connection's domain is assumed. 
          
        OUTPUT:
        
        - the 2-form `\theta^i`, as an instance of 
          :class:`~sage.geometry.manifolds.diffform.DiffForm`
        
        EXAMPLES:
        
        Torsion 2-forms on a 3-dimensional manifold::
        
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.torsion_form(1)
            2-form 'nabla torsion 2-form (1)' on the 3-dimensional manifold 'M'
            sage: nab.torsion_form(1)[:]                               
            [               0             -x^2      (y^2 + y)*z]
            [             x^2                0  x^3 - x^2 + y^2]
            [    -(y^2 + y)*z -x^3 + x^2 - y^2                0]
            
        Torsion 2-forms w.r.t. a non-holonomic frame::
            
            sage: ch_basis = M.automorphism_field()                      
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = M.default_frame().new_frame(ch_basis, 'e')
            sage: e[1][:], e[2][:], e[3][:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: ef = e.coframe()
            sage: ef[1][:], ef[2][:], ef[3][:]
            ([1/y, 0, 0], [0, 1/z, 0], [0, 0, 1/x])
            sage: nab.torsion_form(1, e)
            2-form 'nabla torsion 2-form (1)' on the 3-dimensional manifold 'M'
            sage: nab.torsion_form(1, e).comp(e)[:]
            [                       0                   -x^2*z          (x*y^2 + x*y)*z]
            [                   x^2*z                        0  (x^4 - x^3 + x*y^2)*z/y]
            [        -(x*y^2 + x*y)*z -(x^4 - x^3 + x*y^2)*z/y                        0]
            
        Cartan's first structure equation is
        
        .. MATH::
        
            \theta^i = \mathrm{d} e^i + \omega^i_{\ \, j} \wedge e^j
            
        where the `\omega^i_{\ \, j}`'s are the connection 1-forms (cf. 
        :meth:`connection_form`). Let us check it on the frame e::
        
            sage: for i in M.irange():
            ...       nab.torsion_form(i, e) == ef[i].exterior_der() + sum(nab.connection_form(i,k,e).wedge(ef[k]) for k in M.irange())
            ...
            True
            True
            True

        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._torsion_forms:
            forms = {}
            frame_dom = frame._domain
            for i1 in self._manifold.irange():
                name = self._name + " torsion 2-form (" + str(i1) + ")"
                latex_name = r"\theta^" + str(i1)
                theta = frame_dom.diff_form(2, name=name, 
                                            latex_name=latex_name)
                ctheta = theta.set_comp(frame)
                for k in self._manifold.irange():
                    for l in self._manifold.irange(start=k+1):
                        ctheta[k,l] = self.torsion().comp(frame)[[i1,k,l]]
                forms[i1] = theta
            self._torsion_forms[frame] = forms
        return  self._torsion_forms[frame][i] 
                    
    def curvature_form(self, i, j, frame=None):
        r"""
        Return the curvature 2-form corresponding to the given index and
        vector frame.
        
        The curvature 2-forms with respect to the frame `(e_i)` are the 
        `n^2` 2-forms `\Omega^i_{\ \, j}` defined by 
        
        .. MATH::
        
            \Omega^i_{\ \, j}(u,v) = R(e^i, u, v, e_j)
            
        where `R` is the connection's Riemann curvature tensor (cf. 
        :meth:`riemann`), `(e^i)` is the coframe dual to `(e_i)` and `(u,v)` is a 
        generic pair of vectors.
        
        INPUT:
        
        - ``i``, ``j`` -- indices identifying the 2-form `\Omega^i_{\ \, j}`
        - ``frame`` -- (default: None) vector frame relative to which the 
          curvature 2-forms are defined; if none is provided, the default frame
          of the connection's domain is assumed. 
          
        OUTPUT:
        
        - the 2-form `\Omega^i_{\ \, j}`, as an instance of 
          :class:`~sage.geometry.manifolds.diffform.DiffForm`
        
        EXAMPLES:
        
        Curvature 2-forms on a 3-dimensional manifold::
        
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: nab = M.aff_connection('nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.curvature_form(1,1)
            2-form 'nabla curvature 2-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.curvature_form(1,1).view()
            nabla curvature 2-form (1,1) = (y^2*z^3 + (x*y^3 - x)*z + 2*x) dx/\dy + (x^3*z^2 - x*y) dx/\dz + (x^4*y*z^2 - z) dy/\dz
            
        Curvature 2-forms w.r.t. a non-holonomic frame::
            
            sage: ch_basis = M.automorphism_field()                      
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = M.default_frame().new_frame(ch_basis, 'e')
            sage: e[1].view(), e[2].view(), e[3].view()
            (e_1 = y d/dx, e_2 = z d/dy, e_3 = x d/dz)
            sage: ef = e.coframe()
            sage: ef[1].view(), ef[2].view(), ef[3].view()
            (e^1 = 1/y dx, e^2 = 1/z dy, e^3 = 1/x dz)
            sage: nab.curvature_form(1,1,e)
            2-form 'nabla curvature 2-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.curvature_form(1,1,e).view(e)
            nabla curvature 2-form (1,1) = (y^3*z^4 + 2*x*y*z + (x*y^4 - x*y)*z^2) e^1/\e^2 + (x^4*y*z^2 - x^2*y^2) e^1/\e^3 + (x^5*y*z^3 - x*z^2) e^2/\e^3
            
        Cartan's second structure equation is
         
        .. MATH::
        
            \Omega^i_{\ \, j} = \mathrm{d} \omega^i_{\ \, j} + \omega^i_{\ \, k} \wedge \omega^k_{\ \, j}
            
        where the `\omega^i_{\ \, j}`'s are the connection 1-forms (cf. 
        :meth:`connection_form`). Let us check it on the frame e::
       
            sage: omega = nab.connection_form
            sage: check = []
            sage: for i in M.irange():
            ...       for j in M.irange():
            ...           check.append( nab.curvature_form(i,j,e) == omega(i,j,e).exterior_der() + \
            ...           sum( omega(i,k,e).wedge(omega(k,j,e)) for k in M.irange()) )
            ...
            sage: check
            [True, True, True, True, True, True, True, True, True]
            
        """
        if frame is None:
            frame = self._domain._def_frame
        if frame not in self._curvature_forms:
            forms = {}
            frame_dom = frame._domain
            for i1 in self._manifold.irange():
                for j1 in self._manifold.irange():
                    name = self._name + " curvature 2-form (" + str(i1) + \
                           "," + str(j1) + ")"
                    latex_name = r"\Omega^" + str(i1) + r"_{\ \, " + \
                                str(j1) + "}"
                    omega = frame_dom.diff_form(2, name=name, 
                                                latex_name=latex_name)
                    comega = omega.set_comp(frame)
                    for k in self._manifold.irange():
                        for l in self._manifold.irange(start=k+1):
                            comega[k,l] = \
                                        self.riemann().comp(frame)[[i1,j1,k,l]]
                    forms[(i1,j1)] = omega
            self._curvature_forms[frame] = forms
        return  self._curvature_forms[frame][(i,j)] 

            
#******************************************************************************
            
class LeviCivitaConnection(AffConnection):
    r"""
    Levi-Civita connection on a pseudo-Riemannian manifold.

    Given a differentiable manifold `M` endowed with a pseudo-Riemannian 
    metric `g` and denoting by `\mathcal{X}(M)` 
    the `C^\infty(M)`-module of vector fields on `M`, the 
    *Levi-Civita connection associated with* `g` is the unique operator

    .. MATH::

        \begin{array}{cccc}
        \nabla: & \mathcal{X}(M)\times \mathcal{X}(M) & \longrightarrow & 
                 \mathcal{X}(M) \\
                & (u,v) & \longmapsto & \nabla_u v
        \end{array}
                
    that
    
    - is a bilinear when considering `\mathcal{X}(M)` as a
      vector space over `\RR`
    - is `C^\infty(M)`-linear w.r.t. the first argument:
      `\forall f\in C^\infty(M),\ \nabla_{fu} v = f\nabla_u v` 
    - obeys Leibniz rule w.r.t. the second argument:
      `\forall f\in C^\infty(M),\ \nabla_u (f v) = \mathrm{d}f(u)\, v + f  \nabla_u v`
    - is torsion-free
    - is compatible with `g`: 
      `\forall (u,v,w)\in \mathcal{X}(M)^3,\ u(g(v,w)) = g(\nabla_u v, w) + g(v, \nabla_u w)`
      

    The Levi-Civita connection `\nabla` gives birth to the *covariant derivative
    operator* acting on tensor fields, denoted by the same symbol:
    
    .. MATH::

        \begin{array}{cccc}
        \nabla: &  T^{(k,l)}(M) & \longrightarrow & T^{(k,l+1)}(M)\\
                & t & \longmapsto & \nabla t
        \end{array}
        
    where `T^{(k,l)}(M)` stands for the `C^\infty(M)`-module of tensor fields
    of type `(k,l)` on `M` (with the convention `T^{(0,0)}(M):=C^\infty(M)`).
    For a vector field `v`,  the covariant derivative `\nabla v` is a 
    type-(1,1) tensor field such that 
    
    .. MATH::
    
        \forall u \in\mathcal{X}(M), \   \nabla_u v = \nabla v(., u) 
        
    More generally for any tensor field `t\in T^{(k,l)}(M)`, we have

    .. MATH::
    
        \forall u \in\mathcal{X}(M), \   \nabla_u t = \nabla t(\ldots, u) 
        
    
    .. NOTE::
    
        The above convention means that, in terms of index notation, 
        the "derivation index" in `\nabla t` is the *last* one:
        
        .. MATH::
        
            \nabla_c t^{a_1\ldots a_k}_{\quad\quad b_1\ldots b_l} =
                (\nabla t)^{a_1\ldots a_k}_{\quad\quad b_1\ldots b_l c}


    INPUT:
    
    - ``metric`` -- the metric `g` defining the Levi-Civita connection, as an
      instance of class :class:`~sage.geometry.manifolds.metric.Metric`
    - ``name`` -- name given to the connection
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the  
      connection
    - ``init_coef`` -- (default: True) determines whether the Chrsitoffel 
      symbols are initialized (in the top charts on the domain, i.e. 
      disregarding the subcharts)

    EXAMPLES:
    
    Levi-Civita connection associated with the Euclidean metric on `\RR^3`
    expressed in spherical coordinates::
    
        sage: M = Manifold(3, 'R^3', start_index=1)
        sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: g = M.metric('g')
        sage: g[1,1], g[2,2], g[3,3] = 1, r^2 , (r*sin(th))^2
        sage: g.view()
        g = dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
        sage: from sage.geometry.manifolds.connection import LeviCivitaConnection
        sage: nab = LeviCivitaConnection(g, 'nabla', r'\nabla') ; nab
        Levi-Civita connection 'nabla' associated with the pseudo-Riemannian metric 'g' on the 3-dimensional manifold 'R^3'
    
    Let us check that the connection is compatible with the metric::
    
        sage: Dg = nab(g) ; Dg
        tensor field 'nabla g' of type (0,3) on the 3-dimensional manifold 'R^3'
        sage: Dg == 0
        True

    and that it is torsionless::
    
        sage: nab.torsion() == 0
        True
        sage: sage.geometry.manifolds.connection.AffConnection.torsion(nab) == 0  # forces the computation of the torsion
        True
 
    The connection coefficients in the manifold's default frame are Christoffel 
    symbols, since the default frame is a coordinate frame::
    
        sage: M.default_frame()
        coordinate frame (R^3, (d/dr,d/dth,d/dph))
        sage: nab.coef()
        3-indices components w.r.t. coordinate frame (R^3, (d/dr,d/dth,d/dph)), with symmetry on the index positions (1, 2)
        sage: # note that the Christoffel symbols are symmetric with respect to their last two indices (positions (1,2))
        sage: nab.coef()[:]
        [[[0, 0, 0], [0, -r, 0], [0, 0, -r*sin(th)^2]], 
        [[0, 1/r, 0], [1/r, 0, 0], [0, 0, -cos(th)*sin(th)]], 
        [[0, 0, 1/r], [0, 0, cos(th)/sin(th)], [1/r, cos(th)/sin(th), 0]]]

    """
    def __init__(self, metric, name, latex_name=None, init_coef=True):
        AffConnection.__init__(self, metric._domain, name, latex_name)
        self._metric = metric
        # Initialization of the derived quantities:
        LeviCivitaConnection._init_derived(self)
        if init_coef:
            # Initialization of the Christoffel symbols in the top charts on
            # the domain (i.e. disregarding the subcharts)
            for chart in self._domain._atlas:
                for other in self._domain._atlas:
                    if chart in other._subcharts and chart is not other:
                        break 
                else:
                    # the chart is not any subchart and therefore the
                    # Christoffel symbols relative to it are computed:
                    self.coef(chart._frame)
        
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "Levi-Civita connection"
        if self._name is not None:
            description += " '%s'" % self._name
        description += " associated with the " + str(self._metric)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        AffConnection._init_derived(self)

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        AffConnection._del_derived(self)

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some subdomain.
        
        If such restriction has not been defined yet, it is constructed here.

        INPUT:
        
        - ``subdomain`` -- open subset `U` of ``self._domain`` (must be an 
          instance of :class:`~sage.geometry.manifolds.domain.OpenDomain`)
          
        OUTPUT:
        
        - instance of :class:`LeviCivitaConnection` representing the restriction.

        EXAMPLE:
        
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subdomain(self._domain):
                raise ValueError("The provided domain is not a subdomain of " + 
                                 "the current connection's domain.")
            resu = LeviCivitaConnection(self._metric.restrict(subdomain), 
                                        name=self._name, 
                                        latex_name=self._latex_name,
                                        init_coef=False)
            for frame in self._coefficients:
                for sframe in subdomain._top_frames: 
                    if sframe in frame._subframes:
                        comp_store = self._coefficients[frame]._comp
                        scoef = resu._new_coef(sframe)
                        scomp_store = scoef._comp
                        # the coefficients of the restriction are evaluated 
                        # index by index:
                        for ind, value in comp_store.iteritems():
                            scomp_store[ind] = value.restrict(sframe._domain)
                        resu._coefficients[sframe] = scoef
            if self._riemann is not None:
                resu._riemann = self._riemann.restrict(subdomain)
            if self._ricci is not None:
                resu._ricci = self._ricci.restrict(subdomain)
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def _new_coef(self, frame): 
        r"""
        Create the connection coefficients w.r.t. the given frame. 
                
        """
        from sage.tensor.modules.comp import Components, CompWithSym
        from scalarfield import ScalarField
        from vectorframe import CoordFrame
        if isinstance(frame, CoordFrame):
            # the Christoffel symbols are symmetric:
            return CompWithSym(frame._domain.scalar_field_algebra(), frame, 3, 
                               start_index=self._manifold._sindex,
                               output_formatter=ScalarField.function_chart,
                               sym=(1,2))
        else:
            # a priori no symmetry in a generic frame:
            return Components(frame._domain.scalar_field_algebra(), frame, 3, 
                              start_index=self._manifold._sindex,
                              output_formatter=ScalarField.function_chart)


    def coef(self, frame=None):
        r"""
        Return the connection coefficients relative to the given frame.
        
        `n` being the manifold's dimension, the connection coefficients 
        relative to the vector frame `(e_i)` are the `n^3` scalar fields 
        `\Gamma^k_{\ \, ij}` defined by 
        
        .. MATH::
            
            \nabla_{e_j} e_i = \Gamma^k_{\ \, ij} e_k
        
        If the connection coefficients are not known already, they are computed

         * as Christoffel symbols if the frame `(e_i)` is a coordinate frame
         * frome the above formula otherwise 
                
        INPUT:
        
        - ``frame`` -- (default: None) vector frame relative to which the 
          connection coefficients are required; if none is provided, the 
          domain's default frame is assumed
 
        OUTPUT: 
        
        - connection coefficients relative to the frame ``frame``, as an 
          instance of the class :class:`~sage.tensor.modules.comp.Components` 
          with 3 indices ordered as `(k,i,j)`; for Christoffel symbols, 
          an instance of the subclass
          :class:`~sage.tensor.modules.comp.CompWithSym` is returned. 
        
        EXAMPLES:
        
        Christoffel symbols of the Levi-Civita connection associated to 
        the Euclidean metric on `\RR^3` expressed in spherical coordinates::
        
            sage: M = Manifold(3, 'R^3', start_index=1)
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2 , (r*sin(th))^2
            sage: g.view()
            g = dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
            sage: from sage.geometry.manifolds.connection import LeviCivitaConnection
            sage: nab = LeviCivitaConnection(g, 'nabla', r'\nabla')
            sage: gam = nab.coef() ; gam
            3-indices components w.r.t. coordinate frame (R^3, (d/dr,d/dth,d/dph)), with symmetry on the index positions (1, 2)
            sage: gam[:]
            [[[0, 0, 0], [0, -r, 0], [0, 0, -r*sin(th)^2]], 
            [[0, 1/r, 0], [1/r, 0, 0], [0, 0, -cos(th)*sin(th)]], 
            [[0, 0, 1/r], [0, 0, cos(th)/sin(th)], [1/r, cos(th)/sin(th), 0]]]
            sage: # The only non-zero Christoffel symbols:
            sage: gam[1,2,2], gam[1,3,3]
            (-r, -r*sin(th)^2)
            sage: gam[2,1,2], gam[2,3,3]
            (1/r, -cos(th)*sin(th))
            sage: gam[3,1,3], gam[3,2,3]
            (1/r, cos(th)/sin(th))
            
        Connection coefficients of the same connection with respect to the 
        orthonormal frame associated to spherical coordinates::
        
            sage: ch_basis = M.automorphism_field() 
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = 1, 1/r, 1/(r*sin(th))
            sage: e = c_spher.frame().new_frame(ch_basis, 'e')
            sage: gam_e = nab.coef(e) ; gam_e
            3-indices components w.r.t. vector frame (R^3, (e_1,e_2,e_3))
            sage: gam_e[:]
            [[[0, 0, 0], [0, -1/r, 0], [0, 0, -1/r]],
            [[0, 1/r, 0], [0, 0, 0], [0, 0, -cos(th)/(r*sin(th))]],
            [[0, 0, 1/r], [0, 0, cos(th)/(r*sin(th))], [0, 0, 0]]]
            sage: # The only non-zero connection coefficients:
            sage: gam_e[1,2,2], gam_e[2,1,2]
            (-1/r, 1/r)
            sage: gam_e[1,3,3], gam_e[3,1,3]
            (-1/r, 1/r)
            sage: gam_e[2,3,3], gam_e[3,2,3]
            (-cos(th)/(r*sin(th)), cos(th)/(r*sin(th)))

        """
        from scalarfield import ScalarField
        from vectorframe import CoordFrame
        if frame is None: 
            frame = self._domain._def_frame
        if frame not in self._coefficients:
            # the coefficients must be computed
            #
            # Check whether frame is a subframe of a frame in which the 
            # coefficients are already known:
            for oframe in self._coefficients:
                if frame in oframe._subframes:
                    self._coefficients[frame] = self._new_coef(frame)
                    comp_store = self._coefficients[frame]._comp
                    ocomp_store = self._coefficients[oframe]._comp
                    for ind, value in ocomp_store.iteritems():
                        comp_store[ind] = value.restrict(frame._domain)
                    break
            else:
                # If not, the coefficients must be computed from scratch:
                manif = self._manifold
                if isinstance(frame, CoordFrame):
                    # Christoffel symbols
                    chart = frame._chart
                    gam = self._new_coef(frame)
                    gg = self._metric.comp(frame)
                    ginv = self._metric.inverse().comp(frame)
                    for ind in gam.non_redundant_index_generator():
                        i, j, k = ind
                        # The computation is performed at the FunctionChart level:
                        rsum = 0
                        for s in manif.irange():
                            rsum += ginv[i,s, chart] * ( 
                                                gg[s,k, chart].diff(j)
                                              + gg[j,s, chart].diff(k)
                                              - gg[j,k, chart].diff(s) )
                        gam[i,j,k, chart] = rsum / 2
                        self._coefficients[frame] = gam
                else:
                    # Computation from the formula defining the connection coef.
                    return AffConnection.coef(self, frame)
        return self._coefficients[frame]

    def torsion(self):
        r""" 
        Return the connection's torsion tensor (identically zero for a 
        Levi-Civita connection). 
        
        See :meth:`AffConnection.torsion` for the general definition of the 
        torsion tensor. 
        
        OUTPUT:
        
        - the torsion tensor `T`, as a vanishing instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          
        """
        if self._torsion is None:
            manif = self._manifold
            resu = self._domain.tensor_field(1, 2, antisym=(1,2))
            for frame in self._coefficients:
                # Initialization of the frame components to zero: 
                resu.add_comp(frame) 
            self._torsion = resu
        return self._torsion 

    def riemann(self, name=None, latex_name=None):
        r""" 
        Return the Riemann curvature tensor associated with the metric.

        This method redefines :meth:`AffConnection.riemann` to set some name
        and the latex_name to the output.
        
        The Riemann curvature tensor is the tensor field `R` of type (1,3) 
        defined by

        .. MATH::
            
            R(\omega, u, v, w) = \left\langle \omega, \nabla_u \nabla_v w
                - \nabla_v \nabla_u w - \nabla_{[u, v]} w \right\rangle
        
        for any 1-form  `\omega`  and any vector fields `u`, `v` and `w`. 

        INPUT:
        
        - ``name`` -- (default: None) name given to the Riemann tensor; 
          if none, it is set to "Riem(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
          Riemann tensor; if none, it is set to "\\mathrm{Riem}(g)", where "g" 
          is the metric's name

        OUTPUT:
        
        - the Riemann curvature tensor `R`, as an instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          
        """
        if self._riemann is None:
            AffConnection.riemann(self)
            if name is None:
                self._riemann._name = "Riem(" + self._metric._name + ")"
            else:
                self._riemann._name = name
            if latex_name is None:
                self._riemann._latex_name = r"\mathrm{Riem}\left(" + \
                                           self._metric._latex_name + r"\right)"
            else:
                self._riemann._latex_name = latex_name
            for rst in self._riemann._restrictions.itervalues():
                rst._name = self._riemann._name
                rst._latex_name = self._riemann._latex_name
        return self._riemann


    def ricci(self, name=None, latex_name=None):
        r""" 
        Return the connection's Ricci tensor.
        
        This method redefines :meth:`AffConnection.ricci` to take into account
        the symmetry of the Ricci tensor for a Levi-Civita connection. 

        The Ricci tensor is the tensor field `Ric` of type (0,2) 
        defined from the Riemann curvature tensor `R` by 

        .. MATH::
            
            Ric(u, v) = R(e^i, u, e_i, v)
        
        for any vector fields `u` and `v`, `(e_i)` being any vector frame and
        `(e^i)` the dual coframe. 
                
        INPUT:
        
        - ``name`` -- (default: None) name given to the Ricci tensor; 
          if none, it is set to "Ric(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
          Ricci tensor; if none, it is set to "\\mathrm{Ric}(g)", where "g" 
          is the metric's name
          
        OUTPUT:
        
        - the Ricci tensor `Ric`, as an instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField` of tensor
          type (0,2) and symmetric
        
        EXAMPLES:
        
        Ricci tensor of the standard connection on the 2-dimensional sphere::
        
            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: c_spher.<th,ph> = M.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2] = 1, sin(th)^2
            sage: g.view() # standard metric on S^2:
            g = dth*dth + sin(th)^2 dph*dph
            sage: nab = g.connection() ; nab
            Levi-Civita connection 'nabla_g' associated with the pseudo-Riemannian metric 'g' on the 2-dimensional manifold 'S^2'
            sage: ric = nab.ricci() ; ric             
            field of symmetric bilinear forms 'Ric(g)' on the 2-dimensional manifold 'S^2'
            sage: ric.view()
            Ric(g) = dth*dth + sin(th)^2 dph*dph            
        
        Checking that the Ricci tensor of the Levi-Civita connection associated
        to Schwarzschild metric is identically zero::
        
            sage: M = Manifold(4, 'M')
            sage: c_BL.<t,r,th,ph> = M.chart(r't r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi') # Boyer-Linquist coordinates
            sage: g = M.metric('g')
            sage: m = var('m')  # mass in Schwarzschild metric
            sage: g[0,0], g[1,1] = -(1-2*m/r), 1/(1-2*m/r)
            sage: g[2,2], g[3,3] = r^2, (r*sin(th))^2
            sage: g.view()
            g = (2*m - r)/r dt*dt - r/(2*m - r) dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
            sage: nab = g.connection() ; nab
            Levi-Civita connection 'nabla_g' associated with the pseudo-Riemannian metric 'g' on the 4-dimensional manifold 'M'
            sage: ric = nab.ricci() ; ric
            field of symmetric bilinear forms 'Ric(g)' on the 4-dimensional manifold 'M'
            sage: ric == 0
            True

        """
        if self._ricci is None:
            manif = self._manifold
            riem = self.riemann()
            resu = self._domain.tensor_field(0,2, sym=(0,1))
            for frame in self._coefficients:
                cric = resu.add_comp(frame)
                criem = riem.comp(frame)
                for i in manif.irange():
                    # symmetry of the Ricci tensor taken into account by j>=i: 
                    for j in manif.irange(start=i):  
                        rsum = 0
                        for k in manif.irange():
                            rsum += criem[[k,i,k,j]]
                        cric[i,j] = rsum
            if name is None:
                resu._name = "Ric(" + self._metric._name + ")"
            else:
                resu._name = name
            if latex_name is None:
                resu._latex_name = r"\mathrm{Ric}\left(" + \
                                         self._metric._latex_name + r"\right)"
            else:
                resu._latex_name = latex_name
            for rst in resu._restrictions.itervalues():
                rst._name = resu._name
                rst._latex_name = resu._latex_name
            self._ricci = resu
        return self._ricci 


