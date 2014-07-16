r"""
Tensor fields of type (1,1)

Three derived classes of 
:class:`~sage.geometry.manifolds.tensorfield.TensorField` 
devoted to type-(1,1) tensor fields are implemented:


* :class:`EndomorphismField` for fields of endomorphisms 
  (type (1,1) tensor fields)

  * :class:`AutomorphismField` for fields of invertible endomorphisms

    * :class:`TangentIdentityField` for fields of identity maps on tangent 
      spaces


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014): initial version

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

from sage.tensor.modules.free_module_tensor_spec import \
    FreeModuleEndomorphism, FreeModuleAutomorphism, FreeModuleIdentityMap
from tensorfield import TensorField, TensorFieldParal


class EndomorphismField(TensorField):
    r"""
    Field of tangent-space endomorphisms with values in an open
    subset of a differentiable manifold. 
    
    An instance of this class is a field of endomorphisms (i.e. linear 
    operators in each tangent space) along an open subset `U` of some immersed 
    submanifold `S` of a manifold `M` with values in an open 
    subset `V` of `M`. 
    The standard case of a field of endomorphisms *on* a manifold corresponds 
    to `U=V` (and hence `S=M`).

    If `V` is parallelizable, the class :class:`EndomorphismFieldParal` must be 
    used instead.
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,V)` of vector 
      fields along `U` with values on `V`
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:


    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        TensorField.__init__(self, vector_field_module, (1,1), name=name, 
                             latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "field of endomorphisms "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.
        
        """
        return self.__class__(self._vmodule)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`TensorField.__call__` to allow for a single 
        argument (module element). 
        """
        if len(arg) > 1:
            # the endomorphism acting as a type (1,1) tensor on a pair 
            # (1-form, vector field), returning a scalar field:
            return TensorField.__call__(self, *arg) 
        # the endomorphism acting as such, on a vector field, returning a
        # vector field:
        vector = arg[0]
        if vector._tensor_type != (1,0):
            raise TypeError("The argument must be a vector field.")
        dom_resu = self._domain.intersection(vector._domain)
        if dom_resu.is_manifestly_parallelizable():
            # call of the EndomorphismFieldParal version:
            return self.restrict(dom_resu)(vector.restrict(dom_resu))
        if self._name is not None and vector._name is not None:
            name_resu = self._name + "(" + vector._name + ")"
        else:
            name_resu = None
        if self._latex_name is not None and vector._latex_name is not None:
            latex_name_resu = self._latex_name + r"\left(" + \
                              vector._latex_name + r"\right)"
        else:
            latex_name_resu = None
        dest_map = vector._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu)
        resu = dom_resu.vector_field(name=name_resu, 
                                     latex_name=latex_name_resu,
                                     dest_map=dest_map_resu)
        for dom in self._common_subdomains(vector):
            if dom.is_subdomain(dom_resu):
                resu._restrictions[dom] = \
                    self._restrictions[dom](vector._restrictions[dom])
        return resu

#******************************************************************************

class AutomorphismField(EndomorphismField):
    r"""
    Field of tangent-space automorphisms with values on a open 
    subset of a differentiable manifold. 
    
    An instance of this class is a field of linear automorphisms (i.e. 
    invertible linear operators in each tangent space) along an open subset 
    `U` of some immersed submanifold `S` of a manifold `M` with values in an 
    open subset `V` of `M`. 
    The standard case of a field of automorphisms *on* a manifold corresponds 
    to `U=V` (and hence `S=M`).
    
    If `V` is parallelizable, the class :class:`AutomorphismFieldParal` must be 
    used instead.

    INPUT:
    
    - ``vector_field_module`` -- module `\mathcal{X}(U,V)` of vector 
      fields along `U` with values on `V`
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLE: 
    
    Field of tangent-space automorphisms on a non-parallelizable 2-dimensional 
    manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U') ; V = M.open_domain('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: a = M.automorphism_field('a') ; a
        field of tangent-space automorphisms 'a' on the 2-dimensional manifold 'M'
        sage: a.parent()
        module T^(1,1)(M) of type-(1,1) tensors fields on the 2-dimensional manifold 'M'

    We first define the components of `a` w.r.t the coordinate frame on `U`::

        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: a[eU,:] = [[1,x], [0,2]]

    We then set the components w.r.t. the coordinate frame on `V` by extending 
    the expressions of the components in the corresponding subframe on 
    `W = U\cap V`::
    
        sage: W = U.intersection(V)
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        
    At this stage, the automorphims field `a` is fully defined::
    
        sage: a.view(eU)
        a = d/dx*dx + x d/dx*dy + 2 d/dy*dy
        sage: a.view(eV)
        a = (1/4*u + 1/4*v + 3/2) d/du*du + (-1/4*u - 1/4*v - 1/2) d/du*dv + (1/4*u + 1/4*v - 1/2) d/dv*du + (-1/4*u - 1/4*v + 3/2) d/dv*dv
    
    In particular, we may ask for its inverse on the whole manifold `M`::
    
        sage: ia = a.inverse() ; ia
        field of tangent-space automorphisms 'a^(-1)' on the 2-dimensional manifold 'M'
        sage: ia.view(eU)
        a^(-1) = d/dx*dx - 1/2*x d/dx*dy + 1/2 d/dy*dy
        sage: ia.view(eV)
        a^(-1) = (-1/8*u - 1/8*v + 3/4) d/du*du + (1/8*u + 1/8*v + 1/4) d/du*dv + (-1/8*u - 1/8*v + 1/4) d/dv*du + (1/8*u + 1/8*v + 3/4) d/dv*dv

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        EndomorphismField.__init__(self, vector_field_module, name=name, 
                                        latex_name=latex_name)
        self._init_derived() # initialization of derived quantities

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "field of tangent-space automorphisms "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)
        
    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        TensorField._init_derived(self)
        self._inverse = None  # inverse not set yet
        
    def _del_derived(self):
        r"""
        Delete the derived quantities.
        """
        # First delete the derived quantities pertaining to the mother class:
        TensorField._del_derived(self)
        # then deletes the inverse automorphism:
        self._inverse = None
        
    def inverse(self):
        r"""
        Return the inverse automorphism.
        
        EXAMPLE: 
        
        Inverse of a field of tangent-space automorphisms on a 
        non-parallelizable 2-dimensional manifold::
    
            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U') ; V = M.open_domain('V') 
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: a = M.automorphism_field('a')
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a[eU,:] = [[1,x], [0,2]]
            sage: W = U.intersection(V)
            sage: a.add_comp_by_continuation(eV, W, c_uv)
            sage: ia = a.inverse() ; ia
            field of tangent-space automorphisms 'a^(-1)' on the 2-dimensional manifold 'M'
            sage: a[eU,:], ia[eU,:]
            (
            [1 x]  [     1 -1/2*x]
            [0 2], [     0    1/2]
            )
            sage: a[eV,:], ia[eV,:]
            (
            [ 1/4*u + 1/4*v + 3/2 -1/4*u - 1/4*v - 1/2]
            [ 1/4*u + 1/4*v - 1/2 -1/4*u - 1/4*v + 3/2],
            [-1/8*u - 1/8*v + 3/4  1/8*u + 1/8*v + 1/4]
            [-1/8*u - 1/8*v + 1/4  1/8*u + 1/8*v + 3/4]
            )
        
        Let us check that ia is indeed the inverse of a::
        
            sage: s = a.contract(ia)
            sage: s[eU,:], s[eV,:]
            (
            [1 0]  [1 0]
            [0 1], [0 1]
            )
            sage: s = ia.contract(a)
            sage: s[eU,:], s[eV,:]
            (
            [1 0]  [1 0]
            [0 1], [0 1]
            )
            
        """        
        if self._inverse is None:
            if self._name is None:
                inv_name = None
            else:
                inv_name = self._name  + '^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                inv_latex_name = self._latex_name + r'^{-1}'
            self._inverse = AutomorphismField(self._vmodule, name=inv_name, 
                                              latex_name=inv_latex_name)
            for dom, rst in self._restrictions.iteritems():
                self._inverse._restrictions[dom] = rst.inverse()
        return self._inverse


#******************************************************************************

class TangentIdentityField(AutomorphismField):
    r"""
    Field of tangent-space identity maps with values on an open subset of a 
    differentiable manifold. 
    
    An instance of this class is a field of identity maps (i.e. identity 
    operator in each tangent space) along an open subset `U` of some immersed 
    submanifold `S` of a manifold `M` with values in an open 
    subset `V` of `M`. 
    The standard case of a field of identity maps *on* a manifold corresponds 
    to `U=V` (and hence `S=M`).
    
    If `V` is parallelizable, the class :class:`TangentIdentityFieldParal` must
    be used instead.

    INPUT:
    
    - ``vector_field_module`` -- module `\mathcal{X}(U,V)` of vector 
      fields along `U` with values on `V`
    - ``name`` -- (default: None) name given to the identity map; if none
      is provided, the value 'Id' is set. 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the identity
      map; if none is provided, the LaTeX symbol is set to `\mathrm{Id}`

    EXAMPLES:

    Field of tangent-space identity maps on a non-parallelizable 2-dimensional 
    manifold::
    
        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U') ; V = M.open_domain('V') 
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: a = M.tangent_identity_field() ; a
        field of tangent-space identity maps 'Id' on the 2-dimensional manifold 'M'
        sage: a.parent()
        module T^(1,1)(M) of type-(1,1) tensors fields on the 2-dimensional manifold 'M'
        
    The components in any frame on M are Kronecker deltas::
    
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: a[eU,:]
        [1 0]
        [0 1]
        sage: a[eV,:]
        [1 0]
        [0 1]
        sage: a[eU.restrict(W),:]
        [1 0]
        [0 1]
        sage: a[eV.restrict(W),:]
        [1 0]
        [0 1]
    
    The identity is its own inverse::
        
        sage: a.inverse() is a
        True
    
    The identity map acting on a vector::
    
        sage: v = M.vector_field('v')
        sage: v[eU,:] = [1-y, x*y]
        sage: v.add_comp_by_continuation(eV, W, c_uv)
        sage: a(v)
        vector field 'v' on the 2-dimensional manifold 'M'
        sage: a(v) is v
        True
        
    When the domains of the identity field and the vector fields are different, 
    their intersection is used for the result::
    
        sage: a(v.restrict(U))
        vector field 'v' on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: a(v.restrict(U)) is v.restrict(U)
        True
    
    ::
    
        sage: a.restrict(U)(v)
        vector field 'v' on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: a.restrict(U)(v) is v.restrict(U)
        True
   
    """
    def __init__(self, vector_field_module, name='Id', latex_name=None):
        if latex_name is None and name == 'Id':
            latex_name = r'\mathrm{Id}'
        AutomorphismField.__init__(self, vector_field_module, name=name, 
                                       latex_name=latex_name)
        for dom in self._domain._subdomains:
            if dom.is_manifestly_parallelizable():
                fmodule = dom.vector_field_module()
                self._restrictions[dom] = TangentIdentityFieldParal(fmodule, 
                                              name=name, latex_name=latex_name)
        self._inverse = self    # the identity is its own inverse
        
    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "field of tangent-space identity maps "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`EndomorphismField.__call__`.
        """
        if len(arg) == 1:
            # The identity map acting as such, on a vector field:
            vector = arg[0]
            if vector._tensor_type != (1,0):
                raise TypeError("The argument must be a vector field.")
            dom = self._domain.intersection(vector._domain)
            return vector.restrict(dom)
        elif len(arg) == 2:
            # self acting as a type-(1,1) tensor on a pair 
            # (1-form, vector field), returning a scalar field:
            oneform = arg[0]
            vector = arg[1]
            dom = self._domain.intersection(
                            oneform._domain).intersection(vector._domain)
            return oneform.restrict(dom)(vector.restrict(dom))
        else:
            raise TypeError("Wrong number of arguments.")

#******************************************************************************

class EndomorphismFieldParal(FreeModuleEndomorphism, TensorFieldParal):
    r"""
    Field of tangent-space endomorphisms with values in a parallelizable open 
    subset of a differentiable manifold. 
    
    An instance of this class is a field of endomorphisms (i.e. linear 
    operators in each tangent space) along an open subset `U` of some immersed 
    submanifold `S` of a manifold `M` with values in a parallelizable open 
    subset `V` of `M`. 
    The standard case of a field of endomorphisms *on* a manifold corresponds 
    to `U=V` (and hence `S=M`).
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,V)` of vector 
      fields along `U` with values on `V`
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A field of endomorphisms on a 3-dimensional manifold::
    
        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: t = M.endomorphism_field('T') ; t
        field of endomorphisms 'T' on the 3-dimensional manifold 'M'
        
    A field of endomorphisms is a tensor field of rank 2 and of type (1,1)::
    
        sage: t.parent()
        free module T^(1,1)(M) of type-(1,1) tensors fields on the 3-dimensional manifold 'M'
        sage: t._tensor_rank
        2
        sage: t._tensor_type
        (1, 1)
    
    Components with respect to a given frame::

        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: t[:] = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: t[:]
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: change_basis = M.automorphism_field()
        sage: change_basis[:] = [[-1,0,1], [0,1,2], [-2,1,2]]
        sage: f = e.new_frame(change_basis, 'f')
        sage: t.comp(f) # computation of t components in the frame f, from those in the frame e
        2-indices components w.r.t. vector frame (M, (f_1,f_2,f_3))
        sage: t.comp(f)[:]
        [  9/2    -3 -15/2]
        [  -11     7    19]
        [ -5/2     2   7/2]
        sage: t.comp(f)[1,1]
        9/2

    An endomorphism maps a vector to a vector::
    
        sage: v = M.vector_field('v')
        sage: v[:] = (1,2,3)
        sage: w = t(v) ; w
        vector field 'T(v)' on the 3-dimensional manifold 'M'
        sage: w[:]
        [14, 32, 50]
        sage: t[:] * matrix([v[i].expr() for i in range(1,4)]).transpose()  # check:
        [14]
        [32]
        [50]
        sage: latex(t(v))
        T\left(v\right)

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        FreeModuleEndomorphism.__init__(self, vector_field_module, name=name, 
                                        latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self) 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "field of endomorphisms "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)
        
    def _new_instance(self):
        r"""
        Create an instance if the same type as ``self`` on the same domain.
        """
        return self.__class__(self._fmodule)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities
        
        INPUT:
        
        - ``del_restrictions`` -- (default: True) determines whether the
          restrictions of ``self`` to subdomains are deleted. 
        
        """
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)

    def __call__(self, *arg):
        r"""
        Redefinition of 
        :meth:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleEndomorphism.__call__` 
        to allow for domain treatment
        """
        if len(arg) > 1:
            # the endomorphism acting as a type (1,1) tensor on a pair 
            # (linear form, module element), returning a scalar:
            return TensorFieldParal.__call__(self, *arg)
        else:
            vector = arg[0]
            dom = self._domain.intersection(vector._domain)
            return FreeModuleEndomorphism.__call__(self.restrict(dom), 
                                                   vector.restrict(dom))

#******************************************************************************

class AutomorphismFieldParal(FreeModuleAutomorphism, EndomorphismFieldParal):
    r"""
    Field of tangent-space automorphisms with values on a parallelizable open 
    subset of a differentiable manifold. 
    
    An instance of this class is a field of linear automorphisms (i.e. 
    invertible linear operators in each tangent space) along an open subset `U` 
    of some immersed submanifold `S` of a manifold `M` with values in a 
    parallelizable open subset `V` of `M`. 
    The standard case of a field of automorphisms *on* a manifold corresponds 
    to `U=V` (and hence `S=M`).
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,V)` of vector 
      fields along `U` with values on `V`
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A pi/3-rotation in the Euclidean 2-plane::
    
        sage: M = Manifold(2,'R^2')
        sage: c_xy.<x,y> = M.chart()
        sage: rot = M.automorphism_field('R') ; rot
        field of tangent-space automorphisms 'R' on the 2-dimensional manifold 'R^2'
        sage: rot[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
        
    An automorphism is a type-(1,1) tensor::
    
        sage: rot.parent()
        free module T^(1,1)(R^2) of type-(1,1) tensors fields on the 2-dimensional manifold 'R^2'
    
    The inverse automorphism is obtained via the method :meth:`inverse`::
    
        sage: inv = rot.inverse() ; inv
        field of tangent-space automorphisms 'R^(-1)' on the 2-dimensional manifold 'R^2'
        sage: latex(inv)
        R^{-1}
        sage: inv[:]
        [1/2*sqrt(3)         1/2]
        [       -1/2 1/2*sqrt(3)]
        sage: rot[:]
        [1/2*sqrt(3)        -1/2]
        [        1/2 1/2*sqrt(3)]
        sage: inv[:] * rot[:]  # check
        [1 0]
        [0 1]

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        FreeModuleAutomorphism.__init__(self, vector_field_module, name=name, 
                                        latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self) 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "field of tangent-space automorphisms "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)
        
    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities
        
        INPUT:
        
        - ``del_restrictions`` -- (default: True) determines whether the
          restrictions of ``self`` to subdomains are deleted. 
        
        """
        # First delete the derived quantities pertaining to the mother class:
        EndomorphismFieldParal._del_derived(self, 
                                            del_restrictions=del_restrictions)
        # then deletes the inverse automorphism:
        self._inverse = None
        
    def inverse(self):
        r"""
        Return the inverse automorphism.
        """        
        from sage.matrix.constructor import matrix
        from sage.tensor.modules.comp import Components
        from vectorframe import CoordFrame
        from utilities import simplify_chain
        if self._inverse is None:
            if self._name is None:
                inv_name = None
            else:
                inv_name = self._name  + '^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                inv_latex_name = self._latex_name + r'^{-1}'
            fmodule = self._fmodule
            si = fmodule._sindex ; nsi = fmodule._rank + si
            self._inverse = AutomorphismFieldParal(fmodule, name=inv_name, 
                                                   latex_name=inv_latex_name)
            for frame in self._components:
                if isinstance(frame, CoordFrame):
                    chart = frame._chart
                else:
                    chart = self._domain._def_chart #!# to be improved
                try:
                    mat_self = matrix(
                              [[self.comp(frame)[i, j, chart]._express
                              for j in range(si, nsi)] for i in range(si, nsi)])
                except (KeyError, ValueError):
                    continue
                mat_inv = mat_self.inverse()
                cinv = Components(fmodule._ring, frame, 2, start_index=si,
                                  output_formatter=fmodule._output_formatter)
                for i in range(si, nsi):
                    for j in range(si, nsi):
                        cinv[i, j] = {chart: simplify_chain(mat_inv[i-si,j-si])}
                self._inverse._components[frame] = cinv
        return self._inverse


#******************************************************************************

class TangentIdentityFieldParal(FreeModuleIdentityMap, AutomorphismFieldParal):
    r"""
    Field of tangent-space identity maps with values on a parallelizable open 
    subset of a differentiable manifold. 
    
    An instance of this class is a field of identity maps (i.e. identity 
    operator in each tangent space) along an open subset `U` of some immersed 
    submanifold `S` of a manifold `M` with values in a parallelizable open 
    subset `V` of `M`. 
    The standard case of a field of identity maps *on* a manifold corresponds 
    to `U=V` (and hence `S=M`).
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,V)` of vector 
      fields along `U` with values on `V`
    - ``name`` -- (default: None) name given to the identity map; if none
      is provided, the value 'Id' is set. 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the identity
      map; if none is provided, the LaTeX symbol is set to `\mathrm{Id}`

    EXAMPLES:

    Field of tangent-space identity maps on a 3-dimensional manifold::
    
        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: a = M.tangent_identity_field(); a
        field of tangent-space identity maps 'Id' on the 3-dimensional manifold 'M'
        sage: latex(a)
        \mathrm{Id}
        
    The tangent-space identity map is a type-(1,1) tensor::
    
        sage: a.parent()
        free module T^(1,1)(M) of type-(1,1) tensors fields on the 3-dimensional manifold 'M'
        sage: a[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: a.comp()
        Kronecker delta of size 3x3
        
    The components are automatically defined in any frame::
    
        sage: e = M.vector_frame('e')
        sage: a.comp(e) 
        Kronecker delta of size 3x3
        sage: a.comp(e)[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The components can be read, but cannot be set::
    
        sage: a[1,1]
        1
        sage: a[1,1] = 2
        Traceback (most recent call last):
        ...
        NotImplementedError: The components of the identity map cannot be changed.

    The tangent-space identity map applied to a vector field::
    
        sage: v = M.vector_field()
        sage: v[:] = (2*x, -3, y+z)
        sage: w = a(v) ; w
        vector field on the 3-dimensional manifold 'M'
        sage: w[:]
        [2*x, -3, y + z]
        sage: w is v  # the output is actually the vector v itself
        True

    The tangent-space identity map acting as a type (1,1) tensor on a pair (1-form, vector)::
    
        sage: om = M.one_form()
        sage: om[:] = (0, x*y, 2)
        sage: s = a(om, v) ; s
        scalar field on the 3-dimensional manifold 'M'
        sage: s == om(v)
        True
        
    The tangent-space identity map is its own inverse::
    
        sage: a.inverse() == a
        True
        sage: a.inverse() is a
        True
        
    """
    def __init__(self, vector_field_module, name='Id', latex_name=None):
        if latex_name is None and name == 'Id':
            latex_name = r'\mathrm{Id}'
        FreeModuleIdentityMap.__init__(self, vector_field_module, name=name, 
                                       latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self) 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "field of tangent-space identity maps "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities
        
        INPUT:
        
        - ``del_restrictions`` -- (default: True) determines whether the
          restrictions of ``self`` to subdomains are deleted. 
        
        """
        # AutomorphismFieldParal._del_derived is bypassed:
        EndomorphismFieldParal._del_derived(self, 
                                            del_restrictions=del_restrictions)

    def __call__(self, *arg):
        r"""
        Redefinition of 
        :meth:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleIdentityMap.__call__` 
        to allow for domain treatment
        """
        if len(arg) == 1:
            # The identity map acting as such, on a vector field:
            vector = arg[0]
            if vector._tensor_type != (1,0):
                raise TypeError("The argument must be a vector field.")
            dom = self._domain.intersection(vector._domain)
            return vector.restrict(dom)
        elif len(arg) == 2:
            # self acting as a type-(1,1) tensor on a pair 
            # (1-form, vector field), returning a scalar field:
            oneform = arg[0]
            vector = arg[1]
            dom = self._domain.intersection(
                            oneform._domain).intersection(vector._domain)
            return oneform.restrict(dom)(vector.restrict(dom))
        else:
            raise TypeError("Wrong number of arguments.")
