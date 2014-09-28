r"""
Tensor field modules

The set of tensor fields along an open subset `U` of some manifold `S`
with values in a open subset `V` of a manifold `M` (possibly `S=M` and `U=V`)
is a module over the algebra `C^\infty(U)` of differentiable scalar fields 
on `U`. It is a free module iff `V` is parallelizable.
Accordingly, two classes are devoted to tensor field modules:

- :class:`TensorFieldModule` for tensor fields with values in a generic (in 
  practice, not parallelizable) open set `V`
- :class:`TensorFieldFreeModule` for tensor fields with values in a 
  parallelizable open set `V`

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version
"""

#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.modules.module import Module
from sage.structure.unique_representation import UniqueRepresentation
from sage.tensor.modules.tensor_free_module import TensorFreeModule
from tensorfield import TensorField, TensorFieldParal

class TensorFieldModule(UniqueRepresentation, Module):
    r"""
    Module of tensor fields of a given type `(k,l)` along an open subset `U` 
    of some manifold `S` with values in a open subset `V` of 
    a manifold `M`.
    
    This is a module over `C^\infty(U)`, the ring (algebra) of differentiable 
    scalar fields on `U`. 
    
    The standard case of tensor fields *on* a manifold corresponds to 
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an 
    immersion.

    If `V` is parallelizable, the class :class:`TensorFieldFreeModule` should
    be used instead.
    
    INPUT:
    
    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`. 
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank and 
      `l` the covariant rank
    
    EXAMPLE:
    
    Module of type-(2,0) tensor fields on the 2-sphere::
    
        sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_domain('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_domain('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: T20 = M.tensor_field_module((2,0)) ; T20
        module T^(2,0)(M) of type-(2,0) tensors fields on the 2-dimensional manifold 'M'
        
    `T^{(2,0)}(M)` is a module over the algebra `C^\infty(M)`::
    
        sage: T20.category()
        Category of modules over algebra of scalar fields on the 2-dimensional manifold 'M'
        sage: T20.base_ring() is M.scalar_field_algebra()
        True
    
    `T^{(2,0)}(M)` is not a free module::
    
        sage: isinstance(T20, FiniteRankFreeModule)
        False

    because `M = S^2` is not parallelizable::
    
        sage: M.is_manifestly_parallelizable()
        False
        
    On the contrary, the module of type-(2,0) tensor fields on `U` is a free 
    module, since `U` is parallelizable (being a coordinate domain)::
    
        sage: T20U = U.tensor_field_module((2,0))
        sage: isinstance(T20U, FiniteRankFreeModule)
        True
        sage: U.is_manifestly_parallelizable()
        True

    The zero element::
    
        sage: z = T20.zero() ; z
        tensor field 'zero' of type (2,0) on the 2-dimensional manifold 'M'
        sage: z is T20(0)
        True
        sage: z[c_xy.frame(),:]
        [0 0]
        [0 0]
        sage: z[c_uv.frame(),:]
        [0 0]
        [0 0]

    The module `T^{(2,0)}(M)` coerces to any module of type-(2,0) tensor fields 
    defined on some subdomain of `M`, for instance `T^{(2,0)}(U)`::
    
        sage: T20U.has_coerce_map_from(T20)
        True

    The reverse is not true::
    
        sage: T20.has_coerce_map_from(T20U)
        False
        
    The coercion::
    
        sage: T20U.coerce_map_from(T20)
        Conversion map:
          From: module T^(2,0)(M) of type-(2,0) tensors fields on the 2-dimensional manifold 'M'
          To:   free module T^(2,0)(U) of type-(2,0) tensors fields on the open domain 'U' on the 2-dimensional manifold 'M'

    The conversion map is actually the *restriction* of tensor fields defined 
    on `M` to `U`. 

    """
    
    Element = TensorField

    def __init__(self, vector_field_module, tensor_type):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "T^(" + str(kcon) + "," + str(lcov) + ")(" + domain._name
        latex_name = r"\mathcal{T}^{(" + str(kcon) + "," + str(lcov) + r")}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")" 
            latex_name += r"\right)" 
        else:
            name += "," + dest_map._name + ")" 
            latex_name += "," + dest_map._latex_name + r"\right)" 
        self._vmodule = vector_field_module
        self._tensor_type = tensor_type
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to 
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra() 
        Module.__init__(self, self._ring)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no element 
        # can be constructed here, to avoid some infinite recursion. 

    #### Methods required for any Parent 

    def _element_constructor_(self, comp=[], frame=None, name=None, 
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor field
        """
        if comp == 0:
            if not hasattr(self, '_zero_element'):
                self._zero_element = self._element_constructor_(name='zero', 
                                                                latex_name='0')
                for frame in self._domain._frames:
                    if self._dest_map.restrict(frame._domain) == \
                                                               frame._dest_map:
                        self._zero_element.add_comp(frame)
                        # (since new components are initialized to zero)
            return self._zero_element
        if isinstance(comp, TensorField):
            if self._tensor_type == comp._tensor_type and \
               self._domain.is_subdomain(comp._domain) and \
               self._ambient_domain.is_subdomain(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the " + str(comp) +
                                "to a tensor field in " + str(self))
        resu = self.element_class(self._vmodule, self._tensor_type, name=name, 
                                  latex_name=latex_name, sym=sym, 
                                  antisym=antisym)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) tensor field
        """
        resu = self.element_class(self._vmodule, self._tensor_type)
        return resu
            
    #### End of methods required for any Parent 

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent
        """
        if isinstance(other, (TensorFieldModule, TensorFieldFreeModule)):
            return self._tensor_type == other._tensor_type and \
                   self._domain.is_subdomain(other._domain) and \
                   self._ambient_domain.is_subdomain(other._ambient_domain)
        else:
            return False

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "module "
        if self._name is not None:
            description += self._name + " "
        description += "of type-(%s,%s)" % \
                           (str(self._tensor_type[0]), str(self._tensor_type[1]))
        description += " tensors fields "
        if self._dest_map is self._domain._identity_map:
            description += "on the " + str(self._domain)
        else:
            description += "along the " + str(self._domain) + \
                           " mapped into the " + str(self._ambient_domain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def base_module(self):
        r""" 
        Return the vector field module on which ``self`` is constructed.
        
        OUTPUT:
        
        - instance of class 
          :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`
          representing the module on which the tensor module is defined. 
        
        """
        return self._vmodule


#******************************************************************************

class TensorFieldFreeModule(TensorFreeModule):
    r"""
    Module of tensor fields of a given type `(k,l)` along an open subset `U` 
    of some manifold `S` with values in a parallelizable open subset `V` of 
    a manifold `M`.
    
    Since `V` is parallelizable, the module is a free module over `C^\infty(U)`,
    the ring (algebra) of differentiable scalar fields on `U`. 
    
    The standard case of tensor fields *on* a manifold corresponds to 
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an 
    immersion.

    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`. 
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank and 
      `l` the covariant rank
      
    EXAMPLE:
    
    Module of type-(2,0) tensor fields on `\RR^3`::
    
        sage: M = Manifold(3, 'R^3')
        sage: c_xyz.<x,y,z> = M.chart()  # Cartesian coordinates 
        sage: T20 = M.tensor_field_module((2,0)) ; T20
        free module T^(2,0)(R^3) of type-(2,0) tensors fields on the 3-dimensional manifold 'R^3'

    `T^{(2,0)}(\RR^3)` is a module over the algebra `C^\infty(\RR^3)`::

        sage: T20.category()
        Category of modules over algebra of scalar fields on the 3-dimensional manifold 'R^3'
        sage: T20.base_ring() is M.scalar_field_algebra()
        True
 
    `T^{(2,0)}(\RR^3)` is a free module::
    
        sage: isinstance(T20, FiniteRankFreeModule)
        True
        
    because `M = R^3` is parallelizable::
    
        sage: M.is_manifestly_parallelizable()
        True
    
    The zero element::
    
        sage: z = T20.zero() ; z
        tensor field 'zero' of type (2,0) on the 3-dimensional manifold 'R^3'
        sage: z[:]
        [0 0 0]
        [0 0 0]
        [0 0 0]

    A random element::
    
        sage: t = T20.an_element() ; t
        tensor field of type (2,0) on the 3-dimensional manifold 'R^3'
        sage: t[:]
        [2 0 0]
        [0 0 0]
        [0 0 0]

    The module `T^{(2,0)}(\RR^3)` coerces to any module of type-(2,0) tensor fields 
    defined on some subdomain of `\RR^3`::

        sage: U = M.open_domain('U', coord_def={c_xyz: x>0})
        sage: T20U = U.tensor_field_module((2,0))
        sage: T20U.has_coerce_map_from(T20)
        True
        sage: T20.has_coerce_map_from(T20U)  # the reverse is not true
        False
        sage: T20U.coerce_map_from(T20)
        Conversion map:
          From: free module T^(2,0)(R^3) of type-(2,0) tensors fields on the 3-dimensional manifold 'R^3'
          To:   free module T^(2,0)(U) of type-(2,0) tensors fields on the open domain 'U' on the 3-dimensional manifold 'R^3'
        
    The conversion map is actually the *restriction* of tensor fields defined 
    on `\RR^3` to `U`. 
    
    """

    Element = TensorFieldParal

    def __init__(self, vector_field_module, tensor_type):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "T^(" + str(kcon) + "," + str(lcov) + ")(" + domain._name
        latex_name = r"\mathcal{T}^{(" + str(kcon) + "," + str(lcov) + r")}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")" 
            latex_name += r"\right)" 
        else:
            name += "," + dest_map._name + ")" 
            latex_name += "," + dest_map._latex_name + r"\right)" 
        TensorFreeModule.__init__(self, vector_field_module, tensor_type, 
                                  name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    def _element_constructor_(self, comp=[], basis=None, name=None, 
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor
        """
        if comp == 0:
            return self._zero_element
        if isinstance(comp, TensorField):
            if self._tensor_type == comp._tensor_type and \
               self._domain.is_subdomain(comp._domain) and \
               self._ambient_domain.is_subdomain(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the " + str(comp) +
                                "to a tensor field in " + str(self))
        resu = self.element_class(self._fmodule, self._tensor_type, name=name, 
                                  latex_name=latex_name, sym=sym, 
                                  antisym=antisym)
        if comp != []:
            resu.set_comp(basis)[:] = comp
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent
        """
        if isinstance(other, (TensorFieldModule, TensorFieldFreeModule)):
            return self._tensor_type == other._tensor_type and \
                   self._domain.is_subdomain(other._domain) and \
                   self._ambient_domain.is_subdomain(other._ambient_domain)
        else:
            return False

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self._name is not None:
            description += self._name + " "
        description += "of type-(%s,%s)" % \
                           (str(self._tensor_type[0]), str(self._tensor_type[1]))
        description += " tensors fields "
        if self._dest_map is self._domain._identity_map:
            description += "on the " + str(self._domain)
        else:
            description += "along the " + str(self._domain) + \
                           " mapped into the " + str(self._ambient_domain)
        return description

