r"""
Submanifolds

The class :class:`Submanifold` implements submanifolds of differentiable 
manifolds over `\RR`. By *submanifold* it is meant *embedded submanifold*, 
the definition of which follows:

Given a differentiable manifold `M`, an *embedded submanifold* is a subset
`S\subset M` such that `S` is a manifold for the topology induced by `M`
and `S` is endowed with a differentiable structure with respect to which the 
inclusion map `\iota: S\hookrightarrow M` is a differentiable embedding.

Let us recall that a *differentiable embedding* of a differentiable manifold `S` 
into a differentiable manifold `M` is a differentiable mapping 
`\Phi: S \rightarrow M` whose differential is injective at each point 
(i.e. `\Phi` is an *immersion*) and that is a homeomorphism onto its image. 

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

from sage.structure.sage_object import SageObject
from manifold import Manifold
from chart import Chart
from diffmapping import DiffMapping

class Submanifold(Manifold):
    r"""
    Base class for submanifolds, i.e. manifolds embedded in a differentiable 
    manifold.
    
    INPUT:
    
    - ``ambient_manifold`` -- the ambient manifold
    - ``n`` -- dimension of the submanifold
    - ``name`` -- name given to the submanifold 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the submanifold
    - ``start_index`` -- (default: 0) lower bound of the range of indices on the
      submanifold

    EXAMPLES:
 
    The sphere `S^2` as a submanifold of the Euclidean space `\RR^3`::
    
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()   # Cartesian coordinates on R^3
        sage: S = M.submanifold(2, 'S^2', start_index=1) 
        sage: U = S.open_domain('U') # U = S minus two poles
        sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coordinates on U
        sage: emb = U.diff_mapping(M, [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)], name='i', latex_name=r'\iota') # the inclusion as an embedding S --> M
        sage: S.def_embedding(emb)
        sage: S
        2-dimensional submanifold 'S^2' of the 3-dimensional manifold 'R^3'
        sage: S._embedding.view()
        i: U --> R^3
           (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))

    A helix as a submanifold of `\RR^3`::
    
        sage: H = M.submanifold(1, 'H')
        sage: c_param.<t> = H.chart()
        sage: emb = H.diff_mapping(M, [cos(t), sin(t), t], name='i', latex_name=r'\iota')
        sage: H.def_embedding(emb)
        sage: H
        1-dimensional submanifold 'H' of the 3-dimensional manifold 'R^3'
        sage: H._embedding.view()
         i: H --> R^3
            t |--> (x, y, z) = (cos(t), sin(t), t)

    The constructed submanifolds are automatically added to the subdomains of 
    the ambient manifold::
    
        sage: M.domains()
        [3-dimensional manifold 'R^3', domain 'S^2' on the 3-dimensional manifold 'R^3', domain 'H' on the 3-dimensional manifold 'R^3']

    Pullback of 1-forms defined on `\RR^3` to `S^2`::
    
        sage: dX = c_cart.coframe() ; dX
        coordinate coframe (R^3, (dx,dy,dz))
        sage: dX[1]
        1-form 'dx' on the 3-dimensional manifold 'R^3'
        sage: S._embedding.pullback(dX[1])
        1-form 'i_*(dx)' on the open domain 'U' on the 2-dimensional submanifold 'S^2' of the 3-dimensional manifold 'R^3'
        sage: S._embedding.pullback(dX[1]).view()
        i_*(dx) = cos(ph)*cos(th) dth - sin(ph)*sin(th) dph
        sage: S._embedding.pullback(dX[2]).view()
        i_*(dy) = cos(th)*sin(ph) dth + cos(ph)*sin(th) dph
        sage: S._embedding.pullback(dX[3]).view()
        i_*(dz) = -sin(th) dth
        
    Pushforward of vector fields defined on `S^2` to `\RR^3`::
    
        sage: e = S.default_frame() ; e
        coordinate frame (U, (d/dth,d/dph))
        sage: e[1]
        vector field 'd/dth' on the open domain 'U' on the 2-dimensional submanifold 'S^2' of the 3-dimensional manifold 'R^3'
        sage: S.pushforward(e[1])
        vector field 'i_*(d/dth)' along the open domain 'U' on the 2-dimensional submanifold 'S^2' of the 3-dimensional manifold 'R^3' with values on the 3-dimensional manifold 'R^3'
        sage: S.pushforward(e[1]).view()
        i_*(d/dth) = cos(ph)*cos(th) d/dx + cos(th)*sin(ph) d/dy - sin(th) d/dz
        sage: S.pushforward(e[2]).view()
        i_*(d/dph) = -sin(ph)*sin(th) d/dx + cos(ph)*sin(th) d/dy

    """
    def __init__(self, ambient_manifold, n, name, latex_name=None, 
                 start_index=0):
        if not isinstance(ambient_manifold, Manifold):
            raise TypeError("The argument ambient_manifold must be a manifold.")
        self._ambient_manifold = ambient_manifold
        Manifold.__init__(self, n, name, latex_name, start_index)
        # The embedding map:
        self._embedding = None  # not defined yet
        # The submanifold as a subdomain of the ambient manifold:
        self._domain_amb = self._ambient_manifold.domain(name, latex_name)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return str(self._dim) + \
            "-dimensional submanifold '%s'" % self._name + \
            " of the " + str(self._ambient_manifold)

    def def_embedding(self, embedding):
        r"""
        Define the inclusion map of ``self``, `S` say, into the ambient 
        manifold, `M` say, as a differentiable embedding `S\rightarrow M`
        
        INPUT:
        
        - ``embedding`` -- the inclusion map as a differentiable embedding 
          `S\rightarrow M`; must be an instance of class :class:`DiffMapping`. 
         
        """
        if not isinstance(embedding, DiffMapping):
            raise TypeError("The argument must an instance of " + 
                            "class DiffMapping.")
        self._embedding = embedding

    def embedding(self):
        r"""
        Return the embedding defining the submanifold.
        """
        return self._embedding

    def plot(self, coord_ranges, local_chart=None, ambient_chart=None, **kwds):
        r"""
        Plot of a submanifold embedded in `\RR^2` or `\RR^3`

        INPUT:

         - ``coord_ranges`` -- list of pairs (u_min, u_max) for each coordinate
           u on the submanifold
         - ``local_chart`` -- (default: None) chart on the submanifold in which 
           the above coordinates are defined; if none is provided, the 
           submanifold default chart is assumed. 
         - ``ambient_chart`` -- (default: None) chart on the ambient manifold 
           in terms of which the embedding is defined; if none is provided, the 
           ambient manifold default chart is assumed. 
         - ``**kwds`` -- (default: None) keywords passed to Sage graphic 
           routines
           
        OUTPUT:
        
          - Graphics3d object (ambient manifold = `\RR^3`) or Graphics object
            (ambient manifold = `\RR^2`)

        EXAMPLES:

        Plot of a torus embedded in `\RR^3`::
            
            sage: M = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart.<x,y,z> = M.chart()  # Cartesian coordinates on R^3
            sage: T = M.submanifold(2, 'T')
            sage: W = T.open_domain('W') # Domain of the torus covered by the cyclic coordinates (u,v)
            sage: c_uv.<u,v> = W.chart(r'u:(0,2*pi) v:(0,2*pi)') # cyclic coordinates on T
            sage: T.def_embedding(T.diff_mapping(M, [(2+cos(u))*cos(v),(2+cos(u))*sin(v),sin(u)]))
            sage: T.plot([(0,2*pi), (0,2*pi)], aspect_ratio=1)

        Plot of a helix embedded in `\RR^3`::
            
            sage: H = M.submanifold(1, 'H')
            sage: c_param.<t> = H.chart()
            sage: H.def_embedding(H.diff_mapping(M, [cos(t), sin(t), t]))
            sage: H.plot([0,20])
           
        Plot of an Archimedean spiral embedded in `\RR^2`::
        
            sage: M = Manifold(2, 'R^2', r'\RR^2')
            sage: c_cart.<x,y> = M.chart()
            sage: S = M.submanifold(1, 'S')
            sage: c_param.<t> = S.chart()
            sage: S.def_embedding(S.diff_mapping(M, [t*cos(t), t*sin(t)]))
            sage: S.plot([0,40])

        """
        from sage.plot.plot import parametric_plot
        if local_chart is None:
            local_chart = self._def_chart
        if ambient_chart is None:
            ambient_chart = self._ambient_manifold._def_chart
        if self._dim > 2:
            raise ValueError("The dimension must be at most 2 " + 
                             "for plotting.")
        coord_functions = \
          self._embedding._coord_expression[(local_chart, ambient_chart)]._functions
        coord_express = [coord_functions[i]._express for i in 
                                             range(self._ambient_manifold._dim)]
        if self._dim == 1:
            urange = (local_chart._xx[0], coord_ranges[0], coord_ranges[1])
            graph = parametric_plot(coord_express, urange, **kwds)
        else:   # self._dim = 2
            urange = (local_chart._xx[0], coord_ranges[0][0], coord_ranges[0][1])
            vrange = (local_chart._xx[1], coord_ranges[1][0], coord_ranges[1][1])
            graph = parametric_plot(coord_express, urange, vrange, **kwds)
        return graph

    def pushforward(self, tensor):
        r""" 
        Pushforward operator associated with the embedding.

        INPUT:
        
        - ``tensor`` -- instance of :class:`TensorField` representing a fully 
          contravariant tensor field `T` on the submanifold, i.e. a tensor 
          field of type (k,0), with k a positive integer. The case k=0 
          corresponds to a scalar field.
          
        OUTPUT:
        
        - instance of :class:`TensorField` representing a field of fully 
          contravariant tensors of the ambient manifold, field defined on 
          the domain occupied by the submanifold. 
          
        EXAMPLES:

        Pushforward of a vector field defined on `S^2`, submanifold of `\RR^3`::
        
            sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_cart.<x,y,z> = M.chart() # Cartesian coordinates on R^3
            sage: S = M.submanifold(2, 'S^2', start_index=1)
            sage: U = S.open_domain('U') # U = S minus two poles
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coordinates on U
            sage: S.def_embedding( S.diff_mapping(M, [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)], name='i', latex_name=r'\iota') )
            sage: v = U.vector_field(name='v') 
            sage: v[2] = 1 ; v.view()  # azimuthal vector field on S^2
            v = d/dph
            sage: iv = S.pushforward(v) ; iv
            vector field 'i_*(v)' along the open domain 'U' on the 2-dimensional submanifold 'S^2' of the 3-dimensional manifold 'R^3' with values on the 3-dimensional manifold 'R^3'
            sage: iv.view() 
            i_*(v) = -sin(ph)*sin(th) d/dx + cos(ph)*sin(th) d/dy
            
        The components of the pushforward vector are scalar fields on the submanifold::
        
            sage: iv.comp()[[1]]
            scalar field on the open domain 'U' on the 2-dimensional submanifold 'S^2' of the 3-dimensional manifold 'R^3'

        Pushforward of a tangent vector to a helix, submanifold of `\RR^3`::
        
            sage: H = M.submanifold(1, 'H')
            sage: c_t.<t> = H.chart()
            sage: H.def_embedding( H.diff_mapping(M, [cos(t), sin(t), t], name='iH') )
            sage: u = H.vector_field(name='u')
            sage: u[0] = 1 ; u.view() # tangent vector to the helix
            u = d/dt
            sage: iu = H.pushforward(u) ; iu
            vector field 'iH_*(u)' along the 1-dimensional submanifold 'H' of the 3-dimensional manifold 'R^3' with values on the 3-dimensional manifold 'R^3'
            sage: iu.view()
            iH_*(u) = -sin(t) d/dx + cos(t) d/dy + d/dz

        """
        from tensorfield import TensorFieldParal        
        from sage.tensor.modules.comp import Components, CompWithSym, \
                                                 CompFullySym, CompFullyAntiSym
        dom1 = tensor._domain
        if not dom1.is_subdomain(self):
            raise TypeError("The tensor field is not defined on " + str(self))
        (ncon, ncov) = tensor._tensor_type
        if ncov != 0:
            raise TypeError("The pushforward cannot be taken on a tensor " + 
                            "with some covariant part.")
        embed = self._embedding
        resu_name = None ; resu_latex_name = None
        if embed._name is not None and tensor._name is not None:
            resu_name = embed._name + '_*(' + tensor._name + ')'
        if embed._latex_name is not None and tensor._latex_name is not None:
            resu_latex_name = embed._latex_name + '_*' + tensor._latex_name                
        if ncon == 0:
            raise NotImplementedError("The case of a scalar field is not " + 
                                      "implemented yet.")
        if not isinstance(tensor, TensorFieldParal):
            raise NotImplementedError("The case of a non-parallelizable " + 
                                      "domain is not implemented yet.")
        # A pair of charts (chart1, chart2) where the computation
        # is feasable is searched, privileging the default chart of the 
        # start domain for chart1
        chart1 = None; chart2 = None
        def_chart1 = dom1._def_chart
        def_chart2 = embed._codomain._def_chart 
        if def_chart1._frame in tensor._components and \
               (def_chart1, def_chart2) in embed._coord_expression:
            chart1 = def_chart1
            chart2 = def_chart2
        else:
            for (chart1n, chart2n) in embed._coord_expression:
                if chart2n == def_chart2 and \
                                    chart1n._frame in tensor._components:
                    chart1 = chart1n
                    chart2 = def_chart2
                    break
        if chart2 is None:
            # It is not possible to have def_chart2 as chart for 
            # expressing the result; any other chart is then looked for:
            for (chart1n, chart2n) in self._coord_expression:
                if chart1n._frame in tensor._components:
                    chart1 = chart1n
                    chart2 = chart2n
                    break
        if chart1 is None:
            raise ValueError("No common chart could be find to compute " +
                             "the pushforward of the tensor field.")
        fmodule2 = dom1.vector_field_module(dest_map=embed)
        frame2 = dom1.vector_frame(dest_map=embed, from_frame=chart2._frame)
        si1 = dom1._manifold._sindex
        si2 = fmodule2._sindex
        ring2 = fmodule2._ring
        of2 = fmodule2._output_formatter
        # Computation at the component level:
        tcomp = tensor._components[chart1._frame]
        # Construction of the pushforward components (ptcomp):
        if isinstance(tcomp, CompFullySym):
            ptcomp = CompFullySym(ring2, frame2, ncon, start_index=si2, 
                                  output_formatter=of2)
        elif isinstance(tcomp, CompFullyAntiSym):
            ptcomp = CompFullyAntiSym(ring2, frame2, ncon, start_index=si2, 
                                      output_formatter=of2)
        elif isinstance(tcomp, CompWithSym):
            ptcomp = CompWithSym(ring2, frame2, ncon, start_index=si2, 
                                 output_formatter=of2, sym=tcomp._sym, 
                                 antisym=tcomp._antisym)
        else:
            ptcomp = Components(ring2, frame2, ncon, start_index=si2, 
                                output_formatter=of2)
        phi = embed._coord_expression[(chart1, chart2)]
        jacob = phi.jacobian()
        # X2 coordinates expressed in terms of X1 ones via the mapping:
        # coord2_1 = phi(*(chart1._xx)) 
        si1 = self._sindex
        si2 = self._ambient_manifold._sindex
        for ind_new in ptcomp.non_redundant_index_generator(): 
            res = 0 
            for ind_old in self.index_generator(ncon): 
                t = tcomp[[ind_old]].function_chart(chart1)
                for i in range(ncon):
                    t *= jacob[ind_new[i]-si2][ind_old[i]-si1]
                res += t
            ptcomp[ind_new] = res
        resu = fmodule2.tensor_from_comp((ncon, 0), ptcomp, name=resu_name, 
                                         latex_name=resu_latex_name)
        return resu
