r"""
Scalar fields

The class :class:`ScalarField` implements scalar fields on differentiable 
manifolds over `\RR`, i.e. differentiable mappings of the form

.. MATH::

    f: U\subset M \longrightarrow \RR
    
where `U` is an open subset of the differentiable manifold `M`.

The class :class:`ScalarField`  inherits from the class  :class:`~sage.structure.element.CommutativeAlgebraElement` (a scalar field on
`U` being an element of the commutative algebra `C^\infty(U)`).

The subclass :class:`ZeroScalarField` deals with null scalar fields. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013,2014): initial version

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

from sage.structure.element import CommutativeAlgebraElement
from sage.rings.integer import Integer
from domain import Domain
from chart import FunctionChart, ZeroFunctionChart, MultiFunctionChart

class ScalarField(CommutativeAlgebraElement):
    r"""
    Class for scalar fields on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold open subset `U` on which the scalar field is 
      defined (must be an instance of class 
      :class:`~sage.geometry.manifolds.domain.OpenDomain`)
    - ``coord_expression`` -- (default: None) coordinate expression of the 
      scalar field
    - ``name`` -- (default: None) name given to the scalar field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the scalar field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A scalar field on the 2-sphere::
    
        sage: M = Manifold(2, 'S^2')
        sage: f = M.scalar_field() ; f
        scalar field on the 2-dimensional manifold 'S^2'
    
    Scalar fields on `M` belong to the algebra `C^\infty(M)`::
     
        sage: f.parent()
        algebra of scalar fields on the 2-dimensional manifold 'S^2'
        sage: f.parent() is M.scalar_field_algebra()
        True

    Named scalar field::
    
        sage: f = M.scalar_field(name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: latex(f)
        f

    Named scalar field with LaTeX symbol specified::
    
        sage: f = M.scalar_field(name='f', latex_name=r'\mathcal{F}') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: latex(f)
        \mathcal{F}
        
    Scalar field defined by its coordinate expression, for instance in terms
    of spherical coordinates defined on the complement `U` of some origin 
    half meridian::
    
        sage: U = M.open_domain('U')
        sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: f = M.scalar_field(sin(th)*cos(ph), chart=c_spher, name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: f.view(chart=c_spher)
        f: S^2 --> R
        on U: (th, ph) |--> cos(ph)*sin(th)
        sage: f.expr(chart=c_spher)
        cos(ph)*sin(th)

    Since c_spher is the default chart on `M` (being the first defined one), 
    it can be omitted in the argument lists::
    
        sage: f = M.scalar_field(sin(th)*cos(ph),  name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: f.view()
        f: S^2 --> R
        on U: (th, ph) |--> cos(ph)*sin(th)
        sage: f.expr()
        cos(ph)*sin(th)



    The coordinate expression of a scalar field can be read by means of 
    :meth:`expr` and set by means of :meth:`set_expr`; both methods can 
    have a chart as argument (if not, the manifold's default chart is 
    assumed)::
    
        sage: f.set_expr(cos(th))  # changing the value of f
        sage: f.expr()
        cos(th)
        sage: f.set_expr(sin(th)*cos(ph)) # restoring the original value
        
    The function :meth:`view` displays the coordinate expression of the scalar
    field::
    
        sage: f.view()
        f: S^2 --> R
        on U: (th, ph) |--> cos(ph)*sin(th)
        sage: f.view(c_spher) # equivalent to above since c_spher is the default chart
        f: S^2 --> R
        on U: (th, ph) |--> cos(ph)*sin(th)
        sage: latex(f.view(c_spher)) # nice LaTeX formatting for the notebook
        \begin{array}{llcl} f:& S^2 & \longrightarrow & \RR \\ \mbox{on}\ U : & \left(\theta, \phi\right) & \longmapsto & \cos\left(\phi\right) \sin\left(\theta\right) \end{array}

    A scalar field can also be defined by an unspecified function of the 
    coordinates::
    
        sage: g = M.scalar_field(function('G', th, ph), name='g') ; g
        scalar field 'g' on the 2-dimensional manifold 'S^2'
        sage: g.expr()
        G(th, ph)
        sage: s = f+g ; s.expr()                               
        cos(ph)*sin(th) + G(th, ph)
        
    In each chart, the scalar field is represented by a function of the 
    coordinates, which is a an instance of the class 
    :class:`~sage.geometry.manifolds.chart.FunctionChart` 
    and can be accessed by the method :meth:`function_chart`::
    
        sage: f.function_chart(c_spher)
        cos(ph)*sin(th)
        sage: f.function_chart() # equivalent to above since c_spher is the default chart
        cos(ph)*sin(th)
        sage: print type(f.function_chart())
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        
    The value returned by the method :meth:`expr` is actually the coordinate
    expression of the function::

        sage: f.expr() is f.function_chart().expr()
        True
        
    By definition, a scalar field acts on the manifold's points::
    
        sage: p = M.point((pi/2, pi))
        sage: f(p)
        -1

    A scalar field can be compared to another scalar field::
    
        sage: g = M.scalar_field(sin(th)*cos(ph), name='g')
        sage: f == g
        True
        sage: g.set_expr(cos(th))
        sage: f == g
        False
        
    ...to a symbolic expression::
    
        sage: f == sin(th)*cos(ph)
        True
        sage: f == ph + th^2
        False
        
    ...to a number::
        
        sage: f == 2
        False

    ...to zero::
    
        sage: f == 0
        False
        sage: f.set_expr(0)
        sage: f == 0
        True

    ...to anything else::
    
        sage: f == M
        False

    Scalar fields can be added::
    
        sage: f.set_expr(sin(th)*cos(ph))
        sage: g.set_expr(cos(th))
        sage: s = f + g ; s
        scalar field 'f+g' on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) + cos(th)
        sage: s = f + cos(th) ; s # direct addition with a symbolic expression is allowed
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) + cos(th)
        sage: s = 1 + f ; s 
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) + 1
        sage: s = +f ; s  # the unary plus operator
        scalar field '+f' on the 2-dimensional manifold 'S^2'
        sage: s == f
        True
      
    Scalar fields can be subtracted::
    
        sage: s = f - g ; s
        scalar field 'f-g' on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) - cos(th)
        sage: s = f - cos(th) ; s  # direct subtraction of a symbolic expression is allowed
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) - cos(th)
        sage: s = 1 - f ; s
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        -cos(ph)*sin(th) + 1
        sage: s = f - g + (g - f)
        sage: s == 0
        True
        sage: s = f + (-f) # check of the unary minus operator
        sage: s == 0
        True
     
    Scalar fields can be multiplied and divided::
     
        sage: s = f*g ; s
        scalar field 'f*g' on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*cos(th)*sin(th)
        sage: h = s / g ; h 
        scalar field 'f*g/g' on the 2-dimensional manifold 'S^2'
        sage: h.expr()
        cos(ph)*sin(th)
        sage: h == f
        True
        sage: h1 = s / f ; h1 
        scalar field 'f*g/f' on the 2-dimensional manifold 'S^2'
        sage: h1.expr()
        cos(th)
        sage: h1 == g
        True
            
    The multiplication and division can be performed by a symbolic expression::
    
        sage: s = f*cos(th) ; s
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*cos(th)*sin(th)
        sage: h = s/cos(th) ; h
        scalar field on the 2-dimensional manifold 'S^2'
        sage: h.expr()
        cos(ph)*sin(th)
        sage: h == f
        True        

    The in-place operators +=, -=, \*= and /= are implemented::
    
        sage: f.expr()
        cos(ph)*sin(th)
        sage: f += cos(th)
        sage: f.expr()    
        cos(ph)*sin(th) + cos(th)
        sage: f -= cos(th)
        sage: f.expr()    
        cos(ph)*sin(th)
        sage: f *= cos(th)
        sage: f.expr()
        cos(ph)*cos(th)*sin(th)
        sage: f /= cos(th)
        sage: f.expr()    
        cos(ph)*sin(th)
        
    Test of the arithmetics of scalar fields defined on multiple domains::
    
        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U')
        sage: c_xy.<x,y> = U.chart()
        sage: V = M.open_domain('V')
        sage: c_uv.<u,v> = V.chart()
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: f = M.scalar_field(x^2)
        sage: f.add_expr(u, c_uv)
        sage: g = M.scalar_field(2*v, c_uv)
        sage: g.add_expr(y, c_xy)
        sage: f._express  # random (dictionary output)
        {chart (U, (x, y)): x^2, chart (V, (u, v)): u}
        sage: g._express  # random (dictionary output)
        {chart (V, (u, v)): 2*v, chart (U, (x, y)): y}
        sage: s = f + g ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s._express  # random (dictionary output)
        {chart (U, (x, y)): x^2 + y, chart (V, (u, v)): u + 2*v}
        sage: g.set_expr(3*x, c_xy)
        sage: g._express
        {chart (U, (x, y)): 3*x}
        sage: s = f + g ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s._express
        {chart (U, (x, y)): x^2 + 3*x}
        sage: g = U.scalar_field(3*x)
        sage: g._express
        {chart (U, (x, y)): 3*x}
        sage: s = f + g ; s
        scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: s._express
        {chart (U, (x, y)): x^2 + 3*x}
    
    Vanishing result::
    
        sage: g = M.scalar_field(-x^2)
        sage: g.add_expr(-u, c_uv)
        sage: s = f + g ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: print type(s)
        <class 'sage.geometry.manifolds.scalarfield.ZeroScalarField'>

    """
    def __init__(self, domain, coord_expression=None, name=None, 
                 latex_name=None):
        CommutativeAlgebraElement.__init__(self, domain.scalar_field_algebra())
        self._manifold = domain._manifold
        self._domain = domain
        self._tensor_type = (0,0)
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._express = {} # dict of coordinate expressions (FunctionChart
                           # instances) with charts as keys
        if coord_expression is not None:
            if isinstance(coord_expression, FunctionChart):
                self._express[coord_expression._chart] = coord_expression
            elif isinstance(coord_expression, dict):
                for chart, expression in coord_expression.iteritems():
                    if isinstance(expression, FunctionChart):
                        self._express[chart] = expression
                    else:
                        self._express[chart] = FunctionChart(chart, expression)
            elif coord_expression == 0:
                for chart in self._domain._atlas:
                    self._express[chart] = chart._zero_function
            else:
                for chart in self._domain._atlas:
                    self._express[chart] = FunctionChart(chart, 
                                                        coord_expression)
        self._init_derived()   # initialization of derived quantities

    ####### Required methods for an algebra element (beside arithmetic) #######
    
    def __nonzero__(self):
        r"""
        Return True if ``self`` is nonzero and False otherwise. 
        
        This method is called by self.is_zero(). 

        EXAMPLES:
        
        Tests on a 2-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y)
            sage: f.is_zero()
            False
            sage: f.set_expr(0)
            sage: f.is_zero()
            True
            sage: g = M.scalar_field(0)
            sage: g.is_zero()
            True

        """
        res = True
        for funct in self._express.itervalues():
            res = res and funct.is_zero()
        return not res

    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a scalar field
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if not isinstance(other, ScalarField):
            try:
                other = self.parent()(other)    # conversion to a scalar field
            except TypeError:
                return False
        if other._domain != self._domain:
            return False
        if other.is_zero():
            return self.is_zero()
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the comparison.")
        resu = True
        for chart in com_charts:
            resu = resu and (self._express[chart] == other._express[chart])
        return resu

    def __ne__(self, other):
        r"""
        Non-equality operator.
        """
        return not self.__eq__(other)
        
    ####### End of required methods for an algebra element (beside arithmetic) #######

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        self._restrictions = {} # dict. of restrictions of self on subdomains  
                                # of self._domain, with the subdomains as keys
        self._differential = None # differential
        self._lie_derivatives = {} # collection of Lie derivatives of self

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._restrictions.clear()
        self._differential = None 
        # First deletes any reference to self in the vectors' dictionary:
        for vid, val in self._lie_derivatives.iteritems():
            del val[0]._lie_der_along_self[id(self)]
        self._lie_derivatives.clear()

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "scalar field"
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

    def set_name(self, name=None, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:
        
        - ``name`` -- (string; default: None) name given to the scalar field
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote 
          the scalar field; if None while ``name`` is provided, the LaTeX 
          symbol is set to ``name``.

        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` and on the same domain.
        
        """
        return self.__class__(self._domain)        

    def domain(self):
        r"""
        Return the domain on which the scalar field is defined.
        
        OUTPUT:
        
        - instance of class :class:`~sage.geometry.manifolds.domain.OpenDomain` 
          representing the manifold's open subset on which ``self`` is defined. 
        
        EXAMPLES::
        
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x+2*y)
            sage: f.domain()
            2-dimensional manifold 'M'
            sage: U = M.open_domain('U', coord_def={c_xy: x<0})
            sage: g = f.restrict(U)
            sage: g.domain()
            open domain 'U' on the 2-dimensional manifold 'M'
        
        """
        return self._domain

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``, i.e. (0,0), when ``self`` is 
        considered as a tensor field on the manifold.
        
        EXAMPLE::
        
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x+2*y)
            sage: f.tensor_type()
            (0, 0)
        
        """
        return self._tensor_type

    def copy(self):
        r"""
        Return an exact copy of ``self``.
        
        EXAMPLES:
        
        Copy on a 2-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')  
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: g = f.copy()
            sage: print type(g)
            <class 'sage.geometry.manifolds.scalarfield.ScalarFieldAlgebra_with_category.element_class'>
            sage: g.expr()
            x*y^2
            sage: g == f
            True
            sage: g is f
            False
        
        """
        result = self.__class__(self._domain)  #!# what about the name ?
        for chart, funct in self._express.iteritems():
            result._express[chart] = funct.copy()
        return result

    def function_chart(self, chart=None, from_chart=None):
        r""" 
        Return the function of the coordinates representing the scalar field 
        in a given chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) chart with respect to which the
          coordinate expression is to be returned; if None, the 
          domain's default chart will be used
        - ``from_chart`` -- (default: None) chart from which the
          required expression is computed if it is not known already in the 
          chart ``chart``; if None, a chart is picked in ``self._express``

        OUTPUT:
        
        - instance of :class:`~sage.geometry.manifolds.chart.FunctionChart` 
          representing the coordinate function of the scalar field in the 
          given chart.

        EXAMPLES:
        
        Coordinate function on a 2-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')            
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: f.function_chart()
            x*y^2
            sage: f.function_chart(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: print type(f.function_chart())
            <class 'sage.geometry.manifolds.chart.FunctionChart'>

        Expression via a change of coordinates::
        
            sage: c_uv.<u,v> = M.chart()
            sage: c_uv.coord_change(c_xy, u+v, u-v)
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: f._express # at this stage, f is expressed only in terms of (x,y) coordinates
            {chart (M, (x, y)): x*y^2}
            sage: f.function_chart(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: f.function_chart(c_uv) == (u+v)*(u-v)^2  # check
            True
            sage: f._express  # random (dict. output); f has now 2 coordinate expressions:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        Usage in a physical context (simple Lorentz transformation - boost in 
        x direction, with relative velocity v between o1 and o2 frames)::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: o1.<t,x> = M.chart()
            sage: o2.<T,X> = M.chart()
            sage: f = M.scalar_field(x^2 - t^2)
            sage: f.function_chart(o1)
            -t^2 + x^2
            sage: v = var('v'); gam = 1/sqrt(1-v^2)
            sage: o2.coord_change(o1, gam*(T - v*X), gam*(X - v*T))
            coordinate change from chart (M, (T, X)) to chart (M, (t, x))
            sage: f.function_chart(o2)
            -T^2 + X^2

        """
        if isinstance(self, ZeroScalarField):
            # to ensure that the ZeroScalarField version is called in case
            # of a direct call to the unbound function (ScalarField.function_chart)
            return self.function_chart(chart) 
        if chart is None:
            chart = self._domain._def_chart
        else:
            if chart not in self._domain._atlas:
                raise TypeError("The " + str(chart) + " has not " + \
                      " been defined on the domain " + str(self._domain))
        if chart not in self._express:
            # Check whether chart corresponds to a subchart of a chart
            # where the expression of self is known:
            for known_chart in self._express:
                if chart in known_chart._subcharts:
                    new_expr = self._express[known_chart].expr()
                    self._express[chart] = FunctionChart(chart, new_expr)
                    return self._express[chart]
            # If this point is reached, the expression must be computed 
            # from that in the chart from_chart, by means of a 
            # change-of-coordinates formula:
            if from_chart is None:
                # from_chart in searched among the charts of known expressions
                # and subcharts of them
                known_express = self._express.copy()
                found = False
                for kchart in known_express:
                    for skchart in kchart._subcharts:
                        if (chart, skchart) in self._domain._coord_changes:
                            from_chart = skchart
                            found = True
                            if skchart not in self._express:
                                self._express[skchart] = FunctionChart(skchart, 
                                                  self._express[kchart].expr())
                            break
                    if found:
                        break
                if not found:
                    raise ValueError("No starting chart could be found to " + 
                           "compute the expression in the " + str(chart))
            change = self._domain._coord_changes[(chart, from_chart)]
            # old coordinates expressed in terms of the new ones:
            coords = [ change._transf._functions[i]._express 
                       for i in range(self._manifold._dim) ]
            new_expr = self._express[from_chart](*coords)
            self._express[chart] = FunctionChart(chart, new_expr)
            self._del_derived()
        return self._express[chart]


    def expr(self, chart=None, from_chart=None):
        r""" 
        Return the coordinate expression of the scalar field in a given 
        chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) chart with respect to which the 
          coordinate expression is required; if None, the domain's default 
          chart will be used
        - ``from_chart`` -- (default: None) chart from which the
          required expression is computed if it is not known already in the 
          chart ``chart``; if None, a chart is picked in ``self._express``
          
        OUTPUT:
        
        - symbolic expression representing the coordinate 
          expression of the scalar field in the given chart.
        
        EXAMPLES:
        
        Expression of a scalar field on a 2-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')            
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: f.expr()
            x*y^2
            sage: f.expr(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: print type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        Expression via a change of coordinates::
        
            sage: c_uv.<u,v> = M.chart()
            sage: c_uv.coord_change(c_xy, u+v, u-v)
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: f._express # at this stage, f is expressed only in terms of (x,y) coordinates
            {chart (M, (x, y)): x*y^2}
            sage: f.expr(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: bool( f.expr(c_uv) == (u+v)*(u-v)^2 ) # check
            True
            sage: f._express  # random (dict. output); f has now 2 coordinate expressions:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        """
        return self.function_chart(chart, from_chart)._express
        
    def set_expr(self, coord_expression, chart=None):
        r"""
        Set the coordinate expression of the scalar field.
        
        The expressions with respect to other charts are deleted, in order to 
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.
        
        INPUT:
        
        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: None) chart in which ``coord_expression`` is 
          defined; if None, the domain's default chart is assumed
        
        EXAMPLES:
        
        Setting scalar field expressions on a 2-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + 2*x*y +1)
            sage: f._express                         
            {chart (M, (x, y)): x^2 + 2*x*y + 1}
            sage: f.set_expr(3*y)
            sage: f._express  # the (x,y) expression has been changed:
            {chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = M.chart()
            sage: f.set_expr(cos(u)-sin(v), c_uv)  
            sage: f._express # the (x,y) expression has been lost:
            {chart (M, (u, v)): cos(u) - sin(v)}
            sage: f.set_expr(3*y)    
            sage: f._express # the (u,v) expression has been lost:                    
            {chart (M, (x, y)): 3*y}

        """
        if chart is None:
            chart = self._domain._def_chart
        self._express.clear()
        self._express[chart] = FunctionChart(chart, coord_expression)
        self._del_derived()

    def add_expr(self, coord_expression, chart=None):
        r"""
        Add some coordinate expression to the scalar field.
        
        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead. 
        
        INPUT:
        
        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: None) chart in which ``coord_expression``
          is defined; if None, the domain's default chart is assumed
          
        .. WARNING::
        
            If the scalar field has already expressions in other charts, it 
            is the user's responsability to make sure that the expression
            to be added is consistent with them. 
        
        EXAMPLES:
        
        Adding scalar field expressions on a 2-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + 2*x*y +1)
            sage: f._express
            {chart (M, (x, y)): x^2 + 2*x*y + 1}             
            sage: f.add_expr(3*y)
            sage: f._express  # the (x,y) expression has been changed:
            {chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = M.chart()
            sage: f.add_expr(cos(u)-sin(v), c_uv)  
            sage: f._express # random (dict. output); f has now 2 expressions:
            {chart (M, (x, y)): 3*y, chart (M, (u, v)): cos(u) - sin(v)}

        """
        if chart is None:
            chart = self._domain._def_chart
        self._express[chart] = FunctionChart(chart, coord_expression)
        self._del_derived()

    def add_expr_by_continuation(self, chart, subdomain):
        r"""
        Set coordinate expression in a chart by continuation of the
        coordinate expression in a subchart.  
        
        The continuation is performed by demanding that the coordinate
        expression is identical to that in the restriction of the chart to 
        a given subdomain.  
        
        INPUT:
        
        - ``chart`` -- coordinate chart `(U,(x^i))` in which the expression of 
          the scalar field is to set
        - ``subdomain`` -- open domain `V\subset U` in which the expression
          in terms of the restriction of the coordinate chart `(U,(x^i))` to 
          `V` is already known or can be evaluated by a change of coordinates.
          
        EXAMPLE:
        
        Scalar field on the sphere `S^2`::
        
            sage: M = Manifold(2, 'S^2')
            sage: U = M.open_domain('U') ; V = M.open_domain('V') # the complement of resp. N pole and S pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W =  U.intersection(V)  # S^2 minus the two poles
            sage: f = M.scalar_field(atan(x^2+y^2), chart=c_xy, name='f')

        The scalar field has been defined only on the domain covered by the
        chart c_xy, i.e. `U`::
        
            sage: f.view()
            f: S^2 --> R
            on U: (x, y) |--> arctan(x^2 + y^2)

        We note that on `W = U\cap V`, the expression of `f` in terms of 
        coordinates `(u,v)` can be deduced from that in the coordinates 
        `(x,y)` thanks to the transition map between the two charts::
        
            sage: f.view(c_uv.restrict(W))
            f: S^2 --> R
            on W: (u, v) |--> arctan(1/(u^2 + v^2))
            
        We use this fact to extend the definition of `f` to domain `V`, 
        covered by the chart c_uv::
        
            sage: f.add_expr_by_continuation(c_uv, W)
            
        Then, `f` is known on the whole sphere::
        
            sage: f.view()
            f: S^2 --> R
            on U: (x, y) |--> arctan(x^2 + y^2)
            on V: (u, v) |--> arctan(1/(u^2 + v^2))

        """
        if not chart._domain.is_subdomain(self._domain):
            raise ValueError("The chart is not defined on a subdomain of " + 
                             "the scalar field domain.")
        schart = chart.restrict(subdomain)
        self._express[chart] = FunctionChart(chart, self.expr(schart))
        self._del_derived()

    def _display_expression(self, chart, result):
        r"""
        Helper function for :meth:`view`.
        """
        from sage.misc.latex import latex
        try:
            expression = self.expr(chart)
            coords = chart[:]
            if len(coords) == 1:
                coords = coords[0]
            if chart._domain == self._domain:
                if self._name is not None:
                    result.txt += "   " 
                result.latex += " & " 
            else:
                result.txt += "on " + chart._domain._name + ": " 
                result.latex += r"\mbox{on}\ " + latex(chart._domain) + r": & " 
            result.txt += repr(coords) + " |--> " + repr(expression) + "\n"
            result.latex += latex(coords) + r"& \longmapsto & " + \
                            latex(expression) + r"\\"
        except (TypeError, ValueError):
            pass

    def view(self, chart=None):
        r""" 
        Display the expression of the scalar field in a given chart. 
        
        Without any argument, this function displays the expressions of the
        scalar field in all the charts defined on the scalar field's domain
        that are not restrictions of another chart to some subdomain 
        (the "top charts"). 
        
        INPUT:
        
        - ``chart`` -- (default: None) chart with respect to which the 
          coordinate expression is to be displayed; if None, the display is
          performed in all the top charts in which the coordinate expression is 
          known. 
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLES:
        
        Various displays::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(sqrt(x+1), name='f')
            sage: f.view()
            f: M --> R
               (x, y) |--> sqrt(x + 1)
            sage: latex(f.view())
            \begin{array}{llcl} f:& M & \longrightarrow & \RR \\ & \left(x, y\right) & \longmapsto & \sqrt{x + 1} \end{array}
            sage: g = M.scalar_field(function('G', x, y), name='g')
            sage: g.view()
            g: M --> R
               (x, y) |--> G(x, y)
            sage: latex(g.view())
            \begin{array}{llcl} g:& M & \longrightarrow & \RR \\ & \left(x, y\right) & \longmapsto & G\left(x, y\right) \end{array}

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if self._name is None:
            symbol = ""
        else:
            symbol = self._name + ": "
        result.txt = symbol + self._domain._name + " --> R\n"
        if self._latex_name is None:
            symbol = ""
        else:
            symbol = self._latex_name + ":"
        result.latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                       latex(self._domain) + r"& \longrightarrow & \RR \\"
        if chart is None:
            for ch in self._domain._top_charts:
                self._display_expression(ch, result)
        else:
            self._display_expression(chart, result)             
        result.txt = result.txt[:-1]
        result.latex = result.latex[:-2] + r"\end{array}"
        return result

    def restrict(self, subdomain):
        r"""
        Restriction of the scalar field to a subdomain of its domain of 
        definition.
        
        INPUT:
        
        - ``subdomain`` -- the subdomain (instance of
          :class:`~sage.geometry.manifolds.domain.OpenDomain`)
        
        OUTPUT:
        
        - instance of :class:`ScalarField` representing the restriction of 
          ``self`` to ``subdomain``.

        EXAMPLE: 
        
        Restriction of a scalar field defined on `\RR^2` to the unit open 
        disc::
        
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # Cartesian coordinates
            sage: U = M.open_domain('U')
            sage: X_U = X.restrict(U, x^2+y^2 < 1)  # U is the unit open disc
            sage: f = M.scalar_field(cos(x*y), name='f')
            sage: f_U = f.restrict(U) ; f_U
            scalar field 'f' on the open domain 'U' on the 2-dimensional manifold 'M'
            sage: f_U.view()
            f: U --> R
               (x, y) |--> cos(x*y)
            sage: f.parent()
            algebra of scalar fields on the 2-dimensional manifold 'M'
            sage: f_U.parent()
            algebra of scalar fields on the open domain 'U' on the 2-dimensional manifold 'M'
        
        The restriction to the whole domain is the identity::
        
            sage: f.restrict(M) is f
            True
            sage: f_U.restrict(U) is f_U
            True

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subdomain(self._domain):
                raise ValueError("The specified domain is not a subdomain " + 
                                 "of the domain of definition of the scalar " + 
                                 "field.")
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.iteritems():
                if subdomain.is_subdomain(dom):
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    break
            else:
            # If this fails, the restriction is obtained via coercion
                resu = subdomain.scalar_field_algebra()(self)
                resu._name = self._name
                resu._latex_name = self._latex_name
                self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def common_charts(self, other):
        r"""
        Find common charts for the expressions of ``self`` and ``other``. 
        
        INPUT:
        
        - ``other`` -- a scalar field
        
        OUPUT:
        
        - list of common charts; if no common chart is found, None is 
          returned (instead of an empty list). 

        EXAMPLES:
        
        Search for common charts on a 2-dimensional manifold with 2 
        overlapping domains::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U')
            sage: c_xy.<x,y> = U.chart()
            sage: V = M.open_domain('V')
            sage: c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: f = U.scalar_field(x^2)
            sage: g = M.scalar_field(x+y)
            sage: f.common_charts(g)
            [chart (U, (x, y))]
            sage: g.add_expr(u, c_uv)
            sage: f._express
            {chart (U, (x, y)): x^2}
            sage: g._express  # random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [chart (U, (x, y))]

        Common charts found as subcharts: the subcharts are introduced via
        a transition map between charts c_xy and c_uv on the intersecting domain
        `W = U\cap V`::
        
            sage: trans = c_xy.transition_map(c_uv, (x+y, x-y), 'W', x<0, u+v<0)
            sage: M.atlas()
            [chart (U, (x, y)), chart (V, (u, v)), chart (W, (x, y)), chart (W, (u, v))]
            sage: c_xy_W = M.atlas()[2]
            sage: c_uv_W = M.atlas()[3]
            sage: trans.inverse()
            coordinate change from chart (W, (u, v)) to chart (W, (x, y))
            sage: f.common_charts(g)
            [chart (U, (x, y))]
            sage: f.expr(c_xy_W)  
            x^2
            sage: f._express  # random (dictionary output)
            {chart (U, (x, y)): x^2, chart (W, (x, y)): x^2}
            sage: g._express  # random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: g.common_charts(f)  # c_xy_W is not returned because it is subchart of 'xy'
            [chart (U, (x, y))]
            sage: f.expr(c_uv_W)
            1/4*u^2 + 1/2*u*v + 1/4*v^2
            sage: f._express  # random (dictionary output)
            {chart (U, (x, y)): x^2, chart (W, (x, y)): x^2, chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2}
            sage: g._express  # random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [chart (U, (x, y)), chart (W, (u, v))]
            sage: # the expressions have been updated on the subcharts
            sage: g._express #  random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u, chart (W, (u, v)): u}

        Common charts found by computing some coordinate changes::
        
            sage: W = U.intersection(V)
            sage: f = W.scalar_field(x^2, c_xy_W)
            sage: g = W.scalar_field(u+1, c_uv_W)
            sage: f._express
            {chart (W, (x, y)): x^2}
            sage: g._express
            {chart (W, (u, v)): u + 1}
            sage: f.common_charts(g)
            [chart (W, (u, v)), chart (W, (x, y))]
            sage: f._express # random (dictionary output)
            {chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2, chart (W, (x, y)): x^2}
            sage: g._express # random (dictionary output)
            {chart (W, (u, v)): u + 1, chart (W, (x, y)): x + y + 1}

        """
        if not isinstance(other, ScalarField):
            raise TypeError("The second argument must be a scalar field.")
        dom1 = self._domain
        dom2 = other._domain
        coord_changes = self._manifold._coord_changes
        resu = []
        #
        # 1/ Search for common charts among the existing expressions, i.e. 
        #    without performing any expression transformation. 
        #    -------------------------------------------------------------
        for chart1 in self._express:
            if chart1 in other._express:
                resu.append(chart1)
        # Search for a subchart:
        known_expr1 = self._express.copy()  
        known_expr2 = other._express.copy()
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if chart2 in chart1._subcharts:
                            self.expr(chart2)
                            resu.append(chart2)
                        if chart1 in chart2._subcharts:
                            other.expr(chart1)
                            resu.append(chart1)
        #
        # 2/ Search for common charts via one expression transformation
        #    ----------------------------------------------------------
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if (chart1, chart2) in coord_changes:
                            self.function_chart(chart2, from_chart=chart1)
                            resu.append(chart2)
                        if (chart2, chart1) in coord_changes:
                            other.function_chart(chart1, from_chart=chart2)
                            resu.append(chart1)
        if resu == []:
            return None
        else:
            return resu

    def __call__(self, p, chart=None):
        r"""
        Compute the value of the scalar field at a given point.

        INPUT:
    
        - ``p`` -- point in the scalar field's domain (type: 
          :class:`~sage.geometry.manifolds.point.Point`)
        - ``chart`` -- (default: None) chart in which the coordinates of p 
          are to be considered; if none is provided, a chart in which both p's 
          coordinates and the expression of ``self`` are known is searched, 
          starting from the default chart of self._domain
        
        OUTPUT:

        - value at p 

        EXAMPLES:
        
        """
        if p not in self._manifold: 
            raise ValueError("The point " + str(p) +
                             " does not belong to the " + str(self._manifold))
        if chart is None:
            # A common chart is searched:
            def_chart = self._domain._def_chart
            if def_chart in p._coordinates and def_chart in self._express:
                chart = def_chart
            else:
                for chart_p in p._coordinates:
                    if chart_p in self._express:
                        chart = chart_p
                        break
        if chart is None:
            # A change of coordinates is attempted for p:
            for chart_s in self._express:
                try:
                    p.coord(chart_s)
                    chart = chart_s
                    break
                except ValueError:
                    pass
            else:
                # A change of coordinates is attempted on the scalar field
                # expressions:
                for chart_p in p._coordinates:
                    try:
                        self.function_chart(chart_p)
                        chart = chart_p
                        break
                    except (TypeError, ValueError):
                        pass
        if chart is None:
            raise ValueError("No common chart has been found to evaluate " \
                "the action of " + str(self) + " on the " + str(p) + ".")
        return self._express[chart](*(p._coordinates[chart]))

    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - an exact copy of ``self``
    
        """
        result = self._new_instance()
        for chart in self._express:
            result._express[chart] = + self._express[chart]
        if self._name is not None:
            result._name = '+' + self._name 
        if self._latex_name is not None:
            result._latex_name = '+' + self._latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the negative of ``self``
    
        """
        result = self._new_instance()
        for chart in self._express:
            result._express[chart] = - self._express[chart]
        if self._name is not None:
            result._name = '-' + self._name 
        if self._latex_name is not None:
            result._latex_name = '-' + self._latex_name
        return result


    #########  CommutativeAlgebraElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Scalar field addition. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the addition of ``self`` and 
          ``other``
        
        """
        if isinstance(other, ZeroScalarField):
            return self.copy()
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the addition.")
        dom = self._domain
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart addition:
            result._express[chart] = self._express[chart] + other._express[chart]
        if result.is_zero():
            return dom._zero_scalar_field
        if self._name is not None and other._name is not None:
            result._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '+' + other._latex_name
        return result

    def _sub_(self, other):
        r"""
        Scalar field subtraction. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the subtraction of ``other`` from 
          ``self``

        """
        if isinstance(other, ZeroScalarField):
            return self.copy()
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the subtraction.")
        dom = self._domain
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart subtraction:
            result._express[chart] = self._express[chart] - other._express[chart]
        if result.is_zero():
            return dom._zero_scalar_field
        if self._name is not None and other._name is not None:
            result._name = self._name + '-' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '-' + other._latex_name
        return result


    def _mul_(self, other):
        r"""
        Scalar field multiplication. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the multiplication of ``self`` by 
          ``other``
        
        """
        from utilities import format_mul_txt, format_mul_latex
        if isinstance(other, ZeroScalarField):
            return other
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the multiplication.")
        dom = self._domain
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart multiplication:
            result._express[chart] = self._express[chart] * other._express[chart]
        #!# the following 2 lines could be skipped:
        if result.is_zero():
            return dom._zero_scalar_field
        result._name = format_mul_txt(self._name, '*', other._name)
        result._latex_name = format_mul_latex(self._latex_name, ' ', 
                                             other._latex_name)        
        return result

    def _div_(self, other):
        r"""
        Scalar field division. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the division of ``self`` by 
          ``other``
        
        """
        from utilities import format_mul_txt, format_mul_latex
        if isinstance(other, ZeroScalarField):
            raise ZeroDivisionError("Division of a scalar field by zero.")
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the division.")
        dom = self._domain
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart division:
            result._express[chart] = self._express[chart] / other._express[chart]
        #!# the following 2 lines could be skipped:
        if result.is_zero():
            return dom._zero_scalar_field
        result._name = format_mul_txt(self._name, '/', other._name)
        result._latex_name = format_mul_latex(self._latex_name, '/', 
                                             other._latex_name)
        return result

    def _lmul_(self, number):
        r"""
        Multiplication on the left of a scalar field by a real number. 
        
        INPUT:
        
        - ``number`` -- an element of the ring on which the algebra is defined; 
          mathematically, this should be a real number; here it is a member of
          the symbolic ring SR. 
        
        OUPUT:
        
        - the scalar field ``number*self`` 
        
        """
        if number == 0:
            return self._domain._zero_scalar_field
        result = self.__class__(self._domain)
        for chart, expr in self._express.iteritems():
            result._express[chart] = number * expr
        return result

    def _rmul_(self, number):
        r"""
        Multiplication on the right of a scalar field by a real number. 
        
        INPUT:
        
        - ``number`` -- an element of the ring on which the algebra is defined; 
          mathematically, this should be a real number; here it is a member of
          the symbolic ring SR. 
        
        OUPUT:
        
        - the scalar field ``number*self`` 
        
        """
        return self._lmul_(number) # since the algebra is commutative


    #########  End of CommutativeAlgebraElement arithmetic operators ########


    def exterior_der(self):
        r"""
        Return the exterior derivative of the scalar field. 
                
        NB: for a scalar field, the exterior derivative is nothing but the
        differential of the field. Accordingly, this method simply calls
        :meth:`differential`. 
        
        OUTPUT:
        
        - the 1-form exterior derivative of ``self``. 
        
        EXAMPLES:
        
        Exterior derivative on a 3-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: f = M.scalar_field(cos(x)*z^3 + exp(y)*z^2, name='f')
            sage: df = f.exterior_der() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.view()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            
        Exterior derivative computed on a chart that is not the default one::
        
            sage: c_uvw.<u,v,w> = M.chart()
            sage: g = M.scalar_field(u*v^2*w^3, c_uvw, name='g')
            sage: dg = g.exterior_der() ; dg
            1-form 'dg' on the 3-dimensional manifold 'M'
            sage: dg._components
            {coordinate frame (M, (d/du,d/dv,d/dw)): 1-index components w.r.t. coordinate frame (M, (d/du,d/dv,d/dw))}
            sage: dg.comp(c_uvw.frame())[:, c_uvw]
            [v^2*w^3, 2*u*v*w^3, 3*u*v^2*w^2]
            sage: dg.view(c_uvw.frame(), c_uvw)
            dg = v^2*w^3 du + 2*u*v*w^3 dv + 3*u*v^2*w^2 dw
            
        The exterior derivative is nilpotent::
        
            sage: ddf = df.exterior_der() ; ddf
            2-form 'ddf' on the 3-dimensional manifold 'M'
            sage: ddf == 0
            True
            sage: ddf[:]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: ddg = dg.exterior_der() ; ddg
            2-form 'ddg' on the 3-dimensional manifold 'M'
            sage: ddg == 0
            True

        """
        return self.differential()

    def differential(self):
        r"""
        Return the differential of the scalar field. 
        
        This method simply calls the method :meth:`exterior_der`.  
        
        OUTPUT:
        
        - the 1-form that is the differential of ``self``. 
        
        EXAMPLES:
        
        Differential 1-form on a 3-dimensional manifold::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: f = M.scalar_field(cos(x)*z**3 + exp(y)*z**2, name='f')
            sage: df = f.differential() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.view()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            
        """
        from utilities import format_unop_txt, format_unop_latex
        if self._differential is None:
            # A new computation is necessary:
            rname = format_unop_txt('d', self._name)
            rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
            self._differential = self._domain.one_form(name=rname, 
                                                             latex_name=rlname)
            for chart, f in self._express.iteritems():
                for i in self._manifold.irange():
                    self._differential.add_comp(chart._frame)[i, chart] \
                        = f.diff(i)
        return self._differential

    def lie_der(self, vector):
        r"""
        Computes the Lie derivative with respect to a vector field.
        
        The Lie derivative is stored in the dictionary 
        :attr:`_lie_derivatives`, so that there is no need to 
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile. 
        
        In the present case (scalar field), the Lie derivative is equal to
        the scalar field resulting from the action of the vector field on 
        ``self``. 
        
        INPUT:
        
        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken
          
        OUTPUT:
        
        - the scalar field that is the Lie derivative of ``self`` with 
          respect to ``vector``
          
        EXAMPLES:
        
        Lie derivative on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2*cos(y))
            sage: v = M.vector_field(name='v')
            sage: v[:] = (-y, x)
            sage: f.lie_der(v)
            scalar field on the 2-dimensional manifold 'M'
            sage: f.lie_der(v).expr()
            -x^3*sin(y) - 2*x*y*cos(y)

        Alternative expressions of the Lie derivative of a scalar field::
        
            sage: f.lie_der(v) == v(f)  # the vector acting on f
            True
            sage: f.lie_der(v) == f.differential()(v)  # the differential of f acting on the vector
            True

        A vanishing Lie derivative::
        
            sage: f.set_expr(x^2 + y^2)
            sage: f.lie_der(v).view()
            M --> R
            (x, y) |--> 0

        """
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            res = vector(self)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]         


    def hodge_star(self, metric):
        r"""
        Compute the Hodge dual of the scalar field with respect to some
        pseudo-Riemannian metric. 
        
        If `f` is ``self``, the Hodge dual is the `n`-form
        `*f` defined by (`n` being the manifold's dimension)
        
        .. MATH::
            
            *f = f \epsilon
                
        where `\epsilon` is the volume form associated with some 
        pseudo-Riemannian metric `g` on the manifold. 
        
        INPUT:
        
        - ``metric``: the pseudo-Riemannian metric `g` defining the Hodge dual, 
          via the volume form `\epsilon`; must be an instance of 
          :class:`~sage.geometry.manifolds.metric.Metric`
        
        OUTPUT:
        
        - the `n`-form `*f` 
        
        EXAMPLES:

        Hodge star of a scalar field in the Euclidean space `R^3`::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: f = M.scalar_field(function('F',x,y,z), name='f')
            sage: sf = f.hodge_star(g) ; sf
            3-form '*f' on the 3-dimensional manifold 'M'
            sage: sf.view()
            *f = F(x, y, z) dx/\dy/\dz
            sage: ssf = sf.hodge_star(g) ; ssf
            scalar field '**f' on the 3-dimensional manifold 'M'
            sage: ssf.view()
            **f: M --> R
               (x, y, z) |--> F(x, y, z)
            sage: ssf == f # must hold for a Riemannian metric
            True
        
        """
        from utilities import format_unop_txt, format_unop_latex
        eps = metric.volume_form()
        dom_resu = self._domain.intersection(eps._domain) # result domain
        resu = self.restrict(dom_resu) * eps.restrict(dom_resu)
        if self._name is None:
            resu_name = None
        else:
            resu_name = format_unop_txt('*', self._name)
        if self._latex_name is None:
            resu_latex_name = None
        else:
            resu_latex_name = format_unop_latex(r'\star ', self._latex_name)
        resu.set_name(name=resu_name, latex_name=resu_latex_name)
        return resu


#******************************************************************************

class ZeroScalarField(ScalarField):
    r"""
    Null scalar field on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the scalar field is defined
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    Zero scalar field on a 2-dimensional manifold::
    
        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(2, 'M')                  
        sage: c_xy.<x,y> = M.chart()
        sage: from sage.geometry.manifolds.scalarfield import ZeroScalarField
        sage: f = ZeroScalarField(M) ; f
        zero scalar field on the 2-dimensional manifold 'M'
        sage: f.expr()
        0
        sage: f.is_zero()
        True
        sage: p = M.point((1,2))
        sage: f(p)
        0

    Each manifold has a predefined zero scalar field::
    
        sage: M._zero_scalar_field
        zero scalar field on the 2-dimensional manifold 'M'
        sage: M._zero_scalar_field(p)
        0
        sage: f == M._zero_scalar_field
        True

    Arithmetics with another instance of :class:`ZeroScalarField`::
    
        sage: h = ZeroScalarField(M)
        sage: s = f+h ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f-h ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f*h ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f/h ; s ; s.expr()
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a scalar field by zero.
        
    Arithmetics with a non-zero instance of :class:`ScalarField`::
    
        sage: g = M.scalar_field(x+y)
        sage: s = f+g ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = g+f ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = f-g ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        -x - y
        sage: s = g-f ; s ; s.expr()                     
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = f*g ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = g*f ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f/g ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = g/f ; s ; s.expr()
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a scalar field by zero.

    """
    def __init__(self, domain, name=None, latex_name=None):
        ScalarField.__init__(self, domain, name=name, latex_name=latex_name)

    ####### Required methods for an algebra element (beside arithmetic) #######
    
    def __nonzero__(self):
        r"""
        Always return False (since ``self`` is zero!). 
        
        This method is called by self.is_zero(). 

        """
        return False

    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a scalar field
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if not isinstance(other, ScalarField):
            try:
                other = self.parent()(other)    # conversion to a scalar field
            except TypeError:
                return False
        if other._domain != self._domain:
            return False
        return other.is_zero()
    
    def __ne__(self, other):
        r"""
        Non-equality operator.
        """
        return not self.__eq__(other)
        
    ####### End of required methods for an algebra element (beside arithmetic) #######

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "zero scalar field"
        if self._name is not None:
            description += " '%s'" % self._name
        description += " on the " + str(self._domain)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`ZeroScalarField` instance on the same domain.
        
        """
        return ZeroScalarField(self._domain)        

    def copy(self):
        r"""
        Return an exact copy of ``self``.
        """
        return ZeroScalarField(self._domain)
        
    def function_chart(self, chart=None):
        r""" 
        Return the function of the coordinates representing the scalar field 
        in a given chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) chart; if None, the domain's default 
          chart will be used
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.chart.ZeroFunctionChart` defined in 
          the specified chart.
        
        """
        if chart is None:
            chart = self._domain._def_chart
        return chart._zero_function
 
    def expr(self, chart=None, from_chart=None):
        r""" 
        Return the coordinate expression of the scalar field in a given 
        chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) unused here
        - ``from_chart`` -- (default: None) unused here
                  
        OUTPUT:
        
        - number zero
        
        """
        return 0
 
    def set_expr(self, coord_expression, chart=None):
        r"""
        Set some coordinate expression of the scalar field.
        
        Not valid for a :class:`ZeroScalarField` object. 
        """
        raise TypeError("set_expr() has no meaning for a zero scalar field.")

    def add_expr(self, coord_expression, chart=None):
        r"""
        Add some coordinate expression to the scalar field.
        
        Not valid for a :class:`ZeroScalarField` object. 
        """
        raise TypeError("add_expr() has no meaning for a zero scalar field.")

    def __call__(self, p):
        r"""
        Computes the image of a point.

        INPUT:
    
        - ``p`` -- point on the manifold (type: 
          :class:`~sage.geometry.manifolds.point.Point`)
        
        OUTPUT:

        - the number zero. 
        
        """
        from point import Point
        if not isinstance(p, Point):
            return TypeError("The argument must be a point.")
        return 0
                
    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - ``self``
    
        """
        return self

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - ``self`` (since ``self`` is zero)
    
        """
        return self


    #########  CommutativeAlgebraElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Scalar field addition. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the addition of ``self`` and 
          ``other``
        
        """
        return other.copy()    
            
    def _sub_(self, other):
        r"""
        Scalar field subtraction. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the subtraction of ``other`` from 
          ``self``          
        
        """
        return -other    

    def _mul_(self, other):
        r"""
        Scalar field multiplication. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the multiplication of ``self`` by 
          ``other``
        
        """
        return self

    def _div_(self, other):
        r"""
        Scalar field division. 
        
        INPUT:
        
        - ``other`` -- a scalar field (in the same algebra as self)
        
        OUPUT:
        
        - the scalar field resulting from the division of ``self`` by 
          ``other``
        
        """
        if other == 0:
            raise ZeroDivisionError("Division of a scalar field by zero.")
        else:
            return self

    def _lmul_(self, number):
        r"""
        Multiplication on the left of the scalar field by a real number. 
        
        INPUT:
        
        - ``number`` -- an element of the ring on which the algebra is defined; 
          mathematically, this should be a real number; here it is a member of
          the symbolic ring SR. 
        
        OUPUT:
        
        - the scalar field ``number*self`` 
        
        """
        return self


    #########  End of CommutativeAlgebraElement arithmetic operators ########

    def differential(self):
        r"""
        Return the exterior derivative of the scalar field, which is zero in 
        the present case. 
                
        OUTPUT:
        
        - the (vanishing) 1-form differential of ``self``. 
                
        """
        from utilities import format_unop_txt, format_unop_latex
        if self._differential is None:
            # A new computation is necessary:
            rname = format_unop_txt('d', self._name)
            rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
            self._differential = self._domain.one_form(name=rname, 
                                                             latex_name=rlname)
            for chart in self._domain._atlas:
                self._differential.add_comp(chart._frame) # since a newly
                                            # created set of components is zero
        return self._differential

