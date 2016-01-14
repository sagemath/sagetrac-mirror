"""
Regular polygons in the upper half model for hyperbolic plane

AUTHORS:

- Javier Honrubia (2016-01)
"""
#******************************************************************************
#       Copyright (C) 2016 Javier Honrubia Gonzalez <jhonrubia6@alumno.uned.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.plot.hyperbolic_polygon import HyperbolicPolygon,hyperbolic_polygon
from sage.rings.all import CC
from sage.plot.misc import options, rename_keyword
from sage.symbolic.constants import pi,e
from sage.functions.hyperbolic import arccosh
from sage.functions.trig import sin,cos,cot
from sage.misc.functional import is_odd
from sage.matrix.constructor import matrix

class HyperbolicRegularPolygon(HyperbolicPolygon):
    """
    Primitive class for regular hyberbolic polygon type.

    See ``hyperbolic_regular_polygon?`` for information about plotting a hyperbolic
    regular polygon in the upper complex halfplane.

    INPUT:

    - ``sides`` -- number of sides of the polygon

    - ``i_angle`` -- interior angle of the polygon
    
    - ``center``-- center point as a complex number of the polygon
    

    EXAMPLES:

    Note that constructions should use :func:`hyperbolic_regular_polygon`::

         sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
         sage: print HyperbolicRegularPolygon(5,pi/2,I, {})
         Hyperbolic regular polygon (sides=5,i_angle=1/2*pi,center=1.00000000000000*I)
         
    The code verifies is there exists a compact hyperbolic regular polygon with the given 
    data, checking
    
    .. MATH::
    
        A(\mathcal{P})=\pi(sides-2) - sides \cdot i\_angle>0
         
    and giving an error if the i_angle is less than the minimum to generate a compact 
    polygon ::

        sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
        sage: P=HyperbolicRegularPolygon(4,pi/2,I, {})
        Traceback (most recent call last):
        ...
        ValueError: There exists no hyperbolic regular compact polygon, for sides=4 the interior angle must be less than 1/2*pi
    
    It is an error to give a center outside the upper half plane in this model ::
    
         sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
         sage: P=HyperbolicRegularPolygon(4,pi/4,1-I, {})
         Traceback (most recent call last):
         ...
         ValueError: center:1.00000000000000 - 1.00000000000000*I is not a valid point in the upper half plane model of the hyperbolic plane
    
    TESTS::
    
         sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
         sage: P=HyperbolicRegularPolygon(4,-pi/4,I,{})
         Traceback (most recent call last):
         ...
         ValueError: Interior angle -1/4*pi must be in (0,pi) interval
         
         sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
         sage: P=HyperbolicRegularPolygon(16,3*pi/2,I,{})
         Traceback (most recent call last):
         ...
         ValueError: Interior angle 3/2*pi must be in (0,pi) interval   
         
         sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
         sage: P=HyperbolicRegularPolygon(2,pi/10,I,{})
         Traceback (most recent call last):
         ...
         ValueError: Degenerated polygons (sides<=2) are not supported
         
    """
    def __init__(self, sides, i_angle, center, options):
        """
        Initialize HyperbolicRegularPolygon.

        EXAMPLES::

            sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
            sage: print HyperbolicRegularPolygon(5,pi/2,I, {})
            Hyperbolic regular polygon (sides=5,i_angle=1/2*pi,center=1.00000000000000*I)
        """
        self.center = CC(center)
        if self.center.imag() <= 0 :
            raise ValueError("center:%s is not a valid point in the upper half plane model of the hyperbolic plane"%(self.center))
        if sides < 3 :
            raise ValueError("Degenerated polygons (sides<=2) are not supported")
        if i_angle <=0 or i_angle >= pi:
            raise ValueError("Interior angle %s must be in (0,pi) interval"%(i_angle))
        if pi*(sides-2) - sides*i_angle <= 0 :
            raise ValueError("There exists no hyperbolic regular compact polygon, for sides=%s the interior angle must be less than %s"%(sides, pi * (sides-2) / sides))
        self.sides = sides
        self.i_angle = i_angle
        beta = 2 * pi / self.sides # compute the rotation angle to be used ahead
        alpha = self.i_angle / 2
        I = CC(0,1)
        # compute using cosine theorem the radius of the circumscribed circle using 
        # the triangle formed by the radius and the three known angles
        r = arccosh(cot(alpha) * (1+cos(beta)) / sin(beta)) 
    
        z_0 = [I*(e**r).n(digits=8)]  # The first point will be always on the imaginary axis
                                     # limited to 8 digits for efficiency in the subsequent
                                     # calculations. 
    
        scale = self.center.imag()   # Compute the dilation isometry used to move the center
                                     # from I to the imaginary part of the given center.
    
        h_disp = self.center.real()    # Compute the parabolic isometry to move the center to the
                                     # real part of the given center.                                
    
      
        d_z_k = [z_0[0]*scale + h_disp]  #d_k has the points for the polygon in the given center
        z_k = z_0                      #z_k has the Re(z)>0 vertices for the I centered polygon 
        r_z_k = []                     #r_z_k has the Re(z)<0 vertices
        if is_odd(self.sides):
            vert = (self.sides-1) / 2
        else:
            vert = self.sides/2 - 1
        for k in range(0, vert):
            new_z_k = self._i_rotation(z_k[-1], beta).n(digits=8) #compute with 8 digits 
                                                                #to accelerate calculations
            z_k = z_k + [new_z_k]
            d_z_k = d_z_k + [new_z_k*scale+h_disp]
            r_z_k=[-(new_z_k).conjugate()*scale+h_disp] + r_z_k
        if is_odd(self.sides) :
            HyperbolicPolygon.__init__(self, d_z_k+r_z_k, options)
        else:
            z_opo = [I*(e**(-r)).n(digits=8)*scale + h_disp]
            HyperbolicPolygon.__init__(self, d_z_k+z_opo+r_z_k, options)
    def _repr_(self):
        """
        String representation of HyperbolicRegularPolygon.

        TESTS::

            sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
            sage: HyperbolicRegularPolygon(5,pi/2,I, {})._repr_()
            'Hyperbolic regular polygon (sides=5,i_angle=1/2*pi,center=1.00000000000000*I)'
        """
        return "Hyperbolic regular polygon (sides=%s,i_angle=%s,center=%s)" \
                %(self.sides,self.i_angle,self.center)

    def _i_rotation(self, z, alpha):
        """
        Function to calculate the resulting point after applying an hyperbolic rotation
        centered in 0+I and angle alpha.
        
        INPUT:
        
        - ``z``-- point in the upper complex halfplane to which apply the isometry
        
        - ``alpha``-- angle of rotation (radians,counterwise)
        
        OUTPUT:
        
        - rotated point in the upper complex halfplane
        
        TESTS::
        
        sage: from sage.plot.hyperbolic_regular_polygon import HyperbolicRegularPolygon
        sage: P=HyperbolicRegularPolygon(4,pi/4,1+I, {})
        sage: P._i_rotation(2+I,pi/2)
        I - 2
        
        """
        _a = alpha / 2
        _c = cos(_a)
        _s = sin(_a)
        G = matrix([[_c, _s], [-_s, _c]])
        return (G[0][0]*z+G[0][1])/(G[1][0]*z+G[1][1])
        
@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor="blue", zorder=2,
         linestyle='solid')
def hyperbolic_regular_polygon(sides, i_angle, center=CC(0,1), **options):
    """
    Return a hyperbolic regular polygon in the upper half model of Hyperbolic plane
        given the number of sides, interior angle and possibly a center.

    Type ``?hyperbolic_regular_polygon`` to see all options.

    INPUT:

    - ``sides`` -- number of sides of the polygon

    - ``i_angle`` -- interior angle of the polygon
    
    - ``center``  -- (default:I) hyperbolic center point (complex number) of the polygon

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` -- (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic regular polygon with 6 sides and square angles ::

        sage: from sage.plot.hyperbolic_regular_polygon import hyperbolic_regular_polygon
        sage: g = hyperbolic_regular_polygon(6,pi/2)
        sage: g.plot()
        Graphics object consisting of 1 graphics primitive
        
    .. PLOT::
         
         from sage.plot.hyperbolic_regular_polygon import hyperbolic_regular_polygon
         g=hyperbolic_regular_polygon(6,pi/2)
         sphinx_plot(g.plot())
         
    With more options::
    
        sage: from sage.plot.hyperbolic_regular_polygon import hyperbolic_regular_polygon
        sage: g= hyperbolic_regular_polygon(6,pi/2,center=3+2*I, fill=True, color='red')
        sage: g.plot()
        Graphics object consisting of 1 graphics primitive
        
    .. PLOT::
         
         from sage.plot.hyperbolic_regular_polygon import hyperbolic_regular_polygon
         g=hyperbolic_regular_polygon(6,pi/2,center=3+2*I, fill=True, color='red')
         sphinx_plot(g.plot())
    
    The code verifies is there exists a hyperbolic regular polygon with the given data, 
    checking
    
    .. MATH::
    
        A(\mathcal{P})=\pi(sides-2) - sides \cdot i\_angle>0
         
    giving an error if the i_angle is less than the minimum to generate a compact polygon ::

        sage: from sage.plot.hyperbolic_regular_polygon import hyperbolic_regular_polygon
        sage: hyperbolic_regular_polygon(4,pi/2)
        Traceback (most recent call last):
        ...
        ValueError: There exists no hyperbolic regular compact polygon, for sides=4 the interior angle must be less than 1/2*pi
    
    It is an error to give a center outside the upper half plane in this model ::
    
         sage: from sage.plot.hyperbolic_regular_polygon import hyperbolic_regular_polygon
         sage: hyperbolic_regular_polygon(4,pi/4,1-I)
         Traceback (most recent call last):
         ...
         ValueError: center:1.00000000000000 - 1.00000000000000*I is not a valid point in the upper half plane model of the hyperbolic plane

    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    g.add_primitive(HyperbolicRegularPolygon(sides, i_angle, center, options))
    g.set_aspect_ratio(1)
    return g