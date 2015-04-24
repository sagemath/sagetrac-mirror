"""
Sage Interacts

Sage interacts are applications of the `@interact decorator <../../sagenb/notebook/interact.html>`_.
They are conveniently accessible in the Sage Notebook via ``interacts.[TAB].[TAB]()``.
The first ``[TAB]`` lists categories and the second ``[TAB]`` reveals the interact examples.

EXAMPLE:

Invoked in the notebook, the following command will produce the fully formatted
interactive mathlet.  In the command line, it will simply return the underlying
HTML and Sage code which creates the mathlet::

    sage: interacts.calculus.taylor_polynomial()
    <html>...</html>

AUTHORS:

- William Stein

- Harald Schilly, Robert Marik (2011-01-16): added many examples (#9623) partially based on work by Lauri Ruotsalainen

"""

#*****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.all import *
from sagenb.notebook.interact import *

x = SR.var('x')

def library_interact(f):
    """
    This is a decorator for using interacts in the Sage library.

    EXAMPLES::

        sage: import sage.interacts.library as library
        sage: @library.library_interact
        ... def f(n=5): print n
        ...
        sage: f()  # an interact appears, if using the notebook, else code
        <html>...</html>
    """
    @sage_wraps(f)
    def library_wrapper():
       # Maybe program around bug (?) in the notebook:
       html("</pre>")
       # This prints out the relevant html code to make
       # the interact appear:
       interact(f)
    return library_wrapper

@library_interact
def demo(n=tuple(range(10)), m=tuple(range(10))):
    """
    This is a demo interact that sums two numbers.

    INPUT:

        - ``n`` -- integer slider
        - ``m`` -- integer slider

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.demo()
        <html>...</html>
    """
    print n+m

@library_interact
def taylor_polynomial(
    title = text_control('<h2>Taylor polynomial</h2>'),
    f=input_box(sin(x)*exp(-x),label="$f(x)=$"), order=slider(range(1,13))):
    """
    An interact which illustrates the Taylor polynomial approximation
    of various orders around `x=0`.

        - ``f`` -- function expression
        - ```order``` -- integer slider

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.taylor_polynomial()
        <html>...</html>
    """
    x0  = 0
    p   = plot(f,(x,-1,5), thickness=2)
    dot = point((x0,f(x=x0)),pointsize=80,rgbcolor=(1,0,0))
    ft = f.taylor(x,x0,order)
    pt = plot(ft,(-1, 5), color='green', thickness=2)
    html('$f(x)\;=\;%s$'%latex(f))
    html('$\hat{f}(x;%s)\;=\;%s+\mathcal{O}(x^{%s})$'%(x0,latex(ft),order+1))
    show(dot + p + pt, ymin = -.5, ymax = 1)

@library_interact
def definite_integral(
    title = text_control('<h2>Definite integral</h2>'),
    f = input_box(default = "3*x", label = '$f(x)=$'),
    g = input_box(default = "x^2", label = '$g(x)=$'),
    interval = range_slider(-10,10,default=(0,3), label="Interval"),
    x_range = range_slider(-10,10,default=(0,3), label = "plot range (x)"),
    selection = selector(["f", "g", "f and g", "f - g"], default="f and g", label="Select")):
    """
    This is a demo interact for plotting the definite integral of a function
    based on work by Lauri Ruotsalainen, 2010.

    INPUT:

        - ``function`` -- input box, function in x
        - ``interval`` -- interval for the definite integral
        - ``x_range`` -- range slider for plotting range
        - ``selection`` -- selector on how to visualize the integrals

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.definite_integral()
        <html>...</html>
    """
    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    g = symbolic_expression(g).function(x)
    f_plot = Graphics(); g_plot = Graphics(); h_plot = Graphics();
    text = ""

    # Plot function f.
    if selection != "g":
        f_plot = plot(f(x), x, x_range, color="blue", thickness=1.5)

    # Color and calculate the area between f and the horizontal axis.
    if selection == "f" or selection == "f and g":
        f_plot += plot(f(x), x, interval, color="blue", fill=True, fillcolor="blue", fillalpha=0.15)
        text += r"$\int_{%.2f}^{%.2f}(\color{Blue}{f(x)})\,\mathrm{d}x=\int_{%.2f}^{%.2f}(%s)\,\mathrm{d}x=%.2f$" % (
            interval[0], interval[1],
            interval[0], interval[1],
            latex(f(x)),
            f(x).nintegrate(x, interval[0], interval[1])[0]
        )

    if selection == "f and g":
        text += r"<br/>"

    # Plot function g. Also color and calculate the area between g and the horizontal axis.
    if selection == "g" or selection == "f and g":
        g_plot = plot(g(x), x, x_range, color="green", thickness=1.5)
        g_plot += plot(g(x), x, interval, color="green", fill=True, fillcolor="yellow", fillalpha=0.5)
        text += r"$\int_{%.2f}^{%.2f}(\color{Green}{g(x)})\,\mathrm{d}x=\int_{%.2f}^{%.2f}(%s)\,\mathrm{d}x=%.2f$" % (
            interval[0], interval[1],
            interval[0], interval[1],
            latex(g(x)),
            g(x).nintegrate(x, interval[0], interval[1])[0]
        )

    # Plot function f-g. Also color and calculate the area between f-g and the horizontal axis.
    if selection == "f - g":
        g_plot = plot(g(x), x, x_range, color="green", thickness=1.5)
        g_plot += plot(g(x), x, interval, color="green", fill=f(x), fillcolor="red", fillalpha=0.15)
        h_plot = plot(f(x)-g(x), x, interval, color="red", thickness=1.5, fill=True, fillcolor="red", fillalpha=0.15)
        text = r"$\int_{%.2f}^{%.2f}(\color{Red}{f(x)-g(x)})\,\mathrm{d}x=\int_{%.2f}^{%.2f}(%s)\,\mathrm{d}x=%.2f$" % (
            interval[0], interval[1],
            interval[0], interval[1],
            latex(f(x)-g(x)),
            (f(x)-g(x)).nintegrate(x, interval[0], interval[1])[0]
        )

    show(f_plot + g_plot + h_plot, gridlines=True)
    html(text)

@library_interact
def function_derivative(
    title = text_control('<h2>Derivative grapher</h2>'),
    function = input_box(default="x^5-3*x^3+1", label="Function:"),
    x_range  = range_slider(-15,15,0.1, default=(-2,2), label="Range (x)"),
    y_range  = range_slider(-15,15,0.1, default=(-8,6), label="Range (y)")):
    """
    This is a demo interact for plotting derivatives of a function based on work by
    Lauri Ruotsalainen, 2010.

    INPUT:

        - ``function`` -- input box, function in x
        - ``x_range`` -- range slider for plotting range
        - ``y_range`` -- range slider for plotting range

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.function_derivative()
        <html>...</html>
    """
    x = SR.var('x')
    f = symbolic_expression(function).function(x)
    df = derivative(f, x)
    ddf = derivative(df, x)
    plots = plot(f(x), x_range, thickness=1.5) + plot(df(x), x_range, color="green") + plot(ddf(x), x_range, color="red")
    if y_range == (0,0):
        show(plots, xmin=x_range[0], xmax=x_range[1])
    else:
        show(plots, xmin=x_range[0], xmax=x_range[1], ymin=y_range[0], ymax=y_range[1])

    html("<center>$\color{Blue}{f(x) = %s}$</center>"%latex(f(x)))
    html("<center>$\color{Green}{f'(x) = %s}$</center>"%latex(df(x)))
    html("<center>$\color{Red}{f''(x) = %s}$</center>"%latex(ddf(x)))

@library_interact
def difference_quotient(
    title = text_control('<h2>Difference quotient</h2>'),
    f = input_box(default="sin(x)", label='f(x)'),
    interval= range_slider(0, 10, 0.1, default=(0.0,10.0), label="Range"),
    a = slider(0, 10, None, 5.5, label = '$a$'),
    x0 = slider(0, 10, None, 2.5, label = '$x_0$ (start point)')):
    """
    This is a demo interact for difference quotient based on work by
    Lauri Ruotsalainen, 2010.

    INPUT:

        - ``f`` -- input box, function in `x`
        - ``interval`` -- range slider for plotting
        - ``a`` -- slider for `a`
        - ``x0`` -- slider for starting point `x_0`

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.difference_quotient()
        <html>...</html>
    """
    html('<h2>Difference Quotient</h2>')
    html('<div style="white-space: normal;">\
         <a href="http://en.wikipedia.org/wiki/Difference_quotient" target="_blank">\
         Wikipedia article about difference quotient</a></div>'
         )

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    fmax = f.find_maximum_on_interval(interval[0], interval[1])[0]
    fmin = f.find_minimum_on_interval(interval[0], interval[1])[0]
    f_height = fmax - fmin
    measure_y = fmin - 0.1*f_height

    measure_0 = line2d([(x0, measure_y), (a, measure_y)], rgbcolor="black")
    measure_1 = line2d([(x0, measure_y + 0.02*f_height), (x0, measure_y-0.02*f_height)], rgbcolor="black")
    measure_2 = line2d([(a, measure_y + 0.02*f_height), (a, measure_y-0.02*f_height)], rgbcolor="black")
    text_x0 = text("x0", (x0, measure_y - 0.05*f_height), rgbcolor="black")
    text_a = text("a", (a, measure_y - 0.05*f_height), rgbcolor="black")
    measure = measure_0 + measure_1 + measure_2 + text_x0 + text_a

    tanf = symbolic_expression((f(x0)-f(a))*(x-a)/(x0-a)+f(a)).function(x)

    fplot = plot(f(x), x, interval[0], interval[1])
    tanplot = plot(tanf(x), x, interval[0], interval[1], rgbcolor="#FF0000")
    points = point([(x0, f(x0)), (a, f(a))], pointsize=20, rgbcolor="#005500")
    dashline = line2d([(x0, f(x0)), (x0, f(a)), (a, f(a))], rgbcolor="#005500", linestyle="--")
    html('<h2>Difference Quotient</h2>')
    show(fplot + tanplot + points + dashline + measure, xmin=interval[0], xmax=interval[1], ymin=fmin-0.2*f_height, ymax=fmax)
    html(r"<br>$\text{Line's equation:}$")
    html(r"$y = %s$<br>"%tanf(x))
    html(r"$\text{Slope:}$")
    html(r"$k = \frac{f(x_0)-f(a)}{x_0-a} = %s$<br>" % (N(derivative(tanf(x), x), digits=5)))

@library_interact
def quadratic_equation(A = slider(-7, 7, 1, 1), B = slider(-7, 7, 1, 1), C = slider(-7, 7, 1, -2)):
    """
    This is a demo interact for solving quadratic equations based on work by
    Lauri Ruotsalainen, 2010.

    INPUT:

        - ``A`` -- integer slider
        - ``B`` -- integer slider
        - ``C`` -- integer slider

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.quadratic_equation()
        <html>...</html>
    """
    x = SR.var('x')
    f = symbolic_expression(A*x**2 + B*x + C).function(x)
    html('<h2>The Solutions of the Quadratic Equation</h2>')
    html("$%s = 0$" % f(x))

    show(plot(f(x), x, (-10, 10), ymin=-10, ymax=10), aspect_ratio=1, figsize=4)

    d = B**2 - 4*A*C

    if d < 0:
        color = "Red"
        sol = r"\text{solution} \in \mathbb{C}"
    elif d == 0:
        color = "Blue"
        sol = -B/(2*A)
    else:
        color = "Green"
        a = (-B+sqrt(B**2-4*A*C))/(2*A)
        b = (-B-sqrt(B**2-4*A*C))/(2*A)
        sol = r"\begin{cases}%s\\%s\end{cases}" % (latex(a), latex(b))

    if B < 0:
        dis1 = "(%s)^2-4*%s*%s" % (B, A, C)
    else:
        dis1 = "%s^2-4*%s*%s" % (B, A, C)
    dis2 = r"\color{%s}{%s}" % (color, d)

    html("$Ax^2 + Bx + C = 0$")
    calc = r"$x = \frac{-B\pm\sqrt{B^2-4AC}}{2A} = " + \
           r"\frac{-%s\pm\sqrt{%s}}{2*%s} = " + \
           r"\frac{-%s\pm\sqrt{%s}}{%s} = %s$"
    html(calc % (B, dis1, A, B, dis2, (2*A), sol))

@library_interact
def trigonometric_properties_triangle(
    a0 = slider(0, 360, 1, 30, label="A"),
    a1 = slider(0, 360, 1, 180, label="B"),
    a2 = slider(0, 360, 1, 300, label="C")):
    """
    This is an interact for demonstrating trigonometric properties
    in a triangle based on work by Lauri Ruotsalainen, 2010.

    INPUT:

        - ``a0`` -- angle
        - ``a1`` -- angle
        - ``a2`` -- angle

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.geometry.trigonometric_properties_triangle()
        <html>...</html>
    """
    import math

    # Returns the distance between points (x1,y1) and (x2,y2)
    def distance(x1_y1, x2_y2):
        (x1, y1) = x1_y1
        (x2, y2) = x2_y2
        return sqrt((x2-x1)**2 + (y2-y1)**2)

    # Returns an angle (in radians) when sides a and b
    # are adjacent and the side c is opposite to the angle
    def angle(a, b, c):
        a,b,c = map(float,[a,b,c])
        return acos((b**2 + c**2 - a**2)/(2.0*b*c))

    # Returns the area of a triangle when an angle alpha
    # and adjacent sides a and b are known
    def area(alpha, a, b):
        return 1.0/2.0*a*b*sin(alpha)

    xy = [0]*3
    html('<h2>Trigonometric Properties of a Triangle</h2>')
    # Coordinates of the angles
    a = map(lambda x : math.radians(float(x)), [a0, a1, a2])
    for i in range(3):
        xy[i] = (cos(a[i]), sin(a[i]))

    # Side lengths (bc, ca, ab) corresponding to triangle vertices (a, b, c)
    al = [distance(xy[1], xy[2]), distance(xy[2], xy[0]), distance(xy[0], xy[1])]

    # The angles (a, b, c) in radians
    ak = [angle(al[0], al[1], al[2]), angle(al[1], al[2], al[0]), angle(al[2], al[0], al[1])]

    # The area of the triangle
    A = area(ak[0], al[1], al[2])

    unit_circle = circle((0, 0), 1, aspect_ratio=1)

    # Triangle
    triangle = line([xy[0], xy[1], xy[2], xy[0]], rgbcolor="black")
    triangle_points = point(xy, pointsize=30)

    # Labels of the angles drawn in a distance from points
    a_label = text("A", (xy[0][0]*1.07, xy[0][1]*1.07))
    b_label = text("B", (xy[1][0]*1.07, xy[1][1]*1.07))
    c_label = text("C", (xy[2][0]*1.07, xy[2][1]*1.07))
    labels = a_label + b_label + c_label

    show(unit_circle + triangle + triangle_points + labels, figsize=[5, 5], xmin=-1, xmax=1, ymin=-1, ymax=1)
    angl_txt = r"$\angle A = {%s}^{\circ},$ $\angle B = {%s}^{\circ},$ $\angle C = {%s}^{\circ}$" % (
            math.degrees(ak[0]),
            math.degrees(ak[1]),
            math.degrees(ak[2])
        )
    html(angl_txt)
    html(r"$AB = %s,$  $BC = %s,$  $CA = %s$"%(al[2], al[0], al[1]))
    html(r"Area of triangle $ABC = %s$"%A)

@library_interact
def unit_circle(
    function = selector([(0, sin(x)), (1, cos(x)), (2, tan(x))]),
    x = slider(0,2*pi, 0.005*pi, 0)):
    """
    This is an interact for Sin, Cos and Tan in the Unit Circle
    based on work by Lauri Ruotsalainen, 2010.

    INPUT:

        - ``function`` -- select Sin, Cos or Tan
        - ``x`` -- slider to select angle in unit circle

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.geometry.unit_circle()
        <html>...</html>
    """
    xy = (cos(x), sin(x))
    t = SR.var('t')
    html('<div style="white-space: normal;">Lines of the same color have\
         the same length</div>')

    # Unit Circle
    C = circle((0, 0), 1, figsize=[5, 5], aspect_ratio=1)
    C_line = line([(0, 0), (xy[0], xy[1])], rgbcolor="black")
    C_point = point((xy[0], xy[1]), pointsize=40, rgbcolor="green")
    C_inner = parametric_plot((cos(t), sin(t)), (t, 0, x + 0.001), color="green", thickness=3)
    C_outer = parametric_plot((0.1 * cos(t), 0.1 * sin(t)), (t, 0, x + 0.001), color="black")
    C_graph = C + C_line + C_point + C_inner + C_outer

    # Graphics related to the graph of the function
    G_line = line([(0, 0), (x, 0)], rgbcolor="green", thickness=3)
    G_point = point((x, 0), pointsize=30, rgbcolor="green")
    G_graph = G_line + G_point

    # Sine
    if function == 0:
        Gf = plot(sin(t), t, 0, 2*pi, axes_labels=("x", "sin(x)"))
        Gf_point = point((x, sin(x)), pointsize=30, rgbcolor="red")
        Gf_line = line([(x, 0),(x, sin(x))], rgbcolor="red")
        Cf_point = point((0, xy[1]), pointsize=40, rgbcolor="red")
        Cf_line1 = line([(0, 0), (0, xy[1])], rgbcolor="red", thickness=3)
        Cf_line2 = line([(0, xy[1]), (xy[0], xy[1])], rgbcolor="purple", linestyle="--")
    # Cosine
    elif function == 1:
        Gf = plot(cos(t), t, 0, 2*pi, axes_labels=("x", "cos(x)"))
        Gf_point = point((x, cos(x)), pointsize=30, rgbcolor="red")
        Gf_line = line([(x, 0), (x, cos(x))], rgbcolor="red")
        Cf_point = point((xy[0], 0), pointsize=40, rgbcolor="red")
        Cf_line1 = line([(0, 0), (xy[0], 0)], rgbcolor="red", thickness=3)
        Cf_line2 = line([(xy[0], 0), (xy[0], xy[1])], rgbcolor="purple", linestyle="--")
    # Tangent
    else:
        Gf = plot(tan(t), t, 0, 2*pi, ymin=-8, ymax=8, axes_labels=("x", "tan(x)"))
        Gf_point = point((x, tan(x)), pointsize=30, rgbcolor="red")
        Gf_line = line([(x, 0), (x, tan(x))], rgbcolor="red")
        Cf_point = point((1, tan(x)), pointsize=40, rgbcolor="red")
        Cf_line1 = line([(1, 0), (1, tan(x))], rgbcolor="red", thickness=3)
        Cf_line2 = line([(xy[0], xy[1]), (1, tan(x))], rgbcolor="purple", linestyle="--")

    C_graph += Cf_point + Cf_line1 + Cf_line2
    G_graph += Gf + Gf_point + Gf_line

    html.table([[r"$\text{Unit Circle}$",r"$\text{Function}$"], [C_graph, G_graph]], header_row=True)

@library_interact
def special_points(
    title = text_control('<h2>Special points in triangle</h2>'),
    a0 = slider(0, 360, 1, 30, label="A"),
    a1 = slider(0, 360, 1, 180, label="B"),
    a2 = slider(0, 360, 1, 300, label="C"),
    show_median = checkbox(False, label="Medians"),
    show_pb = checkbox(False, label="Perpendicular Bisectors"),
    show_alt = checkbox(False, label="Altitudes"),
    show_ab = checkbox(False, label="Angle Bisectors"),
    show_incircle = checkbox(False, label="Incircle"),
    show_euler = checkbox(False, label="Euler's Line")):
    """
    This interact demo shows special points in a triangle
    based on work by Lauri Ruotsalainen, 2010.

    INPUT:

        - ``a0`` -- angle
        - ``a1`` -- angle
        - ``a2`` -- angle
        - ``show_median`` -- checkbox
        - ``show_pb`` -- checkbox to show perpendicular bisectors
        - ``show_alt`` -- checkbox to show altitudes
        - ``show_ab`` -- checkbox to show angle bisectors
        - ``show_incircle`` -- checkbox to show incircle
        - ``show_euler`` -- checkbox to show euler's line

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.geometry.special_points()
        <html>...</html>
    """
    import math
    # Return the intersection point of the bisector of the angle <(A[a],A[c],A[b]) and the unit circle. Angles given in radians.
    def half(A, a, b, c):
        if (A[a] < A[b] and (A[c] < A[a] or A[c] > A[b])) or (A[a] > A[b] and (A[c] > A[a] or A[c] < A[b])):
            p = A[a] + (A[b] - A[a]) / 2.0
        else:
            p = A[b] + (2*pi - (A[b]-A[a])) / 2.0
        return (math.cos(p), math.sin(p))

    # Returns the distance between points (x1,y1) and (x2,y2)
    def distance(x1_y1, x2_y2):
        (x1, y1) = x1_y1
        (x2, y2) = x2_y2
        return math.sqrt((x2-x1)**2 + (y2-y1)**2)

    # Returns the line (graph) going through points (x1,y1) and (x2,y2)
    def line_to_points(x1_y1, x2_y2, **plot_kwargs):
        (x1, y1) = x1_y1 
        (x2, y2) = x2_y2
        return plot((y2-y1) / (x2-x1) * (x-x1) + y1, (x,-3,3), **plot_kwargs)

    # Coordinates of the angles
    a = map(lambda x : math.radians(float(x)), [a0, a1, a2])
    xy = [(math.cos(a[i]), math.sin(a[i])) for i in range(3)]

    # Labels of the angles drawn in a distance from points
    a_label = text("A", (xy[0][0]*1.07, xy[0][1]*1.07))
    b_label = text("B", (xy[1][0]*1.07, xy[1][1]*1.07))
    c_label = text("C", (xy[2][0]*1.07, xy[2][1]*1.07))
    labels = a_label + b_label + c_label

    C = circle((0, 0), 1, aspect_ratio=1)

    # Triangle
    triangle = line([xy[0], xy[1], xy[2], xy[0]], rgbcolor="black")
    triangle_points = point(xy, pointsize=30)

    # Side lengths (bc, ca, ab) corresponding to triangle vertices (a, b, c)
    ad = [distance(xy[1], xy[2]), distance(xy[2], xy[0]), distance(xy[0], xy[1])]

    # Midpoints of edges (bc, ca, ab)
    a_middle = [
        ((xy[1][0] + xy[2][0])/2.0, (xy[1][1] + xy[2][1])/2.0),
        ((xy[2][0] + xy[0][0])/2.0, (xy[2][1] + xy[0][1])/2.0),
        ((xy[0][0] + xy[1][0])/2.0, (xy[0][1] + xy[1][1])/2.0)
    ]

    # Incircle
    perimeter = float(ad[0] + ad[1] + ad[2])
    incircle_center = (
        (ad[0]*xy[0][0] + ad[1]*xy[1][0] + ad[2]*xy[2][0]) / perimeter,
        (ad[0]*xy[0][1] + ad[1]*xy[1][1] + ad[2]*xy[2][1]) / perimeter
    )

    if show_incircle:
        s = perimeter/2.0
        incircle_r = math.sqrt((s - ad[0]) * (s - ad[1]) * (s - ad[2]) / s)
        incircle_graph = circle(incircle_center, incircle_r) + point(incircle_center)
    else:
        incircle_graph = Graphics()

    # Angle Bisectors
    if show_ab:
        a_ab = line([xy[0], half(a, 1, 2, 0)], rgbcolor="blue", alpha=0.6)
        b_ab = line([xy[1], half(a, 2, 0, 1)], rgbcolor="blue", alpha=0.6)
        c_ab = line([xy[2], half(a, 0, 1, 2)], rgbcolor="blue", alpha=0.6)
        ab_point = point(incircle_center, rgbcolor="blue", pointsize=28)
        ab_graph = a_ab + b_ab + c_ab + ab_point
    else:
        ab_graph = Graphics()

    # Medians
    if show_median:
        a_median = line([xy[0], a_middle[0]], rgbcolor="green", alpha=0.6)
        b_median = line([xy[1], a_middle[1]], rgbcolor="green", alpha=0.6)
        c_median = line([xy[2], a_middle[2]], rgbcolor="green", alpha=0.6)
        median_point = point(
            (
                (xy[0][0]+xy[1][0]+xy[2][0])/3.0,
                (xy[0][1]+xy[1][1]+xy[2][1])/3.0
            ), rgbcolor="green", pointsize=28)
        median_graph = a_median + b_median + c_median + median_point
    else:
        median_graph = Graphics()

    # Perpendicular Bisectors
    if show_pb:
        a_pb = line_to_points(a_middle[0], half(a, 1, 2, 0), rgbcolor="red", alpha=0.6)
        b_pb = line_to_points(a_middle[1], half(a, 2, 0, 1), rgbcolor="red", alpha=0.6)
        c_pb = line_to_points(a_middle[2], half(a, 0, 1, 2), rgbcolor="red", alpha=0.6)
        pb_point = point((0, 0), rgbcolor="red", pointsize=28)
        pb_graph = a_pb + b_pb + c_pb + pb_point
    else:
        pb_graph = Graphics()

    # Altitudes
    if show_alt:
        xA, xB, xC = xy[0][0], xy[1][0], xy[2][0]
        yA, yB, yC = xy[0][1], xy[1][1], xy[2][1]
        a_alt = plot(((xC-xB)*x+(xB-xC)*xA)/(yB-yC)+yA, (x,-3,3), rgbcolor="brown", alpha=0.6)
        b_alt = plot(((xA-xC)*x+(xC-xA)*xB)/(yC-yA)+yB, (x,-3,3), rgbcolor="brown", alpha=0.6)
        c_alt = plot(((xB-xA)*x+(xA-xB)*xC)/(yA-yB)+yC, (x,-3,3), rgbcolor="brown", alpha=0.6)
        alt_lx = (xA*xB*(yA-yB)+xB*xC*(yB-yC)+xC*xA*(yC-yA)-(yA-yB)*(yB-yC)*(yC-yA))/(xC*yB-xB*yC+xA*yC-xC*yA+xB*yA-xA*yB)
        alt_ly = (yA*yB*(xA-xB)+yB*yC*(xB-xC)+yC*yA*(xC-xA)-(xA-xB)*(xB-xC)*(xC-xA))/(yC*xB-yB*xC+yA*xC-yC*xA+yB*xA-yA*xB)
        alt_intersection = point((alt_lx, alt_ly), rgbcolor="brown", pointsize=28)
        alt_graph = a_alt + b_alt + c_alt + alt_intersection
    else:
        alt_graph = Graphics()

    # Euler's Line
    if show_euler:
        euler_graph = line_to_points(
            (0, 0),
            (
                (xy[0][0]+xy[1][0]+xy[2][0])/3.0,
                (xy[0][1]+xy[1][1]+xy[2][1])/3.0
            ),
            rgbcolor="purple",
            thickness=2,
            alpha=0.7
        )
    else:
        euler_graph = Graphics()

    show(
        C + triangle + triangle_points + labels + ab_graph + median_graph +
        pb_graph + alt_graph + incircle_graph + euler_graph,
        figsize=[5,5], xmin=-1, xmax=1, ymin=-1, ymax=1
    )


@library_interact
def coin(n = slider(2,10000, 100, default=1000, label="Number of Tosses"), interval = range_slider(0, 1, default=(0.45, 0.55), label="Plotting range (y)")):
    """
    This interact demo simulates repeated tosses of a coin,
    based on work by Lauri Ruotsalainen, 2010.

    The points give the cumulative percentage of tosses which
    are heads in a given run of the simulation, so that the
    point `(x,y)` gives the percentage of the first `x` tosses
    that were heads; this proportion should approach .5, of
    course, if we are simulating a fair coin.

    INPUT:

      - ``n`` -- number of tosses
      - ``interval`` -- plot range along
        vertical axis

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.statistics.coin()
        <html>...</html>
    """
    from random import random
    c = []
    k = 0.0
    for i in range(1, n + 1):
        k += random()
        c.append((i, k/i))
    show(point(c[1:], gridlines=[None, [0.5]], pointsize=1), ymin=interval[0], ymax=interval[1])


@library_interact
def bisection_method(
    title = text_control('<h2>Bisection method</h2>'),
    f = input_box("x^2-2", label='f(x)'),
    interval = range_slider(-5,5,default=(0, 4), label="range"),
    d = slider(1, 8, 1, 3, label="$10^{-d}$ precision"),
    maxn = slider(0,50,1,10, label="max iterations")):
    """
    Interact explaining the bisection method, based on similar interact
    explaining secant method and Wiliam Stein's example from wiki.

    INPUT:

      - ``f`` -- function
      - ``interval`` -- range slider for the search interval
      - ``d`` -- slider for the precision (`10^{-d}`)
      - ``maxn`` -- max number of iterations

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.secant_method()
        <html>...</html>
    """
    def _bisection_method(f, a, b, maxn, eps):
        intervals = [(a,b)]
        round = 1
        two = float(2)
        while True:
            c = (b+a)/two
            if abs(f(c)) < h or round >= maxn:
                break
            fa = f(a); fb = f(b); fc = f(c)
            if abs(fc) < eps:
                return c, intervals
            if fa*fc < 0:
                a, b = a, c
            elif fc*fb < 0:
                a, b = c, b
            else:
                raise ValueError("f must have a sign change in the interval (%s,%s)"%(a,b))
            intervals.append((a,b))
            round += 1
        return c, intervals

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    a, b = interval
    h = 10**(-d)
    try:
        c, intervals = _bisection_method(f, float(a), float(b), maxn, h)
    except ValueError:
        print "f must have opposite sign at the endpoints of the interval"
        show(plot(f, a, b, color='red'), xmin=a, xmax=b)
    else:
        html(r"$\text{Precision }h = 10^{-d}=10^{-%s}=%.5f$"%(d, float(h)))
        html(r"${c = }%s$"%latex(c))
        html(r"${f(c) = }%s"%latex(f(c)))
        html(r"$%s \text{ iterations}"%len(intervals))
        P = plot(f, a, b, color='red')
        k = (P.ymax() - P.ymin())/ (1.5*len(intervals))
        L = sum(line([(c,k*i), (d,k*i)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(c,k*i-k/4), (c,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(d,k*i-k/4), (d,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        show(P + L, xmin=a, xmax=b)

@library_interact
def secant_method(
    title = text_control('<h2>Secant method for numerical root finding</h2>'),
    f = input_box("x^2-2", label='f(x)'),
    interval = range_slider(-5,5,default=(0, 4), label="range"),
    d = slider(1, 16, 1, 3, label="10^-d precision"),
    maxn = slider(0,15,1,10, label="max iterations")):
    """
    Interact explaining the secant method, based on work by
    Lauri Ruotsalainen, 2010.
    Originally this is based on work by William Stein.

    INPUT:

      - ``f`` -- function
      - ``interval`` -- range slider for the search interval
      - ``d`` -- slider for the precision (10^-d)
      - ``maxn`` -- max number of iterations

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.secant_method()
        <html>...</html>
    """
    def _secant_method(f, a, b, maxn, h):
        intervals = [(a,b)]
        round = 1
        while True:
            c = b-(b-a)*f(b)/(f(b)-f(a))
            if abs(f(c)) < h or round >= maxn:
                break
            a, b = b, c
            intervals.append((a,b))
            round += 1
        return c, intervals

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    a, b = interval
    h = 10**(-d)
    if float(f(a)*f(b)) > 0:
        print "f must have opposite sign at the endpoints of the interval"
        show(plot(f, a, b, color='red'), xmin=a, xmax=b)
    else:
        c, intervals = _secant_method(f, float(a), float(b), maxn, h)
        html(r"$\text{Precision }h = 10^{-d}=10^{-%s}=%.5f$"%(d, float(h)))
        html(r"${c = }%s$"%latex(c))
        html(r"${f(c) = }%s"%latex(f(c)))
        html(r"$%s \text{ iterations}"%len(intervals))
        P = plot(f, a, b, color='red')
        k = (P.ymax() - P.ymin())/ (1.5*len(intervals))
        L = sum(line([(c,k*i), (d,k*i)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(c,k*i-k/4), (c,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        L += sum(line([(d,k*i-k/4), (d,k*i+k/4)]) for i, (c,d) in enumerate(intervals) )
        S = sum(line([(c,f(c)), (d,f(d)), (d-(d-c)*f(d)/(f(d)-f(c)), 0)], color="green") for  (c,d) in intervals)
        show(P + L + S, xmin=a, xmax=b)

@library_interact
def newton_method(
    title = text_control('<h2>Newton method</h2>'),
    f = input_box("x^2 - 2"),
    c = slider(-10,10, default=6, label='Start ($x$)'),
    d = slider(1, 16, 1, 3, label="$10^{-d}$ precision"),
    maxn = slider(0, 15, 1, 10, label="max iterations"),
    interval = range_slider(-10,10, default = (0,6), label="Interval"),
    list_steps = checkbox(default=False, label="List steps")):
    """
    Interact explaining the newton method, based on work by
    Lauri Ruotsalainen, 2010.
    Originally this is based on work by William Stein.

    INPUT:

      - ``f`` -- function
      - ``c`` -- starting position (`x`)
      - ``d`` -- slider for the precision (`10^{-d}`)
      - ``maxn`` -- max number of iterations
      - ``interval`` -- range slider for the search interval
      - ``list_steps`` -- checkbox, if true shows the steps numerically

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.newton_method()
        <html>...</html>
    """
    def _newton_method(f, c, maxn, h):
        midpoints = [c]
        round = 1
        while True:
            c = c-f(c)/f.derivative(x)(x=c)
            midpoints.append(c)
            if f(c-h)*f(c+h) < 0 or round == maxn:
                break
            round += 1
        return c, midpoints

    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    a, b = interval
    h = 10**(-d)
    c, midpoints = _newton_method(f, float(c), maxn, h/2.0)
    html(r"$\text{Precision } 2h = %s$"%latex(float(h)))
    html(r"${c = }%s$"%c)
    html(r"${f(c) = }%s"%latex(f(c)))
    html(r"$%s \text{ iterations}"%len(midpoints))
    if list_steps:
        s = [["$n$","$x_n$","$f(x_n)$", "$f(x_n-h)\,f(x_n+h)$"]]
        for i, c in enumerate(midpoints):
            s.append([i+1, c, f(c), (c-h)*f(c+h)])
        html.table(s,header_row=True)
    else:
        P = plot(f, x, interval, color="blue")
        L = sum(line([(c, 0), (c, f(c))], color="green") for c in midpoints[:-1])
        for i in range(len(midpoints) - 1):
            L += line([(midpoints[i], f(midpoints[i])), (midpoints[i+1], 0)], color="red")
        show(P + L, xmin=interval[0], xmax=interval[1], ymin=P.ymin(), ymax=P.ymax())

@library_interact
def trapezoid_integration(
    title = text_control('<h2>Trapezoid integration</h2>'),
    f = input_box(default = "x^2-5*x + 10", label='$f(x)=$'),
    n = slider(1,100,1,5, label='# divisions'),
    interval_input = selector(['from slider','from keyboard'], label='Integration interval', buttons=True),
    interval_s = range_slider(-10,10,default=(0,8), label="slider: "),
    interval_g = input_grid(1,2,default=[[0,8]], label="keyboard: "),
    output_form = selector(['traditional','table','none'], label='Computations form', buttons=True)
    ):
    """
    Interact explaining the trapezoid method for definite integrals, based on work by
    Lauri Ruotsalainen, 2010 (based on the application "Numerical integrals with various rules"
    by Marshall Hampton and Nick Alexander)

    INPUT:

      - ``f`` -- function of variable x to integrate
      - ``n`` -- number of divisions
      - ``interval_input`` -- swithes the input for interval between slider and keyboard
      - ``interval_s`` -- slider for interval to integrate
      - ``interval_g`` -- input grid for interval to integrate
      - ``output_form`` -- the computation is formatted in a traditional form, in a table or missing

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.trapezoid_integration()
        <html>...</html>
    """
    xs = []
    ys = []
    if interval_input == 'from slider':
        interval = interval_s
    else:
        interval = interval_g[0]
    h = float(interval[1]-interval[0])/n
    x = SR.var('x')
    f = symbolic_expression(f).function(x)

    trapezoids = Graphics()

    for i in range(n):
        xi = interval[0] + i*h
        yi = f(xi)
        trapezoids += line([[xi, 0], [xi, yi], [xi + h, f(xi + h)],[xi + h, 0],[xi, 0]], rgbcolor = (1,0,0))
        xs.append(xi)
        ys.append(yi)
    xs.append(xi + h)
    ys.append(f(xi + h))

    html(r'Function $f(x)=%s$'%latex(f(x)))
    show(plot(f, interval[0], interval[1]) + trapezoids, xmin = interval[0], xmax = interval[1])

    numeric_value = integral_numerical(f, interval[0], interval[1])[0]
    approx = h *(ys[0]/2 + sum([ys[i] for i in range(1,n)]) + ys[n]/2)

    html(r'Integral value to seven decimal places is: $\displaystyle\int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x} = %.6f$'%(
            interval[0], interval[1], N(numeric_value, digits=7))
         )

    if output_form == 'traditional':
        sum_formula_html = r"\frac {d}{2} \cdot \left[f(x_0) + %s + f(x_{%s})\right]" % (
            ' + '.join([ "2 f(x_{%s})"%i for i in range(1,n)]),
            n
            )
        sum_placement_html = r"\frac{%.2f}{2} \cdot \left[f(%.2f) + %s + f(%.2f)\right]" % (
            h,
            N(xs[0], digits=5),
            ' + '.join([ "2 f(%.2f)" %N(i, digits=5) for i in xs[1:-1]]),
            N(xs[n], digits=5)
            )
        sum_values_html = r"\frac{%.2f}{2} \cdot \left[%.2f + %s + %.2f\right]" % (
            h,
            N(ys[0], digits=5),
            ' + '.join([ "2\cdot %.2f" % N(i, digits=5) for i in ys[1:-1]]),
            N(ys[n], digits=5)
            )

        html(r'''
            <div class="math">
            \begin{align*}
            \int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x}
                & \approx %s \\
                & = %s \\
                & = %s \\
                & = %s
            \end{align*}
            </div>
        ''' % (
                interval[0], interval[1],
                sum_formula_html, sum_placement_html, sum_values_html,
                N(approx, digits=7)
        ))
    elif output_form == 'table':
        s = [['$i$','$x_i$','$f(x_i)$','$m$','$m\cdot f(x_i)$']]
        for i in range(0,n+1):
            if i==0 or i==n:
                j = 1
            else:
                j = 2
            s.append([i, xs[i], ys[i],j,N(j*ys[i])])
        html.table(s,header_row=True)

@library_interact
def simpson_integration(
    title = text_control('<h2>Simpson integration</h2>'),
    f = input_box(default = 'x*sin(x)+x+1', label='$f(x)=$'),
    n = slider(2,100,2,6, label='# divisions'),
    interval_input = selector(['from slider','from keyboard'], label='Integration interval', buttons=True),
    interval_s = range_slider(-10,10,default=(0,10), label="slider: "),
    interval_g = input_grid(1,2,default=[[0,10]], label="keyboard: "),
    output_form = selector(['traditional','table','none'], label='Computations form', buttons=True)):
    """
    Interact explaining the simpson method for definite integrals, based on work by
    Lauri Ruotsalainen, 2010 (based on the application "Numerical integrals with various rules"
    by Marshall Hampton and Nick Alexander)

    INPUT:

      - ``f`` -- function of variable x to integrate
      - ``n`` -- number of divisions (mult. of 2)
      - ``interval_input`` -- swithes the input for interval between slider and keyboard
      - ``interval_s`` -- slider for interval to integrate
      - ``interval_g`` -- input grid for interval to integrate
      - ``output_form`` -- the computation is formatted in a traditional form, in a table or missing

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.simpson_integration()
        <html>...</html>
    """
    x = SR.var('x')
    f = symbolic_expression(f).function(x)
    if interval_input == 'from slider':
        interval = interval_s
    else:
        interval = interval_g[0]
    def parabola(a, b, c):
        from sage.all import solve
        A, B, C = SR.var("A, B, C")
        K = solve([A*a[0]**2+B*a[0]+C==a[1], A*b[0]**2+B*b[0]+C==b[1], A*c[0]**2+B*c[0]+C==c[1]], [A, B, C], solution_dict=True)[0]
        f = K[A]*x**2+K[B]*x+K[C]
        return f
    xs = []; ys = []
    dx = float(interval[1]-interval[0])/n

    for i in range(n+1):
        xs.append(interval[0] + i*dx)
        ys.append(f(x=xs[-1]))

    parabolas = Graphics()
    lines = Graphics()

    for i in range(0, n-1, 2):
        p = parabola((xs[i],ys[i]),(xs[i+1],ys[i+1]),(xs[i+2],ys[i+2]))
        parabolas += plot(p(x=x), (x, xs[i], xs[i+2]), color="red")
        lines += line([(xs[i],ys[i]), (xs[i],0), (xs[i+2],0)],color="red")
        lines += line([(xs[i+1],ys[i+1]), (xs[i+1],0)], linestyle="-.", color="red")

    lines += line([(xs[-1],ys[-1]), (xs[-1],0)], color="red")

    html(r'Function $f(x)=%s$'%latex(f(x)))

    show(plot(f(x),x,interval[0],interval[1]) + parabolas + lines, xmin = interval[0], xmax = interval[1])

    numeric_value = integral_numerical(f,interval[0],interval[1])[0]
    approx = dx/3 *(ys[0] + sum([4*ys[i] for i in range(1,n,2)]) + sum([2*ys[i] for i in range(2,n,2)]) + ys[n])

    html(r'Integral value to seven decimal places is: $\displaystyle\int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x} = %.6f$'%
         (interval[0],interval[1],
         N(numeric_value,digits=7)))

    if output_form == 'traditional':
        sum_formula_html = r"\frac{d}{3} \cdot \left[ f(x_0) + %s + f(x_{%s})\right]" % (
            ' + '.join([ r"%s \cdot f(x_{%s})" %(i%2*(-2)+4, i+1) for i in range(0,n-1)]),
            n
            )

        sum_placement_html = r"\frac{%.2f}{3} \cdot \left[ f(%.2f) +  %s + f(%.2f)\right]" % (
            dx,
            N(xs[0],digits=5),
            ' + '.join([ r"%s \cdot f(%.2f)" %(i%2*(-2)+4, N(xk, digits=5)) for i, xk in enumerate(xs[1:-1])]),
            N(xs[n],digits=5)
            )

        sum_values_html = r"\frac{%.2f}{3} \cdot \left[ %s %s %s\right]" %(
            dx,
            "%.2f + "%N(ys[0],digits=5),
            ' + '.join([ r"%s \cdot %.2f" %(i%2*(-2)+4, N(yk, digits=5)) for i, yk in enumerate(ys[1:-1])]),
            " + %.2f"%N(ys[n],digits=5)
            )

        html(r'''
        <div class="math">
        \begin{align*}
        \int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x}
            & \approx %s \\
            & = %s \\
            & = %s \\
            & = %.6f
        \end{align*}
        </div>
        ''' % (
                interval[0], interval[1],
                sum_formula_html, sum_placement_html, sum_values_html,
                N(approx,digits=7)
                ))
    elif output_form == 'table':
        s = [['$i$','$x_i$','$f(x_i)$','$m$','$m\cdot f(x_i)$']]
        for i in range(0,n+1):
            if i==0 or i==n:
                j = 1
            else:
                j = (i+1)%2*(-2)+4
            s.append([i, xs[i], ys[i],j,N(j*ys[i])])
        s.append(['','','','$\sum$','$%s$'%latex(3/dx*approx)])
        html.table(s,header_row=True)
        html(r'$\int_{%.2f}^{%.2f} {f(x) \, \mathrm{d}x}\approx\frac {%.2f}{3}\cdot %s=%s$'%
             (interval[0], interval[1],dx,latex(3/dx*approx),latex(approx)))

@library_interact
def riemann_sum(
    title = text_control('<h2>Riemann integral with random sampling</h2>'),
    f = input_box("x^2+1", label = "$f(x)=$", width=40),
    n = slider(1,30,1,5, label='# divisions'),
    hr1 = text_control('<hr>'),
    interval_input = selector(['from slider','from keyboard'], label='Integration interval', buttons=True),
    interval_s = range_slider(-5,10,default=(0,2), label="slider: "),
    interval_g = input_grid(1,2,default=[[0,2]], label="keyboard: "),
    hr2 = text_control('<hr>'),
    list_table = checkbox(default=False, label="List table"),
    auto_update = False):
    """
    Interact explaining the definition of Riemann integral

    INPUT:

    - ``f`` -- function of variable x to integrate
    - ``n`` -- number of divisions
      - ``interval_input`` -- swithes the input for interval between slider and keyboard
      - ``interval_s`` -- slider for interval to integrate
      - ``interval_g`` -- input grid for interval to integrate
    - ``list_table`` -- print table with values of the function

    EXAMPLES:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.riemann_sum()
        <html>...</html>

    AUTHORS:

    - Robert Marik (08-2010)
    """
    x = SR.var('x')
    from random import random
    if interval_input == 'from slider':
        a = interval_s[0]
        b = interval_s[1]
    else:
        a = interval_g[0][0]
        b = interval_g[0][1]
    func = symbolic_expression(f).function(x)
    division = [a]+[a+random()*(b-a) for i in range(n-1)]+[b]
    division = sorted([i for i in division])
    xs = [division[i]+random()*(division[i+1]-division[i]) for i in range(n)]
    ys = [func(x_val) for x_val in xs]
    rects = Graphics()
    for i in range(n):
        body=[[division[i],0],[division[i],ys[i]],[division[i+1],ys[i]],[division[i+1],0]]
        if ys[i].n()>0:
            color_rect='green'
        else:
            color_rect='red'
        rects = rects +polygon2d(body, rgbcolor = color_rect,alpha=0.1)\
         + point((xs[i],ys[i]), rgbcolor = (1,0,0))\
         + line(body,rgbcolor='black',zorder=-1)
    html('<small>Adjust your data and click Update button. Click repeatedly for another random values.</small>')

    show(plot(func(x),(x,a,b),zorder=5) + rects)
    delka_intervalu=[division[i+1]-division[i] for i in range(n)]
    if list_table:
        html.table([["$i$","$[x_{i-1},x_i]$","$\eta_i$","$f(\eta_i)$","$x_{i}-x_{i-1}$"]]\
        +[[i+1,[division[i],division[i+1]],xs[i],ys[i],delka_intervalu[i]] for i in range(n)],\
        header_row=True)

    html('Riemann sum: $\displaystyle\sum_{i=1}^{%s} f(\eta_i)(x_i-x_{i-1})=%s$ '%
         (latex(n),latex(sum([ys[i]*delka_intervalu[i] for i in range(n)]))))
    html('Exact value of the integral $\displaystyle\int_{%s}^{%s}%s\,\mathrm{d}x=%s$'%
         (latex(a),latex(b),latex(func(x)),latex(integral_numerical(func(x),a,b)[0])))


x = SR.var('x')
@library_interact
def function_tool(f=sin(x), g=cos(x), xrange=range_slider(-3,3,default=(0,1),label='x-range'),
      yrange='auto',
      a=1,
      action=selector(['f', 'df/dx', 'int f', 'num f', 'den f', '1/f', 'finv',
                       'f+a', 'f-a', 'f*a', 'f/a', 'f^a', 'f(x+a)', 'f(x*a)',
                       'f+g', 'f-g', 'f*g', 'f/g', 'f(g)'],
             width=15, nrows=5, label="h = "),
      do_plot = ("Draw Plots", True)):
    """
    `Function Plotting Tool <http://wiki.sagemath.org/interact/calculus#Functiontool>`_
    (by William Stein (?))

    INPUT:

      - ``f`` -- function f(x)
      - ``g`` -- function g(x)
      - ``xrange`` -- range for plotting (x)
      - ``yrange`` -- range for plotting ('auto' is default, otherwise a tuple)
      - ``a`` -- factor ``a``
      - ``action`` -- select given operation on or combination of functions
      - ``do_plot`` -- if true, a plot is drawn

    EXAMPLE:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.calculus.function_tool()
        <html>...</html>
    """
    x = SR.var('x')
    try:
        f = SR(f); g = SR(g); a = SR(a)
    except TypeError as msg:
        print msg[-200:]
        print "Unable to make sense of f,g, or a as symbolic expressions in single variable x."
        return
    if not (isinstance(xrange, tuple) and len(xrange) == 2):
          xrange = (0,1)
    h = 0; lbl = ''
    if action == 'f':
        h = f
        lbl = 'f'
    elif action == 'df/dx':
        h = f.derivative(x)
        lbl = r'\frac{\mathrm{d}f}{\mathrm{d}x}'
    elif action == 'int f':
        h = f.integrate(x)
        lbl = r'\int f \,\mathrm{d}x'
    elif action == 'num f':
        h = f.numerator()
        lbl = r'\text{numer(f)}'
    elif action == 'den f':
        h = f.denominator()
        lbl = r'\text{denom(f)}'
    elif action == '1/f':
        h = 1/f
        lbl = r'\frac{1}{f}'
    elif action == 'finv':
        h = solve(f == var('y'), x)[0].rhs()
        lbl = 'f^{-1}(y)'
    elif action == 'f+a':
        h = f+a
        lbl = 'f + a'
    elif action == 'f-a':
        h = f-a
        lbl = 'f - a'
    elif action == 'f*a':
        h = f*a
        lbl = r'f \times a'
    elif action == 'f/a':
        h = f/a
        lbl = r'\frac{f}{a}'
    elif action == 'f^a':
        h = f**a
        lbl = 'f^a'
    elif action == 'f^a':
        h = f**a
        lbl = 'f^a'
    elif action == 'f(x+a)':
        h = f.subs(x=x+a)
        lbl = 'f(x+a)'
    elif action == 'f(x*a)':
        h = f.subs(x=x*a)
        lbl = 'f(xa)'
    elif action == 'f+g':
        h = f+g
        lbl = 'f + g'
    elif action == 'f-g':
        h = f-g
        lbl = 'f - g'
    elif action == 'f*g':
        h = f*g
        lbl = r'f \times g'
    elif action == 'f/g':
        h = f/g
        lbl = r'\frac{f}{g}'
    elif action == 'f(g)':
        h = f(g)
        lbl = 'f(g)'
    html('<center><font color="red">$f = %s$</font></center>'%latex(f))
    html('<center><font color="green">$g = %s$</font></center>'%latex(g))
    html('<center><font color="blue"><b>$h = %s = %s$</b></font></center>'%(lbl, latex(h)))
    if do_plot:
        P = plot(f, xrange, color='red', thickness=2) +  \
            plot(g, xrange, color='green', thickness=2) + \
            plot(h, xrange, color='blue', thickness=2)
        if yrange == 'auto':
            show(P, xmin=xrange[0], xmax=xrange[1])
        else:
            yrange = sage_eval(yrange)
            show(P, xmin=xrange[0], xmax=xrange[1], ymin=yrange[0], ymax=yrange[1])

@library_interact
def julia(expo = slider(-10,10,0.1,2),
    c_real = slider(-2,2,0.01,0.5, label='real part const.'),
    c_imag = slider(-2,2,0.01,0.5, label='imag part const.'),
    iterations=slider(1,100,1,20, label='# iterations'),
    zoom_x = range_slider(-2,2,0.01,(-1.5,1.5), label='Zoom X'),
    zoom_y = range_slider(-2,2,0.01,(-1.5,1.5), label='Zoom Y'),
    plot_points = slider(20,400,20, default=150, label='plot points'),
    dpi = slider(20, 200, 10, default=80, label='dpi')):
    """
    Julia Fractal, based on
    `Julia by Harald Schilly <http://wiki.sagemath.org/interact/fractal#Julia>`_.

    INPUT:

        - ``exponent`` -- exponent ``e`` in `z^e+c`
        - ``c_real`` -- real part of the constant ``c``
        - ``c_imag`` -- imaginary part of the constant ``c``
        - ``iterations`` -- number of iterations
        - ``zoom_x`` -- range slider for zoom in x direction
        - ``zoom_y`` -- range slider for zoom in y direction
        - ``plot_points`` -- number of points to plot
        - ``dpi`` -- dots-per-inch parameter for the plot

    EXAMPLE:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.fractals.julia()
        <html>...</html>
    """
    z = SR.var('z')
    I = CDF.gen()
    f = symbolic_expression(z**expo + c_real + c_imag*I).function(z)
    ff_j = fast_callable(f, vars=[z], domain=CDF)

    from sage.interacts.library_cython import julia

    html('<h2>Julia Fractal</h2>')
    html(r'Recursive Formula: $z \leftarrow z^{%.2f} + (%.2f+%.2f*\mathbb{I})$' % (expo, c_real, c_imag))
    complex_plot(lambda z: julia(ff_j, z, iterations), zoom_x, zoom_y, plot_points=plot_points, dpi=dpi).show(frame=True, aspect_ratio=1)

@library_interact
def mandelbrot(expo = slider(-10,10,0.1,2),
    iterations=slider(1,100,1,20, label='# iterations'),
    zoom_x = range_slider(-2,2,0.01,(-2,1), label='Zoom X'),
    zoom_y = range_slider(-2,2,0.01,(-1.5,1.5), label='Zoom Y'),
    plot_points = slider(20,400,20, default=150, label='plot points'),
    dpi = slider(20, 200, 10, default=80, label='dpi')):
    """
    Mandelbrot Fractal, based on
    `Mandelbrot by Harald Schilly <http://wiki.sagemath.org/interact/fractal#Mandelbrot>`_.

    INPUT:

        - ``exponent`` -- exponent ``e`` in `z^e+c`
        - ``iterations`` -- number of iterations
        - ``zoom_x`` -- range slider for zoom in x direction
        - ``zoom_y`` -- range slider for zoom in y direction
        - ``plot_points`` -- number of points to plot
        - ``dpi`` -- dots-per-inch parameter for the plot

    EXAMPLE:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.fractals.mandelbrot()
        <html>...</html>
    """
    x, z, c = SR.var('x, z, c')
    f = symbolic_expression(z**expo + c).function(z, c)
    ff_m = fast_callable(f, vars=[z,c], domain=CDF)

    from sage.interacts.library_cython import mandel

    html('<h2>Mandelbrot Fractal</h2>')
    html(r'Recursive Formula: $z \leftarrow z^{%.2f} + c$ for $c \in \mathbb{C}$' % expo)
    complex_plot(lambda z: mandel(ff_m, z, iterations), zoom_x, zoom_y, plot_points=plot_points, dpi=dpi).show(frame=True, aspect_ratio=1)


@library_interact
def cellular_automaton(
    N=slider(1,500,1,label='Number of iterations',default=100),
    rule_number=slider(0, 255, 1, default=110, label='Rule number'),
    size = slider(1, 11, step_size=1, default=6, label='size of graphic')):
    """
    Yields a matrix showing the evolution of a
    `Wolfram's cellular automaton <http://mathworld.wolfram.com/CellularAutomaton.html>`_.

    `Based on work by Pablo Angulo <http://wiki.sagemath.org/interact/misc#CellularAutomata>`_.

    INPUT:

        - ``N`` -- iterations
        - ``rule_number`` -- rule number (0 to 255)
        - ``size`` -- size of the shown picture

    EXAMPLE:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: interacts.fractals.cellular_automaton()
        <html>...</html>
    """
    from sage.all import Integer
    if not 0 <= rule_number <= 255:
        raise ValueError('Invalid rule number')
    binary_digits = Integer(rule_number).digits(base=2)
    rule = binary_digits + [0]*(8-len(binary_digits))

    html('<h2>Cellular Automaton</h2>'+
         '<div style="white-space: normal;">"A cellular automaton is a collection of "colored" cells \
         on a grid of specified shape that evolves through a number of \
         discrete time steps according to a set of rules based on the \
         states of neighboring cells." &mdash; \
         <a target="_blank" href="http://mathworld.wolfram.com/CellularAutomaton.html">Mathworld,\
         Cellular Automaton</a></div>\
         <div>Rule %s expands to %s</div>' % (rule_number, ''.join(map(str,rule)))
        )

    from sage.interacts.library_cython import cellular
    M = cellular(rule, N)
    plot_M = matrix_plot(M, cmap='binary')
    plot_M.show(figsize=[size,size])


@library_interact
def polar_prime_spiral(
    interval = range_slider(1, 4000, 10, default=(1, 1000), label="range"),
    show_factors = True,
    highlight_primes = True,
    show_curves = True,
    n = slider(1,200, 1, default=89, label="number $n$"),
    dpi = slider(10,300, 10, default=100, label="dpi")):
    """
    Polar Prime Spiral interact, based on work by David Runde.

    For more information about the factors in the spiral,
    `visit John Williamson's website <http://www.dcs.gla.ac.uk/~jhw/spirals/index.html>`_.

    INPUT:

        - ``interval`` -- range slider to specify start and end
        - ``show_factors`` -- if true, show factors
        - ``highlight_primes`` -- if true, prime numbers are highlighted
        - ``show_curves`` -- if true, curves are plotted
        - ``n`` -- number `n`
        - ``dpi`` -- dots per inch resolution for plotting

    EXAMPLE:

    Invoked in the notebook, the following command will produce
    the fully formatted interactive mathlet.  In the command line,
    it will simply return the underlying HTML and Sage code which
    creates the mathlet::

        sage: sage.interacts.algebra.polar_prime_spiral()
        <html>...</html>

    """

    html('<h2>Polar Prime Spiral</h2> \
          <div style="white-space: normal;">\
          For more information about the factors in the spiral, visit \
          <a href="http://www.dcs.gla.ac.uk/~jhw/spirals/index.html" target="_blank">\
          Number Spirals by John Williamson</a>.</div>'
        )

    start, end = interval
    from sage.ext.fast_eval import fast_float
    from math import floor, ceil
    from sage.plot.colors import hue

    if start < 1 or end <= start:
        print "invalid start or end value"
        return
    if n > end:
        print "WARNING: n is greater than end value"
        return
    if n < start:
        print "n < start value"
        return
    nn = SR.var('nn')
    f1 = fast_float(sqrt(nn)*cos(2*pi*sqrt(nn)), 'nn')
    f2 = fast_float(sqrt(nn)*sin(2*pi*sqrt(nn)), 'nn')
    f = lambda x: (f1(x), f2(x))

    list =[]
    list2=[]
    if show_factors == False:
        for i in srange(start, end, include_endpoint = True):
            if Integer(i).is_pseudoprime(): list.append(f(i-start+1)) #Primes list
            else: list2.append(f(i-start+1)) #Composites list
        P = points(list)
        R = points(list2, alpha = .1) #Faded Composites
    else:
        for i in srange(start, end, include_endpoint = True):
            list.append(disk((f(i-start+1)),0.05*pow(2,len(factor(i))-1), (0,2*pi))) #resizes each of the dots depending of the number of factors of each number
            if Integer(i).is_pseudoprime() and highlight_primes: list2.append(f(i-start+1))
        P = plot(list)
        p_size = 5 #the orange dot size of the prime markers
        if not highlight_primes: list2 = [(f(n-start+1))]
        R = points(list2, hue = .1, pointsize = p_size)

    if n > 0:
        html('$n = %s$' % factor(n))

        p = 1
    #The X which marks the given n
        W1 = disk((f(n-start+1)), p, (pi/6, 2*pi/6))
        W2 = disk((f(n-start+1)), p, (4*pi/6, 5*pi/6))
        W3 = disk((f(n-start+1)), p, (7*pi/6, 8*pi/6))
        W4 = disk((f(n-start+1)), p, (10*pi/6, 11*pi/6))
        Q = plot(W1 + W2 + W3 + W4, alpha = .1)

        n = n - start +1        #offsets the n for different start values to ensure accurate plotting
        if show_curves:
            begin_curve = 0
            t = SR.var('t')
            a=1.0
            b=0.0
            if n > (floor(sqrt(n)))**2 and n <= (floor(sqrt(n)))**2 + floor(sqrt(n)):
                c = -((floor(sqrt(n)))**2 - n)
                c2= -((floor(sqrt(n)))**2 + floor(sqrt(n)) - n)
            else:
                c = -((ceil(sqrt(n)))**2 - n)
                c2= -((floor(sqrt(n)))**2 + floor(sqrt(n)) - n)
            html('Pink Curve:  $n^2 + %s$' % c)
            html('Green Curve: $n^2 + n + %s$' % c2)
            m = SR.var('m')
            g = symbolic_expression(a*m**2+b*m+c).function(m)
            r = symbolic_expression(sqrt(g(m))).function(m)
            theta = symbolic_expression(r(m)- m*sqrt(a)).function(m)
            S1 = parametric_plot(((r(t))*cos(2*pi*(theta(t))),(r(t))*sin(2*pi*(theta(t)))),
                 (begin_curve, ceil(sqrt(end-start))), color=hue(0.8), thickness = .3) #Pink Line

            b = 1
            c = c2;
            g = symbolic_expression(a*m**2+b*m+c).function(m)
            r = symbolic_expression(sqrt(g(m))).function(m)
            theta = symbolic_expression(r(m)- m*sqrt(a)).function(m)
            S2 = parametric_plot(((r(t))*cos(2*pi*(theta(t))),(r(t))*sin(2*pi*(theta(t)))),
                 (begin_curve, ceil(sqrt(end-start))), color=hue(0.6), thickness = .3) #Green Line

            show(R+P+S1+S2+Q, aspect_ratio = 1, axes = False, dpi = dpi)
        else: show(R+P+Q, aspect_ratio = 1, axes = False, dpi = dpi)
    else: show(R+P, aspect_ratio = 1, axes = False, dpi = dpi)



