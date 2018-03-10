from sage.plot.plot import plot
from sage.plot.line import line
from sage.rings.real_double import RDF
from sage.misc.prandom import random
from sage.misc.misc import srange
from sage.plot.point import point


def cobweb(a_function, start, mask=0, iterations=20, xmin=0, xmax=1):
    """
    Return a graphics object of a plot of the function and a cobweb
    trajectory starting from the value start.

    INPUT:

    - a_function -- a function of one variable

    - start -- the starting value of the iteration

    - mask -- (optional) the number of initial iterates to ignore

    - iterations -- (optional) the number of iterations to draw,
      following the masked iterations

    - xmin -- (optional) the lower end of the plotted interval

    - xmax -- (optional) the upper end of the plotted interval

    EXAMPLES::

        sage: from sage.dynamics.cobweb import cobweb
        sage: f = lambda x: 3.9*x*(1-x)
        sage: show(cobweb(f,.01,iterations=200), xmin=0, xmax=1, ymin=0)

    .. NOTE::

        This is very slow with symbolic functions.
    """
    basic_plot = plot(a_function, xmin=xmin, xmax=xmax)
    id_plot = plot(lambda x: x, xmin=xmin, xmax=xmax)
    iter_list = []
    current = start
    for i in range(mask):
        current = a_function(current)
    for i in range(iterations):
        iter_list.append([current, a_function(current)])
        current = a_function(current)
        iter_list.append([current, current])
    cobweb = line(iter_list)
    return basic_plot + id_plot + cobweb


def orbit_diagram(a_function, parameter_interval, domain=[0, 1],
                  mask=50, iterations=200, param_num=500.0):
    """
    Return a plot of the iterations of a function as a function of a
    parameter value.

    INPUT:

    - a_function -- a function of one variable

    - parameter_interval -- a two-element list of the lowest and highest
      parameters to plot.

    - domain -- (optional) a two-element list of the lowest and highest
      input values to iterate

    - mask -- (optional) the number of initial iterates to ignore

    - iterations -- (optional) the number of iterations to draw,
      following the masked iterations

    EXAMPLES::

        sage: from sage.dynamics.cobweb import orbit_diagram
        sage: f = lambda x,m: m*x*(1-x)
        sage: show(orbit_diagram(f, [3.4,4], mask=100, iterations=500),
        ....:    xmin=3.4, ymin=0)

    .. NOTE::

        This is pretty crude so far.
    """
    point_list = []
    plen = RDF(parameter_interval[1] - parameter_interval[0])
    seed = random() * (domain[1] - domain[0]) + domain[0]
    for i in srange(parameter_interval[0],
                    parameter_interval[1],
                    plen/param_num):
        for x in range(mask):
            seed = a_function(seed, i)
        for x in srange(iterations):
            seed = a_function(seed, i)
            point_list.append((i, seed))
    return point(point_list, pointsize=1, rgbcolor=(0, 0, 0))
