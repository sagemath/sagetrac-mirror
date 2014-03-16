"""
Lattice plot functionality

AUTHORS:

- Jan Poeschko (2012-08-15): initial version

EXAMPLES::

    sage: L = special_lattice('SimpleCubic')
    sage: g = plot_lattice(L)
    
You can also add the Voronoi cell to the plot::
    
    sage: V = L.voronoi_cell()
    sage: g += V.plot()
    sage: g.show(viewer='tachyon')
"""

#*****************************************************************************
#       Copyright (C) 2012 Jan Poeschko <jan@poeschko.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage import plot

def subsets(s):
    """
    Generate all subsets of a list of elements.
    
    INPUT:
    
    - ``s`` -- a list of elements.
    
    OUTPUT:
    
    A generator iterating through all subsets of `s`.
    
    EXAMPLES::
        
        sage: from sage.lattices.plot import subsets
        sage: list(subsets([1, 2, 3]))
        [[], [1], [2], [1, 2], [3], [1, 3], [2, 3], [1, 2, 3]]
    """
    n_elems = len(s)
    n_subsets = 2**len(s)
    for i in range(0,n_subsets):
        sub = []
        for j in range(0,n_elems):
            if (i >> j & 1):
                sub.append(s[j])
        yield sub

def plot_lattice(L, **options):
    """
    Plot a lattice.
    
    The sum of each subset of basis vectors (including their negated counterparts)
    is plotted.
    
    INPUT:
    
    - ``L`` -- a `Lattice` instance;
    
    - ``options`` - options to be passed to `list_plot`.
    
    OUTPUT:
    
    Graphics containing the plot of the lattice.
    
    EXAMPLES::
    
        sage: L = Lattice([[2, 0], [0, 3]])
        sage: plot_lattice(L)
        
    Three-dimensional lattices are supported, too::
    
        sage: L = Lattice([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        sage: plot_lattice(L)
        
    Higher dimensions are not supported, though::
    
        sage: L = random_lattice(4)
        sage: plot_lattice(L)
        Traceback (most recent call last):
        ...
        ValueError: only 2-dimensional and 3-dimensional lattices can be plotted
    """
    dim = L.dimension()
    basis = L.embedded_basis()
    basis = sum(([-v, v] for v in basis), [])
    points = [sum(v, L(0)) for v in subsets(basis)]
    points = [tuple(v.list()) for v in points]
    points = set(points)
    if dim not in (2, 3):
        raise ValueError("only 2-dimensional and 3-dimensional lattices can be plotted")
    return plot.point.points(points)
