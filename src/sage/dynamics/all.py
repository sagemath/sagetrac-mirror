"""
Some code about flat surfaces and interval exchanges has been removed
from Sage. The package ``surface_dynamics`` contains all that code
and much more. For more information, see

    https://www.labri.fr/perso/vdelecro/surface-dynamics/latest/

You can install the ``surface_dynamics`` package via:

    sage -pip install surface_dynamics

If you do not have write access to the Sage installation you can
alternatively do:

    sage -pip install surface_dynamics --user
"""
from sage.misc.lazy_import import lazy_import

from sage.dynamics.arithmetic_dynamics.all import *
lazy_import('sage.dynamics.complex_dynamics.mandel_julia',
            ["mandelbrot_plot", "external_ray", "kneading_sequence", "julia_plot"])
from sage.dynamics.cellular_automata.all import *

# Discrete dynamical systems
lazy_import('sage.dynamics.finite_dynamical_system',
            ['DiscreteDynamicalSystem'])

lazy_import('sage.dynamics', 'finite_dynamical_system_catalog', 'finite_dynamical_systems')
