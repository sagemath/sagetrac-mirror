"""
Differential Geometry of Parametrized Surfaces

AUTHORS:

- Mikhail Malakhaltsev (2010-09-25): initial version
- Joris Vankerschaver  (2010-10-25): implementation, doctests

"""
# ****************************************************************************
#       Copyright (C) 2010  Mikhail Malakhaltsev <mikarm@gmail.com>
#       Copyright (C) 2010  Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import

lazy_import('sage.manifolds.differentiable.examples.parametrized_surface3d',
            'ParametrizedSurface3D',
            deprecation=32228)
