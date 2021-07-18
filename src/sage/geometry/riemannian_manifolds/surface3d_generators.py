r"""
Common parametrized surfaces in 3D.

AUTHORS::

- Joris Vankerschaver (2012-06-16)

"""
#*****************************************************************************
#       Copyright (C) 2010  Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_import import lazy_import

lazy_import('sage.manifolds.differentiable.examples.surface3d_generators',
            ('SurfaceGenerators', 'surfaces'),
            deprecation=32228)
