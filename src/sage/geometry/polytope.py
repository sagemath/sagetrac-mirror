"""
Polytopes

This file is a wrapper for the polymake module, which cannot be
imported by default because it is linked to an optional shared
library

"""

########################################################################
#       Copyright (C) 2012 Timo Kluck <tkluck@infty.nl>
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

class PolymakeProxy(object):
    # TODO: use some wrapper decorator to make sure that this proxy has
    # the right tab completion and docstrings
    def __getattr__(self, attr):
        try:
            import sage.libs.polymake.polymake as pm
        except:
            raise RuntimeError(
            """Polymake requires the optional polymake package. It can be installed
            by typing

                 sage -i polymake

            on the command line""")
        return getattr(pm, attr)
            

polymake = PolymakeProxy()
