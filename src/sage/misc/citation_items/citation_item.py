r"""
An abstract class for citable items.

All items in ``sage/misc/citation_items`` should inherit from
:class:`CitationItem`.  See ``sage/misc/citation_items/all`` for a
list of all of them.

AUTHORS:

- Martin Westerholt-Raum (martin@raum-brothers.eu)

..todo:

At the current stage, functionality is very limited.

- Include bibliographic information.
"""

###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

import re

class CitationItem( SageObject, UniqueRepresentation ):
    r"""
    EXAMPLES::

        sage: from sage.misc.citation_items.axiom import *
        sage: axiom = Axiom_CitationItem()
        sage: axiom is Axiom_CitationItem()
        True
    """

    def re(self):
        r"""
        Regular expressions matching modules which can be cited using
        this item.

        OUTPUT:

        A list of compiled regular expressions.

        EXAMPLES::

            sage: from sage.misc.citation_items.axiom import *
            sage: [r.pattern for r in Axiom_CitationItem().re()]
            [...'^sage.interfaces.axiom'...]
        """
        return [re.compile(s) for s in self._re]

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.misc.citation_items.axiom import *
            sage: Axiom_CitationItem()
            Axiom
        """
        return self._name

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.misc.citation_items.axiom import *
            sage: latex(Axiom_CitationItem())
            \text{Axiom}
        """
        return r"\text{{{}}}".format(self._name)
