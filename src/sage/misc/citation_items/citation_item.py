r"""
An abstract class for citable items.

All items in ``sage/misc/citation_items`` should inherit from
:class:`CitationItem`.  See ``sage/misc/citation_items/all`` for a
list of all of them.

AUTHORS:

- Martin Westerholt-Raum (martin@raum-brothers.eu)
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

    def _latex_citation_label(self):
        r"""
        A LaTeX citation label which can be used in combination with
        BibTeX code given by :meth:`bibtex` to cite the present item.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.misc.citation_items.axiom import *
            sage: Axiom_CitationItem()._latex_citation_label()
            'software-axiom'
        """
        if not self._bibtex:
            return ""

        start = self._bibtex.find('{') + 1
        end = self._bibtex.find(',', start)
        return self._bibtex[start:end]

    def bibtex(self):
        r"""
        Bibtex code that can be used to cite this item in LaTeX/BibTeX.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.misc.citation_items.axiom import *
            sage: print Axiom_CitationItem().bibtex()
            @Book{software-axiom,
                Author = {Jenks, Richard D. and Sutor, Robert S.},
                Title = {{Axiom. The scientific computation system.}},
                Pages = {xxiv + 742},
                Year = {1992},
                Publisher = {Springer-Verlag, Berlin},
                }

        .. TODO:

        Use optional parsing to generate bibitems for those who prefer
        this variant.
        """
        if not self._bibtex:
            return ""

        return "\n    ".join(s.strip() for s in self._bibtex.splitlines() if s.strip() != "")

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
