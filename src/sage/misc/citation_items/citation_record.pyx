r"""
A record of citation items.

This is a list of citation items, which allows the use to access LaTeX
cite statements and BibTeX code.

.. TODO:

Add pybtex support and replace the workaround as soon as pybtex is a
standard package.

AUTHORS:

- Martin Westerholt-Raum (martin@raum-brothers.eu)

EXAMPLES::

    sage: from sage.misc.citation_items.all import CitationRecord
    sage: from sage.misc.citation_items.ecm import ECM_CitationItem
    sage: from sage.misc.citation_items.pari import PARI_CitationItem
    sage: rec = CitationRecord()
    sage: rec.append(ECM_CitationItem())
    sage: rec
    [ECM]
    sage: rec.append(PARI_CitationItem())
    sage: rec.print_latex_citation()
    Print BibTeX code by calling method print_bibtex().
    \cite{software-ecm, software-pari}
    sage: rec.print_bibtex()
    @Misc{software-ecm,
        Author = {Charron, T. and Daminelli, N. and Granlund, Torbjorn and Leyland, P., and Zimmermann, Paul},
        Title = {{P. The ECMNET Project}},
        Howpublished = {{\url{http://www.loria.fr/∼zimmerma/ecmnet/}}},
        }
    <BLANKLINE>
    @misc{software-pari,
        title        = {{Pari/GP}}
        author = {{Cohen, Henri and Belabas, Karim and others}},
        organization = {The Pari~Group},
        note = {{Version~2.5.5}},
        howpublished          = {{\url{http://pari.math.u-bordeaux.fr}}},
        }

"""

###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.latex import latex
from sage.structure.sage_object import SageObject

class CitationRecord( SageObject, list ):

    def bibtex(self):
        r"""
        Bibtex string for this record.

        OUTPUT:

        A string.

        .. SEEALSO::

            :meth:`print_bibtex`

        EXAMPLES::

            sage: from sage.misc.citation_items.all import CitationRecord
            sage: from sage.misc.citation_items.ecm import ECM_CitationItem
            sage: from sage.misc.citation_items.pari import PARI_CitationItem
            sage: print CitationRecord([PARI_CitationItem(), ECM_CitationItem()]).bibtex()
            @misc{software-pari,
                title        = {{Pari/GP}}
                author = {{Cohen, Henri and Belabas, Karim and others}},
                organization = {The Pari~Group},
                note = {{Version~2.5.5}},
                howpublished          = {{\url{http://pari.math.u-bordeaux.fr}}},
                }
            <BLANKLINE>
            @Misc{software-ecm,
                Author = {Charron, T. and Daminelli, N. and Granlund, Torbjorn and Leyland, P., and Zimmermann, Paul},
                Title = {{P. The ECMNET Project}},
                Howpublished = {{\url{http://www.loria.fr/∼zimmerma/ecmnet/}}},
                }
        """
        return "\n\n".join(item.bibtex() for item in self)

    def print_bibtex(self):
        r"""
        Print the BibTeX entries for this record.

        .. SEEALSO::

            :meth:`bibtex`

        EXAMPLES::

            sage: from sage.misc.citation_items.all import CitationRecord
            sage: from sage.misc.citation_items.ecm import ECM_CitationItem
            sage: from sage.misc.citation_items.pari import PARI_CitationItem
            sage: CitationRecord([PARI_CitationItem(), ECM_CitationItem()]).print_bibtex()
            @misc{software-pari,
                title        = {{Pari/GP}}
                author = {{Cohen, Henri and Belabas, Karim and others}},
                organization = {The Pari~Group},
                note = {{Version~2.5.5}},
                howpublished          = {{\url{http://pari.math.u-bordeaux.fr}}},
                }
            <BLANKLINE>
            @Misc{software-ecm,
                Author = {Charron, T. and Daminelli, N. and Granlund, Torbjorn and Leyland, P., and Zimmermann, Paul},
                Title = {{P. The ECMNET Project}},
                Howpublished = {{\url{http://www.loria.fr/∼zimmerma/ecmnet/}}},
                }
        """
        print self.bibtex()

    def latex_citation(self):
        r"""
        A LaTeX citation statement that can be used in conjunction
        with the BibTeX entries given by :meth:`print_bibtex`.

        OUTPUT:

        A string.

        .. SEEALSO::

            :meth:`print_latex_citation`

        EXAMPLES::

            sage: from sage.misc.citation_items.all import CitationRecord
            sage: from sage.misc.citation_items.ecm import ECM_CitationItem
            sage: from sage.misc.citation_items.pari import PARI_CitationItem
            sage: CitationRecord([PARI_CitationItem(), ECM_CitationItem()]).latex_citation()
            '\\cite{software-pari, software-ecm}'
        """
        return "\\cite{{{}}}".format(", ".join(item._latex_citation_label()
                                               for item in self))

    def print_latex_citation(self):
        r"""
        Print a LaTeX citation statement that can be used in conjunction
        with the BibTeX entries given by :meth:`print_bibtex`.

        OUTPUT:

        A string.

        .. SEEALSO::

            :meth:`latex_citation`

        EXAMPLES::

            sage: from sage.misc.citation_items.all import CitationRecord
            sage: from sage.misc.citation_items.ecm import ECM_CitationItem
            sage: from sage.misc.citation_items.pari import PARI_CitationItem
            sage: CitationRecord([PARI_CitationItem(), ECM_CitationItem()]).print_latex_citation()
            Print BibTeX code by calling method print_bibtex().
            \cite{software-pari, software-ecm}
        """
        print "Print BibTeX code by calling method print_bibtex()."
        print self.latex_citation()

    def _repr_(self):
        return repr(list(self))

    def _latex_(self):
        return latex(list(self))
