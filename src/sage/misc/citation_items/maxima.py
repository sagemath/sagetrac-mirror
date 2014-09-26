###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Maxima_CitationItem( CitationItem ):
    _name = "Maxima"

    ## based on http://sourceforge.net/p/maxima/wiki/FAQ/#what-is-the-appropriate-way-to-cite-maxima-in-an-academic-context
    _bibtex = r"""
@misc{software-maxima,
    title = {{Maxima, a Computer Algebra System.}},
    note = {{Version~5.33}},
    howpublished = {{\url{http://maxima.sourceforge.net}}},
}
    """
