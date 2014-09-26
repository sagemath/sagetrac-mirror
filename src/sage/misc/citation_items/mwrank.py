###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class mwrank_CitationItem( CitationItem ):
    _name = "mwrank"

    _bibtex = r"""
@misc{software-mwrank,
    author = {{Cremona, John E.}},
    title = {{mwrank}},
    note = {{Version~12.12.1}},
    howpublished = {{\url{http://homepages.warwick.ac.uk/staff/J.E.Cremona/mwrank}}},
}
    """
