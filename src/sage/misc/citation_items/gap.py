###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class GAP_CitationItem( CitationItem ):
    _name = "GAP"

    ## provided by GAP homepage
    _bibtex = r"""
 @manual{software-gap,
    title        = {{GAP -- Groups, Algorithms, and Programming}}
    organization = {The GAP~Group},
    year         = {2014},
    note  = {{Version~4.7.5, \url{http://www.gap-system.org}}},
}
    """
