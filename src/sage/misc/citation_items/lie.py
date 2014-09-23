###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class LiE_CitationItem( CitationItem ):
    _name = "LiE"

    _bibtex = r"""
@book{software-lie,
    author = {van Leeuwen, M. A. A. and Cohen, A. M. and Lisser, B.},
    title = {{LiE, A Package for Lie Group Computations}},
    Publisher = {Computer Algebra Nederland, Amsterdam},
    year = {1992},
}

    """
