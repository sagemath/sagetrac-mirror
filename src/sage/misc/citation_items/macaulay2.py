###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Macaulay2_CitationItem( CitationItem ):
    _name = "Macaulay2"

    _bibtex = r"""
@misc{software-macaulay2,
    author = {Grayson, Daniel R. and Stillman, Michael E.},
    title = {Macaulay2, a software system for research in algebraic geometry},
    howpublished = {{\ulr{http://www.math.uiuc.edu/Macaulay2}},
}
    """
