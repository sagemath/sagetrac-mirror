###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class scipy_CitationItem( CitationItem ):
    _name = "scipy"

    _bibtex = r"""
 @misc{software-scipy,
    organization = {The SciPy~Group},
    title        = {{SciPy}},
    note = {{Version~0.14}},
    howpublished          = {{\url{http://www.scipy.org}}},
}
    """
