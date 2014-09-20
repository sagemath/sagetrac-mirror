###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class numpy_CitationItem( CitationItem ):
    _name = "numpy"

    _bibtex = r"""
 @misc{software-numpy,
    organization = {The NumPy~Group},
    title        = {{NumPy}}
    note = {{Version~1.8.1}},
    howpublished = {{\url{http://www.numpy.org}}},
}
    """
