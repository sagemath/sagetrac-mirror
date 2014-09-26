###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class MuPAD_CitationItem( CitationItem ):
    _name = "MuPAD"

    _bibtex = r"""
@misc{software-mupad,
    title = {{MuPAD}},
    organization = {{The MathWorks, Inc.}},
    howpublished = {{\url{http://www.mathworks.de/discovery/mupad.html}}},
}
    """
