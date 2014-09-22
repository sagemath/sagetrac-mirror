###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Givaro_CitationItem( CitationItem ):
    _name = "Givaro"

    _bibtex = r"""
@Misc{software-givaro,
    author = {{Gautier, Thierry and Roch, Jean-Louis and Villard, Gilles and others}},
    title = {{Givaro}},
    note = {{Version~3.7.1}},
    howpublished = {{\url{http://ljk.imag.fr/CASYS/LOGICIELS/givaro}}},
}
    """
