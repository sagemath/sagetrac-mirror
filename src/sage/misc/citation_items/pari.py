###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class PARI_CitationItem( CitationItem ):
    _name = "PARI"

    _bibtex = r"""
 @misc{software-pari,
    title        = {{Pari/GP}}
    author = {{Cohen, Henri and Belabas, Karim and others}},
    organization = {The Pari~Group},
    note = {{Version~2.5.5}},
    howpublished          = {{\url{http://pari.math.u-bordeaux.fr}}},
}
    """
