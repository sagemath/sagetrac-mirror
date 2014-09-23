###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class KASH_CitationItem( CitationItem ):
    _name = "KASH"

    _bibtex = r"""
 @manual{software-gap,
    organization = {The KANT~Group},
    title        = {{KASH, Introduction to KASH3}}
    year         = {2006},
    url          = {{\url{http://www.math.tu-berlin.de/~kant}}},
}
    """
