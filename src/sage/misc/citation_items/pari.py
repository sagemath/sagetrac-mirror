###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class PARI_CitationItem( CitationItem ):
    _name = "PARI"

    _re = [r"^sage.libs.pari",
           "^sage.interfaces.gp"]
