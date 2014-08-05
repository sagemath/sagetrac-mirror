###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class M4RI_CitationItem( CitationItem ):
    _name = "M4RI"

    _re = [r"^sage.matrix.matrix_mod2_dense"]
