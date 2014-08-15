###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Tachyon_CitationItem( CitationItem ):
    _name = "Tachyon"

    ## is not mathematical software
    _bibtex = None

    _re = [r"^sage.interfaces.tachyon"]
