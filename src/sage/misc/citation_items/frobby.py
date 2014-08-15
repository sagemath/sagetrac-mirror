###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Frobby_CitationItem( CitationItem ):
    _name = "Frobby"

    _bibtes = r"""
@misc{software-frobby,
    author = {Roune, Bjarke Hammersholt},
    title = {Frobby},
    howpublished = {{\url{http://www.broune.com/frobby}}},
}
    """

    _re = [r"^sage.interfaces.frobby"]
