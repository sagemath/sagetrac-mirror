###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class qsieve_CitationItem( CitationItem ):
    _name = "Qsieve"

    _bibtex = r"""
 @misc{software-qsieve,
    title        = {{Qsieve}}
    author = {{Reinecke, Thorsten}},
    howpublished          = {{\url{http://www.thorstenreinecke.de/qsieve}}},
}
"""

    _re = [r"^sage.interfaces.qsieve"]
