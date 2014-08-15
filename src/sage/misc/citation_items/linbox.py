###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Linbox_CitationItem( CitationItem ):
    _name = "Linbox"

    _bibtex = r"""
@Misc{software-linbox,
    organization = {{The Linbox Group}},
    title = {{Linbox: exact linear algebra with dense and blackbox matrices}},
    note = {{Version~1.3.2}},
    howpublished = {{\url{http://www.linalg.org}}},
}
    """

    _re = [r"^sage.libs.linbox"]
