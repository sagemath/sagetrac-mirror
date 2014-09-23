###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class M4RI_CitationItem( CitationItem ):
    _name = "M4RI"

    _bibtex = r"""
@misc{software-m4ri,
    author       = {{Albrecht, Martin and Bard, Gregory}},
    organization = {{The M4RI~Team}},
    title        = {{The M4RI Library}},
    year         = {2012},
    note         = {{Version~20130416}},
    howpublished = {{\url{http://m4ri.sagemath.org}}},
}
    """
