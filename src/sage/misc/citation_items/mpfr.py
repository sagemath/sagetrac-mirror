###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class MPFR_CitationItem( CitationItem ):
    _name = "MPFR"

    ## based on http://www.mpfr.org/faq.html#cite
    _bibtex = r"""
@misc{software-mpfr,
    title = {{The GNU MPFR Library}},
    note = {{Version~3.1.2}},
    howpublished = {{\url{http://www.mpfr.org}}},
}
    """
