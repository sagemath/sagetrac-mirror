###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class MPFI_CitationItem( CitationItem ):
    _name = "MPFI"

    _bibtex = r"""
@misc{software-mpfi,
    title = {{MPFI, a multiple precision interval arithmetic library based on MPFR}},
    note = {{Version~1.5.1}},
    howpublished = {{\url{http://perso.ens-lyon.fr/nathalie.revol/software.html}}},
}
    """
