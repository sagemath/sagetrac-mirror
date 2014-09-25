###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class Magma_CitationItem( CitationItem ):
    _name = "Magma"

    _bibtex = r"""
@article{software-magma,
    Author = {Bosma, Wieb and Cannon, John and Playoust, Catherine},
    Title = {{The Magma algebra system. I: The user language.}},
    FJournal = {{Journal of Symbolic Computation}},
    Journal = {{J. Symb. Comput.}},
    Volume = {24},
    Number = {3-4},
    Pages = {235--265},
    Year = {1997},
    Publisher = {Elsevier Science, London},
}
    """
