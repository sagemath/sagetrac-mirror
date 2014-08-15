###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class ECM_CitationItem( CitationItem ):
    _name = "ECM"

    _bibtex = r"""
@Misc{software-ecm,
    Author = {Charron, T. and Daminelli, N. and Granlund, Torbjorn and Leyland, P., and Zimmermann, Paul}, 
    Title = {{P. The ECMNET Project}},
    Howpublished = {{\url{http://www.loria.fr/∼zimmerma/ecmnet/}}},
}
    """

    _re = [r"^sage.interfaces.ecm"]
