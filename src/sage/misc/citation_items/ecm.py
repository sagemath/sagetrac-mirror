###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.misc.citation_items.citation_item import CitationItem

class ECM_CitationItem( CitationItem ):
    _name = "ECM"

    ## based on reference given by Paul Zimmermann via email
    _bibtex = r"""
@InProceedings{software-ecm,
  author =       {Zimmermann, Paul and Dodson, Bruce},
  title =        {20 years of {ECM}},
  booktitle =    {Proceedings of the 7th Algorithmic Number Theory Symposium
                  (ANTS VII)},
  pages =        {525--542},
  year =         2006,
  editor =       {Hess, F. and Pauli, S. and Pohst, M.},
  volume =       4076,
  series =       {{LNCS}},
  address =      {Berlin Heidelberg},
  publisher =    {{Springer}}
}
    """

