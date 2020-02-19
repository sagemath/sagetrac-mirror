r"""
Index of decoders

The ``codes.decoders`` object may be used to access the decoders that Sage can build.

It is usually not necessary to access these directly: rather, the ``decoder``
method directly on a code allows you to construct all compatible decoders for
that code (:meth:`sage.coding.linear_code.AbstractLinearCode.decoder`).

**Extended code decoders**

- :class:`extended_code.ExtendedCodeOriginalCodeDecoder <sage.coding.extended_code.ExtendedCodeOriginalCodeDecoder>`

**Subfield subcode decoder**
- :class:`subfield_subcode.SubfieldSubcodeOriginalCodeDecoder <sage.coding.subfield_subcode.SubfieldSubcodeOriginalCodeDecoder>`

**Generalized Reed-Solomon code decoders**

- :class:`grs_code.GRSBerlekampWelchDecoder <sage.coding.grs_code.GRSBerlekampWelchDecoder>`
- :class:`grs_code.GRSErrorErasureDecoder <sage.coding.grs_code.GRSErrorErasureDecoder>`
- :class:`grs_code.GRSGaoDecoder <sage.coding.grs_code.GRSGaoDecoder>`
- :class:`grs_code.GRSKeyEquationSyndromeDecoder <sage.coding.grs_code.GRSKeyEquationSyndromeDecoder>`
- :class:`guruswami_sudan.gs_decoder.GRSGuruswamiSudanDecoder <sage.coding.guruswami_sudan.gs_decoder.GRSGuruswamiSudanDecoder>`

**Generic decoders**

- :class:`linear_code.LinearCodeNearestNeighborDecoder <sage.coding.linear_code.LinearCodeNearestNeighborDecoder>`
- :class:`linear_code.LinearCodeSyndromeDecoder <sage.coding.linear_code.LinearCodeSyndromeDecoder>`
- :class:`information_set_decoder.LinearCodeInformationSetDecoder <sage.coding.information_set_decoder.LinearCodeInformationSetDecoder>`

**Cyclic code decoder**

- :class:`cyclic_code.CyclicCodeSurroundingBCHDecoder <sage.coding.cyclic_code.CyclicCodeSurroundingBCHDecoder>`

**BCH code decoder**

- :class:`bch_code.BCHUnderlyingGRSDecoder <sage.coding.bch_code.BCHUnderlyingGRSDecoder>`

**Punctured codes decoders**

- :class:`punctured_code.PuncturedCodeOriginalCodeDecoder <sage.coding.punctured_code.PuncturedCodeOriginalCodeDecoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.decoders_catalog import *
"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#                     2015 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import('sage.coding.bch_code', 'BCHUnderlyingGRSDecoder')
lazy_import('sage.coding.cyclic_code', 'CyclicCodeSurroundingBCHDecoder')
lazy_import('sage.coding.extended_code', 'ExtendedCodeOriginalCodeDecoder')
lazy_import('sage.coding.grs_code', ['GRSBerlekampWelchDecoder',
                                      'GRSErrorErasureDecoder',
                                      'GRSGaoDecoder',
                                      'GRSKeyEquationSyndromeDecoder'])

from .guruswami_sudan.gs_decoder import GRSGuruswamiSudanDecoder
lazy_import('sage.coding.linear_code', ['LinearCodeNearestNeighborDecoder',
                                         'LinearCodeSyndromeDecoder',
                                         'LinearCodeInformationSetDecoder'])

lazy_import('sage.coding.punctured_code', 'PuncturedCodeOriginalCodeDecoder')
lazy_import('sage.coding.subfield_subcode', 'SubfieldSubcodeOriginalCodeDecoder')
lazy_import('sage.coding.information_set_decoder', 'LinearCodeInformationSetDecoder')

del lazy_import
