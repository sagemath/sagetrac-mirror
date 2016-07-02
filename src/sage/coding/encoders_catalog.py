r"""
Index of encoders

The ``codes.encoders`` object may be used to access the encoders that Sage can build.

**Generic encoders**

:class:`linear_code.LinearCodeGeneratorMatrixEncoder <sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder>`
:class:`linear_code.LinearCodeParityCheckEncoder <sage.coding.linear_code.LinearCodeParityCheckEncoder>`

**Generalized Reed-Solomon code encoders**

- :class:`grs.GRSEvaluationVectorEncoder <sage.coding.grs.GRSEvaluationVectorEncoder>`
- :class:`grs.GRSEvaluationPolynomialEncoder <sage.coding.grs.GRSEvaluationPolynomialEncoder>`

**Extended code encoders**

- :class:`extended_code.ExtendedCodeExtendedMatrixEncoder <sage.coding.extended_code.ExtendedCodeExtendedMatrixEncoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.encoders_catalog import *
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

from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import('sage.coding.linear_code', ['LinearCodeGeneratorMatrixEncoder',
                                         'LinearCodeParityCheckEncoder'])
_lazy_import('sage.coding.grs', ['GRSEvaluationVectorEncoder', 'GRSEvaluationPolynomialEncoder'])
_lazy_import('sage.coding.reed_muller_code', ['ReedMullerVectorEncoder', 'ReedMullerPolynomialEncoder'])
_lazy_import('sage.coding.extended_code', 'ExtendedCodeExtendedMatrixEncoder')
_lazy_import('sage.coding.punctured_code', 'PuncturedCodePuncturedMatrixEncoder')
_lazy_import('sage.coding.subfield_subcode', 'SubfieldSubcodeParityCheckEncoder')
