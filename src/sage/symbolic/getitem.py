r"""
Operands (deprecated module)

This module consists only of deprecated lazy imports from
:mod:`sage.symbolic.expression`.
"""


from sage.misc.lazy_import import lazy_import
lazy_import('sage.symbolic.expression',
            ['normalize_index_for_doctests', 'OperandsWrapper', 'restore_op_wrapper'],
            deprecation=32386)
