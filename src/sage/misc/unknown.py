"""
The Unknown and Undecidable truth values

"Undecidable" refers to truth values that are known to be formally undecidable,
while "unknown" may be either undecidable or indeterminate due to limitations
in knowledge or algorithms.

We use Kleene's three-valued logic for operations between a boolean and either
``Unknown`` or ``Undecidable``. For operations between ``Unknown`` and
``Undecidable``, ``Unknown`` is always returned, since ``Unknown`` alsoincludes
the possibility of undecidability.

The truth tables are thus

+-------------+-------+-------------+-------------+---------+
| AND         | False | True        | Undecidable | Unknown |
+=============+=======+=============+=============+=========+
| False       | False |             |             |         |
+-------------+-------+-------------+-------------+---------+
| True        | False | True        |             |         |
+-------------+-------+-------------+-------------+---------+
| Undecidable | False | Undecidable | Undecidable |         |
+-------------+-------+-------------+-------------+---------+
| Unknown     | False | Unknown     | Unknown     | Unknown |
+-------------+-------+-------------+-------------+---------+

+-------------+-------------+-------+-------------+---------+
| OR          | False       | True  | Undecidable | Unknown |
+=============+=============+=======+=============+=========+
| False       | False       |       |             |         |
+-------------+-------------+-------+-------------+---------+
| True        | True        | True  |             |         |
+-------------+-------------+-------+-------------+---------+
| Undecidable | Undecidable | True  | Undecidable |         |
+-------------+-------------+-------+-------------+---------+
| Unknown     | Unknown     | True  | Unknown     | Unknown |
+-------------+-------------+-------+-------------+---------+

+-------------+-------------+
| NOT         |             |
+=============+=============+
| False       | True        |
+-------------+-------------+
| True        | False       |
+-------------+-------------+
| Undecidable | Undecidable |
+-------------+-------------+
| Unknown     | Unknown     |
+-------------+-------------+

AUTHORS:

- Florent Hivert (2010): initial version
- Eviatar Bach (2016): adding Undecidable value, other changes
"""

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
class UnknownClass(UniqueRepresentation, SageObject):
    """
    TESTS::

        sage: TestSuite(Unknown).run()
    """
    def __init__(self, unknown_type):
        """
        The ``Unknown`` or ``Undecidable`` truth value.

        EXAMPLES::

            sage: l = [False, True, Undecidable, Unknown]
            sage: for a in l: print([a & b for b in l])
            [False, False, False, False]
            [False, True, Undecidable, Unknown]
            [False, Undecidable, Undecidable, Unknown]
            [False, Unknown, Unknown, Unknown]

            sage: for a in l: print([a | b for b in l])
            [False, True, Undecidable, Unknown]
            [True, True, True, True]
            [Undecidable, True, Undecidable, Unknown]
            [Unknown, True, Unknown, Unknown]

        ..warning:: In the following cases, ``and``, ``not`` and ``or`` raise
                    ValueErrors. This is due to the fact that the logical
                    operators cannot be overridden in Python::

            sage: not Unknown         # should return Unknown
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
            sage: Unknown and False   # should return False
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
            sage: Unknown or False    # should return Unknown
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown

        To do logical operations on Unknown, use the bitwise operators::
            sage: Unknown & True
            Unknown
            sage: ~Unknown
            Unknown
            sage: Unknown | False
            Unknown
        """
        self.unknown_type = unknown_type

    def _repr_(self):
        """
        TESTS::

            sage: Unknown
            Unknown
            sage: Undecidable
            Undecidable
        """
        return self.unknown_type.title()

    def __nonzero__(self):
        """
        When evaluated in a boolean context ``Unknown()`` raises a ValueError

        EXAMPLES::

            sage: bool(Unknown)
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
            sage: not Unknown
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
            sage: bool(Undecidable)
            Traceback (most recent call last):
            ...
            ValueError: Truth value undecidable
        """
        raise ValueError("Truth value {utype}".format(utype=self.unknown_type))

    def __and__(self, other):
        """
        The ``and`` logical connector.

        ..warning:: This is not used by ``and``

        EXAMPLES::

            sage: Unknown & False
            False
            sage: Unknown & Unknown
            Unknown
            sage: Unknown & True
            Unknown
            sage: Unknown & Undecidable
            Unknown
            sage: Undecidable & Unknown
            Unknown

        Compare with::

            sage: Unknown and False
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
        """
        if other is False:
            return False
        elif other is True:
            return self
        else:
            if self.unknown_type == 'unknown':
                return self
            elif other.unknown_type == 'undecidable':
                # self Undecidable, other Undecidable
                return self
            else:
                # self Undecidable, other Unknown
                return other

    __rand__ = __and__

    def __or__(self, other):
        """
        The ``or`` logical connector.

        ..warning:: This is not used by ``or``

        EXAMPLES::

            sage: Unknown | False
            Unknown
            sage: Unknown | Unknown
            Unknown
            sage: Unknown | True
            True
            sage: Unknown | Undecidable
            Unknown
            sage: Undecidable | Unknown
            Unknown

        Compare with::

            sage: Unknown or False    # should return Unknown
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
        """
        if other is True:
            return True
        elif other is False:
            return self
        else:
            if self.unknown_type == 'unknown':
                return self
            elif other.unknown_type == 'undecidable':
                # self Undecidable, other Undecidable
                return self
            else:
                # self Undecidable, other Unknown
                return other

    __ror__ = __or__

    def __invert__(self):
        """
        The ``not`` logical connector.

        ..warning:: This is not used by ``not``

        EXAMPLES::

            sage: ~Unknown
            Unknown
            sage: ~Undecidable
            Undecidable

        Compare with::

            sage: not Unknown  # should return Unknown
            Traceback (most recent call last):
            ...
            ValueError: Truth value unknown
        """
        return self

    def __cmp__(self, other):
        """
        Comparison of truth value.

        EXAMPLES::

            sage: l = [False, Unknown, True]
            sage: for a in l: print ([a < b for b in l])
            [False, True, True]
            [False, False, True]
            [False, False, False]

            sage: for a in l: print ([a <= b for b in l])
            [True, True, True]
            [False, True, True]
            [False, False, True]
        """
        if other is self:
            return 0
        if isinstance(other, bool):
            if other:
                return -1
            else:
                return +1
        else:
            raise ValueError("Unable to compare {} with {}".format(self, other))

Unknown = UnknownClass('unknown')
Undecidable = UnknownClass('undecidable')
