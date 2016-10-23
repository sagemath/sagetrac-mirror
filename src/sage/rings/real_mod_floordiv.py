r"""
Transitional module for remainder and floor division of real numbers

The behavior of the remainder and floordiv operators (``%`` and ``//``) used to be inconsistent
for real numbers. As discussed in https://groups.google.com/forum/#!topic/sage-devel/PfMop0nyiL0
this will be changed in a future SageMath release as follows:

The remainder `x % y` for two real numbers is the unique real number in `[0, y)` (if `y` is
positive) or `(y, 0]` (if `y` is negative) of the form `x + n y` with `n` an integer. This
integer `n` is the floor division `x // y`.

This module defines a global variable ``NEW`` so that when set to ``False`` the
old beaviour is preserved but a deprecation warning is raised.  While if set to
``True`` the modulo and floordiv operators follow the above specifications.

EXAMPLES::

    sage: import sage.rings.real_mod_floordiv
    sage: sage.rings.real_mod_floordiv.NEW = True
    sage: 5/3 % 1
    2/3
    sage: (7/5) % (-1/3)
    -4/15

    sage: sage.rings.real_mod_floordiv.NEW = False
    sage: 5/3 % 1
    doctest:...: DeprecationWarning: In a future Sage release the behavior of the modulo
    operator % and floor division operator // involving real numbers will change. If you
    want to use the new behavior and get rid of this message execute the two following
    commands
        import sage.rings.real_mod_floordiv
        sage.rings.real_mod_floordiv.NEW = True
    See http://trac.sagemath.org/21745 for details.
    0
    sage: (7/5) % (-1/3)
    Traceback (most recent call last):
    ...
    TypeError: no conversion of this rational to integer

TESTS:

Check that both implementations of remainders coincide on integer entries::

    sage: import sage.rings.real_mod_floordiv
    sage: input = [(17, 3), (34123, 222), (25, -3), (1234, -555434), (-12, 5),
    ....:     (-1, 4323), (-22, -7), (-11212, -532)]
    sage: a_old = [QQ(x) % QQ(y) for x,y in input]
    sage: sage.rings.real_mod_floordiv.NEW = True
    sage: a_new = [QQ(x) % QQ(y) for x,y in input]
    sage: sage.rings.real_mod_floordiv.NEW = False
    sage: a_old == a_new
    True

Check that initialization of integer mod elements works with both old and new
versions::

    sage: import sage.rings.real_mod_floordiv
    sage: R1 = Zmod(2**10)
    sage: R2 = Zmod(2**40)
    sage: R3 = Zmod(2**100)
    sage: R4 = GF(3**20)
    sage: R5 = GF(18446744073709551629)
    sage: b_old = (R1(5/3), R2(5/3), R3(5/3), R4(5/7), R5(12/35))
    sage: sage.rings.real_mod_floordiv.NEW = True
    sage: b_new = (R1(5/3), R2(5/3), R3(5/3), R4(5/7), R5(12/35))
    sage: sage.rings.real_mod_floordiv.NEW = False
    sage: b_old == b_new
    True
"""
NEW = False    # whether or not the new behavior is used

def mod_floordiv_deprecation():
    r"""
    Specific deprecation for mod/floordiv of real numbers

    See :trac:`21745`.

    TESTS::

        sage: from sage.rings.real_mod_floordiv import mod_floordiv_deprecation
        sage: mod_floordiv_deprecation()
        doctest:...: DeprecationWarning: In a future Sage release the behavior
        of the modulo operator % and floor division operator // involving real
        numbers will change. If you want to use the new behavior and get rid of
        this message execute the two following commands
            import sage.rings.real_mod_floordiv
            sage.rings.real_mod_floordiv.NEW = True
        See http://trac.sagemath.org/21745 for details.
    """
    from sage.misc.superseded import deprecation
    deprecation(21745, "In a future Sage release the behavior of the modulo operator % and floor division operator // involving real numbers will change. If you want to use the new behavior and get rid of this message execute the two following commands\n    import sage.rings.real_mod_floordiv\n    sage.rings.real_mod_floordiv.NEW = True")
