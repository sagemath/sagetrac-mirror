include "cysignals/signals.pxi"

def debugstack():
    r"""
    Print the internal PARI variables ``top`` (top of stack), ``avma``
    (available memory address, think of this as the stack pointer),
    ``bot`` (bottom of stack).

    EXAMPLE::

    sage: pari.debugstack()  # random
    top =  0x60b2c60
    avma = 0x5875c38
    bot =  0x57295e0
    size = 1000000
    """
    # We deliberately use low-level functions to minimize the
    # chances that something goes wrong here (for example, if we
    # are out of memory).
    printf("top =  %p\navma = %p\nbot =  %p\nsize = %lu\n", pari_mainstack.top, avma, pari_mainstack.bot, <unsigned long>pari_mainstack.rsize)
    fflush(stdout)

cdef inline void clear_stack():
    """
    Call ``sig_off()``. If we are leaving the outermost
    ``sig_on() ... sig_off()`` block, then clear the PARI stack.
    """
    #global avma
    if sig_on_count <= 1:
        avma = pari_mainstack.top
    sig_off()

def stacksize():
    r"""
    Return the current size of the PARI stack, which is `10^6`
    by default.  However, the stack size is automatically doubled
    when needed up to some maximum.

    .. SEEALSO::
    
    - :meth:`stacksizemax` to get the maximum stack size
    - :meth:`allocatemem` to change the current or maximum
    stack size
    
    EXAMPLES::
    
    sage: pari.stacksize()
    1000000
    sage: pari.allocatemem(2^18, silent=True)
    sage: pari.stacksize()
    262144
    """
    return pari_mainstack.size

def stacksizemax():
    r"""
    Return the maximum size of the PARI stack, which is determined
    at startup in terms of available memory. Usually, the PARI
    stack size is (much) smaller than this maximum but the stack
    will be increased up to this maximum if needed.

    .. SEEALSO::

    - :meth:`stacksize` to get the current stack size
    - :meth:`allocatemem` to change the current or maximum
    stack size

    EXAMPLES::

    sage: pari.allocatemem(2^18, 2^26, silent=True)
    sage: pari.stacksizemax()
    67108864
    """
    return pari_mainstack.vsize

def allocatemem(size_t s=0, size_t sizemax=0, *, silent=False):
    r"""
    Change the PARI stack space to the given size ``s`` (or double
    the current size if ``s`` is `0`) and change the maximum stack
    size to ``sizemax``.

    PARI tries to use only its current stack (the size which is set
    by ``s``), but it will increase its stack if needed up to the
    maximum size which is set by ``sizemax``.

    The PARI stack is never automatically shrunk.  You can use the
    command ``pari.allocatemem(10^6)`` to reset the size to `10^6`,
    which is the default size at startup.  Note that the results of
    computations using Sage's PARI interface are copied to the
    Python heap, so they take up no space in the PARI stack.
    The PARI stack is cleared after every computation.

    It does no real harm to set this to a small value as the PARI
    stack will be automatically doubled when we run out of memory.

    INPUT:

    - ``s`` - an integer (default: 0).  A non-zero argument is the
    size in bytes of the new PARI stack.  If `s` is zero, double
    the current stack size.

    - ``sizemax`` - an integer (default: 0).  A non-zero argument
    is the maximum size in bytes of the PARI stack.  If
    ``sizemax`` is 0, the maximum of the current maximum and
    ``s`` is taken.

    EXAMPLES::

    sage: pari.allocatemem(10^7)
    PARI stack size set to 10000000 bytes, maximum size set to 67108864
    sage: pari.allocatemem()  # Double the current size
    PARI stack size set to 20000000 bytes, maximum size set to 67108864
    sage: pari.stacksize()
    20000000
    sage: pari.allocatemem(10^6)
    PARI stack size set to 1000000 bytes, maximum size set to 67108864

    The following computation will automatically increase the PARI
    stack size::

    sage: a = pari('2^100000000')

    ``a`` is now a Python variable on the Python heap and does not
    take up any space on the PARI stack.  The PARI stack is still
    large because of the computation of ``a``::

    sage: pari.stacksize()
    16000000

    Setting a small maximum size makes this fail::

    sage: pari.allocatemem(10^6, 2^22)
    PARI stack size set to 1000000 bytes, maximum size set to 4194304
    sage: a = pari('2^100000000')
    Traceback (most recent call last):
    ...
    PariError: _^s: the PARI stack overflows (current size: 1000000; maximum size: 4194304)
    You can use pari.allocatemem() to change the stack size and try again

    TESTS:

    Do the same without using the string interface and starting
    from a very small stack size::

    sage: pari.allocatemem(1, 2^26)
    PARI stack size set to 1024 bytes, maximum size set to 67108864
    sage: a = pari(2)^100000000
    sage: pari.stacksize()
    16777216

    We do not allow ``sizemax`` less than ``s``::

    sage: pari.allocatemem(10^7, 10^6)
    Traceback (most recent call last):
    ...
    ValueError: the maximum size (10000000) should be at least the stack size (1000000)
    """
    if s == 0:
        s = pari_mainstack.size * 2
        if s < pari_mainstack.size:
            raise OverflowError("cannot double stack size")
    elif s < 1024:
        s = 1024  # arbitrary minimum size
    if sizemax == 0:
        # For the default sizemax, use the maximum of current
        # sizemax and the given size s.
        if pari_mainstack.vsize > s:
            sizemax = pari_mainstack.vsize
        else:
            sizemax = s
    elif sizemax < s:
        raise ValueError("the maximum size ({}) should be at least the stack size ({})".format(s, sizemax))
    sig_on()
    paristack_setsize(s, sizemax)
    sig_off()
    if not silent:
        print("PARI stack size set to {} bytes, maximum size set to {}".
                  format(stacksize(), stacksizemax()))
