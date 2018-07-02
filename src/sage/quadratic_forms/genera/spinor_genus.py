r"""
Spinor genus computations.

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Simon Brandhorst (2018-07-1): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


class AdelicSquareClasses(AbelianGroupGap):
    r"""

    INPUT:

    - a tuple of primes

    EXAMPLES::


    """
    def __init__(self, primes):
        r"""
        """
        if primes[0] != 2:
            raise ValueError("first prime must be 2")
        self._primes = primes
        orders = len(self._primes)*[2] + [2]
        # unit, val, unit, val, unit, val
        order = tuple(orders)
        AbelianGroupGap.__init__(self, orders)

    def to_square_class(self, x, p):
        r"""
        Return `(1,...,1,x,1,...,1)` with the square class of `x` at position `p`.

        INPUT:

        - ``p`` -- a prime

        - ``x```-- a non zero rational number

        EXAMPLES::

            sage: AS = AdelicSquareClasses((2,3,7))
            sage: AS.to_square_class(5,7)
            f6
            sage: AS.to_square_class(5,2)
            f2
            sage: AS.to_square_class(-5,2)
            f1
            sage: AS.to_square_class(7,2)
            f1*f2
        """
        x = QQ(x)
        if x == 0:
            raise ValueError("x must be non zero")
        if not p in self._primes:
            raise ValueError("not a coordinate prime")
        v, u = x.val_unit(p)
        if v != 0:
            raise ValueError("x(=%s) must be a p-adic unit" %x)
        y = self.one()
        if p == 2:
            u = u % 8
            if u == 3:
                y *= self.gens()[0]
            if u == 5:
                y *= self.gens()[1]
            if u == 7:
                y *= self.gens()[0] * self.gens()[1]
            return y
        i = 1 + self._primes.index(p)
        if not u.is_padic_square(p):
            y *= self.gens()[i]
        return y

    def delta(self, r, p=None):
        r"""
        Diagonal embedding of rational square classes.

        INPUT:

        - ``r`` -- a non zero rational number

        - ``p`` -- a prime

        EXAMPLES::

            sage: AS = AdelicSquareClasses((2,3,7))
            sage: AS.delta(2,p=3)
            f4
            sage: AS.delta(2)
            f3*f4
        """
        r = QQ(r)
        if p is None:
            return self.prod([self.to_square_class(r, p) for p in self._primes])
        if p == -1:
            r = r.sign()
            return self.prod([self.to_square_class(r, p) for p in self._primes])
        v, u = r.val_unit(p)
        pv = p**v
        y = self.prod([self.to_square_class(pv,q) for q in self._primes if q!=p])
        y *= self.to_square_class(u, p)
        return y
