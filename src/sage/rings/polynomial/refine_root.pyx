"""
Refine polynomial roots using Newton--Raphson

This is an implementation of the Newton--Raphson algorithm to
approximate roots of real or complex polynomials. The implementation
is based on interval arithmetic.

AUTHORS:

- Carl Witty (2007-11-18): initial version

- Jeroen Demeyer (2015-10-06): improved code, also allow real input
  and use bisection if needed (see :trac:`19362`)
"""

#*****************************************************************************
#       Copyright (C) 2007 Carl Witty
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.real_mpfi cimport RealIntervalField_class


def refine_root(poly, deriv, root, field, long steps=10):
    """
    Given a polynomial and an approximation to a root (defined over
    a real or complex interval field), refine the root until we have
    isolated a unique root.

    INPUT:

    - ``poly`` -- a polynomial with interval coefficients.

    - ``deriv`` -- the derivative of ``poly``.

    - ``root`` -- an interval approximating a root of ``poly``. This
      does not need to contain a root yet.

    - ``field`` -- a real or complex interval field in which all
      computations take place. The base ring of ``poly``, ``deriv``
      and the parent of ``root`` all must be equal to ``field``.

    - ``steps`` -- (default: 10) maximum number of Newton-Raphson steps
      to do.

    OUTPUT: an interval with an isolated root. Unless the number of
    ``steps`` is too small, the accuracy of the interval should be
    close to the best possible.

    If we fail to refine the root, we raise ``ArithmeticError``.

    EXAMPLES::

        sage: from sage.rings.polynomial.refine_root import refine_root
        sage: R.<x> = RIF[]
        sage: pol = x^2 - x - 1
        sage: r = refine_root(pol, pol.derivative(), RIF(2), RIF)
        sage: r
        1.618033988749895?
        sage: r.relative_diameter()
        1.37231112862212e-16

    A complex example::

        sage: R.<x> = ZZ[]
        sage: p = x^9 - 1
        sage: ip = CIF['x'](p); ip
        x^9 - 1
        sage: ipd = CIF['x'](p.derivative()); ipd
        9*x^8
        sage: irt = CIF(ComplexField(20)(cos(2*pi/9), sin(2*pi/9))); irt
        0.76604366302490235? + 0.64278793334960938?*I
        sage: ip(irt)
        -3.505875839?e-6 + 6.744354905?e-6*I
        sage: refine_root(ip, ipd, irt, CIF)
        0.766044443118978? + 0.642787609686540?*I

    This example requires bisection since the given interval contains
    a zero of the derivative::

        sage: R.<x> = RIF[]
        sage: pol = 10*x^6 - 10*x^2 + 1
        sage: refine_root(pol, pol.derivative(), RIF(-1/5, 2/3), RIF)
        0.3178541456092885?

    With too few steps, the refining does not work::

        sage: refine_root(ip, ipd, irt, CIF, steps=1)
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot refine polynomial root (not enough steps?)
        sage: refine_root(ip, ipd, irt, CIF, steps=2)
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot refine polynomial root (not enough steps?)
        sage: refine_root(ip, ipd, irt, CIF, steps=3)
        0.766044443118978? + 0.642787609686540?*I

    Giving a very wrong root also gives an error::

        sage: irt = CIF(0)
        sage: refine_root(ip, ipd, irt, CIF)
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot refine polynomial root (precision too low?)

    If we give an imprecise polynomial, we get an imprecise result::

        sage: R.<x> = RIF[]
        sage: pol = x^3 - RIF(4^3, 5^3)
        sage: r = refine_root(pol, pol.derivative(), RIF(10), RIF)
        sage: r.endpoints()
        (3.00612159221905, 6.10936449050236)
    """
    # We start with a basic fact: if we do an interval Newton-Raphson
    # step, and the refined interval is contained in the original interval,
    # then the refined interval contains exactly one root.

    # Unfortunately, our initial estimated root almost certainly does
    # not contain the actual root (typically, our initial interval is a
    # point, which is exactly equal to whatever floating-point estimate
    # we got from the external solver).  So we need to do multiple
    # Newton-Raphson steps, and check this inclusion property on each
    # step.

    # After a few steps of refinement, if the initial root estimate was
    # close to a root, we should have an essentially perfect interval
    # bound on the root (since Newton-Raphson has quadratic convergence),
    # unless either the real or imaginary component of the root is zero.
    # If the real or imaginary component is zero, then we could spend
    # a long time computing closer and closer approximations to that
    # component.  (This doesn't happen for non-zero components, because
    # of the imprecision of floating-point numbers combined with the
    # outward interval rounding; but close to zero, MPFI provides
    # extremely precise numbers.)

    # If the slope in our current interval is not bounded away
    # from zero, we instead try to refine the root using bisection.
    # This might be needed if the original interval was very large
    # (typically the output of real root isolation algorithms).

    # After every refinement step, we check to see if the real or
    # imaginary component of our interval is close to zero.  If so, we
    # try setting it to exactly zero.  This gives us a good chance of
    # detecting real roots.  However, we do this replacement at most
    # once per component.

    cdef bint is_real = isinstance(field, RealIntervalField_class)
    cdef bint smashed_real = False
    cdef bint smashed_imag = False
    cdef bint converging = False

    if is_real:
        real_field = field
    else:
        real_field = field._real_field()

    epsilon = real_field(-1, 1) >> field.prec()

    cdef long i, n
    for i in range(steps):
        # Check for possibly exact zero real/imaginary components
        if is_real:
            if not smashed_real and root in epsilon:
                root = field.zero()
                smashed_real = True
                converging = False
        else:
            if not smashed_imag and root.imag() in epsilon:
                root = field(root.real(), 0)
                smashed_imag = True
                converging = False
            elif not smashed_real and root.real() in epsilon:
                root = field(0, root.imag())
                smashed_real = True
                converging = False

        slope = deriv(root)
        if slope.contains_zero():
            # We cannot use Newton-Raphson => use bisection instead
            pieces = root.bisection()
            pieces = [p for p in pieces if poly(p).contains_zero()]
            n = len(pieces)
            if n == 1:
                nroot = pieces[0]
            elif n == 0:
                raise ArithmeticError("cannot refine polynomial root (precision too low?)")
            else:
                # Still more than 1 piece left: it could be that the
                # zero of the polynomial is very close to the center of
                # the interval. We try an additional bisection. We take
                # the union of all pieces containing a zero.
                nroot = None
                for p in pieces:
                    for q in p.bisection():
                        if poly(q).contains_zero():
                            if nroot is None:
                                nroot = q
                            else:
                                nroot = nroot.union(q)

                if nroot is None:
                    raise ArithmeticError("cannot refine polynomial root (precision too low?)")
                if not nroot._cmp_(root):
                    # nroot and root are exactly the same interval
                    raise ArithmeticError("cannot refine polynomial root (multiple roots?)")
        else:
            # Use interval Newton-Raphson to refine the root
            center = field(root.center())
            nroot = center - poly(center) / slope

            if converging:
                # If we are converging, make sure we never enlarge our
                # interval again. This is in particular needed when the
                # interval gets small to avoid oscillations.
                nroot = nroot.intersection(root)
            elif nroot in root:
                # Once we are converging, assume that we keep converging.
                converging = True

            if converging:
                if i == steps - 1:
                    # This is the last iteration
                    return nroot
                if (nroot.diameter() << 2) >= root.diameter():
                    # If the new diameter is within a factor 4 of the
                    # original diameter, then we have converged reasonably
                    # well.
                    return nroot
            else:
                # If the new interval still isn't contained in the old
                # after a while, try tripling the size of the region
                if 2*i >= steps:
                    # This gives an interval with the same center
                    # but with a larger diameter
                    nroot += nroot - nroot

        root = nroot

    raise ArithmeticError("cannot refine polynomial root (not enough steps?)")
