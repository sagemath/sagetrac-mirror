"""
Refine polynomial roots using Newton--Raphson

This is an implementation of the Newton--Raphson algorithm to
approximate roots of real or complex polynomials. The implementation
is based on interval arithmetic.

AUTHORS:

- Carl Witty (2007-11-18): initial version

- Jeroen Demeyer (2015-10-06): improved code, also allow real input
  and trim intervals if needed (see :trac:`19362`)
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

DEF DEBUG_REFINE_ROOT = 0

cdef inline intervalstr(x):
    return x.str(style="brackets")


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

    The following 2 examples require trimming since the given interval
    contains a zero of the derivative::

        sage: R.<x> = RIF[]
        sage: pol = 10*x^6 - 10*x^2 + 1
        sage: refine_root(pol, pol.derivative(), RIF(-1/5, 2/3), RIF)
        0.3178541456092885?
        sage: pol = (x+2)*x*(x-2)
        sage: refine_root(pol, pol.derivative(), RIF(1/20, 10), RIF)
        2
        sage: refine_root(pol, pol.derivative(), RIF(6, 10), RIF)
        2

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
        (3.80661823145203, 5.34351492272488)
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
    # from zero, we instead try to trim the interval containing the
    # root by removing a small part along the edges. This might be
    # needed if the original interval was very large (typically the
    # output of real root isolation algorithms).

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

    IF DEBUG_REFINE_ROOT:
        print("Calling refine_root() with polynomial %r" % poly)

    cdef long i, j, jtrim
    cdef double t
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

        IF DEBUG_REFINE_ROOT:
            print("step %i: root = %s" % (i, intervalstr(root)))

        slope = deriv(root)
        if slope.contains_zero():
            # We cannot use Newton-Raphson => trim interval instead.
            # If a certain part along the edge of our interval does not
            # contain a zero of the polynomial, we remove that part.
            edges = root.edges()
            jtrim = -1  # Index of last trimmed edge (-1 if none)
            for j in range(len(edges)):
                e = edges[j]
                if poly(e).contains_zero():
                    # Polynomial is zero on the edge => cannot trim
                    continue

                opposite = edges[j^1]  # Opposite edge
                diff = opposite.center() - e.center()
                # t = fraction of the interval to trim
                if j == jtrim | 1:
                    # We just trimmed the opposite edge: set t to the
                    # value which corresponds to the center of the
                    # original interval (before trimming).
                    t = 0.5/(1 - t)
                else:
                    t = 0.5

                while t >= 0.001:
                    trim_edge = e + field(diff*t)
                    remove = e.union(trim_edge)
                    if not poly(remove).contains_zero():
                        # The polynomial is certainly not zero
                        # on the trimmed part => trim root
                        root = trim_edge.union(opposite)
                        edges = root.edges()
                        jtrim = j
                        IF DEBUG_REFINE_ROOT:
                            print("trim edge[%s] by %s"%(j, t))
                        break
                    t *= 0.5
            if jtrim < 0:
                raise ArithmeticError("cannot refine polynomial root (multiple roots?)")
            # Trimming only makes sense if there is a zero in the
            # remaining part
            if not poly(root).contains_zero():
                raise ArithmeticError("cannot refine polynomial root (precision too low?)")
        else:
            # Use interval Newton-Raphson to refine the root
            center = field(root.center())
            centerval = poly(center)
            nroot = center - centerval / slope

            IF DEBUG_REFINE_ROOT:
                print("nroot = %s" % intervalstr(nroot))

            if converging:
                # If we are converging, make sure we never enlarge our
                # interval again. This is in particular needed when the
                # interval gets small to avoid oscillations.
                nroot = nroot.intersection(root)
                if (nroot.diameter() << 2) >= root.diameter():
                    # If the new diameter is within a factor 4 of the
                    # original diameter, then we have converged reasonably
                    # well.
                    return nroot
            elif nroot in root:
                # We are converging => assume that we keep converging.
                IF DEBUG_REFINE_ROOT:
                    print("start converging")
                converging = True
            else:
                # It could be that the original interval is simply too
                # small for convergence. Try a new Newton-Raphson with
                # a starting interval 3 times as large.
                root3 = root + root - root
                slope = deriv(root3)
                if not slope.contains_zero():
                    nroot3 = (center - centerval / slope)
                    if nroot3 in root3:
                        IF DEBUG_REFINE_ROOT:
                            print("start converging (triple size)")
                        converging = True
                        root = nroot3
                        continue

                # Try to take the intersection of the old and new
                # intervals which is likely a better approximation
                # than just the new interval. In particular, if our
                # old interval does contain a zero, the new interval
                # will also contain that zero.
                #
                # If the intersection is either empty or the
                # intersection is the old interval, we don't do this.
                try:
                    if root not in nroot:
                        nroot = nroot.intersection(root)
                except ValueError:
                    pass
#            elif 2*i >= steps:
#                if root in nroot:
#                    # We are enlarging our interval => ignore the new
#                    # interval and try trimming instead.
#                    IF DEBUG_REFINE_ROOT:
#                        print("trim forced")
#                    do_trim = True
#                    continue
#                else:
#                    # Try to take the intersection of the old and new
#                    # intervals which is likely a better approximation than
#                    # just the new interval. In particular, if our old
#                    # interval does contain a zero, the new interval will
#                    # also contain that zero.
#                    try:
#                        nroot = nroot.intersection(root)
#                    except ValueError:
#                        # If the new interval still isn't contained in
#                        # the old, try tripling the size of the region.
#                        # This gives an interval with the same center
#                        # but with a larger diameter.
#                        nroot += nroot - nroot

            root = nroot

    # If we are converging, return the last interval
    if converging:
        return root
    raise ArithmeticError("cannot refine polynomial root (not enough steps?)")
