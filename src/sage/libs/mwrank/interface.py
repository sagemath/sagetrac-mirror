r"""
Sage interface to Cremona's eclib library (also known as mwrank)

This is the Sage interface to John Cremona's ``eclib`` C++ library for
arithmetic on elliptic curves.  The classes defined in this module
give Sage interpreter-level access to some of the functionality of
``eclib``.  For most purposes, it is not necessary to directly use these
classes. Instead, one can create an
:class:`EllipticCurve <sage.schemes.elliptic_curves.constructor.EllipticCurve>`
and call methods that are implemented using this module.

.. note::

   This interface is a direct library-level interface to ``eclib``,
   including the 2-descent program ``mwrank``.
"""

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import IntegerRing

# Need to permit tabs in order to doctest verbose output.
"""
SAGE_DOCTEST_ALLOW_TABS
"""

def get_precision():
    r"""
    Return the global NTL real number precision.

    See also :meth:`set_precision`.

    .. warning::

       The internal precision is binary.  This function multiplies the
       binary precision by 0.3 (`=\log_2(10)` approximately) and
       truncates.

    OUTPUT:

    (int) The current decimal precision.

    EXAMPLES::

        sage: mwrank_get_precision()
        50
    """
    # don't want to load mwrank every time Sage starts up, so we do
    # the import here.
    from sage.libs.mwrank.mwrank import get_precision
    return get_precision()

def set_precision(n):
    r"""
    Set the global NTL real number precision.  This has a massive
    effect on the speed of mwrank calculations.  The default (used if
    this function is not called) is ``n=50``, but it might have to be
    increased if a computation fails.  See also :meth:`get_precision`.

    INPUT:

    - ``n`` (long) -- real precision used for floating point
      computations in the library, in decimal digits.

    .. warning::

       This change is global and affects *all* future calls of eclib
       functions by Sage.

    EXAMPLES::

        sage: mwrank_set_precision(20)
    """
    # don't want to load mwrank every time Sage starts up, so we do
    # the import here.
    from sage.libs.mwrank.mwrank import set_precision
    set_precision(n)

class mwrank_EllipticCurve(SageObject):
    r"""
    The :class:`mwrank_EllipticCurve` class represents an elliptic
    curve using the ``Curvedata`` class from ``eclib``, called here an 'mwrank
    elliptic curve'.

    Create the mwrank elliptic curve with invariants
    ``ainvs``, which is a list of 5 or less *integers* `a_1`,
    `a_2`, `a_3`, `a_4`, and `a_5`.

    If strictly less than 5 invariants are given, then the *first*
    ones are set to 0, so, e.g., ``[3,4]`` means `a_1=a_2=a_3=0` and
    `a_4=3`, `a_5=4`.

    INPUT:

    - ``ainvs`` (list or tuple) -- a list of 5 or less integers, the
      coefficients of a nonsingular Weierstrass equation.

    - ``verbose`` (bool, default ``False``) -- verbosity flag.  If ``True``,
      then all Selmer group computations will be verbose.

    EXAMPLES:

    We create the elliptic curve `y^2 + y = x^3 + x^2 - 2x`::

        sage: e = mwrank_EllipticCurve([0, 1, 1, -2, 0])
        sage: e.ainvs()
        [0, 1, 1, -2, 0]

    This example illustrates that omitted `a`-invariants default to `0`::

        sage: e = mwrank_EllipticCurve([3, -4])
        sage: e
        y^2 = x^3 + 3*x - 4
        sage: e.ainvs()
        [0, 0, 0, 3, -4]

    The entries of the input list are coerced to :class:`int`.
    If this is impossible, then an error is raised::

        sage: e = mwrank_EllipticCurve([3, -4.8]); e
        Traceback (most recent call last):
        ...
        TypeError: ainvs must be a list or tuple of integers.

    When you enter a singular model you get an exception::

        sage: e = mwrank_EllipticCurve([0, 0])
        Traceback (most recent call last):
        ...
        ArithmeticError: Invariants (= [0, 0, 0, 0, 0]) do not describe an elliptic curve.
    """

    def __init__(self, ainvs, verbose=False):
        r"""
        Create the mwrank elliptic curve with invariants
        ``ainvs``, which is a list of 5 or less *integers* `a_1`,
        `a_2`, `a_3`, `a_4`, and `a_5`.

        See the docstring of this class for full documentation.

        EXAMPLES:

        We create the elliptic curve `y^2 + y = x^3 + x^2 - 2x`::

            sage: e = mwrank_EllipticCurve([0, 1, 1, -2, 0])
            sage: e.ainvs()
            [0, 1, 1, -2, 0]
        """
        # import here to save time during startup (mwrank takes a while to init)

        from sage.libs.mwrank.mwrank import _Curvedata

        # if not isinstance(ainvs, list) and len(ainvs) <= 5:
        if not isinstance(ainvs, (list,tuple)) or not len(ainvs) <= 5:
            raise TypeError("ainvs must be a list or tuple of length at most 5.")

        # Pad ainvs on the beginning by 0's, so e.g.
        # [a4,a5] works.
        ainvs = [0]*(5-len(ainvs)) + ainvs

        # Convert each entry to an int
        try:
            a_int = [IntegerRing()(x) for x in ainvs]
        except (TypeError, ValueError):
            raise TypeError("ainvs must be a list or tuple of integers.")
        self.__ainvs = a_int
        self.__curve = _Curvedata(a_int[0], a_int[1], a_int[2],
                                  a_int[3], a_int[4])

        if verbose:
            self.__verbose = True
        else:
            self.__verbose = False

        # place holders
        self.__saturate = -2  # not yet saturated

    def __reduce__(self):
        r"""
        Standard Python function used in pickling.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.__reduce__()
            (<class 'sage.libs.mwrank.interface.mwrank_EllipticCurve'>, ([0, 0, 1, -7, 6], False))


        """
        return mwrank_EllipticCurve, (self.__ainvs, self.__verbose)


    def set_verbose(self, verbose):
        """
        Set the verbosity of printing of output by the :meth:`two_descent()` and
        other functions.

        INPUT:

        - ``verbose`` (int) -- if positive, print lots of output when
          doing 2-descent.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.saturate() # no output
            sage: E.gens()
            [[0, -1, 1]]

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.set_verbose(1)
            sage: E.saturate() # tol 1e-14
            Basic pair: I=48, J=-432
            disc=255744
            2-adic index bound = 2
            By Lemma 5.1(a), 2-adic index = 1
            2-adic index = 1
            One (I,J) pair
            Looking for quartics with I = 48, J = -432
            Looking for Type 2 quartics:
            Trying positive a from 1 up to 1 (square a first...)
            (1,0,-6,4,1)        --trivial
            Trying positive a from 1 up to 1 (...then non-square a)
            Finished looking for Type 2 quartics.
            Looking for Type 1 quartics:
            Trying positive a from 1 up to 2 (square a first...)
            (1,0,0,4,4) --nontrivial...(x:y:z) = (1 : 1 : 0)
            Point = [0:0:1]
                height = 0.0511114082399688402358
            Rank of B=im(eps) increases to 1 (The previous point is on the egg)
            Exiting search for Type 1 quartics after finding one which is globally soluble.
            Mordell rank contribution from B=im(eps) = 1
            Selmer  rank contribution from B=im(eps) = 1
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            Searching for points (bound = 8)...done:
              found points which generate a subgroup of rank 1
              and regulator 0.0511114082399688402358
            Processing points found during 2-descent...done:
              now regulator = 0.0511114082399688402358
            Saturating (with bound = -1)...done:
              points were already saturated.
        """
        self.__verbose = verbose


    def _curve_data(self):
        r"""
        Returns the underlying :class:`_Curvedata` class for this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-1,0])
            sage: E._curve_data()
            [0,0,1,-1,0]
            b2 = 0	 b4 = -2	 b6 = 1	 b8 = -1
            c4 = 48		c6 = -216
            disc = 37	(# real components = 2)
            #torsion not yet computed
        """
        return self.__curve

    def ainvs(self):
        r"""
        Returns the `a`-invariants of this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-1,0])
            sage: E.ainvs()
            [0, 0, 1, -1, 0]
        """
        return self.__ainvs

    def isogeny_class(self, verbose=False):
        r"""
        Returns the isogeny class of this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.isogeny_class()
            ([[0, -1, 1, 0, 0], [0, -1, 1, -10, -20], [0, -1, 1, -7820, -263580]], [[0, 5, 0], [5, 0, 5], [0, 5, 0]])
        """
        return self.__curve.isogeny_class(verbose)

    def __repr__(self):
        r"""
        Returns the string representation of this mwrank elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.__repr__()
            'y^2+ y = x^3 - x^2 '
        """
        # TODO: Is the use (or omission) of spaces here intentional?
        a = self.ainvs()
        s = "y^2"
        if a[0] == -1:
            s += "- x*y "
        elif a[0] == 1:
            s += "+ x*y "
        elif a[0] != 0:
            s += "+ %s*x*y "%a[0]
        if a[2] == -1:
            s += " - y"
        elif a[2] == 1:
            s += "+ y"
        elif a[2] != 0:
            s += "+ %s*y"%a[2]
        s += " = x^3 "
        if a[1] == -1:
            s += "- x^2 "
        elif a[1] == 1:
            s += "+ x^2 "
        elif a[1] != 0:
            s += "+ %s*x^2 "%a[1]
        if a[3] == -1:
            s += "- x "
        elif a[3] == 1:
            s += "+ x "
        elif a[3] != 0:
            s += "+ %s*x "%a[3]
        if a[4] == -1:
            s += "-1"
        elif a[4] == 1:
            s += "+1"
        elif a[4] != 0:
            s += "+ %s"%a[4]
        s = s.replace("+ -","- ")
        return s


    def two_descent(self,
                    verbose = True,
                    selmer_only = False,
                    first_limit = 20,
                    second_limit = 8,
                    n_aux = -1,
                    second_descent = True):
        """
        Compute 2-descent data for this curve.

        INPUT:

        - ``verbose`` (bool, default ``True``) --  print what mwrank is doing.

        - ``selmer_only`` (bool, default ``False``) -- ``selmer_only`` switch.

        - ``first_limit`` (int, default 20) -- bound on `|x|+|z|` in
          quartic point search.

        - ``second_limit`` (int, default 8) -- bound on
          `\log \max(|x|,|z|)`, i.e. logarithmic.

        - ``n_aux`` (int, default -1) -- (only relevant for general
          2-descent when 2-torsion trivial) number of primes used for
          quartic search.  ``n_aux=-1`` causes default (8) to be used.
          Increase for curves of higher rank.

        - ``second_descent`` (bool, default ``True``) -- (only relevant
          for curves with 2-torsion, where mwrank uses descent via
          2-isogeny) flag determining whether or not to do second
          descent.  *Default strongly recommended.*


        OUTPUT:

        Nothing -- nothing is returned.

        TESTS:

        See :trac:`7992`::

            sage: EllipticCurve([0, prod(prime_range(10))]).mwrank_curve().two_descent()
            Basic pair: I=0, J=-5670
            disc=-32148900
            2-adic index bound = 2
            2-adic index = 2
            Two (I,J) pairs
            Looking for quartics with I = 0, J = -5670
            Looking for Type 3 quartics:
            Trying positive a from 1 up to 5 (square a first...)
            Trying positive a from 1 up to 5 (...then non-square a)
            (2,0,-12,19,-6)	--nontrivial...(x:y:z) = (2 : 4 : 1)
            Point = [-2488:-4997:512]
                height = 6.46767239...
            Rank of B=im(eps) increases to 1
            Trying negative a from -1 down to -3
            Finished looking for Type 3 quartics.
            Looking for quartics with I = 0, J = -362880
            Looking for Type 3 quartics:
            Trying positive a from 1 up to 20 (square a first...)
            Trying positive a from 1 up to 20 (...then non-square a)
            Trying negative a from -1 down to -13
            Finished looking for Type 3 quartics.
            Mordell rank contribution from B=im(eps) = 1
            Selmer  rank contribution from B=im(eps) = 1
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: EllipticCurve([0, prod(prime_range(100))]).mwrank_curve().two_descent()
            Traceback (most recent call last):
            ...
            RuntimeError: Aborted

        Calling this method twice does not cause a segmentation fault
        (see :trac:`10665`)::

            sage: E = EllipticCurve([1, 1, 0, 0, 528])
            sage: E.two_descent(verbose=False)
            True
            sage: E.two_descent(verbose=False)
            True

        """
        from sage.libs.mwrank.mwrank import _two_descent # import here to save time
        first_limit = int(first_limit)
        second_limit = int(second_limit)
        n_aux = int(n_aux)
        second_descent = int(second_descent)	# convert from bool to (int) 0 or 1
        # TODO:  Don't allow limits above some value...???
        #   (since otherwise mwrank just sets limit tiny)
        self.__descent = _two_descent()
        self.__descent.do_descent(self.__curve,
                                      verbose,
                                      selmer_only,
                                      first_limit,
                                      second_limit,
                                      n_aux,
                                      second_descent)
        if not self.__descent.ok():
            raise RuntimeError("A 2-descent did not complete successfully.")
        self.__saturate = -2  # not yet saturated

    def __two_descent_data(self):
        r"""
        Returns the 2-descent data for this elliptic curve.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E._mwrank_EllipticCurve__two_descent_data()
            <sage.libs.mwrank.mwrank._two_descent object at ...>
        """
        try:
            return self.__descent
        except AttributeError:
            self.two_descent(self.__verbose)
            return self.__descent

    def conductor(self):
        """
        Return the conductor of this curve, computed using Cremona's
        implementation of Tate's algorithm.

        .. note::

           This is independent of PARI's.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([1, 1, 0, -6958, -224588])
            sage: E.conductor()
            2310
        """
        return self.__curve.conductor()

    def rank(self):
        """
        Returns the rank of this curve, computed using :meth:`two_descent()`.

        In general this may only be a lower bound for the rank; an
        upper bound may be obtained using the function :meth:`rank_bound()`.
        To test whether the value has been proved to be correct, use
        the method :meth:`certain()`.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.rank()
            0
            sage: E.certain()
            True

        ::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.rank()
            0
            sage: E.certain()
            False

        """
        return self.__two_descent_data().getrank()

    def rank_bound(self):
        """
        Returns an upper bound for the rank of this curve, computed
        using :meth:`two_descent()`.

        If the curve has no 2-torsion, this is equal to the 2-Selmer
        rank.  If the curve has 2-torsion, the upper bound may be
        smaller than the bound obtained from the 2-Selmer rank minus
        the 2-rank of the torsion, since more information is gained
        from the 2-isogenous curve or curves.

        EXAMPLES:

        The following is the curve 960D1, which has rank 0,
        but Sha of order 4::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.rank_bound()
            0
            sage: E.rank()
            0

        In this case the rank was computed using a second descent,
        which is able to determine (by considering a 2-isogenous
        curve) that Sha is nontrivial.  If we deliberately stop the
        second descent, the rank bound is larger::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.two_descent(second_descent = False, verbose=False)
            sage: E.rank_bound()
            2

        In contrast, for the curve 571A, also with rank 0 and Sha
        of order 4, we only obtain an upper bound of 2::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.rank_bound()
            2

        In this case the value returned by :meth:`rank()` is only a
        lower bound in general (though this is correct)::

            sage: E.rank()
            0
            sage: E.certain()
            False

        """
        return self.__two_descent_data().getrankbound()

    def selmer_rank(self):
        r"""
        Returns the rank of the 2-Selmer group of the curve.

        EXAMPLES:

        The following is the curve 960D1, which has rank 0, but Sha of
        order 4.  The 2-torsion has rank 2, and the Selmer rank is 3::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.selmer_rank()
            3

        Nevertheless, we can obtain a tight upper bound on the rank
        since a second descent is performed which establishes the
        2-rank of Sha::

            sage: E.rank_bound()
            0

        To show that this was resolved using a second descent, we do
        the computation again but turn off ``second_descent``::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.two_descent(second_descent = False, verbose=False)
            sage: E.rank_bound()
            2

        For the curve 571A, also with rank 0 and Sha of order 4,
        but with no 2-torsion, the Selmer rank is strictly greater
        than the rank::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.selmer_rank()
            2
            sage: E.rank_bound()
            2

        In cases like this with no 2-torsion, the rank upper bound is
        always equal to the 2-Selmer rank.  If we ask for the rank,
        all we get is a lower bound::

            sage: E.rank()
            0
            sage: E.certain()
            False

        """
        return self.__two_descent_data().getselmer()

    def regulator(self):
        r"""
        Return the regulator of the saturated Mordell-Weil group.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.regulator()
            0.05111140823996884
        """
        self.saturate()
        if not self.certain():
            raise RuntimeError("Unable to saturate Mordell-Weil group.")
        R = self.__two_descent_data().regulator()
        return float(R)

    def saturate(self, bound=-1):
        """
        Compute the saturation of the Mordell-Weil group at all
        primes up to ``bound``.

        INPUT:

        - ``bound`` (int, default -1) -- Use `-1` (the default) to
          saturate at *all* primes, `0` for no saturation, or `n` (a
          positive integer) to saturate at all primes up to `n`.

        EXAMPLES:

        Since the 2-descent automatically saturates at primes up to
        20, it is not easy to come up with an example where saturation
        has any effect::

            sage: E = mwrank_EllipticCurve([0, 0, 0, -1002231243161, 0])
            sage: E.gens()
            [[-1001107, -4004428, 1]]
            sage: E.saturate()
            sage: E.gens()
            [[-1001107, -4004428, 1]]

        Check that :trac:`18031` is fixed::

            sage: E = EllipticCurve([0,-1,1,-266,968])
            sage: Q1 = E([-1995,3674,125])
            sage: Q2 = E([157,1950,1])
            sage: E.saturation([Q1,Q2])
            ([(1 : -27 : 1), (157 : 1950 : 1)], 3, 0.801588644684981)
        """
        bound = int(bound)
        if self.__saturate < bound:
            self.__two_descent_data().saturate(bound)
            self.__saturate = bound

    def gens(self):
        """
        Return a list of the generators for the Mordell-Weil group.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.gens()
            [[0, -1, 1]]
        """
        self.saturate()
        from sage.rings.all import Integer
        L = eval(self.__two_descent_data().getbasis().replace(":",","))
        return [[Integer(x), Integer(y), Integer(z)] for (x,y,z) in L]

    def certain(self):
        r"""
        Returns ``True`` if the last :meth:`two_descent()` call provably correctly
        computed the rank.  If :meth:`two_descent()` hasn't been
        called, then it is first called by :meth:`certain()`
        using the default parameters.

        The result is ``True`` if and only if the results of the methods
        :meth:`rank()` and :meth:`rank_bound()` are equal.

        EXAMPLES:

        A 2-descent does not determine `E(\QQ)` with certainty
        for the curve `y^2 + y = x^3 - x^2 - 120x - 2183`::

            sage: E = mwrank_EllipticCurve([0, -1, 1, -120, -2183])
            sage: E.two_descent(False)
            ...
            sage: E.certain()
            False
            sage: E.rank()
            0

        The previous value is only a lower bound; the upper bound is greater::

            sage: E.rank_bound()
            2

        In fact the rank of `E` is actually 0 (as one could see by
        computing the `L`-function), but Sha has order 4 and the
        2-torsion is trivial, so mwrank cannot conclusively
        determine the rank in this case.
        """
        return bool(self.__two_descent_data().getcertain())

    #def fullmw(self):
    #    return self.__two_descent_data().getfullmw()

    def CPS_height_bound(self):
        r"""
        Return the Cremona-Prickett-Siksek height bound.  This is a
        floating point number `B` such that if `P` is a point on the
        curve, then the naive logarithmic height `h(P)` is less than
        `B+\hat{h}(P)`, where `\hat{h}(P)` is the canonical height of
        `P`.

        .. warning::

           We assume the model is minimal!

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 0, -1002231243161, 0])
            sage: E.CPS_height_bound()
            14.163198527061496
            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.CPS_height_bound()
            0.0
        """
        return self.__curve.cps_bound()

    def silverman_bound(self):
        r"""
        Return the Silverman height bound.  This is a floating point
        number `B` such that if `P` is a point on the curve, then the
        naive logarithmic height `h(P)` is less than `B+\hat{h}(P)`,
        where `\hat{h}(P)` is the canonical height of `P`.

        .. warning::

           We assume the model is minimal!

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0, 0, 0, -1002231243161, 0])
            sage: E.silverman_bound()
            18.29545210468247
            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.silverman_bound()
            6.284833369972403
        """
        return self.__curve.silverman_bound()


class mwrank_MordellWeil(SageObject):
    r"""
    The :class:`mwrank_MordellWeil` class represents a subgroup of a
    Mordell-Weil group.  Use this class to saturate a specified list
    of points on an :class:`mwrank_EllipticCurve`, or to search for
    points up to some bound.

    INPUT:

    - ``curve`` (:class:`mwrank_EllipticCurve`) -- the underlying
      elliptic curve.

    - ``verbose`` (bool, default ``False``) -- verbosity flag (controls
      amount of output produced in point searches).

    - ``pp`` (int, default 1) -- process points flag (if nonzero,
      the points found are processed, so that at all times only a
      `\ZZ`-basis for the subgroup generated by the points found
      so far is stored; if zero, no processing is done and all
      points found are stored).

    - ``maxr`` (int, default 999) -- maximum rank (quit point
      searching once the points found generate a subgroup of this
      rank; useful if an upper bound for the rank is already
      known).

    EXAMPLE::

        sage: E = mwrank_EllipticCurve([1,0,1,4,-6])
        sage: EQ = mwrank_MordellWeil(E)
        sage: EQ
        Subgroup of Mordell-Weil group: []
        sage: EQ.search(2)
        P1 = [0:1:0]	 is torsion point, order 1
        P1 = [1:-1:1]	 is torsion point, order 2
        P1 = [2:2:1]	 is torsion point, order 3
        P1 = [9:23:1]	 is torsion point, order 6

        sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
        sage: EQ = mwrank_MordellWeil(E)
        sage: EQ.search(2)
        P1 = [0:1:0]	 is torsion point, order 1
        P1 = [-3:0:1]	  is generator number 1
        ...
        P4 = [-91:804:343]	 = -2*P1 + 2*P2 + 1*P3 (mod torsion)
        sage: EQ
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

    Example to illustrate the verbose parameter::

        sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
        sage: EQ = mwrank_MordellWeil(E, verbose=False)
        sage: EQ.search(1)
        sage: EQ
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

        sage: EQ = mwrank_MordellWeil(E, verbose=True)
        sage: EQ.search(1)
        P1 = [0:1:0]	 is torsion point, order 1
        P1 = [-3:0:1]	  is generator number 1
        saturating up to 20...Checking 2-saturation
        Points have successfully been 2-saturated (max q used = 7)
        Checking 3-saturation
        Points have successfully been 3-saturated (max q used = 7)
        Checking 5-saturation
        Points have successfully been 5-saturated (max q used = 23)
        Checking 7-saturation
        Points have successfully been 7-saturated (max q used = 41)
        Checking 11-saturation
        Points have successfully been 11-saturated (max q used = 17)
        Checking 13-saturation
        Points have successfully been 13-saturated (max q used = 43)
        Checking 17-saturation
        Points have successfully been 17-saturated (max q used = 31)
        Checking 19-saturation
        Points have successfully been 19-saturated (max q used = 37)
        done
        P2 = [-2:3:1]	  is generator number 2
        saturating up to 20...Checking 2-saturation
        possible kernel vector = [1,1]
        This point may be in 2E(Q): [14:-52:1]
        ...and it is!
        Replacing old generator #1 with new generator [1:-1:1]
        Points have successfully been 2-saturated (max q used = 7)
        Index gain = 2^1
        Checking 3-saturation
        Points have successfully been 3-saturated (max q used = 13)
        Checking 5-saturation
        Points have successfully been 5-saturated (max q used = 67)
        Checking 7-saturation
        Points have successfully been 7-saturated (max q used = 53)
        Checking 11-saturation
        Points have successfully been 11-saturated (max q used = 73)
        Checking 13-saturation
        Points have successfully been 13-saturated (max q used = 103)
        Checking 17-saturation
        Points have successfully been 17-saturated (max q used = 113)
        Checking 19-saturation
        Points have successfully been 19-saturated (max q used = 47)
        done (index = 2).
        Gained index 2, new generators = [ [1:-1:1] [-2:3:1] ]
        P3 = [-14:25:8]	  is generator number 3
        saturating up to 20...Checking 2-saturation
        Points have successfully been 2-saturated (max q used = 11)
        Checking 3-saturation
        Points have successfully been 3-saturated (max q used = 13)
        Checking 5-saturation
        Points have successfully been 5-saturated (max q used = 71)
        Checking 7-saturation
        Points have successfully been 7-saturated (max q used = 101)
        Checking 11-saturation
        Points have successfully been 11-saturated (max q used = 127)
        Checking 13-saturation
        Points have successfully been 13-saturated (max q used = 151)
        Checking 17-saturation
        Points have successfully been 17-saturated (max q used = 139)
        Checking 19-saturation
        Points have successfully been 19-saturated (max q used = 179)
        done (index = 1).
        P4 = [-1:3:1]	 = -1*P1 + -1*P2 + -1*P3 (mod torsion)
        P4 = [0:2:1]	 = 2*P1 + 0*P2 + 1*P3 (mod torsion)
        P4 = [2:13:8]	 = -3*P1 + 1*P2 + -1*P3 (mod torsion)
        P4 = [1:0:1]	 = -1*P1 + 0*P2 + 0*P3 (mod torsion)
        P4 = [2:0:1]	 = -1*P1 + 1*P2 + 0*P3 (mod torsion)
        P4 = [18:7:8]	 = -2*P1 + -1*P2 + -1*P3 (mod torsion)
        P4 = [3:3:1]	 = 1*P1 + 0*P2 + 1*P3 (mod torsion)
        P4 = [4:6:1]	 = 0*P1 + -1*P2 + -1*P3 (mod torsion)
        P4 = [36:69:64]	 = 1*P1 + -2*P2 + 0*P3 (mod torsion)
        P4 = [68:-25:64]	 = -2*P1 + -1*P2 + -2*P3 (mod torsion)
        P4 = [12:35:27]	 = 1*P1 + -1*P2 + -1*P3 (mod torsion)
        sage: EQ
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

    Example to illustrate the process points (``pp``) parameter::

        sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
        sage: EQ = mwrank_MordellWeil(E, verbose=False, pp=1)
        sage: EQ.search(1); EQ # generators only
        Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]
        sage: EQ = mwrank_MordellWeil(E, verbose=False, pp=0)
        sage: EQ.search(1); EQ # all points found
        Subgroup of Mordell-Weil group: [[-3:0:1], [-2:3:1], [-14:25:8], [-1:3:1], [0:2:1], [2:13:8], [1:0:1], [2:0:1], [18:7:8], [3:3:1], [4:6:1], [36:69:64], [68:-25:64], [12:35:27]]
    """

    def __init__(self, curve, verbose=True, pp=1, maxr=999):
        r"""
        Constructor for the :class:`mwrank_MordellWeil` class.

        See the docstring of this class for full documentation.

        EXAMPLE::

            sage: E = mwrank_EllipticCurve([1,0,1,4,-6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ
            Subgroup of Mordell-Weil group: []
        """
        if not isinstance(curve, mwrank_EllipticCurve):
            raise TypeError("curve (=%s) must be an mwrank_EllipticCurve"%curve)
        self.__curve = curve
        self.__verbose = verbose
        self.__pp = pp
        self.__maxr = maxr
        if verbose:
            verb = 1
        else:
            verb = 0
        from sage.libs.mwrank.mwrank import _mw # import here to save time
        self.__mw = _mw(curve._curve_data(), verb, pp, maxr)

    def __reduce__(self):
        r"""
        Standard Python function used in pickling.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.__reduce__()
            (<class 'sage.libs.mwrank.interface.mwrank_MordellWeil'>, (y^2+ y = x^3 - 7*x + 6, True, 1, 999))
        """
        return mwrank_MordellWeil, (self.__curve, self.__verbose, self.__pp, self.__maxr)

    def __repr__(self):
        r"""
        String representation of this Mordell-Weil subgroup.

        OUTPUT:

        (string) String representation of this Mordell-Weil subgroup.

        EXAMPLE::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E, verbose=False)
            sage: EQ.__repr__()
            'Subgroup of Mordell-Weil group: []'
            sage: EQ.search(1)
            sage: EQ.__repr__()
            'Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]'
        """
        # FIXED: Changing "Mordell Weil" to "Mordell-Weil" did break doctests...
        return "Subgroup of Mordell-Weil group: %s"%self.__mw

    def process(self, v, sat=0):
        """
        This function allows one to add points to a :class:`mwrank_MordellWeil` object.

        Process points in the list ``v``, with saturation at primes up to
        ``sat``.  If ``sat`` is zero (the default), do no saturation.

        INPUT:

        - ``v`` (list of 3-tuples or lists of ints or Integers) -- a
          list of triples of integers, which define points on the
          curve.

        - ``sat`` (int, default 0) -- saturate at primes up to ``sat``, or at
          *all* primes if ``sat`` is zero.

        OUTPUT:

        None.  But note that if the ``verbose`` flag is set, then there
        will be some output as a side-effect.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.gens()
            [[1, -1, 1], [-2, 3, 1], [-14, 25, 8]]
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1, -1, 1], [-2, 3, 1], [-14, 25, 8]])
            P1 = [1:-1:1]	  is generator number 1
            P2 = [-2:3:1]	  is generator number 2
            P3 = [-14:25:8]	  is generator number 3

        ::

            sage: EQ.points()
            [[1, -1, 1], [-2, 3, 1], [-14, 25, 8]]

        Example to illustrate the saturation parameter ``sat``::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], sat=20)
            P1 = [1547:-2967:343]	  is generator number 1
            ...
            Gained index 5, new generators = [ [-2:3:1] [-14:25:8] [1:-1:1] ]

            sage: EQ.points()
            [[-2, 3, 1], [-14, 25, 8], [1, -1, 1]]

        Here the processing was followed by saturation at primes up to
        20.  Now we prevent this initial saturation::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], sat=0)
            P1 = [1547:-2967:343]	  is generator number 1
            P2 = [2707496766203306:864581029138191:2969715140223272]	  is generator number 2
            P3 = [-13422227300:-49322830557:12167000000]	  is generator number 3
            sage: EQ.points()
            [[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]]
            sage: EQ.regulator()
            375.42919921875
            sage: EQ.saturate(2)  # points were not 2-saturated
            saturating basis...Saturation index bound = 93
            WARNING: saturation at primes p > 2 will not be done;
            ...
            Gained index 2
            New regulator =  93.857300720636393209
            (False, 2, '[ ]')
            sage: EQ.points()
            [[-2, 3, 1], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]]
            sage: EQ.regulator()
            93.8572998046875
            sage: EQ.saturate(3)  # points were not 3-saturated
            saturating basis...Saturation index bound = 46
            WARNING: saturation at primes p > 3 will not be done;
            ...
            Gained index 3
            New regulator =  10.4285889689595992455
            (False, 3, '[ ]')
            sage: EQ.points()
            [[-2, 3, 1], [-14, 25, 8], [-13422227300, -49322830557, 12167000000]]
            sage: EQ.regulator()
            10.4285888671875
            sage: EQ.saturate(5)  # points were not 5-saturated
            saturating basis...Saturation index bound = 15
            WARNING: saturation at primes p > 5 will not be done;
            ...
            Gained index 5
            New regulator =  0.417143558758383969818
            (False, 5, '[ ]')
            sage: EQ.points()
            [[-2, 3, 1], [-14, 25, 8], [1, -1, 1]]
            sage: EQ.regulator()
            0.4171435534954071
            sage: EQ.saturate()   # points are now saturated
            saturating basis...Saturation index bound = 3
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            done
            (True, 1, '[ ]')
        """
        if not isinstance(v, list):
            raise TypeError("v (=%s) must be a list"%v)
        sat = int(sat)
        for P in v:
            if not isinstance(P, (list,tuple)) or len(P) != 3:
                raise TypeError("v (=%s) must be a list of 3-tuples (or 3-element lists) of ints"%v)
            self.__mw.process(P, sat)

    def regulator(self):
        """
        Return the regulator of the points in this subgroup of
        the Mordell-Weil group.

        .. note::

           ``eclib`` can compute the regulator to arbitrary precision,
           but the interface currently returns the output as a ``float``.

        OUTPUT:

        (float) The regulator of the points in this subgroup.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.regulator()
            1.0

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: E.regulator()
            0.417143558758384
        """
        return self.__mw.regulator()

    def rank(self):
        """
        Return the rank of this subgroup of the Mordell-Weil group.

        OUTPUT:

        (int) The rank of this subgroup of the Mordell-Weil group.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,-1,1,0,0])
            sage: E.rank()
            0

        A rank 3 example::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.rank()
            0
            sage: EQ.regulator()
            1.0

        The preceding output is correct, since we have not yet tried
        to find any points on the curve either by searching or
        2-descent::

            sage: EQ
            Subgroup of Mordell-Weil group: []

        Now we do a very small search::

            sage: EQ.search(1)
            P1 = [0:1:0]	 is torsion point, order 1
            P1 = [-3:0:1]	  is generator number 1
            saturating up to 20...Checking 2-saturation
            ...
            P4 = [12:35:27]	 = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ
            Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]
            sage: EQ.rank()
            3
            sage: EQ.regulator()
            0.4171435534954071

        We do in fact now have a full Mordell-Weil basis.

        """
        return self.__mw.rank()

    def saturate(self, max_prime=-1, odd_primes_only=False):
        r"""
        Saturate this subgroup of the Mordell-Weil group.

        INPUT:

        - ``max_prime`` (int, default -1) -- saturation is performed for
          all primes up to ``max_prime``. If `-1` (the default), an
          upper bound is computed for the primes at which the subgroup
          may not be saturated, and this is used; however, if the
          computed bound is greater than a value set by the ``eclib``
          library (currently 97) then no saturation will be attempted
          at primes above this.

        - ``odd_primes_only`` (bool, default ``False``) -- only do
          saturation at odd primes.  (If the points have been found
          via :meth:``two_descent()`` they should already be 2-saturated.)

        OUTPUT:

        (3-tuple) (``ok``, ``index``, ``unsatlist``) where:

        - ``ok`` (bool) -- ``True`` if and only if the saturation was
          provably successful at all primes attempted.  If the default
          was used for ``max_prime`` and no warning was output about
          the computed saturation bound being too high, then ``True``
          indicates that the subgroup is saturated at *all*
          primes.

        - ``index`` (int) -- the index of the group generated by the
          original points in their saturation.

        - ``unsatlist`` (list of ints) -- list of primes at which
          saturation could not be proved or achieved.  Increasing the
          decimal precision should correct this, since it happens when
          a linear combination of the points appears to be a multiple
          of `p` but cannot be divided by `p`.  (Note that ``eclib``
          uses floating point methods based on elliptic logarithms to
          divide points.)

        .. note::

           We emphasize that if this function returns ``True`` as the
           first return argument (``ok``), and if the default was used for the
           parameter ``max_prime``, then the points in the basis after
           calling this function are saturated at *all* primes,
           i.e., saturating at the primes up to ``max_prime`` are
           sufficient to saturate at all primes.  Note that the
           function might not have needed to saturate at all primes up
           to ``max_prime``.  It has worked out what prime you need to
           saturate up to, and that prime might be smaller than ``max_prime``.

        .. note::

           Currently (May 2010), this does not remember the result of
           calling :meth:`search()`.  So calling :meth:`search()` up
           to height 20 then calling :meth:`saturate()` results in
           another search up to height 18.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)

        We initialise with three points which happen to be 2, 3 and 5
        times the generators of this rank 3 curve.  To prevent
        automatic saturation at this stage we set the parameter
        ``sat`` to 0 (which is in fact the default)::

            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], sat=0)
            P1 = [1547:-2967:343]	  is generator number 1
            P2 = [2707496766203306:864581029138191:2969715140223272]	  is generator number 2
            P3 = [-13422227300:-49322830557:12167000000]	  is generator number 3
            sage: EQ
            Subgroup of Mordell-Weil group: [[1547:-2967:343], [2707496766203306:864581029138191:2969715140223272], [-13422227300:-49322830557:12167000000]]
            sage: EQ.regulator()
            375.42919921875

        Now we saturate at `p=2`, and gain index 2::

            sage: EQ.saturate(2)  # points were not 2-saturated
            saturating basis...Saturation index bound = 93
            WARNING: saturation at primes p > 2 will not be done;
            ...
            Gained index 2
            New regulator =  93.857300720636393209
            (False, 2, '[ ]')
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [2707496766203306:864581029138191:2969715140223272], [-13422227300:-49322830557:12167000000]]
            sage: EQ.regulator()
            93.8572998046875

        Now we saturate at `p=3`, and gain index 3::

            sage: EQ.saturate(3)  # points were not 3-saturated
            saturating basis...Saturation index bound = 46
            WARNING: saturation at primes p > 3 will not be done;
            ...
            Gained index 3
            New regulator =  10.4285889689595992455
            (False, 3, '[ ]')
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [-14:25:8], [-13422227300:-49322830557:12167000000]]
            sage: EQ.regulator()
            10.4285888671875

        Now we saturate at `p=5`, and gain index 5::

            sage: EQ.saturate(5)  # points were not 5-saturated
            saturating basis...Saturation index bound = 15
            WARNING: saturation at primes p > 5 will not be done;
            ...
            Gained index 5
            New regulator =  0.417143558758383969818
            (False, 5, '[ ]')
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [-14:25:8], [1:-1:1]]
            sage: EQ.regulator()
            0.4171435534954071

        Finally we finish the saturation.  The output here shows that
        the points are now provably saturated at all primes::

            sage: EQ.saturate()   # points are now saturated
            saturating basis...Saturation index bound = 3
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            done
            (True, 1, '[ ]')

        Of course, the :meth:`process()` function would have done all this
        automatically for us::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.process([[1547, -2967, 343], [2707496766203306, 864581029138191, 2969715140223272], [-13422227300, -49322830557, 12167000000]], sat=5)
            P1 = [1547:-2967:343]	  is generator number 1
            ...
            Gained index 5, new generators = [ [-2:3:1] [-14:25:8] [1:-1:1] ]
            sage: EQ
            Subgroup of Mordell-Weil group: [[-2:3:1], [-14:25:8], [1:-1:1]]
            sage: EQ.regulator()
            0.4171435534954071

        But we would still need to use the :meth:`saturate()` function to
        verify that full saturation has been done::

            sage: EQ.saturate()
            saturating basis...Saturation index bound = 3
            Checking saturation at [ 2 3 ]
            Checking 2-saturation
            Points were proved 2-saturated (max q used = 11)
            Checking 3-saturation
            Points were proved 3-saturated (max q used = 13)
            done
            (True, 1, '[ ]')

        Note the output of the preceding command: it proves that the
        index of the points in their saturation is at most 3, then
        proves saturation at 2 and at 3, by reducing the points modulo
        all primes of good reduction up to 11, respectively 13.
        """
        ok, index, unsat = self.__mw.saturate(int(max_prime), odd_primes_only)
        return bool(ok), int(str(index)), unsat

    def search(self, height_limit=18, verbose=False):
        r"""
        Search for new points, and add them to this subgroup of the
        Mordell-Weil group.

        INPUT:

        - ``height_limit`` (float, default: 18) -- search up to this
          logarithmic height.

        .. note::

          On 32-bit machines, this *must* be < 21.48 else
          `\exp(h_{\text{lim}}) > 2^{31}` and overflows.  On 64-bit machines, it
          must be *at most* 43.668.  However, this bound is a logarithmic
          bound and increasing it by just 1 increases the running time
          by (roughly) `\exp(1.5)=4.5`, so searching up to even 20
          takes a very long time.

        .. note::

           The search is carried out with a quadratic sieve, using
           code adapted from a version of Michael Stoll's
           ``ratpoints`` program.  It would be preferable to use a
           newer version of ``ratpoints``.

        - ``verbose`` (bool, default ``False``) -- turn verbose operation on
          or off.

        EXAMPLES:

        A rank 3 example, where a very small search is sufficient to
        find a Mordell-Weil basis::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.search(1)
            P1 = [0:1:0]	 is torsion point, order 1
            P1 = [-3:0:1]	  is generator number 1
            ...
            P4 = [12:35:27]	 = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ
            Subgroup of Mordell-Weil group: [[1:-1:1], [-2:3:1], [-14:25:8]]

        In the next example, a search bound of 12 is needed to find a
        non-torsion point::

            sage: E = mwrank_EllipticCurve([0, -1, 0, -18392, -1186248]) #1056g4
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.search(11); EQ
            P1 = [0:1:0]	 is torsion point, order 1
            P1 = [161:0:1]	 is torsion point, order 2
            Subgroup of Mordell-Weil group: []
            sage: EQ.search(12); EQ
            P1 = [0:1:0]	 is torsion point, order 1
            P1 = [161:0:1]	 is torsion point, order 2
            P1 = [4413270:10381877:27000]	  is generator number 1
            ...
            Subgroup of Mordell-Weil group: [[4413270:10381877:27000]]
        """
        height_limit = float(height_limit)
        if height_limit >= 21.4:	# TODO: docstring says 21.48 (for 32-bit machines; what about 64-bit...?)
            raise ValueError("The height limit must be < 21.4.")

        moduli_option = 0  # Use Stoll's sieving program... see strategies in ratpoints-1.4.c

        ##            moduli_option -- int (default: 0); if > 0; a flag used to determine
        ##                    the moduli that are used in sieving
        ##                       1 -- first 10 odd primes; the first one you
        ##                            would think of.
        ##                       2 -- three composites; $2^6\cdot 3^4$, ... TODO
        ##                             (from German mathematician J. Gebel;
        ##                               personal conversation about SIMATH)
        ##                       3 -- nine prime powers; $2^5, \ldots$;
        ##                            like 1 but includes powers of small primes
        ##                    TODO: Extract the meaning from mwprocs.cc; line 776 etc.

        verbose = bool(verbose)
        self.__mw.search(height_limit, moduli_option, verbose)

    def points(self):
        """
        Return a list of the generating points in this Mordell-Weil
        group.

        OUTPUT:

        (list) A list of lists of length 3, each holding the
        primitive integer coordinates `[x,y,z]` of a generating
        point.

        EXAMPLES::

            sage: E = mwrank_EllipticCurve([0,0,1,-7,6])
            sage: EQ = mwrank_MordellWeil(E)
            sage: EQ.search(1)
            P1 = [0:1:0]	 is torsion point, order 1
            P1 = [-3:0:1]	  is generator number 1
            ...
            P4 = [12:35:27]	 = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ.points()
            [[1, -1, 1], [-2, 3, 1], [-14, 25, 8]]

        """
        L = eval(self.__mw.getbasis().replace(":",","))
        from sage.rings.all import Integer
        return [[Integer(x), Integer(y), Integer(z)] for (x,y,z) in L]
