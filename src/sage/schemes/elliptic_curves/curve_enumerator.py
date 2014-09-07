r"""
A class to enumerate all elliptic curves over `\QQ` up to isomorphism
in a given Weierstrass family, ordered by height.

AUTHORS:

- Simon Spicer (2012-09): First version
"""

#*****************************************************************************
#       Copyright (C) 2012 William Stein and Simon Spicer
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import time
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil
from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.misc import powerset, srange
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.arith import valuation


class CurveEnumerator_abstract(object):
    r"""
    The CurveEnumerator class will enumerate all elliptic curves
    over `\QQ` in a specified Weierstrass form ordered by height, where
    height is a function of the Weierstrass equation of the curve
    as described in each subclass's __init__() method.
    """

    def __init__(self):
        r"""
        This should never be called, as CurveEnumerater_abstract
        has been designed only to be inherited from its subclasses

        EXAMPLES::

        """
        raise NameError("Abstract class cannot be instantiated.")

    def __repr__(self):
        """
        Representation of self.

        Prints what family of curves is being considered, the model
        description and the height function on that family.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass"); C
            Height iterator for elliptic curves over QQ
            Family:             Short Weierstrass
            Model:              Y^2 = X^3 + A*X + B
            Coefficients:       [A,B]
            Height function:    H = min{|A|^3,|B|^2}

            sage: C = EllipticCurveEnumerator(family="F_2(2)"); C
            Height iterator for elliptic curves over QQ
            Family:             Rank Two + Two Torsion
            Model:              Y^2 = X^3 + A*X^2 + B, where
                                A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                                B = w0*w1*w2*w3
            Coefficients:       [w0,w1,w2,w3]
            Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
        """

        s = "Height iterator for elliptic curves over QQ\n"
        s += "Family:             "+self._name+"\n"
        s += "Model:              "+self._model+"\n"
        s += "Coefficients:       "+self._coeff_names+"\n"
        s += "Height function:    "+self._height_function
        return s

    def _height_increment(self, coeffs):
        r"""
        Given a tuple of coefficients of a Weierstrass equation with
        a certain height, return the next largest permissible height
        and the range of coefficients that achieve that height.

        INPUT:

        - ``coeffs`` -- A list or tuple of coefficients of the
          same length as the number of coeffients in the model.

        OUTPUT:

        A tuple of three entries consisting of:
        1. The smallest permissible height greater than the height
           of the input coefficient list;
        2. A list of coeffients, some of which attain the above
           height;
        3. A list of indices of which of the coefficients in the
           above list achieve this height. The remaining entries
           in the coefficient list indicate the maximum absolute
           value that coefficient can attain without affecting the
           curve's height.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: C._height_increment([2, 2])
            (9, [2, 3], [1])
            sage: C._height_increment([3,7])
            (64, [4, 8], [0, 1])
        """
        I = range(len(coeffs))
        height_candidates = [(coeffs[i]+1)**(self._pows[i]) for i in I]
        next_height = min(height_candidates)
        index = []

        new_coeffs = list(coeffs)
        for i in I:
            if height_candidates[i] == next_height:
                new_coeffs[i] += 1
                index.append(i)
        return (next_height, new_coeffs, index)

    def next_height(self, N):
        r"""
        Return the next permissible height greater than or equal to N for
        curves in self's family.

        .. WARNING::

            This function my return a height for which only singular
            curves exist. For example, in the short Weierstrass case
            height 0 is permissible, as the curve Y^2 = X^3 (uniquely)
            has height zero.

        INPUT:

        - ``N`` -- An integer. Note that N may be negative, even though height
          is always non-negative.

        OUTPUT:

        A tuple consisting of three elements of the form (H, C, I) such that
        1. H is the smallest height >= N (so if N<0, H will be 0)
        2. C is a list of coefficients for curves of this height
        3. I is list of indices indicating which of the above coefficients
           achieve this height. The remaining values in C  indicate the
           max absolute value those coefficients are allowed to obtain
           without altering the height.

           For example, the tuple (4, [1, 2], [1]) for the short Weierstrass
           case denotes set of curves with height 4; these are all of the
           form Y^2 = X^3 + A*X + B, where B=2 and A ranges between -1 and 1.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: C.next_height(4)
            (4, [1, 2], [1])
            sage: C.next_height(60)
            (64, [4, 8], [0, 1])
            sage: C.next_height(-100)
            (0, [0, 0], [0, 1])
        """

        # Negative heights don't exist
        if N < 0:
            N = ZZ(0)

        coeffs = [ceil(N**(1/n))-1 for n in self._pows]
        # height = max([coeffs[i]**(self._pows[i]) for i in range(self._num_coeffs)])
        return self._height_increment(coeffs)

    def heights(self, lowerbound, upperbound):
        r"""
        Return a list of permissible curve heights in the
        specified range (bounds inclusive), and for each height the
        equation coefficients that produce curves of that height.

        .. WARNING::

            This function my return heights for which only singular
            curves exist. For example, in the short Weierstrass case
            height 0 is permissible, as the curve Y^2 = X^3 (uniquely)
            has height zero.

        INPUT:

        - ``lowerbound`` -- Lower bound for the height range;
        - ``upperbound`` -- Upper bound for the height range. Heights
          returned are up to and including both bounds.

        OUTPUT:

        A list of tuples, each consisting of three elements of the form
        (H, C, I) such that
        1. H is the smallest height >= N
        2. C is a list of coefficients for curves of this height
        3. I is a list of indices indicating which of the above coefficients
           achieve this height. The remaining values in C  indicate the
           max absolute value those coefficients are allowed to obtain
           without altering the height.

        For example, the tuple (4, [1, 2], [1]) for the short
        Weierstrass case denotes set of curves with height 4; these
        are all of the form Y^2 = X^3 + A*X + B, where B=2 and A
        ranges between -1 and 1.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: C.heights(100, 150)
            [(100, [4, 10], [1]), (121, [4, 11], [1]), (125, [5, 11], [0]), (144, [5, 12], [1])]

            sage: C.heights(150, 100)
            Traceback (most recent call last):
            ...
            ValueError: Height upper bound must be greater than or equal to lower bound.
            sage: C.heights(-100, 100)
            Traceback (most recent call last):
            ...
            ValueError: Height lower bound must be non-negative.
        """

        if not lowerbound >= 0:
            raise ValueError("Height lower bound must be non-negative.")
        if not upperbound >= lowerbound:
            raise ValueError("Height upper bound must be greater "
                             "than or equal to lower bound.")

        coeffs = [ceil(lowerbound**(1/n))-1 for n in self._pows]
        height = max([coeffs[i]**(self._pows[i]) for i in range(self._num_coeffs)])

        L = []
        while height <= upperbound:
            C = self._height_increment(coeffs)
            if C[0] > upperbound:
                break
            else:
                height = C[0]
                coeffs = C[1]
                L.append(C)
        return L

    def _is_singular(self, C):
        r"""
        Tests if the a-invariants in 5-tuple C specify a singular
        elliptic curve.

        INPUT:

        - ``C`` -- A 5-tuple/list of a-invariants of a potential
          elliptic curve over `\QQ`

        OUTPUT:

        ``True`` or ``False``

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: C._is_singular([0,0,0,3, 2])
            False
            sage: EllipticCurve([0,0,0,3, 2])
            Elliptic Curve defined by y^2 = x^3 + 3*x + 2 over Rational Field
            sage: C._is_singular([0,0,0,-3, 2])
            True
            sage: EllipticCurve([0,0,0,-3, 2])
            Traceback (most recent call last):
            ...
            ArithmeticError: Invariants [0, 0, 0, -3, 2] define a singular curve.
        """
        a1 = C[0]
        a2 = C[1]
        a3 = C[2]
        a4 = C[3]
        a6 = C[4]

        b2 = a1**2 + 4*a2
        b4 = 2*a4 + a1*a3
        b6 = a3**2 + 4*a6
        b8 = (a1**2)*a6 + 4*a2*a6 - a1*a3*a4 + a2*(a3**2) - a4**2

        Delta = -(b2**2)*b8 - 8*(b4**3) - 27*(b6**2) + 9*b2*b4*b6
        return Delta == 0

    def _coeffs_from_height(self, height_tuple):
        r"""
        Returns a list of tuples of a-invariants of all curves
         described by height_tuple.

        INPUT:

        - ``height_tuple`` -- A tuple of the form
          (H, C, I) such that
          H: The smallest height >= N
          C: A list of coefficients for curves of this height
          I: A list of indices indicating which of the above coefficients
          achieve this height. The remaining values in C  indicate the
          max absolute value those coefficients are allowed to obtain
          without altering the height.

          For example, the tuple (4, [1, 2], [1]) for the short Weierstrass
          case denotes set of curves with height 4; these are all of the
          form Y^2 = X^3 + A*X + B, where B=2 and A ranges between -1 and 1.

        OUTPUT:

        A list of 2-tuples, each consisting of the given height,
        followed by a tuple of a-invariants of a curve of that height.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: B = C.next_height(4); B
            (4, [1, 2], [1])
            sage: L = C._coeffs_from_height(B)
            sage: for ell in L: print ell
            (4, [0, 0, 0, -1, -2])
            (4, [0, 0, 0, -1, 2])
            (4, [0, 0, 0, 0, -2])
            (4, [0, 0, 0, 0, 2])
            (4, [0, 0, 0, 1, -2])
            (4, [0, 0, 0, 1, 2])
        """
        height = height_tuple[0]
        coeffs = height_tuple[1]
        index = height_tuple[2]

        # Produce list of all coefficient tuples with given height
        L = []
        for S in list(powerset(index))[1:]:
            B = []
            for j in range(len(coeffs)):
                if j in S:
                    B.append([-coeffs[j], coeffs[j]])
                elif j in index:
                    B.append(srange(-coeffs[j]+1, coeffs[j]))
                else:
                    B.append(srange(-coeffs[j], coeffs[j]+1))
            C = CartesianProduct(*B).list()
            for c in C:
                L.append(c)

        # Convert coefficient tuples to a-invariants
        L2 = []
        for c in L:
            C = (height, self._coeffs_to_a_invariants(c))
            if not self._is_singular(C[1]):
                # Some families can produce duplicate sets of coefficients
                if not self._duplicates or not C in L2:
                    L2.append(C)
        return L2

    def _coeffs_from_height_list(self, coefficient_list):
        r"""
        Return all height/a-invariant tuples of elliptic curves from a
         list of curve height/coefficient/index tuples.

        INPUT:

        - ``coefficient_list`` -- A list of height/coefficient/index
          tuples. See the documentation for _height_increment() for
          a description of this format.

        OUTPUT:

        A list of tuples, each consisting of a height and a tuple of
        a-invariants defining an elliptic curve over `\QQ` of that height.
        The list will be ordered by increasing height.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: B = C.heights(1, 4); B
            [(1, [1, 1], [0, 1]), (4, [1, 2], [1])]
            sage: L = C._coeffs_from_height_list(B)
            sage: for ell in L: print(ell)
            (1, [0, 0, 0, -1, 0])
            (1, [0, 0, 0, 1, 0])
            (1, [0, 0, 0, 0, -1])
            (1, [0, 0, 0, 0, 1])
            (1, [0, 0, 0, -1, -1])
            (1, [0, 0, 0, -1, 1])
            (1, [0, 0, 0, 1, -1])
            (1, [0, 0, 0, 1, 1])
            (4, [0, 0, 0, -1, -2])
            (4, [0, 0, 0, -1, 2])
            (4, [0, 0, 0, 0, -2])
            (4, [0, 0, 0, 0, 2])
            (4, [0, 0, 0, 1, -2])
            (4, [0, 0, 0, 1, 2])
        """
        L2 = []
        for C in coefficient_list:
            L2 += self._coeffs_from_height(C)
        return L2

    def coefficients_over_height_range(self, lowerbound, upperbound,
                                       output_filename=None, return_data=True):
        r"""
        Return all a-invariant tuples of elliptic curves over a given
         height range, bounds inclusive.

        INPUT:

        - ``lowerbound``  -- The lower height bound

        - ``upperbound``  -- The upper height bound

        - ``output_filename`` -- (Default ``None``): If not ``None``, the string
          name of a file to which the output will be saved.

        - ``return_data`` -- (Default: ``True``) If ``False``, the computed
          data will not be returned.

        OUTPUT:

        If specified, the output is written to file. Each line specifies
        a single curve, and consists of six tab-separated integers in the
        format
        H  a1  a2  a3  a4  a6
        H: The height of the curve
        a1..a6: The a-invariants of the curve
        Curves will be ordered by increasing height.

        If return_data==True, the output is returned in the form of a
        list of tuples. Each tuple is of the form
        (H, [a1,a2,a3,a4,a6])
        where H and a1..a6 are as above.
        The list will be ordered by increasing height.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(0, 4)
            sage: for ell in L: print(ell)
            (1, [0, 0, 0, -1, 0])
            (1, [0, 0, 0, 1, 0])
            (1, [0, 0, 0, 0, -1])
            (1, [0, 0, 0, 0, 1])
            (1, [0, 0, 0, -1, -1])
            (1, [0, 0, 0, -1, 1])
            (1, [0, 0, 0, 1, -1])
            (1, [0, 0, 0, 1, 1])
            (4, [0, 0, 0, -1, -2])
            (4, [0, 0, 0, -1, 2])
            (4, [0, 0, 0, 0, -2])
            (4, [0, 0, 0, 0, 2])
            (4, [0, 0, 0, 1, -2])
            (4, [0, 0, 0, 1, 2])
        """
        H = self.heights(lowerbound, upperbound)
        L = self._coeffs_from_height_list(H)

        #WAS:   open(savefile,'w').write('\n'.join('\t'.join([str(a) for a in C])))
        #WAS: maybe leave in, but use consistent naming, e.g., output_filename...

        # Save data to file
        if output_filename is not None:
            out_file = open(output_filename, "w")
#            for C in [flatten[C] for C in L]:
#                out_file.write("\t".join([str(c) for c in C])+"\n")
            for C in L:
                out_file.write(str(C[0])+"\t")
                for a in C[1]:
                    out_file.write(str(a)+"\t")
                out_file.write("\n")

        if return_data:
            return L

    def rank(self, curves, output_filename=None, problems_filename=None,
             return_data=True, print_timing=True, **rank_opts):
        r"""
        Compute the algebraic rank for a list of curves ordered by height.

        INPUT:

        - ``curves``            -- A list of height/a-invariant tuples of
          curves, as returned by the coefficients_over_height_range() method
          Each tuple is of the form
          (H, [a1,a2,a3,a4,a6]) where
          H is the height of the curve, and
          [a1,...,a6] the curve's a-invariants

        - ``output_filename`` -- (Default ``None``): If not ``None``,
          the string name of the file to which the output will be
          saved.

        - ``problems_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which problem
          curves will be written.

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the data is not returned at the end of computation.

        - ``print_timing`` -- (Default ``True``): If set to ``False``,
          wall time of total computation will not be printed.

        Additional arguments are passed to the rank() method of the
        EllipticCurve class.

        OUTPUT:

        If specified, writes computed data to file. Each line of the
        written file consists of seven tab separated entries of the
        form H a1 a2 a3 a4 a6 d 1. H is the curve's height
        2. a1,...,a6 are the curve's a-invariants 3. d is the computed
        datum for that curve

        If specified, writes problem curves to file. Each line of the
        written file consists of seven tab separated entries of the
        form H a1 a2 a3 a4 a6 1. H is the curve's height 2. a1,...,a6
        are the curve's a-invariants

        If return_data==True:  A list consisting of two lists is returned:
        The first is a list of triples of the form
        (H, (a1,a2,a3,a4,a6), d)
        where the entries are as above.
        The second is a list of curve for which the datum could not be provably
        computed; each entry of this list is just a pair consisting of height
        and a-invariants.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(0, 4)
            sage: R, P = C.rank(L, return_data=True,print_timing=False)
            sage: P
            []
            sage: for r in R: print(r)
            (1, [0, 0, 0, -1, 0], 0)
            (1, [0, 0, 0, 1, 0], 0)
            (1, [0, 0, 0, 0, -1], 0)
            (1, [0, 0, 0, 0, 1], 0)
            (1, [0, 0, 0, -1, -1], 0)
            (1, [0, 0, 0, -1, 1], 1)
            (1, [0, 0, 0, 1, -1], 1)
            (1, [0, 0, 0, 1, 1], 1)
            (4, [0, 0, 0, -1, -2], 1)
            (4, [0, 0, 0, -1, 2], 0)
            (4, [0, 0, 0, 0, -2], 1)
            (4, [0, 0, 0, 0, 2], 1)
            (4, [0, 0, 0, 1, -2], 0)
            (4, [0, 0, 0, 1, 2], 0)
        """
        if print_timing:
            t = time.time()

        if output_filename is not None:
            out_file = open(output_filename, "w")
        if problems_filename is not None:
            prob_file = open(problems_filename, "w")
        if return_data:
            output = []
            problems = []

        for C in curves:
            # Attempt to compute rank and write curve+rank to file
            try:
                E = EllipticCurve(C[1])
                d = E.rank(**rank_opts)

                if output_filename is not None:
                    out_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        out_file.write(str(a)+"\t")
                    out_file.write(str(d)+"\n")
                    out_file.flush()

                if return_data:
                    output.append((C[0], C[1], d))

            # Write to problem file or append to problems list if fail
            except:
                if problems_filename is not None:
                    prob_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        prob_file.write(str(a)+"\t")
                    prob_file.write("\n")
                    prob_file.flush()

                if return_data:
                    problems.append(C)

        if output_filename is not None:
            out_file.close()
        if problems_filename is not None:
            prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output, problems

    def two_selmer(self, curves, rank=True, reduced=False,
                   output_filename=None, problems_filename=None,
                   return_data=True, print_timing=True):
        r"""
        Compute rank or size of two-Selmer for a list of curves
        ordered by height.

        INPUT:

        - ``curves``            -- A list of height/a-invariant tuples of
          curves, as returned by the coefficients_over_height_range() method.
          Each tuple is of the form
          (H, [a1,a2,a3,a4,a6]) where
          H is the height of the curve, and
          [a1,...,a6] the curve's a-invariants

        - ``rank`` -- (Default ``True``) Compute the rank versus size
          of the curve's 2-Selmer group. If ``True``, rank is
          computed; if set to ``False`` size (i.e. 2^rank) is computed
          instead.

        - ``reduced`` -- (Default ``False``) Compute full 2-Selmer or
          reduced 2-Selmer. If ``True``, full 2-Selmer is computed; if
          ``False``, the reduced group rank/size (i.e. 2-Selmer rank -
          2-torsion rank or 2^(2-Selmer rank - 2-torsion rank) as per
          whether 'rank' is set to ``True`` or ``False``) is computed.

          - ``output_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which the output
          will be saved.

        - ``problems_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which problem
          curves will be written.

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the data is not returned at the end of computation.

        - ``print_timing`` -- (Default ``True``): If set to ``False``,
          wall time of total computation will not be printed.

        OUTPUT:

         Writes data to file. Each line of the written file consists of seven
         tab separated entries of the form
         H, a1, a2, a3, a4, a6, d
         H: The curve's height
         a1,...,a6: The curve's a-invariants
         d: The computed datum for that curve

         (only if return_data==True) A list consisting of two lists:
         The first is a list of triples of the form
         (H, (a1,a2,a3,a4,a6), d)
         where the entries are as above.
         The second is a list of curve for which the datum could not be provably
         computed; each entry of this list is just a pair consisting of height
         and a-invariants.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(4, 4)
            sage: R, P = C.two_selmer(L, rank=True,return_data=True,print_timing=False)
            sage: P
            []
            sage: for r in R: print(r)
            (4, [0, 0, 0, -1, -2], 1)
            (4, [0, 0, 0, -1, 2], 0)
            (4, [0, 0, 0, 0, -2], 1)
            (4, [0, 0, 0, 0, 2], 1)
            (4, [0, 0, 0, 1, -2], 1)
            (4, [0, 0, 0, 1, 2], 1)
            sage: R, P = C.two_selmer(L, rank=False,return_data=True,print_timing=False)
            sage: for r in R: print(r)
            (4, [0, 0, 0, -1, -2], 2)
            (4, [0, 0, 0, -1, 2], 1)
            (4, [0, 0, 0, 0, -2], 2)
            (4, [0, 0, 0, 0, 2], 2)
            (4, [0, 0, 0, 1, -2], 2)
            (4, [0, 0, 0, 1, 2], 2)
            sage: R, P = C.two_selmer(L, reduced=True,print_timing=False)
            sage: for r in R: print(r)
            (4, [0, 0, 0, -1, -2], 1)
            (4, [0, 0, 0, -1, 2], 0)
            (4, [0, 0, 0, 0, -2], 1)
            (4, [0, 0, 0, 0, 2], 1)
            (4, [0, 0, 0, 1, -2], 0)
            (4, [0, 0, 0, 1, 2], 0)
            sage: R, P = C.two_selmer(L, rank=False,reduced=True,print_timing=False)
            sage: for r in R: print(r)
            (4, [0, 0, 0, -1, -2], 2)
            (4, [0, 0, 0, -1, 2], 1)
            (4, [0, 0, 0, 0, -2], 2)
            (4, [0, 0, 0, 0, 2], 2)
            (4, [0, 0, 0, 1, -2], 1)
            (4, [0, 0, 0, 1, 2], 1)
        """
        if print_timing:
            t = time.time()

        if output_filename is not None:
            out_file = open(output_filename, "w")
        if problems_filename is not None:
            prob_file = open(problems_filename, "w")
        if return_data:
            output = []
            problems = []

        for C in curves:
            # Attempt to compute datum and write curve+datum to file
            try:
                E = EllipticCurve(C[1])

                d = E.selmer_rank()
                if reduced:
                    d -= E.two_torsion_rank()
                if not rank:
                    d = 2**d

                if output_filename is not None:
                    out_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        out_file.write(str(a)+"\t")
                    out_file.write(str(d)+"\n")
                    out_file.flush()

                if return_data:
                    output.append((C[0], C[1], d))

            # Write to problem file if fail
            except:
                if problems_filename is not None:
                    prob_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        prob_file.write(str(a)+"\t")
                    prob_file.write("\n")
                    prob_file.flush()

                if return_data:
                    problems.append(C)

        if output_filename is not None:
            out_file.close()
        if problems_filename is not None:
            prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output, problems

    def three_selmer(self, curves, rank=True, reduced=False,
                     output_filename=None, problems_filename=None,
                     proof=True, return_data=True, print_timing=True):
        r"""
        Compute rank or size of two-Selmer for a list of curves
        ordered by height using Magma.

        .. WARNING::

            This function will only work if Magma is installed.

        INPUT:

        - ``curves``            -- A list of height/a-invariant tuples of
          curves, as returned by the coefficients_over_height_range() method.
          Each tuple is of the form
          (H, [a1,a2,a3,a4,a6]) where
          H is the height of the curve, and
          [a1,...,a6] the curve's a-invariants

        - ``rank`` -- (Default ``True``) Compute the rank versus size
          of the curve's 3-Selmer group. If ``True``, rank is
          computed; if set to ``False`` size (i.e. 3^rank) is computed
          instead.

        - ``reduced`` -- (Default ``False``) Compute full 3-Selmer or
          reduced 3-Selmer. If ``True``, full 3-Selmer is computed; if
          ``False``, the reduced group rank/size (i.e. 3-Selmer rank -
          3-torsion rank or 3^(3-Selmer rank - 3-torsion rank) as per
          whether 'rank' is set to ``True`` or ``False``) is computed.

          - ``output_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which the output
          will be saved.

        - ``problems_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which problem
          curves will be written.

        - ``proof`` -- (Default ``True``): If ``False``, the
          Generalized Riemann Hypothesis will not be assumed in the
          computing of class group bounds in Magma (and thus will be
          slower); if ``True``, GRH will be assumed and computation
          will be quicker.

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the data is not returned at the end of computation.

        - ``print_timing`` -- (Default ``True``): If set to ``False``,
          wall time of total computation will not be printed.

        OUTPUT:

         Writes data to file. Each line of the written file consists of seven
         tab separated entries of the form
         H, a1, a2, a3, a4, a6, d
         H: The curve's height
         a1,...,a6: The curve's a-invariants
         d: The computed datum for that curve

         (only if return_data==True) A list consisting of two lists:
         The first is a list of triples of the form
         (H, (a1,a2,a3,a4,a6), d)
         where the entries are as above.
         The second is a list of curve for which the datum could not be provably
         computed; each entry of this list is just a pair consisting of height
         and a-invariants.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(4, 4)
            sage: R, P = C.three_selmer(L, rank=True,return_data=True,print_timing=False) # optional - magma
            sage: R # optional - magma
            [(4, [0, 0, 0, -1, -2], 1),
            (4, [0, 0, 0, -1, 2], 0),
            (4, [0, 0, 0, 0, -2], 1),
            (4, [0, 0, 0, 0, 2], 1),
            (4, [0, 0, 0, 1, -2], 0),
            (4, [0, 0, 0, 1, 2], 0)]
            sage: P # optional - magma
            []
        """
        from sage.interfaces.all import magma

        if print_timing:
            t = time.time()

        if not proof:
            magma.eval('SetClassGroupBounds("GRH")')
        else:
            magma.eval('SetClassGroupBounds("PARI")')

        if output_filename is not None:
            out_file = open(output_filename, "w")
        if problems_filename is not None:
            prob_file = open(problems_filename, "w")
        if return_data:
            output = []
            problems = []

        for C in curves:
            # Attempt to compute datum and write curve+datum to file
            try:
                E = EllipticCurve(C[1])

                d = E.three_selmer_rank()
                if reduced:
                    d -= valuation(E.torsion_order(), 3)
                if not rank:
                    d = 3**d

                if output_filename is not None:
                    out_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        out_file.write(str(a)+"\t")
                    out_file.write(str(d)+"\n")
                    out_file.flush()

                if return_data:
                    output.append((C[0], C[1], d))

            # Write to problem file if fail
            except:
                if problems_filename is not None:
                    prob_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        prob_file.write(str(a)+"\t")
                    prob_file.write("\n")
                    prob_file.flush()

                if return_data:
                    problems.append(C)

        if output_filename is not None:
            out_file.close()
        if problems_filename is not None:
            prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output, problems

    def two_torsion(self, curves, rank=True,
                    output_filename=None, problems_filename=None,
                    return_data=True, print_timing=True):
        r"""
        Compute rank or size of two-torsion groups for a list of
        curves ordered by height.

        INPUT:

        - ``curves``            -- A list of (height,a-invariant) tuples of
          curves, as returned by the coefficients_over_height_range() method.
          Each tuple is of the form
          (H, [a1,a2,a3,a4,a6]) where
          H is the height of the curve, and
          [a1,...,a6] the curve's a-invariants

        - ``rank`` -- (Default ``True``) Compute the rank versus size
          of the curve's 2-torsion group. If ``True``, rank is
          computed; if set to ``False`` size (i.e. 2^rank) is computed
          instead.

          - ``output_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which the output
          will be saved.

        - ``problems_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which problem
          curves will be written.

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the data is not returned at the end of computation.

        - ``print_timing`` -- (Default ``True``): If set to ``False``,
          wall time of total computation will not be printed.

        OUTPUT:

        Writes data to file. Each line of the written file consists of seven
        tab separated entries of the form
        H, a1, a2, a3, a4, a6, d
        H: The curve's height
        a1,...,a6: The curve's a-invariants
        d: The computed datum for that curve

        (only if return_data==True) A list consisting of two lists:
        The first is a list of triples of the form
        (H, (a1,a2,a3,a4,a6), d)
        where the entries are as above.
        The second is a list of curve for which the datum could not be provably
        computed; each entry of this list is just a pair consisting of height
        and a-invariants.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(4, 4)
            sage: R, P = C.two_torsion(L, rank=True,return_data=True,print_timing=False)
            sage: P
            []
            sage: for r in R: print(r)
            (4, [0, 0, 0, -1, -2], 0)
            (4, [0, 0, 0, -1, 2], 0)
            (4, [0, 0, 0, 0, -2], 0)
            (4, [0, 0, 0, 0, 2], 0)
            (4, [0, 0, 0, 1, -2], 1)
            (4, [0, 0, 0, 1, 2], 1)
            sage: R, P = C.two_torsion(L, rank=False,return_data=True,print_timing=False)
            sage: for r in R: print(r)
            (4, [0, 0, 0, -1, -2], 1)
            (4, [0, 0, 0, -1, 2], 1)
            (4, [0, 0, 0, 0, -2], 1)
            (4, [0, 0, 0, 0, 2], 1)
            (4, [0, 0, 0, 1, -2], 2)
            (4, [0, 0, 0, 1, 2], 2)
        """
        if print_timing:
            t = time.time()

        if output_filename is not None:
            out_file = open(output_filename, "w")
        if problems_filename is not None:
            prob_file = open(problems_filename, "w")
        if return_data:
            output = []
            problems = []

        for C in curves:
            # Attempt to compute datum and write curve+datum to file
            try:
                E = EllipticCurve(C[1])

                d = E.two_torsion_rank()
                if not rank:
                    d = 2**d

                if output_filename is not None:
                    out_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        out_file.write(str(a)+"\t")
                    out_file.write(str(d)+"\n")
                    out_file.flush()

                if return_data:
                    output.append((C[0], C[1], d))

            # Write to problem file if fail
            except:
                if problems_filename is not None:
                    prob_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        prob_file.write(str(a)+"\t")
                    prob_file.write("\n")
                    prob_file.flush()

                if return_data:
                    problems.append(C)

        if output_filename is not None:
            out_file.close()
        if problems_filename is not None:
            prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output, problems

    def n_torsion(self, curves, n, rank=True,
                  output_filename=None, problems_filename=None,
                  return_data=True, print_timing=True):
        r"""
        Compute size of torsion groups for a list of curves ordered by height;
        or compute the size of n-torsion of the curves;
        or compute the rank of n-torsion of the curves.

        INPUT:

        - ``curves``            -- A list of (height,a-invariant) tuples of
          curves, as returned by the coefficients_over_height_range() method.
          Each tuple is of the form
          (H, [a1,a2,a3,a4,a6]) where
          H is the height of the curve, and
          [a1,...,a6] the curve's a-invariants

        - ``n`` -- The n for which the rank/size of the n-torsion
          subgroup will be computed

        - ``rank`` -- (Default ``True``): Compute the rank versus size
          of the curve's torsion/n-torsion subgroup. If ``True``, the
          n-torsion rank is computed (i.e. 0, 1 or 2). If set to
          ``False``, the size of the n-torsion subgroup is computed
          instead (i.e. n^rank)

          - ``output_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which the output
          will be saved.

        - ``problems_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which problem
          curves will be written.

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the data is not returned at the end of computation.

        - ``print_timing`` -- (Default ``True``): If set to ``False``,
          wall time of total computation will not be printed.

        OUTPUT:

        Writes data to file. Each line of the written file consists of seven
        tab separated entries of the form
        H, a1, a2, a3, a4, a6, d
        H: The curve's height
        a1,...,a6: The curve's a-invariants
        d: The computed datum for that curve

        (only if return_data==True) A list consisting of two lists:
        The first is a list of triples of the form
        (H, (a1,a2,a3,a4,a6), d)
        where the entries are as above.
        The second is a list of curve for which the datum could not be provably
        computed; each entry of this list is just a pair consisting of height
        and a-invariants.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(1, 1)
            sage: R, P = C.n_torsion(L, n=3,rank=True,return_data=True,print_timing=False)
            sage: P
            []
            sage: for r in R: print(r)
            (1, [0, 0, 0, -1, 0], 0)
            (1, [0, 0, 0, 1, 0], 0)
            (1, [0, 0, 0, 0, -1], 0)
            (1, [0, 0, 0, 0, 1], 1)
            (1, [0, 0, 0, -1, -1], 0)
            (1, [0, 0, 0, -1, 1], 0)
            (1, [0, 0, 0, 1, -1], 0)
            (1, [0, 0, 0, 1, 1], 0)
            sage: R, P = C.n_torsion(L, n=3,rank=False,return_data=True,print_timing=False)
            sage: for r in R: print(r)
            (1, [0, 0, 0, -1, 0], 1)
            (1, [0, 0, 0, 1, 0], 1)
            (1, [0, 0, 0, 0, -1], 1)
            (1, [0, 0, 0, 0, 1], 3)
            (1, [0, 0, 0, -1, -1], 1)
            (1, [0, 0, 0, -1, 1], 1)
            (1, [0, 0, 0, 1, -1], 1)
            (1, [0, 0, 0, 1, 1], 1)
            sage: R, P = C.n_torsion(L, n=6,rank=False,return_data=True,print_timing=False)
            sage: for r in R: print(r)
            (1, [0, 0, 0, -1, 0], 1)
            (1, [0, 0, 0, 1, 0], 1)
            (1, [0, 0, 0, 0, -1], 1)
            (1, [0, 0, 0, 0, 1], 6)
            (1, [0, 0, 0, -1, -1], 1)
            (1, [0, 0, 0, -1, 1], 1)
            (1, [0, 0, 0, 1, -1], 1)
            (1, [0, 0, 0, 1, 1], 1)
        """
        if not n > 0:
            raise ValueError("n must be a positive integer.")

        # Use self.two_torsion() for n=2
        if n == 2:
            return self.two_torsion(curves, rank=rank,
                                    output_filename=output_filename,
                                    problems_filename=problems_filename,
                                    return_data=return_data,
                                    print_timing=print_timing)

        if print_timing:
            t = time.time()

        if output_filename is not None:
            out_file = open(output_filename, "w")
        if problems_filename is not None:
            prob_file = open(problems_filename, "w")
        if return_data:
            output = []
            problems = []

        for C in curves:
            # Attempt to compute datum and write curve+datum to file
            try:
                E = EllipticCurve(C[1])

#####THIS REQUIRES CHECKING#####
                # Compute datum
                # By Mazur's Torsion Theorem, the possible torsion subgroups
                # for curves over QQ are C1 thru C10, and C2xC2, C2xC4, C2xC6
                # and C2xC8, where Cm is the cyclic group of order m.
                d = E.torsion_order()
                t = E.two_torsion_rank()
                if t == 2:
                    d = d/2
                if n.divides(d):
                    # The only n for which E may have full n-torsion is 2,
                    # which is covered by self.two_torsion()
                    d = 1
                    if not rank:
                        d = n
                        # If n is even and the curve has full 2-torsion,
                        # there are twice as many points of order n
                        if t == 2 and n % 2 == 0:
                            d = 2*d
                # Boring case: no n-torsion
                else:
                    d = 0
                    if not rank:
                        d = 1
#####

                # Write to file and/or append to return list
                if output_filename is not None:
                    out_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        out_file.write(str(a)+"\t")
                    out_file.write(str(d)+"\n")
                    out_file.flush()

                if return_data:
                    output.append((C[0], C[1], d))

            # Write to problem file if fail
            except:
                if problems_filename is not None:
                    prob_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        prob_file.write(str(a)+"\t")
                    prob_file.write("\n")
                    prob_file.flush()

                if return_data:
                    problems.append(C)

        if output_filename is not None:
            out_file.close()
        if problems_filename is not None:
            prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output, problems

    def torsion_order(self, curves, output_filename=None, problems_filename=None,
                      return_data=True, print_timing=True):
        r"""
        Compute size of the torsion subgroups of a list of curves
        ordered by height.

        INPUT:

        - ``curves``            -- A list of (height,a-invariant) tuples of
          curves, as returned by the coefficients_over_height_range() method.
          Each tuple is of the form
          (H, [a1,a2,a3,a4,a6]) where
          H is the height of the curve, and
          [a1,...,a6] the curve's a-invariants

          - ``output_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which the output
          will be saved.

        - ``problems_filename`` -- (Default ``None``): If not
          ``None``, the string name of the file to which problem
          curves will be written.

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the data is not returned at the end of computation.

        - ``print_timing`` -- (Default ``True``): If set to ``False``,
          wall time of total computation will not be printed.

        OUTPUT:

        Writes data to file. Each line of the written file consists of seven
        tab separated entries of the form
        H, a1, a2, a3, a4, a6, d
        H: The curve's height
        a1,...,a6: The curve's a-invariants
        d: The computed datum for that curve

        (only if return_data==True) A list consisting of two lists:
        The first is a list of triples of the form
        (H, (a1,a2,a3,a4,a6), d)
        where the entries are as above.
        The second is a list of curve for which the datum could not be provably
        computed; each entry of this list is just a pair consisting of height
        and a-invariants.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(1, 1)
            sage: R, P = C.torsion_order(L, return_data=True,print_timing=False)
            sage: P
            []
            sage: for r in R: print(r)
            (1, [0, 0, 0, -1, 0], 4)
            (1, [0, 0, 0, 1, 0], 2)
            (1, [0, 0, 0, 0, -1], 2)
            (1, [0, 0, 0, 0, 1], 6)
            (1, [0, 0, 0, -1, -1], 1)
            (1, [0, 0, 0, -1, 1], 1)
            (1, [0, 0, 0, 1, -1], 1)
            (1, [0, 0, 0, 1, 1], 1)
        """
        if print_timing:
            t = time.time()

        if output_filename is not None:
            out_file = open(output_filename, "w")
        if problems_filename is not None:
            prob_file = open(problems_filename, "w")
        if return_data:
            output = []
            problems = []

        for C in curves:
            # Attempt to compute datum and write curve+datum to file
            try:
                E = EllipticCurve(C[1])
                d = E.torsion_order()

                # Write to file and/or append to return list
                if output_filename is not None:
                    out_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        out_file.write(str(a)+"\t")
                    out_file.write(str(d)+"\n")
                    out_file.flush()

                if return_data:
                    output.append((C[0], C[1], d))

            # Write to problem file if fail
            except:
                if problems_filename is not None:
                    prob_file.write(str(C[0])+"\t")
                    for a in C[1]:
                        prob_file.write(str(a)+"\t")
                    prob_file.write("\n")
                    prob_file.flush()

                if return_data:
                    problems.append(C)

        if output_filename is not None:
            out_file.close()
        if problems_filename is not None:
            prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output, problems

    def averaged_data(self, input, output_filename=None, return_data=True):
        """
        INPUT:

        - ``input`` -- Either: A list where each element is
          (H, (a1,a2,a3,a4,a6), d), and
          H: curve height
          a1,...,a6: list of a-invariants of elliptic curve giving that height
          d: the type of data that's being averaged, e.g., rank, 2-Selmer,
          Or: String name of a file of array of data, where each line is
          H, a1, a2, a3, a4, a6, d
          where the elements are as above

        - ``output_filename`` -- (Default ``None``): If not ``None``,
          the string name of the file to which computed averaged data
          will be written

        - ``return_data`` -- (Default ``True``): If set to ``False``,
          the computed data is not returned.

        OUTPUT:

        If specified, writes the averaged data to file, where the two
        tab-separated elements of each line are of the form
        H  a
        where H is height, and a the average datum up to that height

        If return_data==True, returns a list of tuples of the form
        (H, a)
        where H is height, and a the average datum up to that height

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: L = C.coefficients_over_height_range(0, 4)
            sage: R = C.rank(L, return_data=true,print_timing=False)
            sage: A = C.averaged_data(R[0],return_data=True); A
            [(1.0, 0.375), (4.0, 0.42857142857142855)]
        """
        import numpy as np

        # Load, format data
        if isinstance(input, str):
            data = np.loadtxt(input)
            X = data[:, 0]
            V = data[:, -1]
        elif isinstance(input, list):
            X = np.array([C[0] for C in input])
            V = np.array([C[-1] for C in input])
        else:
            raise IOError("Input must either be a string or a list")

        #Compute running average
        N = np.arange(X.shape[0], dtype=np.float64) + 1
        Y = np.cumsum(V)/N

        # Retain only lines with new heights
        I = X[:-1] != X[1:]
        I = np.append(I, True)

        # Save, return output
        Z = np.vstack([X[I], Y[I, :].T]).T
        if output_filename is not None:
            np.savetxt(output_filename, Z)
        if return_data:
            return [(C[0], C[1]) for C in Z]


class CurveEnumeratorShortWeierstrass(CurveEnumerator_abstract):
    """
    Height iterator for curves in short Weierstrass form.
    Family:             Short Weierstrass
    Model:              Y^2 = X^3 + A*X + B
    Coefficients:       [A,B]
    Height function:    H = min{|A|^3,|B|^2}
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorShortWeierstrass
            sage: C = CurveEnumeratorShortWeierstrass(); C
            Height iterator for elliptic curves over QQ
            Family:             Short Weierstrass
            Model:              Y^2 = X^3 + A*X + B
            Coefficients:       [A,B]
            Height function:    H = min{|A|^3,|B|^2}
        """
        # For __repr__()
        self._name = "Short Weierstrass"
        self._model = "Y^2 = X^3 + A*X + B"
        self._coeff_names = "[A,B]"
        self._height_function = "H = min{|A|^3,|B|^2}"

        # The following constants must be Sage Integers; if not, some
        # methods won't work
        self._num_coeffs = ZZ(2)
        self._pows = (ZZ(3), ZZ(2))

        # Each coefficient tuple corresponds to a unique curve
        self._duplicates = False

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator(family="short_weierstrass")
            sage: C._coeffs_to_a_invariants([4,9])
            [0, 0, 0, 4, 9]
        """
        return [0, 0, 0, c[0], c[1]]


# class CurveEnumeratorFullWeierstrass(CurveEnumerator_abstract):
#     """
#     Height iterator for elliptic curves in full Weierstrass form.
#     Not yet implemented.

#     EXAMPLES::

#         sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorFullWeierstrass
#         sage: C = CurveEnumeratorFullWeierstrass()
#         Traceback (most recent call last):
#         ...
#         NotImplementedError: Family not yet implemented.

#     """
#     def __init__(self):
#         raise NotImplementedError("Family not yet implemented.")

#     def __repr__(self):
#         raise NotImplementedError("Family not yet implemented.")

#     def _coeffs_to_a_invariants(self, c):
#         raise NotImplementedError("Family not yet implemented.")


class CurveEnumeratorRankOne(CurveEnumerator_abstract):
    """
    Height iterator for family with a high incidence of rank one curves.

    Family:             Rank One
    Model:              Y^2 + A*Y = X^3 + B*X^2 + C*X
    Coefficients:       [A,B,C]
    Height function:    H = min{|A|^6,|B|^4,|C|^3}

    The point (0, 0) is often free.
    """
    def __init__(self):
        """
        Instantiates self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorRankOne
            sage: C = CurveEnumeratorRankOne(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank One
            Model:              Y^2 + A*Y = X^3 + B*X^2 + C*X
            Coefficients:       [A,B,C]
            Height function:    H = min{|A|^6,|B|^4,|C|^3}
        """
        # For __init__()
        self._name = "Rank One"
        self._model = "Y^2 + A*Y = X^3 + B*X^2 + C*X"
        self._coeff_names = "[A,B,C]"
        self._height_function = "H = min{|A|^6,|B|^4,|C|^3}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(3)
        self._pows = (ZZ(6), ZZ(4), ZZ(3))

        # Each coefficient tuple corresponds to a unique curve
        self._duplicates = False

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("rank_one")
            sage: C._coeffs_to_a_invariants([1, 2, 3])
            [0, 1, 2, 3, 0]
        """
        return [0, c[0], c[1], c[2], 0]


class CurveEnumeratorRankTwo(CurveEnumerator_abstract):
    """
    Height iterator for family with a high incidence of rank two curves.

    Family:             Rank Two
    Model:              Y^2 = X^3 - 27*I*X - 27*J, where
                        I = 3*a4^2 + b4^2 - 3*a2*a6, and
                        J = -27/4*a2^2*a4^2 + 18*a4^2*b4 - 2*b4^3 + 9*a2*b4*a6 - 27*a6^2
    Coefficients:       [a2,a4,b4,c6]
    Height function:    H = min{|a2|^6,|a4|^3,|b4|^3,|c6|^2}
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorRankTwo
            sage: C = CurveEnumeratorRankTwo(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank Two
            Model:              Y^2 = X^3 - 27*I*X - 27*J, where
                                I = 3*a4^2 + b4^2 - 3*a2*a6, and
                                J = -27/4*a2^2*a4^2 + 18*a4^2*b4 - 2*b4^3 + 9*a2*b4*a6 - 27*a6^2
            Coefficients:       [a2,a4,b4,c6]
            Height function:    H = min{|a2|^6,|a4|^3,|b4|^3,|c6|^2}
        """
        # For __init__()
        self._name = "Rank Two"
        self._model = "Y^2 = X^3 - 27*I*X - 27*J, where \n"\
            + "                    I = 3*a4^2 + b4^2 - 3*a2*a6, and \n"\
            + "                    J = -27/4*a2^2*a4^2 + 18*a4^2*b4 - 2*b4^3 + 9*a2*b4*a6 - 27*a6^2"
        self._coeff_names = "[a2,a4,b4,c6]"
        self._height_function = "H = min{|a2|^6,|a4|^3,|b4|^3,|c6|^2}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(4)
        self._pows = (ZZ(6), ZZ(3), ZZ(3), ZZ(2))

        # Each coefficient tuple corresponds to a unique curve
        self._duplicates = False

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::
            sage: C = EllipticCurveEnumerator("rank_two")
            sage: C._coeffs_to_a_invariants([1, 2, 3, 4])
            [0, 0, 0, -243, 5130]

            sage: C._coeffs_to_a_invariants([1, 3,5,7])
            [0, 0, 0, -837, 13797]
        """
        a2, a4, b4, a6 = c[0], c[1], c[2], c[3]
        I = 3*(a4**2)+b4**2-3*a2*a6
        J = -27/4*(a2**2)*a4**2+18*(a4**2)*b4 - 2*(b4**3) \
            + 9*a2*b4*a6 - 27*(a6**2)

        # J may have denominator 2 or 4. If so, the following produces
        # an isomorphic curve with integral coefficients
        if J.denominator() > 1:
            I = I*16
            J = J*64

        return [0, 0, 0, -27*I, -27*J]


class CurveEnumeratorTwoTorsion(CurveEnumerator_abstract):
    """
    Height iterator for family of curves with guaranteed two torsion.

    Family:             Two Torsion
    Model:              Y^2  = X^3 + A*X^2 + B*X
    Coefficients:       [A,B]
    Height function:    H = min{|A|^6,|C|^3}

    The point (0,0) is two torsion.
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorTwoTorsion
            sage: C = CurveEnumeratorTwoTorsion(); C
            Height iterator for elliptic curves over QQ
            Family:             Two Torsion
            Model:              Y^2  = X^3 + A*X^2 + B*X
            Coefficients:       [A,B]
            Height function:    H = min{|A|^6,|C|^3}
        """
        # For __repr__()
        self._name = "Two Torsion"
        self._model = "Y^2  = X^3 + A*X^2 + B*X"
        self._coeff_names = "[A,B]"
        self._height_function = "H = min{|A|^6,|C|^3}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(2)
        self._pows = (ZZ(6), ZZ(3))

        # Each coefficient tuple corresponds to a unique curve
        self._duplicates = False

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("two_torsion")
            sage: C._coeffs_to_a_invariants([1, 2])
            [0, 1, 0, 2, 0]
        """
        return [0, c[0], 0, c[1], 0]


class CurveEnumeratorThreeTorsion(CurveEnumerator_abstract):
    """
    Height Iterator for family of curves with guaranteed three torsion.

    Family:             Three Torsion
    Model:              Y^2 + A*X*Y + B*Y = X^3
    Coefficients:       [A,B]
    Height function:    H = min{|A|^12,|B|^4}

    The point (0, 0) is three torsion.
    """
    def __init__(self):
        """
        Crates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorThreeTorsion
            sage: C = CurveEnumeratorThreeTorsion(); C
            Height iterator for elliptic curves over QQ
            Family:             Three Torsion
            Model:              Y^2 + A*X*Y + B*Y = X^3
            Coefficients:       [A,B]
            Height function:    H = min{|A|^12,|B|^4}
        """
        # For __repr__()
        self._name = "Three Torsion"
        self._model = "Y^2 + A*X*Y + B*Y = X^3"
        self._coeff_names = "[A,B]"
        self._height_function = "H = min{|A|^12,|B|^4}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(2)
        self._pows = (ZZ(12), ZZ(4))

        # Each coefficient tuple corresponds to a unique curve
        self._duplicates = False

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("three_torsion")
            sage: C._coeffs_to_a_invariants([1, 2])
            [1, 0, 2, 0, 0]
        """
        return [c[0], 0, c[1], 0, 0]


class CurveEnumeratorF_12(CurveEnumerator_abstract):
    """
    Height iterator for family of with high incidence of rank one curves with
    two torsion.

    Family:             Rank One + Two Torsion
    Model:              Y^2 = X^3 + (A^2-B-C)*X^2 + (B*C)*X
    Coefficients:       [A,B,C]
    Height function:    H = min{|A|,|B|^2,|C|^2}
    """
    def __init__(self):
        """
        Create an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorF_12
            sage: C = CurveEnumeratorF_12(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank One + Two Torsion
            Model:              Y^2 = X^3 + (A^2-B-C)*X^2 + (B*C)*X
            Coefficients:       [A,B,C]
            Height function:    H = min{|A|,|B|^2,|C|^2}
        """
        # for __repr__()
        self._name = "Rank One + Two Torsion"
        self._model = "Y^2 = X^3 + (A^2-B-C)*X^2 + (B*C)*X"
        self._coeff_names = "[A,B,C]"
        self._height_function = "H = min{|A|,|B|^2,|C|^2}"

        # The following constants must be Sage Integers; if not, kittens die
        self._num_coeffs = ZZ(3)
        self._pows = (ZZ(1), ZZ(2), ZZ(2))

        # Coefficients, converted to a-invariants, will produce duplicate curves
        self._duplicates = True

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("F_1(2)")
            sage: C._coeffs_to_a_invariants([1, 2, 3])
            [0, -4, 0, 6, 0]
        """
        return [0, c[0]**2-c[1]-c[2], 0, c[1]*c[2], 0]


class CurveEnumeratorF_22(CurveEnumerator_abstract):
    """
    Height iterator for family with high incidence of rank two curves with
    two torsion.

    Family:             Rank Two + Two Torsion
    Model:              Y^2 = X^3 + A*X^2 + B, where
                        A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                        B = w0*w1*w2*w3
    Coefficients:       [w0,w1,w2,w3]
    Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorF_22
            sage: C = CurveEnumeratorF_22(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank Two + Two Torsion
            Model:              Y^2 = X^3 + A*X^2 + B, where
                                A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                                B = w0*w1*w2*w3
            Coefficients:       [w0,w1,w2,w3]
            Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
        """
        # For __repr__()
        self._name = "Rank Two + Two Torsion"
        self._model = "Y^2 = X^3 + A*X^2 + B, where \n" + (" "*10) +\
            "A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and \n"\
            + (" "*10) + "B = w0*w1*w2*w3"
        self._coeff_names = "[w0,w1,w2,w3]"
        self._height_function = "H = min{|w0|,|w1|,|w2|,|w3|}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(4)
        self._pows = (ZZ(1), ZZ(1), ZZ(1), ZZ(1))

        # Coefficients, when converted to a-invariants, will produce
        # duplicate curves
        self._duplicates = True

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("F_2(2)")
            sage: C._coeffs_to_a_invariants([1, 2, 3, 4])
            [0, -10, 0, 0, 24]
        """
        w0, w1, w2, w3 = c[0], c[1], c[2], c[3]
        A = (2*(w0**2+w1**2+w2**2+w3**2) - (w0+w1+w2+w3)**2)/4
        B = w0*w1*w2*w3

        # A may be be a rational with denominator 2 or 4.
        # If so, the following converts to an integral model:
        if A.denominator() > 1:
            A = A*4
            B = B*64

        return [0, A, 0, 0, B]


class CurveEnumeratorF_13(CurveEnumerator_abstract):
    """
    Height iterator for family with high incidence of rank one curves with
    three torsion.

    Model:              Y^2 + A*X*Y + B*Y = X^3, where
                        A = w0 + w1 + w2, and
                        B = w0*w1*w2
    Coefficients:       [w0,w1,w2]
    Height function:    H = min{|w0|,|w1|,|w2|}
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorF_13
            sage: C = CurveEnumeratorF_13(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank One + Three Torsion
            Model:              Y^2 + A*X*Y + B*Y = X^3, where
                                A = w0 + w1 + w2, and
                                B = w0*w1*w2
            Coefficients:       [w0,w1,w2]
            Height function:    H = min{|w0|,|w1|,|w2|}
        """
        # for __repr__()
        self._name = "Rank One + Three Torsion"
        self._model = "Y^2 + A*X*Y + B*Y = X^3, where \n" + \
            (" " * 10) + "A = w0 + w1 + w2, and \n" + \
            (" " * 10) + "B = w0*w1*w2"
        self._coeff_names = "[w0,w1,w2]"
        self._height_function = "H = min{|w0|,|w1|,|w2|}"

        # The following constants must be Sage Integers;
        self._num_coeffs = ZZ(3)
        self._pows = (ZZ(1), ZZ(1), ZZ(1))

        # Coefficients, when converted to a-invariants, will produce
        # duplicate curves
        self._duplicates = True

    def _coeffs_to_a_invariants(self, c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

        - ``c`` -- The list of coefficients of the equation of a curve
          as per the family model description. See the __init__() method
          of this class for more info.

        OUTPUT:

        A list of five integers corresponding to the a-invariants of
        the curve.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("F_1(3)")
            sage: C._coeffs_to_a_invariants([1, 2, 4])
            [7, 0, 8, 0, 0]
        """
        A = c[0]+c[1]+c[2]
        B = c[0]*c[1]*c[2]
        return [A, 0, B, 0, 0]


class CurveEnumeratorF_14(CurveEnumeratorF_22):
    """
    Height iterator for family with high incidence of rank one curves with
    four torsion.

    Family:             Rank One + Four Torsion
    Model:              Y^2 = X^3 + A*X^2 + B, where
                        A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                        B = w0*w1*w2*w3
    Coefficients:       [w0,w1,w2,w3], with w0*w1 = w2*w3
    Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorF_14
            sage: C = CurveEnumeratorF_14(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank One + Four Torsion
            Model:              Y^2 = X^3 + A*X^2 + B, where
                                A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                                B = w0*w1*w2*w3
            Coefficients:       [w0,w1,w2,w3], with w0*w1 = w2*w3
            Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
        """
        # For repr__()
        self._name = "Rank One + Four Torsion"
        self._model = "Y^2 = X^3 + A*X^2 + B, where \n" +\
                      (" " * 10) + "A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and \n" +\
                      (" " * 10) + "B = w0*w1*w2*w3"
        self._coeff_names = "[w0,w1,w2,w3], with w0*w1 = w2*w3"
        self._height_function = "H = min{|w0|,|w1|,|w2|,|w3|}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(4)
        self._pows = (ZZ(1), ZZ(1), ZZ(1), ZZ(1))

        # Coefficients, when converted to a-invariants, will produce
        # duplicate curves
        self._duplicates = True

    def _coeffs_from_height(self, height_tuple):
        """
        Returns a list of tuples of a-invariants of all curves
         described by height_tuple.

        INPUT:

        - ``height_tuple`` -- A tuple of the form
          (H, C, I) such that
          H: The smallest height >= N
          C: A list of coefficients for curves of this height
          I: A list of indices indicating which of the above coefficients
          achieve this height. The remaining values in C  indicate the
          max absolute value those coefficients are allowed to obtain
          without altering the height.

          For example, the tuple (4, [1, 2], [1]) for the short Weierstrass
          case denotes set of curves with height 4; these are all of the
          form Y^2 = X^3 + A*X + B, where B=2 and A ranges between -1 and 1.

        OUTPUT:

        A list of 2-tuples, each consisting of the given height,
        followed by a tuple of a-invariants of a curve of that height.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("F_1(4)")
            sage: B = C.next_height(3); B
            (3, [3, 3, 3, 3], [0, 1, 2, 3])
            sage: L = C._coeffs_from_height(B)
            sage: for ell in L: print(ell)
            (3, [0, -12, 0, 0, 36])
            (3, [0, 13, 0, 0, 36])
            (3, [0, -6, 0, 0, 9])
            (3, [0, 10, 0, 0, 9])
            (3, [0, 6, 0, 0, 9])
            (3, [0, 12, 0, 0, 36])
            (3, [0, -18, 0, 0, 81])
            (3, [0, 18, 0, 0, 81])
        """
        height = height_tuple[0]
        coeffs = height_tuple[1]
        index = height_tuple[2]

        # Produce list of all coefficient tuples with given height
        L = []
        for S in list(powerset(index))[1:]:
            B = []
            for j in range(len(coeffs)):
                if j in S:
                    B.append([-coeffs[j], coeffs[j]])
                elif j in index:
                    B.append(srange(-coeffs[j]+1, coeffs[j]))
                else:
                    B.append(srange(-coeffs[j], coeffs[j]+1))
            C = CartesianProduct(*B).list()
            for c in C:
                # This family has the additional constraint that
                # c[0]*c[1] = c[2]*c[3]
                if c[0]*c[1] == c[2]*c[3]:
                    L.append(c)

        # Convert coefficient tuples to a-invariants
        L2 = []
        for c in L:
            C = (height, self._coeffs_to_a_invariants(c))
            if not self._is_singular(C[1]):
                # This family produces duplicate curves
                if not C in L2:
                    L2.append(C)
        return L2


class CurveEnumeratorF_12x2(CurveEnumeratorF_22):
    """
    Height iterator for family with high incidence of rank one curves with
    full two torsion.

    Family:             Rank One + Four Torsion
    Model:              Y^2 = X^3 + A*X^2 + B, where
                        A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                        B = w0*w1*w2*w3
    Coefficients:       [w0,w1,w2,w3], with w0 + w1 = w2 + w3
    Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
    """
    def __init__(self):
        """
        Creates an instance of self.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.curve_enumerator import CurveEnumeratorF_12x2
            sage: C = CurveEnumeratorF_12x2(); C
            Height iterator for elliptic curves over QQ
            Family:             Rank One + Full Two Torsion
            Model:              Y^2 = X^3 + A*X^2 + B, where
                                A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                                B = w0*w1*w2*w3
            Coefficients:       [w0,w1,w2,w3], with w0 + w1 = w2 + w3
            Height function:    H = min{|w0|,|w1|,|w2|,|w3|}
        """
        # for __repr__()
        self._name = "Rank One + Full Two Torsion"
        self._model = "Y^2 = X^3 + A*X^2 + B, where \n" +\
                      (" " * 10) + "A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and \n" +\
                      (" " * 10) + "B = w0*w1*w2*w3"
        self._coeff_names = "[w0,w1,w2,w3], with w0 + w1 = w2 + w3"
        self._height_function = "H = min{|w0|,|w1|,|w2|,|w3|}"

        # The following constants must be Sage Integers
        self._num_coeffs = ZZ(4)
        self._pows = (ZZ(1), ZZ(1), ZZ(1), ZZ(1))

        # Coefficients, when converted to a-invariants, will produce
        # duplicate curves
        self._duplicates = True

    def _coeffs_from_height(self, height_tuple):
        """
        Returns a list of tuples of a-invariants of all curves
         described by height_tuple.

        INPUT:

        - ``height_tuple`` -- A tuple of the form
          (H, C, I) such that
          H: The smallest height >= N
          C: A list of coefficients for curves of this height
          I: A list of indices indicating which of the above coefficients
          achieve this height. The remaining values in C  indicate the
          max absolute value those coefficients are allowed to obtain
          without altering the height.

          For example, the tuple (4, [1, 2], [1]) for the short Weierstrass
          case denotes set of curves with height 4; these are all of the
          form Y^2 = X^3 + A*X + B, where B=2 and A ranges between -1 and 1.

        OUTPUT:

        A list of 2-tuples, each consisting of the given height,
        followed by a tuple of a-invariants of a curve of that height.

        EXAMPLES::

            sage: C = EllipticCurveEnumerator("F_1(2x2)")
            sage: B = C.next_height(3); B
            (3, [3, 3, 3, 3], [0, 1, 2, 3])
            sage: L = C._coeffs_from_height(B)
            sage: for ell in L: print(ell)
            (3, [0, -7, 0, 0, 12])
            (3, [0, 2, 0, 0, -3])
            (3, [0, 8, 0, 0, 12])
            (3, [0, 13, 0, 0, 36])
            (3, [0, 10, 0, 0, 9])
            (3, [0, -12, 0, 0, 36])
            (3, [0, -6, 0, 0, 9])
            (3, [0, 6, 0, 0, 9])
            (3, [0, 12, 0, 0, 36])
            (3, [0, -18, 0, 0, 81])
            (3, [0, 18, 0, 0, 81])
        """
        height = height_tuple[0]
        coeffs = height_tuple[1]
        index = height_tuple[2]

        # Produce list of all coefficient tuples with given height
        L = []
        for S in list(powerset(index))[1:]:
            B = []
            for j in range(len(coeffs)):
                if j in S:
                    B.append([-coeffs[j], coeffs[j]])
                elif j in index:
                    B.append(srange(-coeffs[j]+1, coeffs[j]))
                else:
                    B.append(srange(-coeffs[j], coeffs[j]+1))
            C = CartesianProduct(*B).list()
            for c in C:
                # This family has the additional constraint that
                # c[0]+c[1] = c[2]+c[3]
                if c[0]+c[1] == c[2]+c[3]:
                    L.append(c)

        # Convert coefficient tuples to a-invariants
        L2 = []
        for c in L:
            C = (height, self._coeffs_to_a_invariants(c))
            if not self._is_singular(C[1]):
                # This family produces duplicate curves
                if not C in L2:
                    L2.append(C)
        return L2


def EllipticCurveEnumerator(family):
    r"""
    Return the correct CurveEnumerator family. This instance will allow
    enumeration of all elliptic curves with a given height range, so that
    values of associated invariants (and averages thereof) can be quickly
    computed.

    INPUT:

    - ``family`` -- string; the family of curves being considered
      Current options are

      * 'short_weierstrass' ::

          Family:             Short Weierstrass
          Model:              Y^2 = X^3 + A*X + B
          Coefficients:       [A,B]
          Height function:    H = min{|A|^3,|B|^2}

      * 'rank_one' ::

          Family:             Rank One
          Model:              "Y^2 + A*Y = X^3 + B*X^2 + C*X"
          Coefficients:       [A,B,C]
          Height function:    H = min{|A|^6,|B|^4,|C|^3}

      * 'rank_two' ::

          Family:             Rank Two
          Model:              Y^2 = X^3 - 27*I*X - 27*J, where
                              I = 3*a4^2 + b4^2 - 3*a2*a6, and
                              J = -27/4*a2^2*a4^2 + 18*a4^2*b4 - 2*b4^3 + 9*a2*b4*a6 - 27*a6^2
          Coefficients:       [a2,a4,b4,c6]
          Height function:    H = min{|a2|^6,|a4|^3,|b4|^3,|c6|^2}

      * 'two_torsion' ::

          Family:             Two Torsion
          Model:              Y^2  = X^3 + A*X^2 + B*X
          Coefficients:       [A,B]
          Height function:    H = min{|A|^6,|C|^3}

      * 'three_torsion' ::

          Family:             Three Torsion
          Model:              Y^2 + A*X*Y + B*Y = X^3
          Coeddicients:       [A,B]
          Height function:    H = min{|A|^12,|B|^4}

      * 'F_1(2)' ::

          Family:             Rank One + Two Torsion
          Model:              Y^2 = X^3 + (A^2-B-C)*X^2 + (B*C)*X
          Coefficients:       [A,B,C]
          Height function:    H = min{|A|,|B|^2,|C|^2}

      * 'F_1(3)' ::

          Family:             Rank One + Three Torsion
          Model:              Y^2 + A*X*Y + B*Y = X^3
                              A = w0 + w1 + w2, and
                              B = w0*w1*w2
          Coefficients:       [w0,w1,w2]
          Height function:    H = min{|w0|,|w1|,|w2|}

      * 'F_1(4)' ::

          Family:             Rank One + Four Torsion
          Model:              Y^2 = X^3 + A*X^2 + B, where
                              A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                              B = w0*w1*w2*w3
          Coefficients:       [w0,w1,w2,w3], with w0*w1 = w2*w3
          Height function:    H = min{|w0|,|w1|,|w2|,|w3|}

      * 'F_1(2x2)' ::

          Family:             Rank One + Full Two Torsion
          Model:              Y^2 = X^3 + A*X^2 + B, where
                              A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                              B = w0*w1*w2*w3
          Coefficients:       [w0,w1,w2,w3], with w0 + w1 = w2 + w3
          Height function:    H = min{|w0|,|w1|,|w2|,|w3|}

      * 'F_2(2)' ::

          Family:             Rank Two + Two Torsion
          Model:              Y^2 = X^3 + A*X^2 + B, where
                              A = 1/4*(2*(w0^2+w1^2+w2^2+w3^2) - (w0+w1+w2+w3)^2), and
                              B = w0*w1*w2*w3
          Coefficients:       [w0,w1,w2,w3] (no constraints on w0,w1,w2,w3)
          Height function:    H = min{|w0|,|w1|,|w2|,|w3|}


    - To be implemented in future:

      * 'full_weierstrass' ::

          Family:             Full Weierstrass
          Model:              Y^2 + a1*X*Y + a3*Y = X^3 + a2*X^2 + a4*X + a6
          Coefficients:       [a1,a2,a3,a4,a6]

    OUTPUT:

    An instance of the relevant child class of CurveEnumeratorAbstract.

    EXAMPLES::

        sage: C = EllipticCurveEnumerator(family="short_weierstrass"); C
        Height iterator for elliptic curves over QQ
        Family:             Short Weierstrass
        Model:              Y^2 = X^3 + A*X + B
        Coefficients:       [A,B]
        Height function:    H = min{|A|^3,|B|^2}

        sage: C = EllipticCurveEnumerator(family='full_weierstrass')
        Traceback (most recent call last):
        ...
        NotImplementedError: Family not yet implemented.

        sage: C = EllipticCurveEnumerator(family="condorcet")
        Traceback (most recent call last):
        ...
        ValueError: 'family' must be a recognized Weierstrass family of elliptic curves.
    """
# This code will work once all classes are implemented
#    D = dict({"short_weierstrass": CurveEnumeratorShortWeierstrass(), \
#              "full_weierstrass": CurveEnumeratorFullWeierstrass(), \
#              "rank_one": CurveEnumeratorRankOne(), \
#              "rank_two": CurveEnumeratorRankTwo(), \
#              "two_torsion": CurveEnumeratorTwoTorsion(), \
#              "three_torsion": CurveEnumeratorThreeTorsion(), \
#              "F_1(2)": CurveEnumeratorF_12(), \
#              "F_1(3)": CurveEnumeratorF_13(), \
#              "F_1(4)": CurveEnumeratorF_14(), \
#              "F_1(2x2)": CurveEnumeratorF_12x2(), \
#              "F_2(2)": CurveEnumeratorF_22()})
#
#    try:
#        return D[family]
#    except:
#        raise ValueError("'family' must be a recognized Weierstrass "
#                         "family of elliptic curves.")

    # This can be rewritten to be more elegant
    if family == "short_weierstrass":
        return CurveEnumeratorShortWeierstrass()
    elif family == "full_weierstrass":
        raise NotImplementedError("Family not yet implemented.")
        # return CurveEnumeratorFullWeierstrass()
    elif family == "rank_one":
        return CurveEnumeratorRankOne()
    elif family == "rank_two":
        return CurveEnumeratorRankTwo()
    elif family == "two_torsion":
        return CurveEnumeratorTwoTorsion()
    elif family == "three_torsion":
        return CurveEnumeratorThreeTorsion()
    elif family == "F_1(2)":
        return CurveEnumeratorF_12()
    elif family == "F_1(3)":
        return CurveEnumeratorF_13()
    elif family == "F_1(4)":
        return CurveEnumeratorF_14()
    elif family == "F_1(2x2)":
        return CurveEnumeratorF_12x2()
    elif family == "F_2(2)":
        return CurveEnumeratorF_22()
    else:
        raise ValueError("'family' must be a recognized Weierstrass "
                         "family of elliptic curves.")
