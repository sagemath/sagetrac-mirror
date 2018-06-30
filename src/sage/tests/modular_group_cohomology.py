r"""
Tests for the optional ``p_group_cohomology`` package.

AUTHOR:

- Simon King (July 2018, see :trac:`18514`)

TESTS::

    sage: from pGroupCohomology import CohomologyRing   # optional - p_group_cohomology
    sage: CohomologyRing.set_workspace(tmp_dir())       # optional - p_group_cohomology
    sage: H0 = CohomologyRing(64,14)                    # optional - p_group_cohomology

Computation of a modular cohomology ring of a prime power group in characteristic 2,
and comparison with stored result in a database::

    sage: H = CohomologyRing(64,14,from_scratch=True)   # optional - p_group_cohomology
    sage: H.make()                                      # optional - p_group_cohomology
    sage: print(H)                                      # optional - p_group_cohomology
    <BLANKLINE>
    Cohomology ring of Small Group number 14 of order 64 with coefficients in GF(2)
    <BLANKLINE>
    Computation complete
    Minimal list of generators:
    [a_2_1: 2-Cocycle in H^*(SmallGroup(64,14); GF(2)),
     c_2_2: 2-Cocycle in H^*(SmallGroup(64,14); GF(2)),
     c_4_4: 4-Cocycle in H^*(SmallGroup(64,14); GF(2)),
     a_1_0: 1-Cocycle in H^*(SmallGroup(64,14); GF(2)),
     a_1_1: 1-Cocycle in H^*(SmallGroup(64,14); GF(2)),
     a_3_3: 3-Cocycle in H^*(SmallGroup(64,14); GF(2))]
    Minimal list of algebraic relations:
    [a_1_0^2,
     a_1_0*a_1_1,
     a_1_1^3,
     a_2_1*a_1_0,
     a_2_1^2+a_2_1*a_1_1^2,
     a_1_1*a_3_3+a_2_1^2,
     a_1_0*a_3_3+a_2_1^2,
     a_2_1*a_3_3,
     a_3_3^2]
    sage: H == H0                                       # optional - p_group_cohomology
    True

Some algebraic constructions in the cohomology ring::

    sage: H.essential_ideal()                           # optional - p_group_cohomology
    a_2_1*a_1_1^2
    sage: ascii_art(H.bar_code('LowerCentralSeries')[2])# optional - p_group_cohomology
            *
          *
          *
        *-*-*
        *-*-*
      *
    *-*-*
    *-*
    *

Computation of a modular cohomology ring of a non prime power group in
characteristic 2::

    sage: H = CohomologyRing(gap(AlternatingGroup(7)), GroupName="A(7)", prime=2, from_scratch=True) # optional - p_group_cohomology
    sage: H.make()                                      # optional - p_group_cohomology
    sage: print(H)                                      # optional - p_group_cohomology
    <BLANKLINE>
    Cohomology ring of A(7) with coefficients in GF(2)
    <BLANKLINE>
    Computation complete
    Minimal list of generators:
    [c_2_0: 2-Cocycle in H^*(A(7); GF(2)),
     b_3_0: 3-Cocycle in H^*(A(7); GF(2)),
     b_3_1: 3-Cocycle in H^*(A(7); GF(2))]
    Minimal list of algebraic relations:
    [b_3_0*b_3_1]

Computing a modular cohomology ring of a non prime power group in odd
characteristic::

    sage: H = CohomologyRing(gap(AlternatingGroup(6)), GroupName="A(6)", prime=3, from_scratch=True) # optional - p_group_cohomology
    sage: H.make()                                      # optional - p_group_cohomology
    sage: print(H)                                      # optional - p_group_cohomology
    <BLANKLINE>
    Cohomology ring of A(6) with coefficients in GF(3)
    <BLANKLINE>
    Computation complete
    Minimal list of generators:
    [a_2_0: 2-Cocycle in H^*(A(6); GF(3)),
     c_4_0: 4-Cocycle in H^*(A(6); GF(3)),
     c_8_1: 8-Cocycle in H^*(A(6); GF(3)),
     c_8_2: 8-Cocycle in H^*(A(6); GF(3)),
     a_3_0: 3-Cocycle in H^*(A(6); GF(3)),
     a_3_1: 3-Cocycle in H^*(A(6); GF(3)),
     a_7_2: 7-Cocycle in H^*(A(6); GF(3)),
     a_7_3: 7-Cocycle in H^*(A(6); GF(3))]
    Minimal list of algebraic relations:
    [a_2_0^2,
     a_2_0*a_3_0,
     a_2_0*a_3_1,
     a_3_0*a_3_1+a_2_0*c_4_0,
     a_2_0*a_7_2,
     a_2_0*a_7_3,
     a_3_0*a_7_2+a_2_0*c_8_2,
     a_3_0*a_7_3+a_2_0*c_8_1,
     a_3_1*a_7_2-a_2_0*c_8_1,
     a_3_1*a_7_3+a_2_0*c_8_2,
     c_8_2*a_3_0-c_8_1*a_3_1+c_4_0*a_7_3,
     c_8_2*a_3_1+c_8_1*a_3_0-c_4_0*a_7_2,
     a_7_2*a_7_3+a_2_0*c_4_0*c_8_2,
     c_8_2*a_7_2+c_8_1*a_7_3+c_4_0*c_8_1*a_3_0-c_4_0^2*a_7_2,
     c_8_2*a_7_3-c_8_1*a_7_2+c_4_0*c_8_1*a_3_1-c_4_0^2*a_7_3,
     c_8_2^2+c_8_1^2-c_4_0^2*c_8_2]

Computing higher algebraic structures in the cohomology ring::

    sage: H.5.massey_power()                            # optional - p_group_cohomology
    -c_8_2+c_4_0^2: 8-Cocycle in H^*(A(6); GF(3))

"""