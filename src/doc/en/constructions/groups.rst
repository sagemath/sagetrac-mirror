.. _chapter-groups:

******
Groups
******

.. index::
   pair: group; permutation

.. _section-permutation:

Permutation groups
==================

A permutation group is a subgroup of some symmetric group
:math:`S_n`. Sage has a Python class ``PermutationGroup``, so you
can work with such groups directly::

    sage: G = PermutationGroup(['(1,2,3)(4,5)'])
    sage: G
    Permutation Group with generators [(1,2,3)(4,5)]
    sage: g = G.gens()[0]; g
    (1,2,3)(4,5)
    sage: g*g
    (1,3,2)
    sage: G = PermutationGroup(['(1,2,3)'])
    sage: g = G.gens()[0]; g
    (1,2,3)
    sage: g.order()
    3

For the example of the Rubik's cube group (a permutation subgroup
of :math:`S_{48}`, where the non-center facets of the Rubik's
cube are labeled :math:`1,2,...,48` in some fixed way), you can
use the GAP-Sage interface as follows.

.. index::
   pair: group; Rubik's cube

.. skip

::

    sage: cube = "cubegp := Group(
    ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
    ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
    (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
    (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
    (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
    (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"
    sage: gap(cube)
    'permutation group with 6 generators'
    sage: gap("Size(cubegp)")
    43252003274489856000'

Another way you can choose to do this:

-  Create a file ``cubegroup.py`` containing the
   lines::

       cube = "cubegp := Group(
       ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
       ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
       (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
       (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
       (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
       (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"

   Then place the file in the subdirectory
   ``$SAGE_ROOT/local/lib/python2.4/site-packages/sage`` of your Sage home
   directory. Last, read (i.e., ``import``) it into Sage:

   .. skip

   ::

       sage: import sage.cubegroup
       sage: sage.cubegroup.cube
       'cubegp := Group(( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)
       (11,35,27,19),( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)
       ( 6,22,46,35),(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)
       ( 8,30,41,11),(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)
       ( 8,33,48,24),(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)
       ( 1,14,48,27),(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)
       (16,24,32,40) )'
       sage: gap(sage.cubegroup.cube)
       'permutation group with 6 generators'
       sage: gap("Size(cubegp)")
       '43252003274489856000'

   (You will have line wrap instead of the above carriage returns in
   your Sage output.)

-  Use the ``CubeGroup`` class::

       sage: rubik = CubeGroup()
       sage: rubik
       The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).
       sage: rubik.order()
       43252003274489856000

   (1) has implemented classical groups (such as :math:`GU(3,\GF{5})`)
   and matrix groups over a finite field with user-defined generators.

   (2) also has implemented finite and infinite (but finitely
   generated) abelian groups.

.. index::
   pair: group; conjugacy classes

.. _section-conjugacy:

Conjugacy classes
=================

You can compute conjugacy classes of a finite group using "natively"::

    sage: G = PermutationGroup(['(1,2,3)', '(1,2)(3,4)', '(1,7)'])
    sage: CG = G.conjugacy_classes_representatives()
    sage: gamma = CG[2]
    sage: CG; gamma
    [(), (4,7), (3,4,7), (2,3)(4,7), (2,3,4,7), (1,2)(3,4,7), (1,2,3,4,7)]
    (3,4,7)

You can use the Sage-GAP interface::

    sage: gap.eval("G := Group((1,2)(3,4),(1,2,3))")
    'Group([ (1,2)(3,4), (1,2,3) ])'
    sage: gap.eval("CG := ConjugacyClasses(G)")
    '[ ()^G, (2,3,4)^G, (2,4,3)^G, (1,2)(3,4)^G ]'
    sage: gap.eval("gamma := CG[3]")
    '(2,4,3)^G'
    sage: gap.eval("g := Representative(gamma)")
    '(2,4,3)'

Or, here's another (more "pythonic") way to do this type of computation::

    sage: G = gap.Group('[(1,2,3), (1,2)(3,4), (1,7)]')
    sage: CG = G.ConjugacyClasses()
    sage: gamma = CG[2]
    sage: g = gamma.Representative()
    sage: CG; gamma; g
    [ ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), () ), 
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (4,7) ), 
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (3,4,7) ), 
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (2,3)(4,7) ), 
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (2,3,4,7) ), 
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2)(3,4,7) ), 
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2,3,4,7) ) ]
    ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (4,7) )
    (4,7)

.. index::
   pair: group; normal subgroups

.. _section-normal:

Normal subgroups
================

If you want to find all the normal subgroups of a permutation group
:math:`G` (up to conjugacy), you can use Sage's interface to GAP::

    sage: G = AlternatingGroup( 5 )
    sage: gap(G).NormalSubgroups()
    [ Group( () ), AlternatingGroup( [ 1 .. 5 ] ) ]

or

::

    sage: G = gap("AlternatingGroup( 5 )")
    sage: G.NormalSubgroups()
    [ Group( () ), AlternatingGroup( [ 1 .. 5 ] ) ]

Here's another way, working more directly with GAP::

    sage: print(gap.eval("G := AlternatingGroup( 5 )"))
    Alt( [ 1 .. 5 ] )
    sage: print(gap.eval("normal := NormalSubgroups( G )"))
    [ Group(()), Alt( [ 1 .. 5 ] ) ]
    sage: G = gap.new("DihedralGroup( 10 )")
    sage: G.NormalSubgroups()
    [ Group( <identity> of ... ), Group( [ f2 ] ), Group( [ f1, f2 ] ) ]
    sage: print(gap.eval("G := SymmetricGroup( 4 )"))
    Sym( [ 1 .. 4 ] )
    sage: print(gap.eval("normal := NormalSubgroups( G );"))
    [ Group(()), Group([ (1,4)(2,3), (1,3)(2,4) ]), Group([ (2,4,3), (1,4)
      (2,3), (1,3)(2,4) ]), Sym( [ 1 .. 4 ] ) ]

.. index::
   pair: groups; center

.. _section-center:

Centers
=======

How do you compute the center of a group in Sage?

Although Sage calls GAP to do the computation of the group center,
``center`` is "wrapped" (i.e., Sage has a class PermutationGroup with
associated class method "center"), so the user does not need to use
the ``gap`` command. Here's an example::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.center()
    Subgroup of (Permutation Group with generators [(3,4), (1,2,3)(4,5)]) generated by [()]

A similar syntax for matrix groups also works::

    sage: G = SL(2, GF(5) )
    sage: G.center()
    Matrix group over Finite Field of size 5 with 1 generators (
    [4 0]
    [0 4]
    )
    sage: G = PSL(2, 5 )
    sage: G.center()
    Subgroup of (The projective special linear group of degree 2 over Finite Field of size 5) generated by [()]

.. NOTE:: ``center`` can be spelled either way in GAP, not so in Sage.

The group id database
=====================

The function ``group_id`` requires that the Small Groups Library of
E. A. O'Brien, B. Eick, and H. U. Besche be installed.  You can do
this by typing ``sage -i database_gap`` in the shell.

::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.order()
    120
    sage: G.group_id()      # optional - database_gap
    [120, 34]

Another example of using the small groups database: ``group_id``

.. skip

::

    sage: gap_console()
    GAP4, Version: 4.4.6 of 02-Sep-2005, x86_64-unknown-linux-gnu-gcc
    gap> G:=Group((4,6,5)(7,8,9),(1,7,2,4,6,9,5,3));
    Group([ (4,6,5)(7,8,9), (1,7,2,4,6,9,5,3) ])
    gap> StructureDescription(G);
    "(((C3 x C3) : Q8) : C3) : C2"

Construction instructions for every group of order less than 32
===============================================================

AUTHORS:

* Davis Shurbert

Every group of order less than 32 is implemented in Sage as a permutation
group. They can all be created easily. We will first show how to build direct
products and semidirect products, then give the commands necessary to build
all of these small groups. 

Let ``G1``, ``G2``, ..., ``Gn`` be permutation groups already initialized in
Sage. The following command can be used to take their direct product (where,
of course, the ellipses are simply being used here as a notation, and you
actually must enter every factor in your desired product explicitly).

.. skip

::

    sage: G = direct_product_permgroups([G1, G2, ..., Gn])

The semidirect product operation can be thought of as a generalization of the
direct product operation. Given two groups, `H` and `K`, their semidirect
product, `H \ltimes_{\phi} K`, (where `\phi : H \rightarrow Aut(K)` is a
homomorphism) is a group whose underlying set is the cartersian product of
`H` and `K`, but with the operation:

.. MATH::

    (h_1, k_1) (h_2, k_2) = (h_1 h_2, k_1^{\phi(h_2)} k_2).

The output is not the group explicity described in the definition of the
operation, but rather an isomorphic group of permutations. In the routine
below, assume ``H`` and ``K`` already have been defined and initialized in
Sage. Also, ``phi`` is a list containing two sublists that define the
underlying homomorphism by giving the images of a set of generators of ``H``.
For each semidirect product in the table below we will show you how to build
``phi``, then assume you have read this passage and understand how to go
from there.

.. skip

::

    sage: G = H.semidirect_product(K, phi)

To avoid unnecessary repitition, we will now give commands that one can use to
create the cyclic group of order `n`, `C_n`, and the dihedral group on `n`
letters, `D_n`. We will present one more example of each to ensure that the
reader understands the command, then it will be withheld.

.. skip

::

    sage: G = CyclicPermutationGroup(n)

    sage: G = DihedralGroup(n)

Note that exponential notation will be used for the direct product operation.
For example, `{C_2}^2 = C_2 \times C_2`. This table was crafted with the help
of *Group Tables*, by AD Thomas and GV Wood (1980, Shiva Publishing).


===== =============================================== =============================================================================================== ===========================
Order Group Description                                Command(s)                                                                                     GAP ID
===== =============================================== =============================================================================================== ===========================
1     The Trivial Group                               ::                                                                                              [1,1]

                                                        sage: G = SymmetricGroup(1)
2     `C_2`                                           ::                                                                                              [2,1]

                                                        sage: G = SymmetricGroup(2)
3     `C_3`                                           ::                                                                                              [3,1]

                                                        sage: G = CyclicPermutationGroup(3)
4     `C_4`                                                                                                                                           [4,1]
4     `C_2 \times C_2`                                ::                                                                                              [4,2]

                                                        sage: G = KleinFourGroup()
5     `C_5`                                                                                                                                           [5,1]
6     `C_6`                                                                                                                                           [6,2]
6     `S_3` (Symmetric Group on 3 letters)            ::                                                                                              [6,1]

                                                        sage: G = SymmetricGroup(3)
7     `C_7`                                                                                                                                           [7,1]
8     `C_8`                                                                                                                                           [8,1]
8     `C_4 \times C_2`                                                                                                                                [8,2]
8     `C_2\times C_2\times C_2`                                                                                                                       [8,5]
8     `D_4`                                           ::                                                                                              [8,3]

                                                        sage: G = DihedralGroup(4)
8     The Quaternion Group (Q)                        ::                                                                                              [8,4]

                                                        sage: G = QuaternionGroup()
9     `C_9`                                                                                                                                           [9,1]
9     `C_3 \times C_3`                                                                                                                                [9,2]
10    `C_{10}`                                                                                                                                        [10,2]
10    `D_5`                                                                                                                                           [10,1]
11    `C_{11}`                                                                                                                                        [11,1]
12    `C_{12}`                                                                                                                                        [12,2]
12    `C_6 \times C_2`                                                                                                                                [12,5]
12    `D_6`                                                                                                                                           [12,4]
12    `A_4` (Alternating Group on 4 letters)          ::                                                                                              [12,3]

                                                        sage: G = AlternatingGroup(4)
12    `Q_6` (DiCyclic group of order 12)              ::                                                                                              [12,1]

                                                        sage: G = DiCyclicGroup(3)
13    `C_{13}`                                                                                                                                        [13,1]
14    `C_{14}`                                                                                                                                        [14,2]
14    `D_{7}`                                                                                                                                         [14,1]
15    `C_{15}`                                                                                                                                        [15,1]
16    `C_{16}`                                                                                                                                        [16,1]
16    `C_8 \times C_2`                                                                                                                                [16,5]
16    `C_4 \times C_4`                                                                                                                                [16,2]
16    `C_4\times C_2\times C_2`                                                                                                                       [16,10]
16    `{C_2}^4`                                                                                                                                       [16,14]
16    `D_4 \times C_2`                                                                                                                                [16,11]
16    `Q \times C_2`                                                                                                                                  [16,12]
16    `D_8`                                                                                                                                           [16,7]
16    `Q_{8}` (Dicyclic group of order 16)            ::                                                                                              [16,9]

                                                        sage: G = DiCyclicGroup(4)
16    Semidihedral Group of order `2^4`               ::                                                                                              [16,8]

                                                        sage: G = SemidihedralGroup(4)
16    Split Metacyclic Group of order `2^4`           ::                                                                                              [16,6]

                                                        sage: G = SplitMetacyclicGroup(2,4)
16    `(C_4 \times C_2) \rtimes_{\phi} C_2`           ::                                                                                              [16,13]

                                                        sage: C2 = SymmetricGroup(2); C4 = CyclicPermutationGroup(4)
                                                        sage: A = direct_product_permgroups([C2,C4])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0],A.gens()[0]^2*A.gens()[1]])
                                                        sage: phi = [[(1,2)],[alpha]]
16    `(C_4 \times C_2) \rtimes_{\phi} C_2`           ::                                                                                              [16,3]

                                                        sage: C2 = SymmetricGroup(2); C4 = CyclicPermutationGroup(4)
                                                        sage: A = direct_product_permgroups([C2,C4])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0]^3*A.gens()[1],A.gens()[1]])
                                                        sage: phi = [[(1,2)],[alpha]]
16    `C_4 \rtimes_{\phi} C_4`                        ::                                                                                              [16,4]

                                                        sage: C4 = CyclicPermutationGroup(4)
                                                        sage: alpha = PermutationGroupMorphism(C4,C4,[C4.gen().inverse()])
                                                        sage: phi = [[(1,2,3,4)],[alpha]]
17    `C_{17}`                                                                                                                                        [17,1]
18    `C_{18}`                                                                                                                                        [18,2]
18    `C_6 \times C_3`                                                                                                                                [18,5]
18    `D_9`                                                                                                                                           [18,1]
18    `S_3 \times C_3`                                                                                                                                [18,3]
18    `Dih(C_3 \times C_3)`                           ::                                                                                              [18,4]

                                                        sage: G = GeneralDihedralGroup([3,3])
19    `C_{19}`                                                                                                                                        [19,1]
20    `C_{20}`                                                                                                                                        [20,2]
20    `C_{10} \times C_2`                                                                                                                             [20,5]
20    `D_{10}`                                                                                                                                        [20,4]
20    `Q_{10}` (Dicyclic Group of order 20)                                                                                                           [20,1]
20    `Hol(C_5)`                                      ::                                                                                              [20,3]

                                                        sage: C5 = CyclicPermutationGroup(5)
                                                        sage: G = C5.holomorph()
21    `C_{21}`                                                                                                                                        [21,2]
21    `C_7 \rtimes_{\phi} C_3`                        ::                                                                                              [21,1]

                                                        sage: C7 = CyclicPermutationGroup(7)
                                                        sage: alpha = PermutationGroupMorphism(C7,C7,[C7.gen()**4])
                                                        sage: phi = [[(1,2,3)],[alpha]]
22    `C_{22}`                                                                                                                                        [22,2]
22    `D_{11}`                                                                                                                                        [22,1]
23    `C_{23}`                                                                                                                                        [23,1]
24    `C_{24}`                                                                                                                                        [24,2]
24    `D_{12}`                                                                                                                                        [24,6]
24    `Q_{12}` (DiCyclic Group of order 24)                                                                                                           [24,4]
24    `C_{12} \times C_2`                                                                                                                             [24,9]
24    `C_6 \times C_2 \times C_2`                                                                                                                     [24,15]
24    `S_4` (Symmetric Group on 4 letters)            ::                                                                                              [24,12]

                                                        sage: G = SymmetricGroup(4)
24    `S_3 \times C_4`                                                                                                                                [24,5]
24    `S_3 \times C_2 \times C_2`                                                                                                                     [24,14]
24    `D_4 \times C_3`                                                                                                                                [24,10]
24    `Q \times C_3`                                                                                                                                  [24,11]
24    `A_4 \times C_2`                                                                                                                                [24,13]
24    `Q_6 \times C_2`                                                                                                                                [24,7]
24    `Q \rtimes_{\phi} C_3`                          ::                                                                                              [24,3]

                                                        sage: Q = QuaternionGroup()
                                                        sage: alpha = PermutationGroupMorphism(Q,Q,[Q.gens()[0]*Q.gens()[1],Q.gens()[0].inverse()])
                                                        sage: phi = [[(1,2,3)],[alpha]]
24    `C_3 \rtimes_{\phi} C_8`                        ::                                                                                              [24,1]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: alpha = PermutationGroupMorphism(C3,C3,[C3.gen().inverse()])
                                                        sage: phi = [[(1,2,3,4,5,6,7,8)],[alpha]]
24    `C_3 \rtimes_{\phi} D_4`                        ::                                                                                              [24,8]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: alpha1 = PermutationGroupMorphism(C3,C3,[C3.gen().inverse()])
                                                        sage: alpha2 = PermutationGroupMorphism(C3,C3,[C3.gen()])
                                                        sage: phi = [[(1,2,3,4),(1,3)],[alpha1,alpha2]]
25    `C_{25}`                                                                                                                                        [25,1]
25    `C_5 \times C_5`                                                                                                                                [25,2]
26    `C_{26}`                                                                                                                                        [26,2]
26    `D_{13}`                                                                                                                                        [26,1]
27    `C_{27}`                                                                                                                                        [27,1]
27    `C_9 \times C_3`                                                                                                                                [27,2]
27    `C_3 \times C_3 \times C_3`                                                                                                                     [27,5]
27    Split Metacyclic Group of order `3^3`           ::                                                                                              [27,4]

                                                        sage: G = SplitMetacyclicGroup(3,3)
27    `(C_3 \times C_3) \rtimes_{\phi} C_3`           ::                                                                                              [27,3]

                                                        sage: C3 = CyclicPermutationGroup(3)
                                                        sage: A = direct_product_permgroups([C3,C3])
                                                        sage: alpha = PermutationGroupMorphism(A,A,[A.gens()[0]*A.gens()[1].inverse(),A.gens()[1]])
                                                        sage: phi = [[(1,2,3)],[alpha]]
28    `C_{28}`                                                                                                                                        [28,2]
28    `C_{14} \times C_2`                                                                                                                             [28,4]
28    `D_{14}`                                                                                                                                        [28,3]
28    `Q_{14}` (DiCyclic Group of order 28)                                                                                                           [28,1]
29    `C_{29}`                                                                                                                                        [29,1]
30    `C_{30}`                                                                                                                                        [30,4]
30    `D_{15}`                                                                                                                                        [30,3]
30    `D_5 \times C_3`                                                                                                                                [30,2]
30    `D_3 \times C_5`                                                                                                                                [30,1]
31    `C_{31}`                                                                                                                                        [31,1]
===== =============================================== =============================================================================================== ===========================


Table By Kevin Halasz

Construction instructions for every finitely presented group of order 15 or less
================================================================================

Sage has the capability to easily construct every group of order 15 or less
as a finitely presented group. We will begin with some discussion on creating
finitely generated abelian groups, as well as direct and semidirect products
of finitely presented groups.

All finitely generated abelian groups can be created using the
``groups.presentation.FGAbelian(ls)`` command, where ``ls`` is a list of
non-negative integers which gets reduced to invariants defining the group
to be returned. For example, to construct
`C_4 \times C_2 \times C_2 \times C_2` we can simply use::

    sage: A = groups.presentation.FGAbelian([4,2,2,2])

The output for a given group is the same regardless of the input list of
integers.  The following example yeilds identical presentations for the
cyclic group of order 30.
::

    sage: A = groups.presentation.FGAbelian([2,3,5])
    sage: B = groups.presentation.FGAbelian([30])

If ``G`` and ``H`` are finitely presented groups, we can use the following
code to create the direct product of ``G`` and ``H``, `G \times H`.

.. skip

::

    sage: D = G.direct_product(H)

Suppose there exists a homomorphism `\phi` from a group `G` to the
automorphism group of a group `H`. Define the semidirect product of `G`
with `H` via `\phi`, as the Cartesian product of `G` and `H`, with the
operation `(g_1, h_1)(g_2, h_2) = (g_1 g_2, \phi_{h_1}(g_2) h_2)` where
`\phi_h = \phi(h)`. To construct this product in Sage for two finitely
presented groups, we must define `\phi` manually using a pair of lists. The
first list consists of generators of the group `G`, while the second list
consists of images of the corresponding generators in the first list. These
automorphisms are similarly defined as a pair of lists, generators in one
and images in the other. As an example, we construct the dihedral group of
order 16 as a semidirect product of cyclic groups.
::

    sage: C2 = groups.presentation.Cyclic(2)
    sage: C8 = groups.presentation.Cyclic(8)
    sage: hom = (C2.gens(), [ ([C8([1])], [C8([-1])]) ])
    sage: D = C2.semidirect_product(C8, hom)

The following table shows the groups of order 15 or less, and how to construct
them in Sage. Repeated commands have been omitted but instead are described
by the following exmples.

The cyclic group of order `n` can be crated with a single command:

.. skip

::

    sage: C = groups.presentation.Cyclic(n)

Similarly for the dihedral group of order `2n`:

.. skip

::

    sage: D = groups.presentation.Dihedral(n)
 
This table was modeled after the preceding table created by Kevin Halasz. 


===== =============================================== =============================================================================================== =========================== 
Order Group Description                                Command(s)                                                                                     GAP ID 
===== =============================================== =============================================================================================== =========================== 
1     The Trivial Group                               ::                                                                                              [1,1] 

                                                        sage: G = groups.presentation.Symmetric(1) 

2     `C_2`                                           ::                                                                                              [2,1] 

                                                        sage: G = groups.presentation.Symmetric(2)

3     `C_3`                                           ::                                                                                              [3,1] 

                                                        sage: G = groups.presentation.Cyclic(3) 

4     `C_4`                                                                                                                                           [4,1] 

4     `C_2 \times C_2`                                ::                                                                                              [4,2] 

                                                        sage: G = groups.presentation.Klein() 

5     `C_5`                                                                                                                                           [5,1] 
6     `C_6`                                                                                                                                           [6,2] 

6     `S_3` (Symmetric Group on 3 letters)            ::                                                                                              [6,1] 

                                                        sage: G = groups.presentation.Symmetric(3) 

7     `C_7`                                                                                                                                           [7,1] 
8     `C_8`                                                                                                                                           [8,1] 

8     `C_4 \times C_2`                                ::                                                                                              [8,2]

                                                        sage: G = groups.presentation.FGAbelian([4,2])

8     `C_2\times C_2\times C_2`                       ::                                                                                              [8,5] 

                                                        sage: G = groups.presentation.FGAbelian([2,2,2])

8     `D_4`                                           ::                                                                                              [8,3] 

                                                        sage: G = groups.presentation.Dihedral(4)
 
8     The Quaternion Group (Q)                        ::                                                                                              [8,4] 

                                                        sage: G = groups.presentation.Quaternion() 

9     `C_9`                                                                                                                                           [9,1] 
9     `C_3 \times C_3`                                                                                                                                [9,2] 
10    `C_{10}`                                                                                                                                        [10,2] 
10    `D_5`                                                                                                                                           [10,1] 
11    `C_{11}`                                                                                                                                        [11,1] 
12    `C_{12}`                                                                                                                                        [12,2] 
12    `C_6 \times C_2`                                                                                                                                [12,5] 
12    `D_6`                                                                                                                                           [12,4] 
12    `A_4` (Alternating Group on 4 letters)          ::                                                                                              [12,3] 

                                                        sage: G = groups.presentation.Alternating(4) 

12    `Q_6` (DiCyclic group of order 12)              ::                                                                                              [12,1] 
       
                                                        sage: G = groups.presentation.DiCyclic(3)
 
13    `C_{13}`                                                                                                                                        [13,1] 
14    `C_{14}`                                                                                                                                        [14,2] 
14    `D_{7}`                                                                                                                                         [14,1] 
15    `C_{15}`                                                                                                                                        [15,1]
===== =============================================== =============================================================================================== ===========================


.. index::
   pair: group; free groups

.. _section-free groups:

Free Groups
===========
The Train-track package was first written by Thierry Coulbois and
received contributions by Matt Clay, Brian Mann and others.

It is primarily intended to implement the computation of a train-track
representative for automorphisms of free groups as introduced by
M.Bestvina and M.Handel [1].

-  Free groups and automorphisms::

   To create ``FreeGroup`` on generators:

   -  Creating free groups::

      The first need to create a ``FreeGroup``. It can be specified by
      its rank or a list of letters::

          sage: F = FreeGroup('a, b');  F
          Free Group on generators {a, b}
          sage: F.rank()
          2
          sage: H = FreeGroup(3, 'x')
          sage: H
          Free Group on generators {x0, x1, x2}
          sage: H.rank()
          3

      Note that they are reduced by default.
      Words can be multiplied and inverted easily::

          sage: w=F('abA')
          sage: w*w
          a*b^2*a^-1
          sage: w.inverse()
          a*b^-1*a^-1
          sage: w**5
          a*b^5*a^-1

   -  Free group automorphisms::

      The creation (and the parsing) of free group automorphisms relies on that of
      substitutions. Most of what you might expect should correctly create a
      free group automorphism::

          sage: phi = FreeGroupAutomorphism('a->ab,b->ac,c->a')
          sage: psi = FreeGroupAutomorphism('a->c,b->ba,c->bcc')
          sage: print(phi*psi)
          a->a,b->a*c*a*b,c->a*c*a^2
          sage: print(phi.inverse())
          a->c,b->c^-1*a,c->c^-1*b
          sage: print(phi**3)
          a->a*b*a*c*a*b*a,b->a*b*a*c*a*b,c->a*b*a*c
          sage: phi('aBc')
          a*b*c^-1

      There is a list of pre-defined automorphisms of free groups taken from the litterature:
      ::
          sage: print free_group_automorphisms.Handel_Mosher_inverse_with_same_lambda()
          Automorphism of the Free Group on generators {a, b, c}: a->b^-1*a*(c*a^-1)^2*b*c^-1*a*(b^-1*a*c)^2*a^-1*c*a^-1*b,b->b^-1*a*(c*a^-1)^2*b*c^-1*a*b^-1*a*c,c->b^-1*a*(c*a^-1)^2*b*c^-1*a

      Also Free group automorphisms can be obtained as composition of
      elementary Nielsen automorphisms  (of the form $a->ab$). Up to now they
      are rather called Dehn twists.

      If the free group as even rank $N=2g$, then it is the fundamental
      group of an oriented surface of genus $g$ with one boundary
      component. In this case the mapping class group\index{mapping class group}
      of $S_{g,1}$ is a
      subgroup of the outer automorphism group of $F_N$ and it is generated
      by a collection of $3g-1$ Dehn twists along curves. Those Dehn
      twists are accessed through:
      ::
          sage: F = FreeGroup(4, 'a,b,c,d')
          sage: FreeGroupAutomorphism.surface_dehn_twist(F, k=2)
          Automorphism of the Free Group on generators {a, b, c, d}: a->a,b->a*b,c->a*c*a^-1,d->a*d*a^-1

      Similarly the braid group $B_N$ is a subgroup of $\Aut(\FN)$ and its
      usual generators are obtained by:
      ::

          sage: F = FreeGroup(4,'a,b,c,d')
          sage:FreeGroupAutomorphism.braid_automorphism(F, 2)
          Automorphism of the Free Group on generators {a, b, c, d}: a->a,b->b*c*b^-1,c->b,d->d

-  Graphs and maps::

   Graphs and maps are used to represent free group automorphisms. A
   graph here is a ``GraphWithInverses``: it has a set of vertices and a set
   of edges in one-to-one correspondance with the letters of an
   ``AlphabetWithInverses``: each non-oriented edge is a pair $\{e,\bar e\}$
   of a letter of the alphabet and its inverse. This is complient with
   Serre's view [5,6]. As the alphabet has a set of positive letters there is a
   default choice of orientation for edges.  first
   The easiest graph is the rose::

       sage: from sage.groups.free_groups.inverse_graph import GraphWithInverses
       sage: A = AlphabetWithInverses(3)
       sage: G = GraphWithInverses.rose_graph(A)
       sage: print G
       a: 0->0, b: 0->0, c: 0->0
       sage: G.plot()

   .. image:: tmp_S25CbR.*


   Otherwise a graph can be given by a variety of inputs like a list of
   edges, etc. Graphs can easily be plotted. Note that
   ``plot()`` tries to lower the number of accidental crossing of
   edges, using some thermodynamics and randomness, thus two calls of
   ``plot()`` may output two different figures.

   A number of operations on graphs are defined: subdividing, folding,
   collapsing edges, etc. But, as of now, not all
   Stallings [7] moves are implemented.

   Graphs come with maps between them: a map is a continuous map from a
   graph to another which maps vertices to vertices and edges to
   edge-paths. Again they can be given by a variety of means. As Graph
   maps are intended to represent free group automorphisms a simple way
   to create a graph map is from a free group automorphism::

       sage: phi = free_group_automorphisms.tribonacci()
       sage: print(phi.rose_representative())
       GraphSelfMap:
       Marked graph: a: 0->0, b: 0->0, c: 0->0
       Marking: a->a, b->b, c->c
       Edge map: a->ab, b->ac, c->a

   Remark that by default the rose graph is marked: it comes
   with a marking from the rose (itself, but you should think of that one
   as fixed) to the graph. Here the graph map is a graph self map as the
   source and the target are the same.

   Graph maps can also be folded, subdivided, etc. If the graphs are
   marked then those operations will carry on the marking.

   Note that to associate an automorphism to a graph self map that is a
   homotopy equivalence we need to fix a base point to
   compute the fundamental group. Thus if we do not fix the base point we
   only get an outer automorphism of the free
   group. However, the program do not handle directly outer automorphism,
   rather ``f.automophism()`` returns an automorphism but with no
   guarantee on how the base is chosen, thus this automorphism is an
   arbitrary representative of the graph self map $f$.
   Moreover, if the base graph is not marked, then the automorphism is
   only defined up to conjugacy in $\Out(\FN)$.
   In this case ``f.automorphism()`` returns an arbitrary
   automorphism in the conjugacy class. We provide a
   ``phi.simple_outer_representative()`` which return an
   automorphism in the outer class of $\phi$ with the smallest possible length
   of images.

-  Train-tracks::

   The main feature and the main achievement of the program is to compute
   train-track representative for (outer) automorphisms of free groups.
   ``phi.train_track()`` computes a train-track representative for
   the (outer) automorphism phi. This train-track can be either an
   absolute train-track or a relative train-track. The celebrated theorem
   of M. Bestvina and M. Handel [bh-traintrack] assures that if $\Phi$ is
   fully irreducible (iwip) then there exists an absolute train-track
   representing $\Phi$.

   The ``{train\_track(relative=False)`` method will terminate with
   either an absolute train-track or with a topological representative
   with a reduction: an invariant strict subgraph with non-trivial
   fundamental group.

   One more feature of train-tracks (absolute or relative) is to lower
   the number of Nielsen paths. Setting the stable=True
   option will return a train-track with at most one indivisible Nielsen path
   (per exponential stratum if it is a relative
   train-track).

   -  Examples::

      Let's start with building absolute train-tracks::

          sage: phi = free_group_automorphisms.tribonacci()
          sage: phi.train_track()
          Train-track map:
          Marked graph: a: 0->0, b: 0->0, c: 0->0
          Marking: a->a, b->b, c->c
          Edge map: a->ab, b->ac, c->a
          Irreducible representative

      Indeed Tribonacci automorphism
      $\phi:a->ab,\ b->ac,\ c->a$
      is a positive automorphism (also called
      substitution), and thus it defines a map from the
      rose to itself which is a train-track map. Note that here the output
      is a ``TrainTrackMap``::

          sage: phi = FreeGroupAutomorphism("a->ab,b->ac,c->c")
          sage: phi.train_track(relative=False)
          Marked graph: a: 0->0, b: 0->0, c: 0->0
          Marking: a->a, b->b, c->c
          Edge map: a->ab, b->ac, c->c
          Strata: [set(['c']), set(['a', 'b'])]

      Here the automorphism is not irreducible (it fixes the free group
      element $c$).
      And the algorithm correctly detect that by returning a stratified
      graph map. Although the rose representative is reducible, it is a
      train track map (because the automorphism is positive). But this is
      not detected by the ``train\_track()`` method. We provide a
      ``is\_train\_track()`` method to test that.
      ::

          sage: phi = FreeGroupAutomorphism("a->ab,b->ac,c->c")
          sage: f = phi.rose_representative()
          sage: f.is_train_track()
          True

      You can promote this map to become a ``TrainTrackMap`` by using
      ``TrainTrackMap(f)``. This can be useful to compute Nielsen paths
      of such reducible train-track maps (but this may cause infinite loops
      in the program).

      Reducible automorphisms always have a relative train-track
      representative.
      ::

          sage: phi = FreeGroupAutomorphism("a->ab,b->ac,c->c")
          sage: phi.train_track()
          Graph self map:
          Marked graph: a: 2->0, b: 1->0, c: 0->0, d: 1->0, e: 2->0
          Marking: a->Ea, b->Db, c->c
          Edge map: a->b, b->ac, c->c, d->e, e->dAe
          Strata: [set(['c']), set(['a', 'b']), set(['e', 'd'])]

      (compare with the above example and note that the default option is
      relative=True). Ask for details of the computation by setting
      option verbose=1 (or 2, or more).

      The default option for this ``train_track`` method is to set
      stable=True, meaning that it looks for a stable
      train-track.

   -  Train-tracks and graph maps::

      In the previous section we computed train-track representatives for
      automorphisms of free group. The process goes by building a graph self
      map on the rose to represent the automorphism (this is called a
      topological representative\index{topological representative} and then
      perform operations on this graph self map.

      The graph on which the topological representative is built can be any
      kind of our graphs: ``GraphWithInverses``, ``MarkedGraph``,
      ``MetricGraph``, ``MarkedMetricGraph``. If the graph is not
      marked, then one give up the possibility to recover the original outer
      automorphism from the train-track. Indeed, all outer automorphisms in
      a conjugacy class in $\Out(\FN)$
      can be represented as the same homotopy equivalence on a graph.

      The train-track algorithm can be called directly on a graph self map
      ``f.train\_track()`` with the same options as for automorpism but
      $f$
      will not be promoted to become a ``TrainTrackMap`` even if it
      could. One can access intermediate operations like
      ``f.stabilize()``, ``f.reduce()``, etc.

   -  Nielsen paths::

      Nielsen paths are a main tool to refine the understanding of
      train-tracks and of automorphisms of free groups. A Nielsen
      path\index{Nielsen path} for a graph self map $f$
      is a path homotopic to its image relative to its endpoints. In our
      context, we only compute and use Nielsen paths in the case of
      train-track maps (or relative train-track maps).

      Nielsen paths of a graph self map ``f`` can be computed
      ``f.indivisible_nielsen_paths()``. The output is a list of pairs
      ``(u,v)`` of paths in the domain of ``f``. The paths
      ``u`` and ``v`` starts at the same vertex and the ends of
      the Nielsen path are inside the last edges of ``u`` and
      ``v``. We also provide the computation of periodic Nielsen
      paths, that-is-to-say Nielsen paths of
      iterates of ``f``. In this case a Nielsen path is coded by
      ``((u,v),period)``. To build longer Nielsen paths we need to
      concatenate the indivisible ones and for that we need to encode the
      endpoints of periodic Nielsen paths. This normal form for points
      inside edges is a little tricky and can be obtained using
      ``TrainTrackMap.periodic_point_normal_form()``.

-  More on free group automorphisms::

   The computation of other invariants for iwip
   automorphisms of free groups. Beware, that Python and Sage let you
   check the requirements: computing the index of a reducible
   automorphism may cause errors or infinite loops by the program.

   A graph self map as Whitehead graphs at each
   vertex and thanks to Brian Mann, they can be computed. The Whitehead
   graph of a graph self map $f:G\to G$
   at a vertex $v$
   as the set of germs of edges outgoing from $v$
   as vertices and as an edge for a germ from another if and only this
   turn is taken by the iterated image of an edge. Stable Whitehead
   graphs are also available: they only keep germs of edges which are periodic.

   Finally the ideal Whitehead graph is an
   invariant of iwip automorphisms. And we can compute them. From the
   ideal Whitehead graph one can compute the index list and the index of
   an iwip automorphism of a free group.

   Using the ``train_track()`` method our program can decide wether
   an automorphism is fully irreducible or not.
   If it is iwip, one can compute the ``index``, ``index-list`` or ``ideal
   Whitehead graphs``. Not that these computations are done using an
   absolute expanding train-track representative: they can be use for a
   broader class than just iwip automorphisms.

-  Convex cores, curve complex and more::

   The programm is also designed to handle trees in Outer space as well
   as simplicial trees in the boundary of Outer space.

   -  Metric simplicial trees and Outer space::

      Recall that M.Culler and K.Vogtmann [2] introduced the Outer
      space of a free group $\FN$
      which we denote by $\CVN$.
      Outer space is made of simplicial metric trees $T$
      with a free minimal action of the free group $\FN$
      by isometries. Alternatively a point in Outer space is a marked metric
      graph $T/\FN$.

      Our classes ``MetricGraph`` and ``MarkedMetricGraph`` allow
      us to handle points in Outer space. In a metric graph edges of length
      $0$
      can be used as an artefact to code simplicial trees with non-free
      action. For instance::

          sage: from sage.groups.free_groups.marked_graph import MarkedMetricGraph
          sage: A = AlphabetWithInverses(3)
          sage: G = MarkedMetricGraph.splitting(2,A)
          sage: print(G)
          Marked graph: a: 0->0, b: 0->0, c: 1->1, d: 0->1
          Marking: a->a, b->b, c->dcD
          Length: a:0, b:0, c:0, d:1

      This graphs codes the splitting\index{splitting} of the free group
      $F_3=F_2*\mathbb Z$.
      HNN-splittings are also available:
      ``MarkedMetricGraph.HNN_splitting()``. Thus the metric graphs
      (with edges of length $0$)
      are a convenient tool to work with the splitting complex.

      Let us emphasize that the splitting complex of a free group is
      becoming a popular tool after being proved hyperbolic by M.Handel and
      L.Mosher [4].

      In the geometric situation, these non-free trees can be used to encode
      arcs in the arc complex
      for a surface $S_{g,1}$
      of genus $g$
      with one puncture, ideal arcs (a curve from the puncture to itself)
      are in one-to-one correspondance with splittings of the free group
      $\pi_1(S_{g,1})$. They can also be used to study curve diagram in the
      context of braid groups.

   -  Convex cores::

      Are also implemented the computation of V.Guirardel [6]
      convex core of two simplicial trees in outer space and its
      boundary. The convex core is a square complex inside the cartesian
      product $T_0\times T_1$
      of two trees with action of the free group. Here it is encoded by its
      quotient $C(T_0\times T_1)/\FN$
      which is a finite square complex.  We have the convention that the
      convex core is connected and thus we give up unicity: instead we
      include twice light squares inside the core.


      The first way to create a convex core is by using a free group
      automorphism ``phi``. Then ``ConvexCore(phi)`` returns the
      convex core of the Cayley graph $T_0$
      of the free group with the tree $T_1$
      which is the same as $T_0$
      but with the action twisted by ``phi``.
      ::
          sage: from sage.groups.free_groups.convex_core import ConvexCore
          sage: phi = FreeGroupAutomorphism("a->ab,b->ac,c->a")**2
          sage: C = ConvexCore(phi)
          sage: C.squares()
          [[3, 0, 2, 1, 'c', 'a']]
          sage: C.one_squeleton(side=1)
          Looped multi-digraph on 4 vertices

      The second way involves creating the two trees. This requires the
      creation of two marked graphs, which can be a little teddious, but
      some methods shorten the typesetting.
      ::

          sage: from sage.groups.free_groups.marked_graph import MarkedGraph
          sage: A = AlphabetWithInverses(2)
          sage: G1 = MarkedGraph(GraphWithInverses.rose_graph(A))
          sage: G2 = MarkedGraph(GraphWithInverses.rose_graph(A))
          sage: C = ConvexCore(G1, G2)
          sage: C.volume()
          0
      Remark that if the automorphism is a mapping class and the trees are
      transverse to ideal curves then the convex core (as a CW-complex) is
      homeomorphic to the surface.

-  References::

    [1] M. Bestvina, and M. Handel, Train tracks and
        automorphisms of free groups. Ann.  of Math.  (2) 135
        (1992), no.  1, 1--51

    [2]  M.Culler, K.Vogtmann, Moduli of graphs and
          automorphisms of free groups. Invent.  Math.  84 (1986), no.
          1, 91--119

    [3] Vincent Guirardel, Cœur et nombre d’intersection pour les actions
        de groupes sur les arbres, Ann. Sci. Ecole Norm. Sup. (4) 38 (2005),
        no. 6, 847–888. MR MR2216833 (2007e:20055)

    [4] Michael Handel and Lee Mosher, The free splitting complex of a free group,
        I: hyperbolicity, Geom. Topol. 17 (2013), no. 3, 1581–1672. MR 3073931

    [5] Jean-Pierre Serre, Arbres, amalgames, SL 2 , Société
        Mathématique de
        France, Paris, 1977, Avec un sommaire anglais, Réedigé avec la collabora-
        tion de Hyman Bass, Astérisque, No. 46. MR 0476875

    [6] __ Trees, Springer Monographs in Mathematics, Springer-Verlag,
        Berlin, 2003, Translated from the French original by John Stillwell,
        Corrected 2nd printing of the 1980 English translation. MR 1954121

    [7] John R. Stallings, Topology of finite graphs, Invent. Math. 71 (1983), no. 3,
        551–565. MR MR695906 (85m:05037a)