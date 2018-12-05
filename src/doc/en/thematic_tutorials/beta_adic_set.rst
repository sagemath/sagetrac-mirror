.. -*- coding: utf-8 -*-
.. _beta_adic_set:

===================
BetaAdicSet in Sage
===================

.. MODULEAUTHOR:: Paul Mercat, and Dominique Benielli
                  Labex Archimede - I2M -
                  AMU Aix-Marseille Universite


The present tool BetaAdicSet is a convenient tool to describe sets like Rauzy fractals and quasicrystals associated to substitutions.


Introduction and Definitions
----------------------------

Subtitution
~~~~~~~~~~~

Given a substitution, let's say for example the Fibonnacci substitution :

.. MATH::
    \left\{
    \begin{array}{rcl}
    a & \mapsto & ab \\
    b & \mapsto & a
    \end{array}
    \right.

A fixed point always exists, let see get the following fixed point :

:math:`abaababaabaababaababaabaababaabaababaaba...`

This infinite sequence has a lot of interesting properties in general.
If we replace letters by intervals of convenient lengths, we get a self-similar tiling of :math:`\mathbf R_+`.
Convenient lengths are given by the non-negative coefficients of a Perron eigenvector of the incidence matrix :math:`M` of the substitution.
Here this matrix is

.. MATH::
    M = \begin{pmatrix}
    1 & 1 \\
    1 & 0
    \end{pmatrix}
  
and :math:`\begin{pmatrix} 1 \\ \beta-1 \end{pmatrix}` is an eigenvector for the Perron eigenvalue :math:`\beta`, which is the golden number.

Self-similar tiling of :math:`\mathbf R_+`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The set of points that appears at the boundaries of intervals in this self-similar tiling are elements of :math:`\mathbf Q(\beta)`, and this set of points have very strong properties, because it is a non-periodic self-similar Meyer set.
Here, we get the points
  
:math:`0, 1,\ \beta,\ \beta + 1,\ \beta + 2,\ 2\beta + 1,\ 2\beta + 2,\ 3\beta + 1,\ 3\beta + 2,\ 3\beta + 3, \ 4\beta + 2,\ 4\beta + 3, ...`

To get this set of points, we start from :math:`0` and read the fixed point : 
we add :math:`1` each time we read a letter :math:`a` and we add :math:`\beta-1` each time we read a letter :math:`b`.
  
Automaton
~~~~~~~~~

This set of points of :math:`\mathbf Q(\beta)` is described by the automaton followed figure describing the quasicrystal 
associated to the Fibonnacci substitution:

.. PLOT::
   :width: 30%

    a = DetAutomaton([('a', 'b', '1'), ('a', 'a', '0'),('b', 'a', '0')], i='a')
    sphinx_plot(a)

The initial state of this automaton is the one labelled by letter 'a', and the two states 'a' and 'b' are final states.
The language recognized by this automaton, defined as the set of sequences of labels of paths from the initial state to a final state,
is here exactly the set of words over the alphabet :math:`\{0,1\}` that does not contain the subword :math:`11`.

The self-similar tiling of :math:`\mathbf R_+` is obtained by taking the embedding of :math:`\mathbf Q(\beta)` in :math:`\mathbf R` corresponding to the Perron eigenvalue of the incidence matrix :math:`M`.
If we look at the others embeddings, corresponding to eigenvalues less than :math:`1` we get a bounded set whose adherence is called Rauzy fractal of the substitution.
Here there is a unic such embedding, which is a real one, corresponding to the root of :math:`x^2-x-1` between :math:`-1` and $:math:`1`.
The Rauzy fractal is here the interval :math:`[-1, \varphi]` of :math:`\mathbf R` where :math:`\varphi` is the golden number (i.e. greatest root of :math:`x^2-x-1`).


Main Definitions
----------------

Broken line
~~~~~~~~~~~

Consider :math:`s` substitution, or in others words, a word morphism over a finite alphabet :math:`A = \{a_1, ..., a_n\}`.
Up to replace :math:`s` by a power, we can assume that :math:`s` has a fixed point :math:`\omega`.
    
We defined the \defi{broken line} associated to :math:`\omega` as the subset of :math:`\mathbf Z^n` defined by

.. MATH::
    \{ {\begin{pmatrix}
    \text{number of occurences of } a_1 \text{ in } \omega_k \\
    \text{number of occurences of } a_2 \text{ in } \omega_k \\
    \vdots \\
    \text{number of occurences of } a_n \text{ in } \omega_k
    \end{pmatrix}
    \in \mathbf Z^n
    } 
    {k \in \mathbf N} \}
  
where :math:`\omega_k` is the prefix of length :math:`k` of the infinite word :math:`\omega`.

This broken line is very interesting since it is a geometrical object which completely encode the substitution and is stable by multiplication by the incidence matrix.

Rauzy fractal
~~~~~~~~~~~~~

The Rauzy fractal is the closure of the projection of the broken line to the contracting space along the expanding line.

Expanding line and contracting space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The \og expanding line \fg\ has dimension :math:`1` for Pisot numbers, but it can have greater dimension for other Perron numbers.

Let :math:`M_s` be the incidence matrix of the substitution :math:`s`.
By definition the coefficient :math:`(i,j)` of this matrix is the number of occurrences of the letter :math:`a_j` in the word :math:`s(a_i)`.
By Perron-Frobenius theorem, there exists an eigenvector :math:`v \in (\mathbf R_+)^n`, unic if the matrix is irreducible, for an eigenvalue :math:`\lambda` which is the spectral radius of :math:`M_s`,
and moreover we can assume that :math:`v \in (\mathbf Q(\lambda))^{n}`.
    
We can define a sort of broken line in $:math:`\mathbf Q(\lambda)`, by the following. 

.. MATH::

    \{ Q_\omega = { \\sum_{k=1}^N v_{a_k} }{ N \in \mathbf N,\ a_1 a_2 ... a_N \text{ prefix of } \omega \text{ of length } N }.
    

This is a projection of the broken line on :math:`\mathbf Q(\lambda)`.
This set is invariant by multiplication by the Perron eigenvalue :math:`\lambda` and gives a self-similar tiling of :math:`\mathbf R_+`.
The definition of :math:`Q_\omega` depends of the choice of an eigenvector.
We prefer to choose an eigenvector whose coefficients belongs to the integer ring :math:`\mathcal O_\lambda`, in order to have :math:`Q_\omega \subset \mathcal O_\lambda`.
    
For :math:`\mathbf Q(\lambda)`, there are natural contracting and expanding spaces for the multiplication by :math:`\lambda`.
Indeed, consider the bigest sets :math:`P_+` and :math:`P_-` of places (i.e. equivalence class of absolute values) 
of :math:`\mathbf Q(\lambda)` such that

.. MATH::

    	\forall v \in P_+,\ |{\lambda}|_v > 1 \quad \text{ and } \quad \forall v \in P_-,\ |{\lambda}|_v < 1.
  
If :math:`\lambda` is an algebraic unit, the set :math:`P_+` corresponds to roots of the minimal polynomial of :math:`\lambda` greater than :math:`1` in absolute value, counting two conjugate complexes only once,
and it is the same for :math:`P_-` with the roots of modulus less than :math:`1`.
    
For each place :math:`v`, we define a space :math:`E_v` as the completion of :math:`\mathbf Q(\lambda)` for the absolute value :math:`v`.
If :math:`v` is a real place (i.e. corresponding to a real root or the minimal polynomial of :math:`\lambda`), then :math:`E_v = \mathbf R`.
If :math:`v` is a complex place (i.e. corresponding to two conjugated complex roots or the minimal polynomial of :math:`\lambda`), then :math:`E_v = \mathbf C`.
Otherwise, :math:`E_v` is a :math:`p`-adic space, which is a finite extension of the :math:`p`-adic field :math:`\mathbf Q_p` (which is the completion of :math:`\mathbf Q` for the :math:`p`-adic absolute value).
    
    
We can define the expanding space

.. MATH::

    	E_\lambda^+ := \prod_{v \in P_+} E_v,

and the contracting one

.. MATH::

    	E_\lambda^- := \prod_{v \in P_-} E_v.


Let's take :math:`\sigma_+` and :math:`\sigma_-` some embeddings of :math:`\mathbf Q(\lambda)`
into the spaces :math:`E_+` and :math:`E_-` respectively.
We will also denote by :math:`\sigma_\beta` the maximal real embedding when :math:`\beta` is a Perron number.
      
So Rauzy fractal of the substitution :math:`s` can be define as the adherence of :math:`\sigma_-(Q_\omega)` in :math:`E_{\lambda}^-`.

Set :math:`P`
^^^^^^^^^^^^^
Let :math:`\beta` be a Pisot number (not necessarly unit), and let :math:`P \subseteq E_\beta^-`.
The set :math:`P` is arbitrarily approximated by Rauzy fractals, for the Hausdorff distance, associated to :math:`\beta^n`, 
if and only if :math:`P` is bounded and :math:`0 \in \overline{P}`.


g-:math:`\beta-sets`
~~~~~~~~~~~~~~~~~~~~

A g-:math:`\beta`-set, for an algebraic number :math:`\beta`, is a subset of :math:`\mathbf Q(\beta)` of the form

.. MATH::

        { \mathbf Q_{\beta,L} := \{ \sum_{i=0}^n a_i \beta^i} { n \in \mathbf N,\ a_0 a_1 ... a_n \in L } \}.

where :math:`L` is a regular language over a finite alphabet :math:`\Sigma \subset \mathbf Q(\beta)`.

Some Properties
^^^^^^^^^^^^^^^

For a fixed algebraic number :math:`\beta` with no conjugate of modulus one,
the set of g-:math:`\beta`-sets is stable by

* intersection
* union
* complementary (in another g-:math:`\beta`-set)
* Minkowski sum (i.e. the sum of two g-:math:`\beta`-sets is a g-:math:`\beta`-set)
* multiplication by an element of :math:`\mathbf Q(\beta)`
* translation by an element of :math:`\mathbf Q(\beta)`
* adherence, interior, boundary, for the topology of :math:`\mathcal O_\beta` induced by :math:`E_-`. 



The fact that g-:math:`\beta`-sets come naturally to describe quasicrystals arising from substitutions
and has a lot of nice properties show that it is an interesting fundamental object.
    

Remarks: on any Shape
^^^^^^^^^^^^^^^^^^^^^

We see from theses properties that we can construct g-:math:`\beta`-sets with any shape in the contracting space :math:`E^-`.
This allows us to construct Rauzy fractals of any shape.



Construction of a domain exchange
---------------------------------
The first step, to construct a substitution from a quasicrystal, is to construct a domain exchange which describe the shift on the quasicrystal.

Let :math:`\beta` be a Pisot number (eventually non unit), and let :math:`Q \subseteq \mathbf Q(\beta)` such that :math:`\sigma_+(Q)` is a quasicrystal of :math:`\mathbf R` or :math:`\mathbf R^+`.
Then there exists a domain exchange with a finite number of pieces such that the union of the pieces is :math:`Q`.
Moreover, this domain exchange is conjugated to the shift on :math:`\sigma_+(Q)`. %defined by the window $:math:`\Omega`. 

.. figure::echange_rond2.pdf, echange_rond1.pdf
  :scale: 40 %


  .. image:: echange_rond2.pdf
  .. image:: echange_rond1.pdf
  Construction of a domain exchange in the unit disk, for the integer ring :math:`\mathcal O_\beta`,
  where :math:`\beta` is the Tribonnacci number. 
  \textcolor{red}{:math:`-2\beta^2+2\beta`}, \quad \textcolor{orange}{:math:`\beta^2-\beta-1`}, \quad \textcolor{lime}{:math:`\beta-1`}, \quad \textcolor{green}{:math:`1`}, \quad \textcolor{cyan}{:math:`-\beta^2+2\beta+1`, \quad \textcolor{bleu}{:math:`\beta^2-\beta`, \quad \textcolor{magenta}{:math:`\beta`

The domain exchange described in the figure for the open unit disk gives exactly the list of Pisot numbers (including non-unit ones) of degree :math:`3` in :math:`\mathbf Q(\beta)`,
where :math:`\beta` is the Tribonnacci number (i.e. greatest root of $x^3-x^2-x-1$).
Indeed if :math:`x` is a Pisot number of degree three in :math:`\mathbf Q(\beta)`, the next Pisot number is obtained by looking in which piece is the conjugate :math:`\overline{x}`,
and adding the corresponding translation to :math:`x`.

Construction of a substitution
------------------------------

If we know that a quasicrystal :math:`\sigma_+(Q)` of `\mathbf R` or :math:`\mathbf R_+` comes from the fixed point of a substitution for a Pisot number :math:`\lambda`,
it is not difficult to guess what is the substitution.
Indeed, it is enough to take intervals between two consecutive points, multiply it by :math:`\lambda`,
and see how the result is covered by others intervals.

.. image:: media/subtitution.png
  :scale: 20 %

Construction of a domain exchange in the disk of radius :math:`1` and center :math:`0`,
for the Tribonnacci number :math:`\beta`.

But we have to take care of the fact that one interval can have several substitutions rules,
corresponding to the fact that several letters of a substitution can give intervals of same lengths.

If we look at what happens in the contracting space :math:`E^-`, we have to do a sort of induction on :math:`\lambda Q`
for the domain exchange on :math:`Q`, and we have to iterate it up to stabilization.
But it's not really an induction : we have to distinguish between different possible 
trajectories for points in :math:`\lambda Q` before they come back to :math:`\lambda Q`,
otherwise the induction only give the same domain exchange on :math:`\lambda Q` than in :math:`Q`.


Examples of Usage of BetaAdicSet
--------------------------------




A Sierpinsky gasket
~~~~~~~~~~~~~~~~~~~

Take the Tribonnacci Pisot number β, root of x 3 − x 2 − x − 1,
and take L the regular language defined by the followed automaton.


This automaton describing the regular language describing a g- :math:`\beta`-set which is a Sierpiński
gasket union a set of non-empty interior for :math:`\beta` the Tribonnacci number.

.. PLOT::
   :width: 80%

   # automaton that describe a Sierpinsky gasket
   a = DetAutomaton([(0,2,0),(0,6,1),(2,3,1),(2,12,0),(6,7,1),(6,9,0),(3,4,1),(3,5,0),(12,13,1),(12,14,0),(7,8,0),(7,15,1),(9,10,0),(9,11,1),(4,0,0),(5,0,0),(5,0,1),(13,0,0),(13,0,1),(14,0,0),(8,0,0),(8,0,1),(15,0,1),(10,0,1),(11,0,1),(11,0,0)], i=0)

   # automaton recognizing a set of non-empty interior
   a2 = DetAutomaton([(0,1,0),(1,2,0),(2,2,0),(2,2,1)],i=0, final_states=[2])
   # multiply by b^2
   a3 = a.unshift1(0, final=True).unshift1(1)
   a = a2.union(a3)
   sphinx_plot(a)

Obtained by the code:

.. code-block:: Python

   # automaton that describe a Sierpinsky gasket
   a = DetAutomaton([(0,2,0),(0,6,1),(2,3,1),(2,12,0),(6,7,1),(6,9,0),(3,4,1),(3,5,0),(12,13,1),(12,14,0),(7,8,0),(7,15,1),(9,10,0),(9,11,1),(4,0,0),(5,0,0),(5,0,1),(13,0,0),(13,0,1),(14,0,0),(8,0,0),(8,0,1),(15,0,1),(10,0,1),(11,0,1),(11,0,0)], i=0)

   # automaton recognizing a set of non-empty interior
   a2 = DetAutomaton([(0,1,0),(1,2,0),(2,2,0),(2,2,1)],i=0, final_states=[2])
   # multiply by b^2
   a3 = a.unshift1(0, final=True).unshift1(1)
   a = a2.union(a3)
   a.plot()



The Domain exchange with :math:`6` pieces, describing the shift on :math:`\sigma_+(Q_L)` for the regular language :math:`L` can be computed.

.. code-block:: Python

   m = BetaAdicSet(x^3-x^2-x-1, a) #choose to work with the alphabet {0,1} and with the Tribonnacci polynomial
   pp = m.b.parent().places()[0] #expanding place
   print pp
   Ring morphism:
     From: Number Field in b with defining polynomial x^3 - x^2 - x - 1
     To:   Real Field with 106 bits of precision
     Defn: b |--> 1.839286755214161132551852564671
   m.plot(nprec=6)

.. image:: media/beta_adic_image1.png
  :scale: 80 %

Now the code to plot the list after the exchange

.. code-block:: Python

   # compute a domain exchange
   l = m.domain_exchange()
   print("Exchange with %s pieces."%len(l))
   Exchange with 6 pieces.
   # plot it
   m.plot_list([a for t,a in l], nprec=6)



.. image:: media/domain1.png
  :scale: 80 %

And plot the domain after exchange

.. code-block:: Python

   # plot it after exchange
   m.plot_list([a.proj(m, t) for t,a in l], nprec=6)

.. image:: media/domain2.png
  :scale: 80 %


Compute the subtitution

.. code-block:: Python

   # compute a substitution whose Rauzy fractal is this BetaAdicSet
   %time d = m.substitution()
   d
   CPU times: user 24 s, sys: 156 ms, total: 24.1 s
   Wall time: 24.1 s
   
   {1: [60, 6],
    2: [19],
    3: [19, 54],
    4: [50, 42],
    5: [57, 9, 58, 3],
    6: [60, 6, 40, 48],
    7: [60, 6, 53],
    8: [21, 35, 48, 60, 1],
    9: [19, 55, 5],
    10: [21, 66, 49, 60, 1],
    11: [64, 6, 15, 5],
    12: [60, 6, 63, 49, 60, 1],
    13: [53, 64, 7, 25, 4],
    14: [54, 20, 33, 4],
    15: [60, 18, 38, 3, 37, 46, 58, 2],
    16: [36, 17, 45, 41, 46, 58, 2],
    17: [64, 6, 53, 5],
    18: [60, 6, 53, 64, 1],
    19: [57, 9, 58, 3, 37, 46],
    20: [57, 9, 58, 3, 52],
    21: [34, 11, 58, 3, 37, 46],
    22: [34, 11, 58, 3, 52],
    23: [52, 41, 3, 52, 4],
    24: [64, 18, 43, 41, 46, 58, 2],
    25: [64, 18, 43, 50, 4],
    26: [57, 9, 58, 3, 37, 46, 58, 2],
    27: [57, 9, 58, 3, 52, 41, 2],
    28: [40, 48, 60, 7, 65, 47, 58, 2],
    29: [35, 48, 22, 61, 47, 58, 2],
    30: [34, 11, 58, 3, 37, 46, 58, 2],
    31: [34, 11, 58, 3, 52, 41, 2],
    32: [41, 46, 45, 41, 46, 58, 2],
    33: [41, 46, 45, 50, 4],
    34: [15],
    35: [16],
    36: [24],
    37: [26],
    38: [28],
    39: [29],
    40: [30],
    41: [32],
    42: [50, 42, 50],
    43: [13, 42, 50],
    44: [14, 42, 50],
    45: [23, 42, 50],
    46: [19, 54, 5],
    47: [21, 54, 5],
    48: [64, 6, 40, 48, 60, 1],
    49: [60, 6, 40, 48, 60, 1],
    50: [50, 42, 50, 4],
    51: [23, 42, 50, 4],
    52: [27, 44, 50, 4],
    53: [31, 39, 3, 37, 46, 58, 2],
    54: [51, 41, 3, 37, 46, 58, 2],
    55: [58, 46, 58, 3, 37, 46, 58, 2],
    56: [37, 46, 58, 3, 37, 46, 58, 2],
    57: [55],
    58: [56],
    59: [59, 12],
    60: [62, 12],
    61: [61, 8],
    62: [63, 49],
    63: [60, 49],
    64: [68, 10],
    65: [64, 49],
    66: [69, 8],
    67: [65, 8],
    68: [66, 49],
    69: [67, 49]}

 The g-:math:`\beta`-set :math:`Q_{]-1,1[}` can be computed, for any quadratic Pisot number :math:`\beta`, and then compute a substitution describing the quasicrystal.

And directly with the WordMorphism of the subtitution and it's rauzy_fractal_plot.

.. code-block:: Python

    #plot the Rauzy fractal from the substitution
    s = WordMorphism(d)
    s.rauzy_fractal_plot()

.. image:: media/domain3.png
  :scale: 100 %


The Dragon Fractal
~~~~~~~~~~~~~~~~~~

.. code-block:: Python

    ################################################
    # The dragon fractal
    ################################################
    m = BetaAdicSet(1/(1+I), [0,1])
    m
    b-adic set with b root of x^2 - x + 1/2, and an automaton of 1 states and 2 letters


.. code-block:: Python

    a = m.relations_automaton(ext=True)
    a.plot()

.. PLOT::
   :width: 60%

    m = BetaAdicSet(1/(1+I), [0,1])
    a = m.relations_automaton(ext=True)
    sphinx_plot(a)

.. code-block:: Python

    mi = m.intersection_words([0], [1])
    m.plot_list([mi])

.. image:: media/dragon1.png
  :scale: 70 %


.. code-block:: Python

    mi.plot(nprec=6)

.. image:: media/dragon2.png
  :scale: 70 %

Compute the Hausdorff dimension.

.. code-block:: Python

    # compute the Hausdorff dimension
    mi.critical_exponent()
    log(y)/log(1.414213562373095?) where y is the max root of x^3 - x^2 - 2, and 1.414213562373095? is root of x^2 - 2.
    1.523627086202492


Any Shape
~~~~~~~~~

Disk
----

that permit to draw a Rauzy fractal of any shape with the mouse, like in a drawing software,
and to compute the corresponding substitution.
The following example has been obtain by drawing randomly using this tool.

Definition of the beta-Adic-Set:

.. code-block:: Python

    ######################################
    # BetaAdicSet approximating a disk
    ######################################
    #. BetaAdicSet approximating a square
    m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
    m
    b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 3 states and 2 letters


The relation automaton associated

.. PLOT::
   :width: 60%

    ######################################
    # BetaAdicSet approximating a disk
    ######################################
    #. BetaAdicSet approximating a square
    m = WordMorphism('a->ab,b->ac,c->a').DumontThomas().mirror()
    a = m.relations_automaton()
    sphinx_plot(a)


.. code-block:: Python

    pm = m.b.parent().places()[1]
    pm
    Ring morphism:
      From: Number Field in b with defining polynomial x^3 - x^2 - x - 1
      To:   Complex Field with 53 bits of precision
      Defn: b |--> -0.419643377607080 + 0.606290729207199*I

The disk definition:

.. code-block:: Python

    md = m.approx(14, lambda x: (pm(x).real())^2 + (pm(x).imag())^2 < .4)
    print(md)
    b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 265 states and 2 letters

.. code-block:: Python

    m.plot_list([md])


.. image:: media/shap1.png
  :scale: 70 %

.. code-block:: Python

    md1 = md.proj(m)
    md1
    b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 253 states and 2 letters

    # domain exchange for this set
    l = md1.domain_exchange()
    print(l)
    [(1, b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 99 states and 2 letters), (b^2 - b, b-adic set with b root of x^3 - x^2
    - x - 1, and an automaton of 70 states and 2 letters), (b, b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 134 states and 2
    letters), (b + 1, b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 99 states and 2 letters), (b^2, b-adic set with b root of
    x^3 - x^2 - x - 1, and an automaton of 164 states and 2 letters), (b^2 + 1, b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of
    61 states and 2 letters), (b^2 + b, b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 37 states and 2 letters), 
    (b^2 + b + 1, b-adic set with b root of x^3 - x^2 - x - 1, and an automaton of 56 states and 2 letters)]

    md1.plot_list([a for t,a in l], nprec=6)

.. image:: media/shap2.png
  :scale: 70 %

And the domain exchange after exchange

.. code-block:: Python

    # plot the domain exchange after exchange
    md1.plot_list([a.proj(md, t) for t,a in l], nprec=6)

.. image:: media/shap3.png
  :scale: 70 %

 compute a substitution whose Rauzy fractal is this approximation of disk

.. code-block:: Python

    %time d = md.substitution()
    d
    CPU times: user 48.6 s, sys: 783 ms, total: 49.4 s
    Wall time: 49.2 s

    {1: [248, 318, 288, 324],
     2: [168, 272, 92],
     3: [264, 274],
     4: [407, 2],
     5: [117, 99],
     6: [352, 218],
     7: [226, 235, 372, 323],
     8: [415, 134, 309, 232, 380, 291, 93],
     9: [411, 6, 232, 288, 273, 208],
     10: [374, 310, 123, 168, 292, 92],
     11: [232, 169, 273, 208],
     12: [232, 288, 281],
     13: [411, 329, 232, 288, 273, 208],
     14: [415, 339, 232, 380, 291, 93],
     15: [237],
     16: [4],
     17: [5],
     18: [152],
     19: [8],
     20: [9],
     21: [10],
     22: [13],
     23: [14],
     24: [317, 341, 368, 322],
     25: [15, 172, 385],
     26: [391, 400, 309],
     27: [168, 292, 92],
     28: [235, 217, 323],
     29: [3, 367, 165],
     30: [107, 334],
     31: [302, 227],
     32: [26, 11],
     33: [56, 258],
     34: [226, 235, 29, 323],
     35: [226, 128, 256, 94],
     36: [312, 226, 128, 256, 94],
     37: [234, 17, 193],
     38: [318, 288, 325],
     39: [411, 399, 309, 232, 288, 273, 208],
     40: [391, 400, 309, 232, 169, 273, 208],
     41: [414, 328, 232, 288, 33],
     42: [16, 248, 318, 288, 325],
     43: [19],
     44: [20],
     45: [21],
     46: [22],
     47: [23],
     48: [299],
     49: [30],
     50: [32],
     51: [24, 361],
     52: [36],
     53: [41],
     54: [42],
     55: [380, 291, 93],
     56: [224, 5],
     57: [110, 152],
     58: [317, 250, 368],
     59: [248, 318, 288],
     60: [58, 367, 361],
     61: [356, 228],
     62: [357, 55],
     63: [404, 25],
     64: [407, 27],
     65: [254, 228],
     66: [255, 55],
     67: [140, 12],
     68: [139, 38],
     69: [270, 11],
     70: [259, 257],
     71: [271, 78],
     72: [226, 235, 300, 323],
     73: [396, 163, 126, 299, 201],
     74: [332, 226, 126, 296, 50],
     75: [31, 162, 122, 269, 209],
     76: [31, 162, 128, 268, 96],
     77: [318, 288, 198],
     78: [232, 288, 37],
     79: [316, 236, 261],
     80: [205, 120, 105],
     81: [89, 286, 318, 288, 70],
     82: [100],
     83: [52],
     84: [53],
     85: [54],
     86: [61],
     87: [62],
     88: [63],
     89: [64],
     90: [65],
     91: [66],
     92: [67],
     93: [68],
     94: [69],
     95: [70],
     96: [71],
     97: [178, 323],
     98: [177, 361],
     99: [176, 95],
     100: [73],
     101: [74],
     102: [75],
     103: [76],
     104: [81],
     105: [128, 216, 96],
     106: [124, 59, 207],
     107: [88, 330, 127],
     108: [89, 248, 318, 288, 324],
     109: [89, 248, 318, 288, 70],
     110: [102],
     111: [103],
     112: [115],
     113: [104],
     114: [108],
     115: [109],
     116: [224, 171, 161],
     117: [224, 171, 113],
     118: [110, 171, 112],
     119: [111, 171, 85],
     120: [162],
     121: [289],
     122: [43],
     123: [44],
     124: [45],
     125: [46],
     126: [47],
     127: [50],
     128: [87],
     129: [89],
     130: [91],
     131: [94],
     132: [110],
     133: [49, 126],
     134: [226, 126],
     135: [137, 231],
     136: [138, 127],
     137: [395, 49, 126, 299],
     138: [332, 226, 126, 296],
     139: [392, 157],
     140: [414, 328],
     141: [124, 1],
     142: [129, 1],
     143: [360, 314],
     144: [86, 72],
     145: [88, 315],
     146: [90, 34],
     147: [183, 180],
     148: [107, 180],
     149: [230, 179],
     150: [117, 54],
     151: [117, 114],
     152: [117, 115],
     153: [117, 106],
     154: [390, 51],
     155: [389, 98],
     156: [345, 99],
     157: [206, 126, 98],
     158: [391, 400, 309, 232, 185, 208],
     159: [82, 398, 162, 122, 269, 209],
     160: [51],
     161: [99],
     162: [147],
     163: [148],
     164: [149],
     165: [153],
     166: [155],
     167: [156],
     168: [79],
     169: [80],
     170: [116],
     171: [117],
     172: [118],
     173: [275, 323],
     174: [39],
     175: [40],
     176: [89, 248, 318, 288],
     177: [317, 250, 368, 322],
     178: [265, 249, 367, 165],
     179: [289, 359],
     180: [379, 150],
     181: [289, 151],
     182: [289, 152],
     183: [360, 340, 127],
     184: [79, 132, 18],
     185: [80, 234, 17],
     186: [135, 228],
     187: [136, 119],
     188: [211, 117],
     189: [52, 118],
     190: [100, 227],
     191: [101, 119],
     192: [276, 38],
     193: [276, 77],
     194: [278, 78],
     195: [364, 310, 123, 184, 92],
     196: [363, 162, 122, 269, 209],
     197: [364, 164, 123, 298, 385],
     198: [234, 167, 258],
     199: [203, 236, 28],
     200: [90, 236, 28],
     201: [174],
     202: [175],
     203: [186],
     204: [187],
     205: [190],
     206: [191],
     207: [324],
     208: [193],
     209: [194],
     210: [195],
     211: [196],
     212: [197],
     213: [82, 398, 162],
     214: [122, 269, 209],
     215: [232, 185, 208],
     216: [247, 367, 165],
     217: [249, 367, 165],
     218: [229, 367, 361],
     219: [230, 121, 362],
     220: [183, 334],
     221: [226, 235, 301, 94],
     222: [226, 130, 297, 94],
     223: [210],
     224: [211],
     225: [212],
     226: [220],
     227: [83, 172, 385],
     228: [83, 172, 84],
     229: [317, 341, 368],
     230: [317, 342, 131],
     231: [201],
     232: [202],
     233: [223],
     234: [224],
     235: [265],
     236: [226],
     237: [316, 35],
     238: [317, 221],
     239: [316, 222],
     240: [264, 7],
     241: [264, 34],
     242: [396, 163, 126],
     243: [238, 179],
     244: [239, 181],
     245: [213, 214],
     246: [270, 215],
     247: [203, 236, 97],
     248: [204, 130, 98],
     249: [264, 236, 173],
     250: [226, 235, 371],
     251: [267],
     252: [243],
     253: [244],
     254: [405, 285, 231],
     255: [415, 339, 232],
     256: [200, 367, 165],
     257: [233, 166, 38],
     258: [233, 166, 77],
     259: [224, 397, 161],
     260: [128, 268, 96],
     261: [128, 256, 94],
     262: [254, 311],
     263: [255, 333],
     264: [262],
     265: [263],
     266: [237, 57],
     267: [237, 182],
     268: [199, 367, 165],
     269: [199, 319, 165],
     270: [225, 400, 309],
     271: [225, 400, 160],
     272: [132, 18],
     273: [234, 17],
     274: [226, 275, 323],
     275: [265, 327, 165],
     276: [223, 155],
     277: [224, 156],
     278: [225, 154],
     279: [56, 68],
     280: [56, 192],
     281: [56, 193],
     282: [347, 105],
     283: [343, 214],
     284: [355, 261],
     285: [49, 126, 299],
     286: [204, 130, 60],
     287: [415, 339, 232, 380, 279],
     288: [282],
     289: [283],
     290: [284],
     291: [56],
     292: [57],
     293: [158],
     294: [159],
     295: [287],
     296: [144, 367, 165],
     297: [146, 367, 165],
     298: [237, 132, 18],
     299: [238, 121, 362],
     300: [240, 367, 165],
     301: [241, 367, 165],
     302: [242, 351],
     303: [143, 119],
     304: [145, 334],
     305: [237, 118],
     306: [293],
     307: [294],
     308: [295],
     309: [218],
     310: [219],
     311: [290, 172, 84],
     312: [321, 227],
     313: [403, 228],
     314: [226, 126, 296, 50],
     315: [226, 130, 296, 50],
     316: [312],
     317: [313],
     318: [306],
     319: [307],
     320: [308],
     321: [242, 48, 231],
     322: [245],
     323: [246],
     324: [277, 257],
     325: [277, 258],
     326: [264, 7, 322],
     327: [264, 274, 322],
     328: [365, 126, 218],
     329: [365, 320, 218],
     330: [226, 130, 296],
     331: [337, 349],
     332: [404, 354],
     333: [380, 170, 93],
     334: [379, 171, 85],
     335: [336, 366, 126],
     336: [337, 251, 387],
     337: [335, 252, 125],
     338: [338, 253, 388],
     339: [226, 126, 218],
     340: [226, 126, 296],
     341: [226, 235, 300],
     342: [226, 235, 301],
     343: [31, 162],
     344: [184, 92],
     345: [188, 113],
     346: [189, 385],
     347: [190, 162],
     348: [266, 385],
     349: [267, 385],
     350: [298, 385],
     351: [299, 201],
     352: [303, 308],
     353: [304, 47],
     354: [305, 385],
     355: [312, 226],
     356: [405, 133, 48, 231],
     357: [415, 134, 309, 232],
     358: [117, 141],
     359: [117, 142],
     360: [332],
     361: [358],
     362: [359],
     363: [302, 346],
     364: [396, 353],
     365: [303],
     366: [304],
     367: [322],
     368: [323],
     369: [364, 164, 123],
     370: [364, 310, 123],
     371: [326, 165],
     372: [327, 165],
     373: [331, 366, 126],
     374: [396, 366, 126],
     375: [343, 260],
     376: [347, 260],
     377: [414, 328, 232, 288, 281],
     378: [410, 329, 232, 288, 280],
     379: [375],
     380: [376],
     381: [377],
     382: [378],
     383: [369, 350],
     384: [370, 344],
     385: [381],
     386: [382],
     387: [385],
     388: [386],
     389: [191, 47],
     390: [191, 308],
     391: [383],
     392: [384],
     393: [373, 48, 125],
     394: [373, 252, 125],
     395: [393, 349],
     396: [394, 349],
     397: [345],
     398: [346],
     399: [352],
     400: [390],
     401: [395, 49, 126],
     402: [396, 49, 126],
     403: [401, 48, 231],
     404: [402, 48, 231],
     405: [395],
     406: [374, 164, 123],
     407: [374, 310, 123],
     408: [406, 348],
     409: [406, 350],
     410: [408],
     411: [409],
     412: [407, 344],
     413: [407, 348],
     414: [412],
     415: [413]}

.. code-block:: Python

    s = WordMorphism(d)
    s.rauzy_fractal_plot()


.. image:: media/shap3.png
  :scale: 70 %



