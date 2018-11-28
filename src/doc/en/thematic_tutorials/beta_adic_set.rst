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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   
    
Let's take :math:`\sigma_+` and :math:`\sigma_-` some embeddings of :math:`  Q(\lambda)` into the spaces :math:`E_+` and :math:`E_-` respectively.
We will also denote by :math:`\sigma_\beta` the maximal real embedding when :math:`\beta is a Perron number.
      
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



.. tikz::
     \draw (0,0) -- (12, 0);
     \draw (0.000, -.1) -- (0.000, .1);
     \draw (1.000, -.1) -- (1.000, .1);
     \draw (1.839, -.1) -- (1.839, .1);
     \draw (2.839, -.1) -- (2.839, .1);
     \draw (3.383, -.1) -- (3.383, .1);
     \draw (4.383, -.1) -- (4.383, .1);
     \draw (5.222, -.1) -- (5.222, .1);
     \draw (6.222, -.1) -- (6.222, .1);
     \draw (7.222, -.1) -- (7.222, .1);
     \draw (8.062, -.1) -- (8.062, .1);
     \draw (9.062, -.1) -- (9.062, .1);
     \draw (9.605, -.1) -- (9.605, .1);
     \draw (10.605, -.1) -- (10.605, .1);
     \draw (11.445, -.1) -- (11.445, .1);

     \draw (0,3) -- (12, 3);
     \draw (0.000, 2.9) -- (0.000, 3.1);
     \draw (1.000, 2.9) -- (1.000, 3.1);
     \draw (1.839, 2.9) -- (1.839, 3.1);
     \draw (2.839, 2.9) -- (2.839, 3.1);
     \draw (3.383, 2.9) -- (3.383, 3.1);
     \draw (4.383, 2.9) -- (4.383, 3.1);
     \draw (5.222, 2.9) -- (5.222, 3.1);
     \draw (6.222, 2.9) -- (6.222, 3.1);
     \draw (7.222, 2.9) -- (7.222, 3.1);
     \draw (8.062, 2.9) -- (8.062, 3.1);
     \draw (9.062, 2.9) -- (9.062, 3.1);
     \draw (9.605, 2.9) -- (9.605, 3.1);
     \draw (10.605, 2.9) -- (10.605, 3.1);
     \draw (11.445, 2.9) -- (11.445, 3.1);
	
     \draw (0.000, 3) -- (0.000, 0);
     \draw (1.000, 3) -- (1.839, 0);
     \draw (1.839, 3) -- (3.383, 0);
     \draw (2.839, 3) -- (5.222, 0);
     \draw (3.383, 3) -- (6.222, 0);
     \draw (4.383, 3) -- (8.062, 0);
     \draw (5.222, 3) -- (9.605, 0);
     \draw (6.222, 3) -- (11.445, 0); 
     \draw (0.500, -.2) node {a};
     \draw (1.420, -.2) node {b};
     \draw (2.339, -.2) node {a};
     \draw (3.111, -.2) node {c};
     \draw (3.883, -.2) node {a};
     \draw (4.803, -.2) node {b};
     \draw (5.722, -.2) node {a};
     \draw (6.722, -.2) node {a};
     \draw (7.642, -.2) node {b};
     \draw (8.562, -.2) node {a};
     \draw (9.333, -.2) node {c};
     \draw (10.105, -.2) node {a};
     \draw (11.025, -.2) node {b};
     \draw (0.500, 3.2) node {a};
     \draw (1.420, 3.2) node {b};
     \draw (2.339, 3.2) node {a};
     \draw (3.111, 3.2) node {c};
     \draw (3.883, 3.2) node {a};
     \draw (4.803, 3.2) node {b};
     \draw (5.722, 3.2) node {a};
     \draw (6.722, 3.2) node {a};
     \draw (7.642, 3.2) node {b};
     \draw (8.562, 3.2) node {a};
     \draw (9.333, 3.2) node {c};
     \draw (10.105, 3.2) node {a};
     \draw (11.025, 3.2) node {b};
     \draw[->] (-.3, 3) arc (150:210:3);
     \draw (-.7 ,1.5) node[left] {:math:`\times \lambda`};

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


