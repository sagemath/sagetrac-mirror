.. -*- coding: utf-8 -*-
.. _language_automaton:

======================================
Automata and Rational language in Sage
======================================

.. MODULEAUTHOR:: Paul Mercat, and Dominique Benielli
                  Labex Archimede - I2M -
                  AMU Aix-Marseille Universite 
                   

Automata and Rational language
------------------------------

Automata are in a way machines that can realize linear time calculation only requiring a fine memory. For more details see [Ca].

Automata
~~~~~~~~


Definition of an automaton
^^^^^^^^^^^^^^^^^^^^^^^^^^

We call automaton a quintuplet :math:`A := (\Sigma,\mathrm{Q},\mathrm{T},\mathrm{I},\mathrm{F})`, where

    - :math:`\Sigma` is a finite set called alphabet,
    - :math:`\mathrm{Q}` is a finite set of states,
    - :math:`\mathrm{T} \subseteq \mathrm{Q} \times \Sigma \times \mathrm{Q}` is a finite set of transitions,
    - :math:`\mathrm{I} \subseteq \mathrm{Q}` is a finite set of initial states,
    - :math:`\mathrm{F} \subseteq \mathrm{Q}` is a finite set of final states.

The automaton is deterministic if 
    - :math:`\sharp \, \mathrm{I} = 1` and 
    - :math:`\left[ \left( p, a, q \right) \in \mathrm{T} \quad and  \quad \left(p, a, r \right) \in \mathrm{T} \right] \Rightarrow q = r`


So, when the automaton :math:`A` is deterministic, :math:`\mathrm{T}` is the graph of a partial function
of transition :math:`\mathrm{Q} \times \Sigma \rightarrow \mathrm{Q}`, and there is only one initial state.
We will sometimes consider infinite automata, that is, automata for
which set of states :math:`\mathrm{Q}` is infinite.

.. NOTE::

    We often denote :math:`p \overset{a}{\rightarrow} q  \quad` if
    :math:`\quad \left( p, a, q \right) \in \mathrm{T}`.

.. NOTE::

    For an alphabet :math:`\Sigma`, we note :math:`\Sigma^* := \Sigma^{(\mathbb N)}` the set of finish words. 
    For a word :math:`u \in \Sigma^{*}`, we note :math:`u^* := \cup_{n \in \mathbb N} = \{ u^n \}^*`.

Graphical representation
^^^^^^^^^^^^^^^^^^^^^^^^

Automata are represented as graphs whose edges are labeled by letters of the alphabet. 
On the drawings in this section, the initial state is in bold, and the final states are the circles drawn with a double line
Determinist Automaton can be created in sage by the use of :class:`sage.combinat.words.FastAutomaton` as follow::

    sage: a = FastAutomaton([(0,0,'(0,0)'),(0,0,'(1,1)'),(0,3,'(1,0)'),(1,2,'(0,1)'),(2,0,'(0,1)'),(2,1,'(1,1)'),(2,1,'(0,0)'),(3,4,'(0,1)'),(4,3,'(0,0)'),(4,0,'(1,0)')])
    sage: a.set_final_states([0])
    sage: a.set_initial_state(0)
    sage: a.add_edge(0,'(1,0)',1)
    sage: a.plot().show()

.. PLOT::

    a = FastAutomaton([(0,0,'(0,0)'),(0,0,'(1,1)'),(0,3,'(1,0)'),(1,2,'(0,1)'),(2,0,'(0,1)'),(2,1,'(1,1)'),(2,1,'(0,0)'),(3,4,'(0,1)'),(4,3,'(0,0)'),(4,0,'(1,0)')])
    a.set_final_states([0])
    a.set_initial_state(0)
    a.add_edge(0,'(1,0)',1)
    sphinx_plot(a)

Automaton with states \{0, 1, 2, 3, 4\}, alphabet \{(0,0), (0,1), (1,0), (1,1)\}, set of inital states \{0\}, and set of final states \{0\}.

.. PLOT::
   :width: 50%

    a = FastAutomaton([(0,0,'*'),(0,1,'0'),(0,3,'1'),(1,2,'1'),(2,0,'1'),(2,1,'*'),(4,0,'0'),(4,3,'*'),(3,4,'0')])
    a.set_final_states([0])
    a.set_initial_state(0)
    sphinx_plot(a)

Automaton with states  \{0, 1, 2, 3, 4\},  alphabet \{0, 1, *\}, set of inital states \{0\} and set of final states \{0\}.

.. PLOT::

    a = FastAutomaton([(0,0,'(0,0)'),(0,1,'(1,1)'),(0,3,'(0,1)'),(0,5,'(1,0)'),(3,4,'(0,1)'),(4,2,'(1,0)'),(2,1,'(1,1)'),(1,5,'(1,0)'),(5,6,'(0,1)'),(6,5,'(0,0)'),(6,5,'(1,1)')])
    a.add_edge(1,'(1,1)',1)
    a.add_edge(1,'(0,0)',2)
    a.add_edge(4,'(0,0)',3)
    a.add_edge(4,'(1,1)',3)
    a.set_final_states([0,1,2])
    a.set_initial_state(0)
    sphinx_plot(a)

Automaton of states \{0, 1, 2, 3, 4, 5, 6\},  alphabet \{(0,0), (0,1), (1,0), (1,1)\}, for inital state \{0\} and finals states \{0, 1, 2\}.

Language
~~~~~~~~

Definition: rational language
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A language is a set of words over a given alphabet.
The language recognized by an automaton :math:`A = (\Sigma, Q, T, I, F)` is the set :math:`L_A` of words :math:`a_1 \dots a_n \in \Sigma^*` such that there  exists a path
:math:`\mathrm{I}  \ni q_0 \xrightarrow{a_1} q_1 \xrightarrow{a_2} \dots \dots \xrightarrow{a_{n-1}} q_{n-1} \xrightarrow{a_n} q_n \in \mathrm{F}`
in the automaton :math:`A` from an initial state to an end state. 

A word :math:`u \in \Sigma^*` is recognized  by the automaton  :math:`A` if we have :math:`u \in L_A`.

A word  $a_1 \dots a_n$ is therefore recognized by the automaton :math:`A` if there exists a path in the graph, labeled by  $a_1, a_2, \dots, a_n$, starting from an initial state and ending to a final state.

.. note::

    If the automaton is deterministic, the path is determined by the sequence of labels.

Examples
^^^^^^^^
some examples of automaton.

.. PLOT::
   :width: 50%

    a = FastAutomaton([(0, 0,'0'),(0, 1, '1'),(1, 0, '1'), (1, 2, '0'), (2, 1, '0'), (2, 2, '1')])
    a.set_final_states([0])
    a.set_initial_state(0)
    sphinx_plot(a)

The above automaton recognize all the numbers written in binaries that are divisible by 3.

.. PLOT::
   :width: 50%

    a = FastAutomaton([(0,1,'a'),(1,2,'b'),(2,0,'a')])
    a.set_final_states([1])
    a.set_initial_state(0)
    sphinx_plot(a)

The above automaton recognize the set of words of the form :math:`a(baa)^n`.

.. PLOT::

    a = FastAutomaton([(0,1,'l'),(1,2,'a'),(2,3,'p') ,(3,4,'i'),(4,10,'n'),(0,5,'l'),(5,6,'a'),(6,7,'i'),(7,8,'t'),(8,9,'u'),(9,11,'e') ])
    a.set_final_states([10,11])
    a.set_initial_state(0)
    b= NFastAutomaton(a)
    b.add_edge(0,'l',1)
    sphinx_plot(b)

The above non deterministic automaton recognize the set of words
\{lapin, laitue\}. Obtained with the followed code and the class :class:`sage.combinat.words.NFastAutomaton`::

    sage: a = FastAutomaton([(0,1,'l'),(1,2,'a'),(2,3,'p') ,(3,4,'i'),(4,10,'n'),(0,5,'l'),(5,6,'a'),(6,7,'i'),(7,8,'t'),(8,9,'u'),(9,11,'e')])
    sage: a.set_final_states([10,11])
    sage: a.set_initial_state(0)
    sage: b = NFastAutomaton(a)
    sage: b.add_edge(0,'l',1)
    sage: b.plot().show()

Equivalent automata
^^^^^^^^^^^^^^^^^^^

Two automata :math:`A` and :math:`A'` are equivalent if they recognize the same language  $L_A = L_{A'}$.
Automaton equivalent to the previous one is::

    sage: c = b.determinise()
    sage: c.plot().show()

.. PLOT::

    a = FastAutomaton([(0,1,'l'),(1,2,'a'),(2,3,'i') ,(3,5,'t'),(5,7,'u'),(7,9,'e'),(2,4,'p'),(4,6,'i'),(6,8,'n') ])
    a.set_final_states([8,9])
    a.set_initial_state(0)
    sphinx_plot(a)

.. NOTE::

    Any automaton is equivalent to a deterministic automaton.

Minimal automata
^^^^^^^^^^^^^^^^

   A minimal automaton of an automaton :math:`A` (or the minimal automaton of the corresponding language) is a deterministic automaton :math:`A '`, equivalent to :math:`A`,
   and having a minimal number of vertices for these properties.

.. NOTE::

   The minimal automaton is unique. Moreover, if the automaton :math:`A` is deterministic,
   then the minimal automaton is obtained like the quotient of the automaton :math:`A` by an equivalence
   relation consisting of identifying vertices between them.

The minimal automaton of the language \{lapin, laitue\} is the following::

    sage: d = c.minimise()
    sage: c.plot().show()

.. PLOT::

    a = FastAutomaton([(7,6,'l'),(6,5,'a'),(5,1,'i') ,(1,8,'t'),(8,2,'u'),(2,0,'e'),(5,4,'p'),(4,3,'i'),(3,0,'n') ])
    a.set_final_states([0])
    a.set_initial_state(7)
    sphinx_plot(a)


Transpose automaton
^^^^^^^^^^^^^^^^^^^

The transposed (or the mirror) automaton of an automaton :math:`A := (\Sigma,\mathrm{Q},\mathrm{T},\mathrm{I},\mathrm{F})` is the automaton

.. MATH::
    A^t := (\Sigma, \mathrm{Q}, \mathrm{T}^t, \mathrm{F}, \mathrm{I})
    \text{ where } \mathrm{T}^t := \{ (p, a, q) \in \mathrm{Q} \times \Sigma \times \mathrm{Q}  |  (q, a, p) \in \mathrm{T} \}

.. NOTE::

   The language recognized by the transposed automaton :math:`A^t` is the transpose of the recognized language by the
   initial automaton :math:`A`.

The transposed of the minimal automaton of the language \{lapin, laitue\} is::

    sage: b = a.transpose()
    sage: b.plot().show()

.. PLOT::

    a = FastAutomaton([(7,6,'l'),(6,5,'a'),(5,1,'i') ,(1,8,'t'),(8,2,'u'),(2,0,'e'),(5,4,'p'),(4,3,'i'),(3,0,'n') ])
    a.set_final_states([0])
    a.set_initial_state(7)
    b = a.transpose()
    sphinx_plot(b)

Emonded automaton
^^^^^^^^^^^^^^^^^

The emonded automaton is the automaton restricted to
states that are reachable from an initial state, and from which we can go to a final state.
An automaton is emonded if it is equal to its emonded.

.. NOTE::

    An automaton (possibly infinite) deterministic emonded, and with a deterministic transposed is minimal.
    In particular, if it is infinite, the language that it recognizes is not rational.

Example of non-emonded automaton::

    sage: a = FastAutomaton([(0,0,'(0,0)'),(0,0,'(1,1)'),(0,3,'(1,0)'),(1,2,'(0,1)'),(2,0,'(0,1)'),(2,1,'(1,1)'),(2,1,'(0,0)'),(3,4,'(0,1)'),(4,3,'(0,0)'),(4,0,'(1,0)')])
    sage: a.set_final_states([0])
    sage: a.set_initial_state(0)
    sage: a.add_edge(0,'(1,0)',1)
    sage: a.plot().show()

.. PLOT::

    a = FastAutomaton([(0,0,'(0,0)'),(0,0,'(1,1)'),(0,3,'(1,0)'),(1,2,'(0,1)'),(2,0,'(0,1)'),(2,1,'(1,1)'),(2,1,'(0,0)'),(3,4,'(0,1)'),(4,3,'(0,0)'),(4,0,'(1,0)')])
    a.set_final_states([0])
    a.set_initial_state(0)
    a.add_edge(0,'(1,0)',1)
    sphinx_plot(a)

And the corresponding emonded automaton::

    sage: b = a.emonde()
    sage: b.plot().show()

This automaton can be saw below:

.. PLOT::
   :width: 50%
      
    a = FastAutomaton([(0,0,'(0,0)'),(0,0,'(1,1)'),(0,3,'(1,0)'),(1,2,'(0,1)'),(2,0,'(0,1)'),(2,1,'(1,1)'),(2,1,'(0,0)'),(3,4,'(0,1)'),(4,3,'(0,0)'),(4,0,'(1,0)')])
    a.set_final_states([0])
    a.set_initial_state(0)
    a.add_edge(0,'(1,0)',1)
    b = a.emonde()
    sphinx_plot(b)

The emonded example automaton.
                   