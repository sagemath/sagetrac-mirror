"""
Extensive Form Games with N players.

This module implements a class for extensive form games
[NN2007]_. Graphical representations of the games are implemented and solution
algorithms are being developed (with an interface to gambit [MMAT2014]_).

A well known example that can be implemented as an extensive form game is the
battle of the sexes. Consider two players, Celine and Bob. The two are deciding
on how to spend their evening, they can either watch Sports or go and see a
Comedy. Bob would prefer to see a Comedy, Celine would prefer to watch a Sports
movie. Depending on who choses what, there are different payoffs, this can be
demonstrated in tree form.

.. PLOT::
    :width: 500 px

    player_1 = EFG_Player('Bob')
    player_2 = EFG_Player('Celine')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    node_3 = EFG_Node(player_2, {'Sports': leaf_1, 'Comedy': leaf_2}, 'c')
    node_2 = EFG_Node(player_2, {'Sports': leaf_3, 'Comedy': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'Sports': node_3, 'Comedy': node_2}, 'a')
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    battle_of_the_sexes.set_info_set([node_2, node_3])
    p = battle_of_the_sexes.plot()
    sphinx_plot(p)

We can see there are three nodes, one for Bob, two for Celine. Connecting those
nodes are actions. These actions represent choices made by one player, and the
actions then lead on to a node of another player.  So Bob either chooses Sports
or Comedy, and Celine chooses Sports or Comedy. The location of the payoffs
correspond to the leaf of the underlying tree and they show the outcome for
each player.  So if Bob chooses Sports, and Celine chooses Sports, we see the
payoff is (2, 3), which represents Bob getting a payoff of 2, and Celine
getting a payoff of 3. Note that the line between Celine's two nodes indicate
that they are in the same 'information set', in other words, Celine does not
know what Bob has picked. Thus the following game corresponds to a different
situation which is easier for both players to coordinate:

.. PLOT::
    :width: 500 px

    player_1 = EFG_Player('Bob')
    player_2 = EFG_Player('Celine')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    node_3 = EFG_Node(player_2, {'Sports': leaf_1, 'Comedy': leaf_2}, 'c')
    node_2 = EFG_Node(player_2, {'Sports': leaf_3, 'Comedy': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'Sports': node_3, 'Comedy': node_2}, 'a')
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    p = battle_of_the_sexes.plot()
    sphinx_plot(p)

The first game (with the information set) corresponds to the following normal
form game (which are also implemented in Sage)::

    sage: A = matrix([[3, 1], [0, 2]])
    sage: B = matrix([[2, 1], [0, 3]])
    sage: battle_of_the_sexes = NormalFormGame([A, B])
    sage: battle_of_the_sexes
    Normal Form Game with the following utilities: {(0, 1): [1, 1],
    (1, 0): [0, 0], (0, 0): [3, 2], (1, 1): [2, 3]}

To generate an extensive form game we need to generate the nodes and assign
them to players as well as describing the actions they have and to which node
each action goes. As such it makes sense to start with the terminal nodes of
the tree, but the initial step is to create players as we will need to create
leafs which will map players to utilities::

    sage: player_1, player_2 = EFG_Player('Bob'), EFG_Player('Celine')
    sage: player_1, player_2
    (EFG Player "Bob", EFG Player "Celine")

Once we have done this, we are ready to create our leafs::

    sage: leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    sage: leaf_1
    EFG Leaf with utilities - (2, 3)
    sage: leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    sage: leaf_2
    EFG Leaf with utilities - (0, 0)
    sage: leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    sage: leaf_3
    EFG Leaf with utilities - (1, 1)
    sage: leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    sage: leaf_4
    EFG Leaf with utilities - (3, 2)

We can then create the parents of these leafs, the general :code:`EFG_Node`
class takes 3 arguments: the player who makes the decision at this node,
a dictionary mapping actions to other nodes, and finally a name for the node::

    sage: node_1 = EFG_Node(player_2, {'Sports': leaf_3,
    ....:                              'Comedy': leaf_4}, 'b')
    sage: node_1
    EFG Node "b"
    sage: node_1.player
    EFG Player "Celine"
    sage: node_1.name
    'b'
    sage: node_2 = EFG_Node(player_2, {'Sports': leaf_1,
    ....:                              'Comedy': leaf_2}, 'c')
    sage: node_2
    EFG Node "c"
    sage: node_2.player
    EFG Player "Celine"
    sage: node_2.name
    'c'

Finally, we create the root of the tree::

    sage: root = EFG_Node(player_1, {'Sports': node_2, 'Comedy': node_1}, 'a')
    sage: root
    EFG Node "a"
    sage: root.player
    EFG Player "Bob"
    sage: root.name
    'a'

The extensive form game can then be created by passing this root (which
recursively has all required information)::

    sage: battle_of_the_sexes = ExtensiveFormGame(root)
    sage: battle_of_the_sexes
    Extensive Form Game with the following underlying tree: {...

By default all nodes are in their own information set. If we plot the tree we
see this::

    sage: battle_of_the_sexes.plot()
    Graphics object consisting of 20 graphics primitives

Here is the output (this is the same tree as above):

.. PLOT::
    :width: 500 px

    player_1 = EFG_Player('Bob')
    player_2 = EFG_Player('Celine')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    node_3 = EFG_Node(player_2, {'Sports': leaf_1, 'Comedy': leaf_2}, 'c')
    node_2 = EFG_Node(player_2, {'Sports': leaf_3, 'Comedy': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'Sports': node_3, 'Comedy': node_2}, 'a')
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    p = battle_of_the_sexes.plot()
    sphinx_plot(p)

An extensive form game where all nodes are in their own information set is said
to have 'perfect information'. The game we have so far still has perfect
information::

    sage: battle_of_the_sexes.info_sets
    [[EFG Node "a"],
     [EFG Node "b"],
     [EFG Node "c"]]
    sage: battle_of_the_sexes.perfect_info()
    True

To set the information sets as described above we use the :code:`set_info_set`
method::

    sage: battle_of_the_sexes.set_info_set([node_1, node_2])
    sage: battle_of_the_sexes.info_sets
    [[EFG Node "a"],
     [EFG Node "b",
      EFG Node "c"]]

Now the game does not have perfect information::

    sage: battle_of_the_sexes.perfect_info()
    False

Information sets are demonstrated visually on the graph we plot by default::

    sage: battle_of_the_sexes.plot()
    Graphics object consisting of 21 graphics primitives

Which will be plotted as follows:

.. PLOT::
    :width: 500 px

    player_1 = EFG_Player('Bob')
    player_2 = EFG_Player('Celine')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    node_3 = EFG_Node(player_1, {'Sports': leaf_1, 'Comedy': leaf_2}, 'c')
    node_2 = EFG_Node(player_1, {'Sports': leaf_3, 'Comedy': leaf_4}, 'b')
    node_1 = EFG_Node(player_2, {'Sports': node_3, 'Comedy': node_2}, 'a')
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    battle_of_the_sexes.set_info_set([node_2, node_3])
    p = battle_of_the_sexes.plot()
    sphinx_plot(p)

If we wish to remove the information sets from the plot, we simply set
``view_info_set_lines`` to be ``False``::

    sage: battle_of_the_sexes.plot(view_info_set_lines = False)
    Graphics object consisting of 20 graphics primitives

Which will be plotted as follows:

.. PLOT::
    :width: 500 px

    player_1 = EFG_Player('Bob')
    player_2 = EFG_Player('Celine')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    node_3 = EFG_Node(player_1, {'Sports': leaf_1, 'Comedy': leaf_2}, 'c')
    node_2 = EFG_Node(player_1, {'Sports': leaf_3, 'Comedy': leaf_4}, 'b')
    node_1 = EFG_Node(player_2, {'Sports': node_3, 'Comedy': node_2}, 'a')
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    battle_of_the_sexes.set_info_set([node_2, node_3])
    p = battle_of_the_sexes.plot(view_info_set_lines = False)
    sphinx_plot(p)

In some situations, it is helpful to view circles around nodes of information
sets, where the color of those circles are dependent on what information set
they are in, we can view them by using the plot option
``view_info_set_circles`` which is set to ``False`` be default::

    sage: battle_of_the_sexes.plot(view_info_set_circles = True)
    Graphics object consisting of 23 graphics primitives

Which will be plotted as follows:

.. PLOT::
    :width: 500 px

    player_1 = EFG_Player('Bob')
    player_2 = EFG_Player('Celine')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 3})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 3, player_2: 2})
    node_3 = EFG_Node(player_1, {'Sports': leaf_1, 'Comedy': leaf_2}, 'c')
    node_2 = EFG_Node(player_1, {'Sports': leaf_3, 'Comedy': leaf_4}, 'b')
    node_1 = EFG_Node(player_2, {'Sports': node_3, 'Comedy': node_2}, 'a')
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    battle_of_the_sexes.set_info_set([node_2, node_3])
    p = battle_of_the_sexes.plot(view_info_set_circles = True)
    sphinx_plot(p)

We can remove information sets which automatically sets all the nodes in that
information set to be in their own information set::

    sage: battle_of_the_sexes.remove_info_set([node_1, node_2])
    sage: battle_of_the_sexes.info_sets
    [[EFG Node "a"],
     [EFG Node "b"],
     [EFG Node "c"]]

With the game set up with the information sets we want, we can then obtain Nash
Equilbria for this game. We do this with the ``obtain_nash`` method. This
currently uses Gambit and the LCP algorithm, and hence needs gambit installed.
To install it type ``sage -i gambit`` at the command line.

We can then obtain the Nash Equilbria for the game, the following shows three
Nash Equilbria profiles obtained. We can do this with the ``obtain_nash``
method. This method returns the Nash Equilibria in a format which will show a
list of tuples containing Nash Equilibria, where in each tuple are tuples with
information sets and the corresponding probability distributions mapping
actions to their probabilities.::

    sage: battle_of_the_sexes.obtain_nash()  # optional - gambit
    [((('a',), {'Comedy': 0.0, 'Sports': 1.0}),
      (('b',), {'Comedy': 0.0, 'Sports': 1.0}),
      (('c',), {'Comedy': 0.0, 'Sports': 1.0})),
     ((('a',), {'Comedy': 0.0, 'Sports': 1.0}),
      (('b',), {'Comedy': 0.5, 'Sports': 0.5}),
      (('c',), {'Comedy': 0.0, 'Sports': 1.0})),
     ((('a',), {'Comedy': 1.0, 'Sports': 0.0}),
      (('b',), {'Comedy': 1.0, 'Sports': 0.0}),
      (('c',), {'Comedy': 0.0, 'Sports': 1.0}))]

We can also check if the game is a constant sum game, with the ``is_constant_sum``
method, which simply checks that the sum of utilities for a Leaf within the game
is the same for every Leaf within that game::

    sage: battle_of_the_sexes.is_constant_sum()
    False

Here is a different example with more players. This is a three player game
where the player Jack has the options to either choose A or B, Jane has the
options to chose C or D, and then Josh has the option to chose E or F.
Initally everyone knows what the others have chosen::

    sage: jack = EFG_Player('Jack')
    sage: jane = EFG_Player('Jane')
    sage: josh = EFG_Player('Josh')
    sage: leaf_8 = EFG_Leaf({jack : 0, jane: 1, josh: 3}, 'Leaf 8')
    sage: leaf_7 = EFG_Leaf({jack : 1, jane: 0, josh: 1}, 'Leaf 7')
    sage: leaf_6 = EFG_Leaf({jack : 2, jane: 4, josh: 2}, 'Leaf 6')
    sage: leaf_5 = EFG_Leaf({jack : 2, jane: 1, josh: 3}, 'Leaf 5')
    sage: leaf_4 = EFG_Leaf({jack : 0, jane: 1, josh: 0}, 'Leaf 4')
    sage: leaf_3 = EFG_Leaf({jack : 1, jane: 0, josh: 2}, 'Leaf 3')
    sage: leaf_2 = EFG_Leaf({jack : 2, jane: 4, josh: 1}, 'Leaf 2')
    sage: leaf_1 = EFG_Leaf({jack : 2, jane: 1, josh: 2}, 'Leaf 1')
    sage: node_6 = EFG_Node(josh, {'E': leaf_7, 'F': leaf_8}, 'Node 6')
    sage: node_5 = EFG_Node(josh, {'E': leaf_5, 'F': leaf_6}, 'Node 5')
    sage: node_4 = EFG_Node(josh, {'E': leaf_3, 'F': leaf_4}, 'Node 4')
    sage: node_3 = EFG_Node(josh, {'E': leaf_1, 'F': leaf_2}, 'Node 3')
    sage: node_2 = EFG_Node(jane, {'C': node_4, 'D': node_5}, 'Node 2')
    sage: node_1 = EFG_Node(jane, {'C': node_3, 'D': node_6}, 'Node 1')
    sage: root = EFG_Node(jack, {'A': node_1, 'B': node_2}, 'Root')
    sage: jjj_game = ExtensiveFormGame(root)

We can see how the information sets are connected by plotting the information
set graph, which will show a DiGraph where the vertices are the information
sets. We can do this by using ``plot_info_sets``::

    sage: jjj_game.plot_info_sets()
    Graphics object consisting of 20 graphics primitives

Which will be plotted as follows:

.. PLOT::
    :width: 500 px

    jack = EFG_Player('Jack')
    jane = EFG_Player('Jane')
    josh = EFG_Player('Josh')
    leaf_8 = EFG_Leaf({jack : 0, jane: 1, josh: 3}, 'Leaf 8')
    leaf_7 = EFG_Leaf({jack : 1, jane: 0, josh: 1}, 'Leaf 7')
    leaf_6 = EFG_Leaf({jack : 2, jane: 4, josh: 2}, 'Leaf 6')
    leaf_5 = EFG_Leaf({jack : 2, jane: 1, josh: 3}, 'Leaf 5')
    leaf_4 = EFG_Leaf({jack : 0, jane: 1, josh: 0}, 'Leaf 4')
    leaf_3 = EFG_Leaf({jack : 1, jane: 0, josh: 2}, 'Leaf 3')
    leaf_2 = EFG_Leaf({jack : 2, jane: 4, josh: 1}, 'Leaf 2')
    leaf_1 = EFG_Leaf({jack : 2, jane: 1, josh: 2}, 'Leaf 1')
    node_6 = EFG_Node(josh, {'E': leaf_7, 'F': leaf_8}, 'Node 6')
    node_5 = EFG_Node(josh, {'E': leaf_5, 'F': leaf_6}, 'Node 5')
    node_4 = EFG_Node(josh, {'E': leaf_3, 'F': leaf_4}, 'Node 4')
    node_3 = EFG_Node(josh, {'E': leaf_1, 'F': leaf_2}, 'Node 3')
    node_2 = EFG_Node(jane, {'C': node_4, 'D': node_5}, 'Node 2')
    node_1 = EFG_Node(jane, {'C': node_3, 'D': node_6}, 'Node 1')
    root = EFG_Node(jack, {'A': node_1, 'B': node_2}, 'Root')
    jjj_game = ExtensiveFormGame(root)
    p = jjj_game.plot_info_sets()
    sphinx_plot(p)

Here is the corresponding Extensive Form Game tree:

.. PLOT::
    :width: 500 px

    jack = EFG_Player('Jack')
    jane = EFG_Player('Jane')
    josh = EFG_Player('Josh')
    leaf_8 = EFG_Leaf({jack : 0, jane: 1, josh: 3}, 'Leaf 8')
    leaf_7 = EFG_Leaf({jack : 1, jane: 0, josh: 1}, 'Leaf 7')
    leaf_6 = EFG_Leaf({jack : 2, jane: 4, josh: 2}, 'Leaf 6')
    leaf_5 = EFG_Leaf({jack : 2, jane: 1, josh: 3}, 'Leaf 5')
    leaf_4 = EFG_Leaf({jack : 0, jane: 1, josh: 0}, 'Leaf 4')
    leaf_3 = EFG_Leaf({jack : 1, jane: 0, josh: 2}, 'Leaf 3')
    leaf_2 = EFG_Leaf({jack : 2, jane: 4, josh: 1}, 'Leaf 2')
    leaf_1 = EFG_Leaf({jack : 2, jane: 1, josh: 2}, 'Leaf 1')
    node_6 = EFG_Node(josh, {'E': leaf_7, 'F': leaf_8}, 'Node 6')
    node_5 = EFG_Node(josh, {'E': leaf_5, 'F': leaf_6}, 'Node 5')
    node_4 = EFG_Node(josh, {'E': leaf_3, 'F': leaf_4}, 'Node 4')
    node_3 = EFG_Node(josh, {'E': leaf_1, 'F': leaf_2}, 'Node 3')
    node_2 = EFG_Node(jane, {'C': node_4, 'D': node_5}, 'Node 2')
    node_1 = EFG_Node(jane, {'C': node_3, 'D': node_6}, 'Node 1')
    root = EFG_Node(jack, {'A': node_1, 'B': node_2}, 'Root')
    jjj_game = ExtensiveFormGame(root)
    p = jjj_game.plot()
    sphinx_plot(p)

If we change the game so now Josh doesn't know what Jane has chosen,
the information set graph changes to show the changes in relationships between
the information sets::

    sage: jjj_game.set_info_set([node_1, node_2])
    sage: jjj_game.plot_info_sets()
    Graphics object consisting of 17 graphics primitives

Which will be plotted as follows:

.. PLOT::
    :width: 500 px

    jack = EFG_Player('Jack')
    jane = EFG_Player('Jane')
    josh = EFG_Player('Josh')
    leaf_8 = EFG_Leaf({jack : 0, jane: 1, josh: 3}, 'Leaf 8')
    leaf_7 = EFG_Leaf({jack : 1, jane: 0, josh: 1}, 'Leaf 7')
    leaf_6 = EFG_Leaf({jack : 2, jane: 4, josh: 2}, 'Leaf 6')
    leaf_5 = EFG_Leaf({jack : 2, jane: 1, josh: 3}, 'Leaf 5')
    leaf_4 = EFG_Leaf({jack : 0, jane: 1, josh: 0}, 'Leaf 4')
    leaf_3 = EFG_Leaf({jack : 1, jane: 0, josh: 2}, 'Leaf 3')
    leaf_2 = EFG_Leaf({jack : 2, jane: 4, josh: 1}, 'Leaf 2')
    leaf_1 = EFG_Leaf({jack : 2, jane: 1, josh: 2}, 'Leaf 1')
    node_6 = EFG_Node(josh, {'E': leaf_7, 'F': leaf_8}, 'Node 6')
    node_5 = EFG_Node(josh, {'E': leaf_5, 'F': leaf_6}, 'Node 5')
    node_4 = EFG_Node(josh, {'E': leaf_3, 'F': leaf_4}, 'Node 4')
    node_3 = EFG_Node(josh, {'E': leaf_1, 'F': leaf_2}, 'Node 3')
    node_2 = EFG_Node(jane, {'C': node_4, 'D': node_5}, 'Node 2')
    node_1 = EFG_Node(jane, {'C': node_3, 'D': node_6}, 'Node 1')
    root = EFG_Node(jack, {'A': node_1, 'B': node_2}, 'Root')
    jjj_game = ExtensiveFormGame(root)
    jjj_game.set_info_set([node_1, node_2])
    p = jjj_game.plot_info_sets()
    sphinx_plot(p)

Here is the corresponding Extensive Form Game tree:

.. PLOT::
    :width: 500 px

    jack = EFG_Player('Jack')
    jane = EFG_Player('Jane')
    josh = EFG_Player('Josh')
    leaf_8 = EFG_Leaf({jack : 0, jane: 1, josh: 3}, 'Leaf 8')
    leaf_7 = EFG_Leaf({jack : 1, jane: 0, josh: 1}, 'Leaf 7')
    leaf_6 = EFG_Leaf({jack : 2, jane: 4, josh: 2}, 'Leaf 6')
    leaf_5 = EFG_Leaf({jack : 2, jane: 1, josh: 3}, 'Leaf 5')
    leaf_4 = EFG_Leaf({jack : 0, jane: 1, josh: 0}, 'Leaf 4')
    leaf_3 = EFG_Leaf({jack : 1, jane: 0, josh: 2}, 'Leaf 3')
    leaf_2 = EFG_Leaf({jack : 2, jane: 4, josh: 1}, 'Leaf 2')
    leaf_1 = EFG_Leaf({jack : 2, jane: 1, josh: 2}, 'Leaf 1')
    node_6 = EFG_Node(josh, {'E': leaf_7, 'F': leaf_8}, 'Node 6')
    node_5 = EFG_Node(josh, {'E': leaf_5, 'F': leaf_6}, 'Node 5')
    node_4 = EFG_Node(josh, {'E': leaf_3, 'F': leaf_4}, 'Node 4')
    node_3 = EFG_Node(josh, {'E': leaf_1, 'F': leaf_2}, 'Node 3')
    node_2 = EFG_Node(jane, {'C': node_4, 'D': node_5}, 'Node 2')
    node_1 = EFG_Node(jane, {'C': node_3, 'D': node_6}, 'Node 1')
    root = EFG_Node(jack, {'A': node_1, 'B': node_2}, 'Root')
    jjj_game = ExtensiveFormGame(root)
    jjj_game.set_info_set([node_1, node_2])
    p = jjj_game.plot()
    sphinx_plot(p)

Here is an example of a larger game: a 3 player rock paper scissors variant
where all players play at the same time.

Here is a method to compute the utilities of a particular strategy profile::

    sage: def strategy_to_utility(P1, P2, P3):
    ....:     what_beats = {'R': 'P', 'P': 'S', 'S': 'R'}
    ....:     what_loses = {'R': 'S', 'P': 'R', 'S': 'P'}
    ....:     p1_score = [P2, P3].count(what_loses[P1])- [P2, P3].count(what_beats[P1])
    ....:     p2_score = [P1, P3].count(what_loses[P2])- [P1, P3].count(what_beats[P2])
    ....:     p3_score = [P2, P1].count(what_loses[P3])- [P2, P1].count(what_beats[P3])
    ....:     return p1_score, p2_score, p3_score

Here we get all possible strategy profiles at each level of the tree::

    sage: import itertools
    sage: all_possible_profiles = list(itertools.product(['R','P','S'], repeat=3))
    sage: all_second_level_profiles = list(itertools.product(['R','P','S'], repeat=2))
    sage: all_first_level_profiles = [('R',), ('P',), ('S',)]

Creating our players::

    sage: p1 = EFG_Player('P1')
    sage: p2 = EFG_Player('P2')
    sage: p3 = EFG_Player('P3')

Putting nodes in a dictionary by level so as to recursively build up the tree::

    sage: third_level_dictionary = {profile:EFG_Leaf({p1:strategy_to_utility(*profile)[0],
    ....: p2:strategy_to_utility(*profile)[1],
    ....: p3:strategy_to_utility(*profile)[2]}) for profile in all_possible_profiles}
    sage: second_level_dictionary = {profile:EFG_Node(p3, {'R':third_level_dictionary[profile+('R',)],
    ....: 'P':third_level_dictionary[profile+('P',)],
    ....: 'S':third_level_dictionary[profile+('S',)]}) for profile in all_second_level_profiles}
    sage: first_level_dictionary = {profile:EFG_Node(p2, {'R':second_level_dictionary[profile+('R',)],
    ....: 'P':second_level_dictionary[profile+('P',)],
    ....: 'S':second_level_dictionary[profile+('S',)]}) for profile in all_first_level_profiles}

Building the node of the tree::

    sage: root = EFG_Node(p1, {'P':first_level_dictionary[('P',)],
    ....: 'R':first_level_dictionary[('R',)],
    ....: 'S':first_level_dictionary[('S',)]})

Building the tree::

    sage: g = ExtensiveFormGame(root)

Setting the information sets::

    sage: info_set_for_p3 = second_level_dictionary.values()
    sage: info_set_for_p2 = first_level_dictionary.values()
    sage: g.set_info_set(info_set_for_p2)
    sage: g.set_info_set(info_set_for_p3)

    sage: p = g.plot()
    sage: p
    Graphics object consisting of 121 graphics primitives

If you would like to see this tree run the follow but beware it's a large
plot!::

    sage: p.show(figsize=[20, 50])  # modifies the size of the plot

We can create a sage ExtensiveFormGame from a gambit game. Here we create a
gambit game with two players, Paul and Mary, where both Paul and Mary have two
actions "Go" and "Stop".

To setup this game we must first import gambit and setup a new ``Game`` in the
form of a tree::

    sage: import gambit  # optional - gambit
    sage: gambit_game = gambit.Game.new_tree()  # optional - gambit

We can then assign players to the game::

    sage: gambit_game.players.add("Paul")  # optional - gambit
    <Player [0] 'Paul' in game ''>
    sage: gambit_game.players.add("Mary")  # optional - gambit
    <Player [1] 'Mary' in game ''>

We can then set up the gambit moves, assigning players to nodes and labeling
any actions::

    sage: paul_move = gambit_game.root.append_move(gambit_game.players["Paul"], int(2))  # optional - gambit
    sage: paul_move.actions[int(0)].label = "Go"  # optional - gambit
    sage: paul_move.actions[int(1)].label = "Stop"  # optional - gambit
    sage: mary_move_1 = gambit_game.root.children[int(0)].append_move(gambit_game.players["Mary"], int(2))  # optional - gambit
    sage: mary_move_1.actions[int(0)].label = "Go"  # optional - gambit
    sage: mary_move_1.actions[int(1)].label = "Stop"  # optional - gambit
    sage: mary_move_2 = gambit_game.root.children[int(1)].append_move(gambit_game.players["Mary"], int(2))  # optional - gambit
    sage: mary_move_2.actions[int(0)].label = "Go"  # optional - gambit
    sage: mary_move_2.actions[int(1)].label = "Stop"  # optional - gambit

Then we can set up outcomes and assign them to terminal nodes::

    sage: outcome = gambit_game.outcomes.add()  # optional - gambit
    sage: outcome[int(0)] = int(1)  # optional - gambit
    sage: outcome[int(1)] = int(6)  # optional - gambit
    sage: gambit_game.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
    sage: outcome = gambit_game.outcomes.add()  # optional - gambit
    sage: outcome[int(0)] = int(6)  # optional - gambit
    sage: outcome[int(1)] = int(2)  # optional - gambit
    sage: gambit_game.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
    sage: outcome = gambit_game.outcomes.add()  # optional - gambit
    sage: outcome[int(0)] = int(4)  # optional - gambit
    sage: outcome[int(1)] = int(0)  # optional - gambit
    sage: gambit_game.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
    sage: outcome = gambit_game.outcomes.add()  # optional - gambit
    sage: outcome[int(0)] = int(6)  # optional - gambit
    sage: outcome[int(1)] = int(2)  # optional - gambit
    sage: gambit_game.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

We can then use this game to generate the ``ExtensiveFormGame`` by simply using
``gambit_game`` as the ``Generator``::

    sage: sage_game = ExtensiveFormGame(gambit_game)  # optional - gambit

Then we can see the game by using the ``plot`` method for
``ExtensiveFormGame``::

    sage: sage_game.plot()  # optional - gambit
    Graphics object consisting of 20 graphics primitives

Here is the output (this is the same tree as above):

.. PLOT::
    :width: 500 px

    import gambit  # optional - gambit
    gambit_game = gambit.Game.new_tree()  # optional - gambit
    gambit_game.players.add("Paul")  # optional - gambit
    gambit_game.players.add("Mary")  # optional - gambit
    paul_move = gambit_game.root.append_move(gambit_game.players["Paul"], int(2))  # optional - gambit
    paul_move.actions[int(0)].label = "Go"  # optional - gambit
    paul_move.actions[int(1)].label = "Stop"  # optional - gambit
    mary_move_1 = gambit_game.root.children[int(0)].append_move(gambit_game.players["Mary"], int(2))  # optional - gambit
    mary_move_1.actions[int(0)].label = "Go"  # optional - gambit
    mary_move_1.actions[int(1)].label = "Stop"  # optional - gambit
    mary_move_2 = gambit_game.root.children[int(1)].append_move(gambit_game.players["Mary"], int(2))  # optional - gambit
    mary_move_2.actions[int(0)].label = "Go"  # optional - gambit
    mary_move_2.actions[int(1)].label = "Stop"  # optional - gambit
    outcome = gambit_game.outcomes.add()  # optional - gambit
    outcome[int(0)] = int(1)  # optional - gambit
    outcome[int(1)] = int(6)  # optional - gambit
    gambit_game.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
    outcome = gambit_game.outcomes.add()  # optional - gambit
    outcome[int(0)] = int(6)  # optional - gambit
    outcome[int(1)] = int(2)  # optional - gambit
    gambit_game.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
    outcome = gambit_game.outcomes.add()  # optional - gambit
    outcome[int(0)] = int(4)  # optional - gambit
    outcome[int(1)] = int(0)  # optional - gambit
    gambit_game.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
    outcome = gambit_game.outcomes.add()  # optional - gambit
    outcome[int(0)] = int(6)  # optional - gambit
    outcome[int(1)] = int(2)  # optional - gambit
    gambit_game.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit
    sage_game = ExtensiveFormGame(gambit_game)  # optional - gambit
    p = sage_game.plot()  # optional - gambit
    sphinx_plot(p)  # optional - gambit

We can convert any sage ``ExtensiveFormGame`` to a Gambit ``Game``, by using
the ``sage_to_gambit`` method::

    sage: gambit_game = sage_game.sage_to_gambit()  # optional - gambit
    sage: type(gambit_game)  # optional - gambit
    <class 'gambit.Game'>

AUTHOR:

- Hannah Lorrimore (07-2015): Original version
"""

from operator import attrgetter
from copy import copy
from parser import Parser
from decimal import Decimal
from fractions import Fraction

from sage.rings.rational import Rational
from sage.rings.real_mpfr import RealLiteral
from sage.rings.integer import Integer
from sage.graphs.all import Graph
from sage.graphs.digraph import DiGraph
from sage.plot.line import line2d
from sage.plot.circle import circle
from sage.graphs.generic_graph import GenericGraph
from sage.plot.colors import rainbow

try:
    from gambit import Game
except ImportError:
    Game = None

class ExtensiveFormGame():
    def __init__(self, generator):
        r"""
        An object representing an Extensive Form Game. Primarily used to
        compute the Nash Equilibria.

        INPUT:

        - ``generator`` - Can be an instance of the Node class which serves as
          the root of the tree. Or it can be a gambit Game which is then
          converted.

        Currently there are two ways of creating an ``ExtensiveFormGame``,
        either through the EFG classes created for it, or through a gambit
        game.

        Here are some examples of a game created with the EFG classes:

        A game with 2 players and 7 nodes::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: leaf_a5 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a6 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a7 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a8 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a1, {'A': leaf_a1, 'B': leaf_a2})
            sage: node_a2 = EFG_Node(player_a1, {'A': leaf_a3, 'B': leaf_a4})
            sage: node_a3 = EFG_Node(player_a1, {'A': leaf_a5, 'B': leaf_a6})
            sage: node_a4 = EFG_Node(player_a1, {'A': leaf_a7, 'B': leaf_a8})
            sage: node_a5 = EFG_Node(player_a2, {'C': node_a1, 'D': node_a2})
            sage: node_a6 = EFG_Node(player_a2, {'C': node_a3, 'D': node_a4})
            sage: root_a = EFG_Node(player_a1, {'A': node_a5, 'B': node_a6})
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.tree
            Graph on 15 vertices

        The generated tree has a variety of attributes::

            sage: egame_a.players
            [EFG Player "Player 1", EFG Player "Player 2"]
            sage: egame_a.nodes
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Node 3",
             EFG Node "Node 4",
             EFG Node "Node 5",
             EFG Node "Node 6",
             EFG Node "Tree Root"]
            sage: egame_a.info_sets
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Node 3"],
             [EFG Node "Node 4"],
             [EFG Node "Node 5"],
             [EFG Node "Node 6"],
             [EFG Node "Tree Root"]]
             sage: egame_a.tree
             Graph on 15 vertices

        We can plot the above tree through the ``plot`` method::

            sage: egame_a.plot()
            Graphics object consisting of 44 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_a1 = EFG_Player('Player 1')
            player_a2 = EFG_Player('Player 2')
            leaf_a1 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            leaf_a2 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            leaf_a3 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            leaf_a4 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            leaf_a5 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            leaf_a6 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            leaf_a7 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            leaf_a8 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            node_a1 = EFG_Node(player_a1, {'A': leaf_a1, 'B': leaf_a2})
            node_a2 = EFG_Node(player_a1, {'A': leaf_a3, 'B': leaf_a4})
            node_a3 = EFG_Node(player_a1, {'A': leaf_a5, 'B': leaf_a6})
            node_a4 = EFG_Node(player_a1, {'A': leaf_a7, 'B': leaf_a8})
            node_a5 = EFG_Node(player_a2, {'C': node_a1, 'D': node_a2})
            node_a6 = EFG_Node(player_a2, {'C': node_a3, 'D': node_a4})
            root_a = EFG_Node(player_a1, {'A': node_a5, 'B': node_a6})
            egame_a = ExtensiveFormGame(root_a)
            p = egame_a.plot()
            sphinx_plot(p)

        It is possible to have games with more than two players::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: player_b3 = EFG_Player('Player 3')
            sage: leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1,
            ....:                     player_b3: -5}, 'Leaf 1')
            sage: leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0,
            ....:                     player_b3: -4}, 'Leaf 2')
            sage: leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4,
            ....:                     player_b3: -3}, 'Leaf 3')
            sage: leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1,
            ....:                     player_b3: -2}, 'Leaf 4')
            sage: node_b1 = EFG_Node(player_b3, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, 'Node 1')
            sage: node_b2 = EFG_Node(player_b2, {'A': leaf_b3,
            ....:                                'B': leaf_b4}, 'Node 2')
            sage: root_b = EFG_Node(player_b1, {'C': node_b1,
            ....:                               'D': node_b2}, 'Root 1')
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: egame_b.players
            [EFG Player "Player 1",
             EFG Player "Player 2",
             EFG Player "Player 3"]
            sage: egame_b.nodes
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Root 1"]
            sage: egame_b.leafs
            [EFG Leaf with utilities - (0, 1, -5),
             EFG Leaf with utilities - (1, 0, -4),
             EFG Leaf with utilities - (2, 4, -3),
             EFG Leaf with utilities - (2, 1, -2)]
            sage: egame_b.tree
            Graph on 7 vertices

        We can plot the above tree through the ``plot`` method::

            sage: egame_b.plot()
            Graphics object consisting of 20 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_b1 = EFG_Player('Player 1')
            player_b2 = EFG_Player('Player 2')
            player_b3 = EFG_Player('Player 3')
            leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1,
                                player_b3: -5}, 'Leaf 1')
            leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0,
                                player_b3: -4}, 'Leaf 2')
            leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4,
                                player_b3: -3}, 'Leaf 3')
            leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1,
                                player_b3: -2}, 'Leaf 4')
            node_b1 = EFG_Node(player_b3, {'A': leaf_b1,
                                           'B': leaf_b2}, 'Node 1')
            node_b2 = EFG_Node(player_b2, {'A': leaf_b3,
                                           'B': leaf_b4}, 'Node 2')
            root_b = EFG_Node(player_b1, {'C': node_b1,
                                          'D': node_b2}, 'Root 1')
            egame_b = ExtensiveFormGame(root_b)
            p = egame_b.plot()
            sphinx_plot(p)

        If we do not name our nodes, unique names are automatically set during
        the initialisation::

            sage: player_c1 = EFG_Player('Player 1')
            sage: player_c2 = EFG_Player('Player 2')
            sage: leaf_c1 = EFG_Leaf({player_c1 : 0, player_c2: 1})
            sage: leaf_c2 = EFG_Leaf({player_c1 : 1, player_c2: 0})
            sage: leaf_c3 = EFG_Leaf({player_c1 : 2, player_c2: 4})
            sage: leaf_c4 = EFG_Leaf({player_c1 : 2, player_c2: 1})
            sage: node_c1 = EFG_Node(player_c2, {'A': leaf_c1, 'B': leaf_c2})
            sage: node_c2 = EFG_Node(player_c2, {'A': leaf_c3, 'B': leaf_c4})
            sage: root_c = EFG_Node(player_c1, {'C': node_c1, 'D': node_c2})
            sage: node_c1.name is node_c2.name is root_c.name
            True
            sage: egame_c = ExtensiveFormGame(root_c)
            sage: node_c1.name is node_c2.name is root_c.name
            False
            sage: sorted([node_c1.name, node_c2.name])
            ['Node 1', 'Node 2']
            sage: root_c.name
            'Tree Root'

        If the generator is an ``EFG_Node``, it needs children, actions and
        players::

            sage: false_root_1 = EFG_Node(args = {'C': node_a1,
            ....:                                 'D': node_a2})
            sage: false_egame_1 = ExtensiveFormGame(false_root_1)
            Traceback (most recent call last):
            ...
            AttributeError: Root node has no player.

            sage: false_root_2 = EFG_Node(player_a1)
            sage: false_egame_2 = ExtensiveFormGame(false_root_2)
            Traceback (most recent call last):
            ...
            AttributeError: Root node has no children.

        Similarly, we cannot create a tree with a Node with a player attribute
        which isn't an instance of the ``Player`` class::

            sage: player_d1 = EFG_Player('Player 1')
            sage: player_d2 = EFG_Player('Player 2')
            sage: leaf_d1 = EFG_Leaf({player_d1 : 0, player_d2: 1}, 'Leaf 1')
            sage: leaf_d2 = EFG_Leaf({player_d1 : 1, player_d2: 0}, 'Leaf 2')
            sage: leaf_d3 = EFG_Leaf({player_d1 : 2, player_d2: 4}, 'Leaf 3')
            sage: leaf_d4 = EFG_Leaf({player_d1 : 2, player_d2: 1}, 'Leaf 4')
            sage: node_d1 = EFG_Node(player_d2, {'A': leaf_d1,
            ....:                                'B': leaf_d2}, 'Node 1')
            sage: node_d2 = EFG_Node(leaf_d2, {'A': leaf_d3,
            ....:                              'B': leaf_d4}, 'Node 2')
            sage: root_d = EFG_Node(player_d1, {'C': node_d1,
            ....:                               'D': node_d2}, 'Root 1')
            sage: egame_d = ExtensiveFormGame(root_d)
            Traceback (most recent call last):
            ...
            TypeError: Cannot assign a non Player as a player to a Node.

        We can also create an ExtensiveFormGame with a gambit game:

        We first need to set up the gambit game as normal (note: we can also
        import a gambit game from a gambit .efg format)::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_2.actions[int(0)].label = "A"  # optional - gambit
            sage: move_2.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2.actions[int(2)].label = "C"  # optional - gambit
            sage: g.root.children[int(0)].label = "Node 1"  # optional - gambit
            sage: g.root.children[int(1)].label = "Node 2"  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "A"  # optional - gambit
            sage: move_3.actions[int(1)].label = "B"  # optional - gambit
            sage: outcome_1 = g.outcomes.add()  # optional - gambit
            sage: outcome_1[int(0)] = int(1)  # optional - gambit
            sage: outcome_1[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome_1  # optional - gambit
            sage: outcome_2 = g.outcomes.add()  # optional - gambit
            sage: outcome_2[int(0)] = int(5)  # optional - gambit
            sage: outcome_2[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome_2  # optional - gambit
            sage: outcome_3 = g.outcomes.add()  # optional - gambit
            sage: outcome_3[int(0)] = int(9)  # optional - gambit
            sage: outcome_3[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome_3  # optional - gambit
            sage: outcome_4 = g.outcomes.add()  # optional - gambit
            sage: outcome_4[int(0)] = int(3)  # optional - gambit
            sage: outcome_4[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome_4  # optional - gambit
            sage: outcome_5 = g.outcomes.add()  # optional - gambit
            sage: outcome_5[int(0)] = int(2)  # optional - gambit
            sage: outcome_5[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome_5  # optional - gambit

        Then we generate the ``ExtensiveFormGame`` with the Gambit ``Game``::

            sage: sage_game = ExtensiveFormGame(g)  # optional - gambit

        Then we can look at all the outputs and attributes of the
        ``ExtensiveFormGame``::

            sage: sage_game.nodes  # optional - gambit
            [EFG Node "Node 1",
            EFG Node "Node 2",
            EFG Node "Tree Root"]
            sage: sage_game.leafs  # optional - gambit
            [EFG Leaf with utilities - (Fraction(1, 1), Fraction(5, 1)),
            EFG Leaf with utilities - (Fraction(2, 1), Fraction(7, 1)),
            EFG Leaf with utilities - (Fraction(3, 1), Fraction(0, 1)),
            EFG Leaf with utilities - (Fraction(5, 1), Fraction(2, 1)),
            EFG Leaf with utilities - (Fraction(9, 1), Fraction(1, 1))]
            sage: sage_game.tree_root  # optional - gambit
            EFG Node "Tree Root"
            sage: sage_game.tree  # optional - gambit
            Graph on 8 vertices
            sage: sage_game.players  # optional - gambit
            [EFG Player "1", EFG Player "2"]
            sage: sage_game.info_sets  # optional - gambit
            [[EFG Node "Node 1"],
            [EFG Node "Node 2"],
            [EFG Node "Tree Root"]]

        We can also see the plot of this game::

            sage: sage_game.plot()  # optional - gambit
            Graphics object consisting of 23 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            import gambit  # optional - gambit
            g = gambit.Game.new_tree()  # optional - gambit
            g.players.add("1")  # optional - gambit
            g.players.add("2")  # optional - gambit
            move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            move_1.actions[int(0)].label = "A"  # optional - gambit
            move_1.actions[int(1)].label = "B"  # optional - gambit
            move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            move_2.actions[int(0)].label = "A"  # optional - gambit
            move_2.actions[int(1)].label = "B"  # optional - gambit
            move_2.actions[int(2)].label = "C"  # optional - gambit
            g.root.children[int(0)].label = "Node 1"  # optional - gambit
            g.root.children[int(1)].label = "Node 2"  # optional - gambit
            move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            move_3.actions[int(0)].label = "A"  # optional - gambit
            move_3.actions[int(1)].label = "B"  # optional - gambit
            outcome_1 = g.outcomes.add()  # optional - gambit
            outcome_1[int(0)] = int(1)  # optional - gambit
            outcome_1[int(1)] = int(5)  # optional - gambit
            g.root.children[int(0)].children[int(0)].outcome = outcome_1  # optional - gambit
            outcome_2 = g.outcomes.add()  # optional - gambit
            outcome_2[int(0)] = int(5)  # optional - gambit
            outcome_2[int(1)] = int(2)  # optional - gambit
            g.root.children[int(0)].children[int(1)].outcome = outcome_2  # optional - gambit
            outcome_3 = g.outcomes.add()  # optional - gambit
            outcome_3[int(0)] = int(9)  # optional - gambit
            outcome_3[int(1)] = int(1)  # optional - gambit
            g.root.children[int(0)].children[int(2)].outcome = outcome_3  # optional - gambit
            outcome_4 = g.outcomes.add()  # optional - gambit
            outcome_4[int(0)] = int(3)  # optional - gambit
            outcome_4[int(1)] = int(0)  # optional - gambit
            g.root.children[int(1)].children[int(0)].outcome = outcome_4  # optional - gambit
            outcome_5 = g.outcomes.add()  # optional - gambit
            outcome_5[int(0)] = int(2)  # optional - gambit
            outcome_5[int(1)] = int(7)  # optional - gambit
            g.root.children[int(1)].children[int(1)].outcome = outcome_5  # optional - gambit
            sage_game = ExtensiveFormGame(g)
            p = sage_game.plot()
            sphinx_plot(p)

        If we try to generate a game with something that isn't an ``EFG_Node``
        or a Gambit game, we'll also return an error::

            sage: false_egame_4 = ExtensiveFormGame(player_a1)
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an input in the form
            of an EFG_Node or a Gambit Game object.

            sage: false_egame_5 = ExtensiveFormGame([node_a1, node_a2])
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an input in the form
            of an EFG_Node or a Gambit Game object.

            sage: false_egame_6 = ExtensiveFormGame(leaf_a1)
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an input in the form
            of an EFG_Node or a Gambit Game object.

        """

        self.nodes = []
        self._sage_to_gambit_info_set_list = []
        self._gambit_to_sage_node_dict = {}
        self._gambit_to_sage_player_dict = {}

        if isinstance(generator, EFG_Node):
            if generator.actions is False:
                raise AttributeError("Root node has no actions.")
            elif generator.children == []:
                raise AttributeError("Root node has no children.")
            elif generator.player is False:
                raise AttributeError("Root node has no player.")
            else:
                self.tree_root = generator
                self.tree, self.tree_dictionary, self.nodes = self._grow_tree()
                self.players = []
                self.info_sets = [[node] for node in self.nodes]
                self.leafs = []

                self._check_node_names_and_find_players(generator)

                self.players.sort(key=attrgetter('name'))
                self.nodes.sort(key=attrgetter('name', 'parent'))
                self.leafs.sort(key=attrgetter('name', 'payoffs'))
                self.info_sets.sort(key=lambda x: x[0].name)

        elif type(generator) is Game:
            self._gambit_to_sage(generator)
        else:
            raise TypeError("Extensive form game must be passed an input in "
                            "the form of an EFG_Node or a Gambit Game object.")


    def _gambit_to_sage(self, gambit_game):
        r"""
        A private method which takes a gambit extensive form game and converts
        it to a sage extensive form game.

        To use this, we first set up a Gambit Game tree::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.label = "a"  # optional - gambit
            sage: move_1.actions[int(0)].label = "X"  # optional - gambit
            sage: move_1.actions[int(1)].label = "W"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_2.label = "b"  # optional - gambit
            sage: move_2.actions[int(0)].label = "B"  # optional - gambit
            sage: move_2.actions[int(1)].label = "A"  # optional - gambit
            sage: g.root.children[int(1)].append_move(move_2)  # optional - gambit
            <Infoset [0] 'b' for player '2' in game ''>
            sage: move_3 = g.root.children[int(0)].children[int(1)].append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_3.label = "d"  # optional - gambit
            sage: move_3.actions[int(0)].label = "Z"  # optional - gambit
            sage: move_3.actions[int(1)].label = "Y"  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].label = 'Node 3'  # optional - gambit
            sage: outcome_1 = g.outcomes.add("1")  # optional - gambit
            sage: outcome_1[int(0)] = int(2)  # optional - gambit
            sage: outcome_1[int(1)] = int(0)  # optional - gambit
            sage: outcome_2 = g.outcomes.add("outcome_2")  # optional - gambit
            sage: outcome_2[int(0)] = int(3)  # optional - gambit
            sage: outcome_2[int(1)] = int(1)  # optional - gambit
            sage: outcome_3 = g.outcomes.add("outcome_3")  # optional - gambit
            sage: outcome_3[int(0)] = int(4)  # optional - gambit
            sage: outcome_3[int(1)] = int(2)  # optional - gambit
            sage: outcome_4 = g.outcomes.add("outcome_4")  # optional - gambit
            sage: outcome_4[int(0)] = int(3)  # optional - gambit
            sage: outcome_4[int(1)] = int(5)  # optional - gambit
            sage: outcome_5 = g.outcomes.add("outcome_5")  # optional - gambit
            sage: outcome_5[int(0)] = int(4)  # optional - gambit
            sage: outcome_5[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome_1 # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(0)].outcome = outcome_2  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(1)].outcome = outcome_3  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome_4  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome_5  # optional - gambit

        Then we generate the ``ExtensiveFormGame`` with the Gambit ``Game``::

            sage: sage_game = ExtensiveFormGame(g)  # optional - gambit

        Then we can look at all the outputs and attributes::

            sage: sage_game.nodes  # optional - gambit
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Node 3",
             EFG Node "Root"]
            sage: sage_game.leafs  # optional - gambit
            [EFG Leaf with utilities - (Fraction(2, 1), Fraction(0, 1)),
             EFG Leaf with utilities - (Fraction(3, 1), Fraction(1, 1)),
             EFG Leaf with utilities - (Fraction(3, 1), Fraction(5, 1)),
             EFG Leaf with utilities - (Fraction(4, 1), Fraction(1, 1)),
             EFG Leaf with utilities - (Fraction(4, 1), Fraction(2, 1))]
            sage: sage_game.tree_root  # optional - gambit
            EFG Node "Root"
            sage: sage_game.tree  # optional - gambit
            Graph on 9 vertices
            sage: sage_game.players  # optional - gambit
            [EFG Player "1", EFG Player "2"]
            sage: sage_game.info_sets  # optional - gambit
            [[EFG Node "Node 1",
              EFG Node "Node 2"],
             [EFG Node "Node 3"],
             [EFG Node "Root"]]

        Another test with a different Gambit tree (with only one labeled node)::

            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(3))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_1.actions[int(2)].label = "C"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_2.actions[int(0)].label = "D"  # optional - gambit
            sage: move_2.actions[int(1)].label = "E"  # optional - gambit
            sage: move_3 = g.root.children[int(0)].children[int(1)].append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "F"  # optional - gambit
            sage: move_3.actions[int(1)].label = "G"  # optional - gambit
            sage: move_4 = g.root.children[int(2)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_4.actions[int(0)].label = "H"  # optional - gambit
            sage: move_4.actions[int(1)].label = "I"  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].label = 'Node 3' # optional - gambit
            sage: outcome_1 = g.outcomes.add()  # optional - gambit
            sage: outcome_1[int(0)] = int(1)  # optional - gambit
            sage: outcome_1[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome_1  # optional - gambit
            sage: outcome_2 = g.outcomes.add()  # optional - gambit
            sage: outcome_2[int(0)] = int(5)  # optional - gambit
            sage: outcome_2[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(0)].outcome = outcome_2  # optional - gambit
            sage: outcome_3 = g.outcomes.add()  # optional - gambit
            sage: outcome_3[int(0)] = int(9)  # optional - gambit
            sage: outcome_3[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(1)].outcome = outcome_3  # optional - gambit
            sage: outcome_4 = g.outcomes.add()  # optional - gambit
            sage: outcome_4[int(0)] = int(3)  # optional - gambit
            sage: outcome_4[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].outcome = outcome_4  # optional - gambit
            sage: outcome_5 = g.outcomes.add()  # optional - gambit
            sage: outcome_5[int(0)] = int(2)  # optional - gambit
            sage: outcome_5[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(2)].children[int(0)].outcome = outcome_5  # optional - gambit
            sage: outcome_6 = g.outcomes.add()  # optional - gambit
            sage: outcome_6[int(0)] = int(1)  # optional - gambit
            sage: outcome_6[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(2)].children[int(1)].outcome = outcome_6  # optional - gambit

        Then we generate the ``ExtensiveFormGame`` with the Gambit ``Game``::

            sage: sage_game = ExtensiveFormGame(g)  # optional - gambit

        Then we can look at all the outputs and attributes::

            sage: sage_game.nodes  # optional - gambit
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Node 3",
             EFG Node "Tree Root"]
            sage: sage_game.leafs  # optional - gambit
            [EFG Leaf with utilities - (Fraction(1, 1), Fraction(5, 1)),
             EFG Leaf with utilities - (Fraction(1, 1), Fraction(5, 1)),
             EFG Leaf with utilities - (Fraction(2, 1), Fraction(7, 1)),
             EFG Leaf with utilities - (Fraction(3, 1), Fraction(0, 1)),
             EFG Leaf with utilities - (Fraction(5, 1), Fraction(2, 1)),
             EFG Leaf with utilities - (Fraction(9, 1), Fraction(1, 1))]
            sage: sage_game.tree_root  # optional - gambit
            EFG Node "Tree Root"
            sage: sage_game.tree  # optional - gambit
            Graph on 10 vertices
            sage: sage_game.players  # optional - gambit
            [EFG Player "1", EFG Player "2"]
            sage: sage_game.info_sets  # optional - gambit
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Node 3"],
             [EFG Node "Tree Root"]]

        Another test with a different tree::

            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: move_2.actions[int(2)].label = "E"  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "F"  # optional - gambit
            sage: move_3.actions[int(1)].label = "G"  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: outcome_1 = g.outcomes.add()  # optional - gambit
            sage: outcome_1[int(0)] = int(1)  # optional - gambit
            sage: outcome_1[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome_1  # optional - gambit
            sage: outcome_2 = g.outcomes.add()  # optional - gambit
            sage: outcome_2[int(0)] = int(5)  # optional - gambit
            sage: outcome_2[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome_2  # optional - gambit
            sage: outcome_3 = g.outcomes.add()  # optional - gambit
            sage: outcome_3[int(0)] = int(9)  # optional - gambit
            sage: outcome_3[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome_3  # optional - gambit
            sage: outcome_4 = g.outcomes.add()  # optional - gambit
            sage: outcome_4[int(0)] = int(3)  # optional - gambit
            sage: outcome_4[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome_4  # optional - gambit
            sage: outcome_5 = g.outcomes.add()  # optional - gambit
            sage: outcome_5[int(0)] = int(2)  # optional - gambit
            sage: outcome_5[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome_5  # optional - gambit

        Then we generate the ``ExtensiveFormGame`` with the Gambit ``Game``::

            sage: sage_game = ExtensiveFormGame(g)  # optional - gambit

        Then we can look at all the outputs and attributes::

            sage: sage_game.nodes  # optional - gambit
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Root"]
            sage: sage_game.leafs  # optional - gambit
            [EFG Leaf with utilities - (Fraction(1, 1), Fraction(5, 1)),
             EFG Leaf with utilities - (Fraction(2, 1), Fraction(7, 1)),
             EFG Leaf with utilities - (Fraction(3, 1), Fraction(0, 1)),
             EFG Leaf with utilities - (Fraction(5, 1), Fraction(2, 1)),
             EFG Leaf with utilities - (Fraction(9, 1), Fraction(1, 1))]
            sage: sage_game.tree_root  # optional - gambit
            EFG Node "Root"
            sage: sage_game.tree  # optional - gambit
            Graph on 8 vertices
            sage: sage_game.players  # optional - gambit
            [EFG Player "1", EFG Player "2"]
            sage: sage_game.info_sets  # optional - gambit
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Root"]]

        Another test::

            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: g.root.children[int(1)].append_move(move_2)  # optional - gambit
            <Infoset [0] '' for player '2' in game ''>
            sage: outcome_1 = g.outcomes.add()  # optional - gambit
            sage: outcome_1[int(0)] = int(1)  # optional - gambit
            sage: outcome_1[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome_1  # optional - gambit
            sage: outcome_2 = g.outcomes.add()  # optional - gambit
            sage: outcome_2[int(0)] = int(5)  # optional - gambit
            sage: outcome_2[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome_2  # optional - gambit
            sage: outcome_3 = g.outcomes.add()  # optional - gambit
            sage: outcome_3[int(0)] = int(3)  # optional - gambit
            sage: outcome_3[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome_3  # optional - gambit
            sage: outcome_4 = g.outcomes.add()  # optional - gambit
            sage: outcome_4[int(0)] = int(2)  # optional - gambit
            sage: outcome_4[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome_4  # optional - gambit

        Then we can use it to set up an ``ExtensiveFormGame``::

            sage: sage_game = ExtensiveFormGame(g)  # optional - gambit

        Then we can look at all the outputs and attributes::

            sage: sage_game.nodes  # optional - gambit
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Tree Root"]
            sage: sage_game.leafs  # optional - gambit
            [EFG Leaf with utilities - (Fraction(1, 1), Fraction(5, 1)),
             EFG Leaf with utilities - (Fraction(2, 1), Fraction(7, 1)),
             EFG Leaf with utilities - (Fraction(3, 1), Fraction(0, 1)),
             EFG Leaf with utilities - (Fraction(5, 1), Fraction(2, 1))]
            sage: sage_game.tree_root  # optional - gambit
            EFG Node "Tree Root"
            sage: sage_game.tree  # optional - gambit
            Graph on 7 vertices
            sage: sage_game.players  # optional - gambit
            [EFG Player "1", EFG Player "2"]
            sage: sage_game.info_sets  # optional - gambit
            [[EFG Node "Node 1",
              EFG Node "Node 2"],
             [EFG Node "Tree Root"]]

        """
        for player in gambit_game.players:
            self._gambit_to_sage_player_dict[player] = EFG_Player(player.label)

        self._gambit_to_sage_organize_nodes_for_looping(gambit_game)

        for i in reversed(range(len(self._gambit_to_sage_node_indexed_dict))):
                for gambit_node in self._gambit_to_sage_node_indexed_dict[i]:
                    self._gambit_to_sage_rename_blank_actions(gambit_node)
                    self._gambit_to_sage_create_node_or_leaf(gambit_node,
                                                             gambit_game)

        self._gambit_to_sage_check_node_names_for_renaming()
        converted_root = self._gambit_to_sage_node_dict[gambit_game.root]
        converted_gambit_game = ExtensiveFormGame(converted_root)
        self._gambit_to_sage_setup_info_sets(gambit_game,
                                             converted_gambit_game)
        self._gambit_to_sage_setup_attributes(converted_root,
                                              converted_gambit_game)

    def _gambit_to_sage_organize_nodes_for_looping(self, gambit_game):
        """
        A method for ``_gambit_to_sage`` which goes through nodes of
        a gambit game, and stores them in a dictionary with corresponding keys
        depending on what 'generation' they are. i.e. the tree root is at 0,
        its children are at 1.

        If we wish to test the method, we need a gambit game set up (this one
        has some unlabeled actions)::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = 'A' # optional - gambit
            sage: move_3.actions[int(1)].label = 'B' # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(5)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(9)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(3)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        We can then see the dictionary::

            sage: sage_game = ExtensiveFormGame(g)
            sage: sage_game._gambit_to_sage_organize_nodes_for_looping(g)
            sage: expected_dict = {0: [g.root],
            ....:                  1: [g.root.children[int(0)],
            ....:                      g.root.children[int(1)]],
            ....:                  2: [g.root.children[int(0)].children[int(0)],
            ....:                      g.root.children[int(0)].children[int(1)],
            ....:                      g.root.children[int(0)].children[int(2)],
            ....:                      g.root.children[int(1)].children[int(0)],
            ....:                      g.root.children[int(1)].children[int(1)]]}
            sage: sage_game._gambit_to_sage_node_indexed_dict == expected_dict
            True
        """
        index = 1
        self._gambit_to_sage_node_indexed_dict = {0: [gambit_game.root]}
        node_list = [gambit_game.root]
        while len(node_list) != 0:
            new_node_list = []
            for node in node_list:
                for child in node.children:
                    new_node_list.append(child)
            if len(new_node_list) != 0:
                self._gambit_to_sage_node_indexed_dict[index] = new_node_list
                index += 1
            node_list = new_node_list

    def _gambit_to_sage_rename_blank_actions(self, gambit_node):
        """
        A method for ``_gambit_to_sage`` which goes through the actions of a
        Gambit Node, and if they are blank they will be renamed via a default
        naming system

        If we wish to test the method, we need a gambit game set up (this one
        has some unlabeled actions)::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = 'A' # optional - gambit
            sage: move_3.actions[int(1)].label = 'B' # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(5)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(9)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(3)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        If we pass the game to an ``ExtensiveFormGame``, this method will be
        used automatically, (i.e. all unnamed actions will be named by running
        this)::

            sage: unnamed_action_game = ExtensiveFormGame(g)  # optional - gambit

        We can then use check the actions of each node::

            sage: unnamed_action_game._gambit_to_sage_node_dict[g.root].actions  # optional - gambit
            ['Action 2', 'Action 1']
            sage: unnamed_action_game._gambit_to_sage_node_dict[g.root.children[int(0)]].actions  # optional - gambit
            ['Action 2', 'Action 3', 'Action 1']

        If we have already labeled the actions before converting the game,
        the actions do not get renamed::

            sage: unnamed_action_game._gambit_to_sage_node_dict[g.root.children[int(1)]].actions  # optional - gambit
            ['A', 'B']

        """
        act_string_index = 1
        gambit_act_index = 0
        if gambit_node.infoset:
            for action in gambit_node.infoset.actions:
                if action.label is '':
                    gambit_node.infoset.actions[gambit_act_index].label = "Action %s" %act_string_index
                    act_string_index += 1
                gambit_act_index += 1

    def _gambit_to_sage_create_node_or_leaf(self, gambit_node, gambit_game):
        r"""
        A method for ``_gambit_to_sage`` which takes a gambit_node and either
        creates an ExtensiveFormGame Leaf or Node, it then stores the Sage Leaf
        or Node in a dictionary where the key is the original Gambit Node. It
        does not allow users to use gambit nodes where they are terminal but
        have no assigned outcome.

        If we wish to test the method, we need a gambit game set up::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "F"  # optional - gambit
            sage: move_3.actions[int(1)].label = "G"  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(5)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(9)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(3)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        If we pass the game to an ``ExtensiveFormGame``, the method will be
        used automatically::

            sage: egame = ExtensiveFormGame(g)  # optional - gambit

        We can then check the corresponding ``EFG_Leaf`` or ``EFG_Node`` for
        each Gambit Node::

            sage: egame._gambit_to_sage_node_dict[g.root.children[int(1)].children[int(1)]]  # optional - gambit
            EFG Leaf with utilities - (Fraction(2, 1), Fraction(7, 1))
            sage: egame._gambit_to_sage_node_dict[g.root.children[int(1)]]  # optional - gambit
            EFG Node "Node 2"

        The code will also stop users from inputting gambit games with terminal
        nodes without outcomes::

            sage: bad_game = gambit.Game.new_tree()  # optional - gambit
            sage: bad_game.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: bad_game.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move = bad_game.root.append_move(bad_game.players["1"], int(2))  # optional - gambit
            sage: ExtensiveFormGame(bad_game)  # optional - gambit
            Traceback (most recent call last):
            ...
            AttributeError: One or more terminal nodes in the Gambit Game have
            no outcome.
        """

        if gambit_node.is_terminal:
            if not gambit_node.outcome :
                raise AttributeError("One or more terminal nodes in the Gambit Game have no outcome.")
            leaf_dictionary = {self._gambit_to_sage_player_dict[gambit_game.players[outcome_index]]:gambit_node.outcome[outcome_index]
                               for outcome_index in range(len(list(gambit_node.outcome)))}
            leaf = EFG_Leaf(leaf_dictionary, gambit_node.label)
            self._gambit_to_sage_node_dict[gambit_node] = leaf
        else:
            node_dictionary = {gambit_node.infoset.actions[gambit_index].label:self._gambit_to_sage_node_dict[gambit_node.children[gambit_index]]
                               for gambit_index in range(len(list(gambit_node.children)))}
            node = EFG_Node(args = node_dictionary,
                            player = self._gambit_to_sage_player_dict[gambit_node.infoset.player],
                            name = gambit_node.label)
            self._gambit_to_sage_node_dict[gambit_node] = node

    def _gambit_to_sage_setup_info_sets(self, gambit_game, converted_gambit_game):
        r"""
        A method for ``_gambit_to_sage`` which sets up the information sets of
        the ``ExtensiveFormGame`` created with a Gambit Game. If we wish to
        test the method, we need a gambit game set up::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "F"  # optional - gambit
            sage: move_3.actions[int(1)].label = "G"  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node A'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node B'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(6)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(4)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        If we pass the game to an ``ExtensiveFormGame``, the method will be
        used automatically::

            sage: egame = ExtensiveFormGame(g)  # optional - gambit

        We can then see the information sets::

            sage: egame.info_sets  # optional - gambit
            [[EFG Node "Node A"],
            [EFG Node "Node B"],
            [EFG Node "Root"]]

        If we then recreate the game by only changing it so that Node A and
        Node B are now the same information set, this will show::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: g.root.children[int(1)].append_move(move_2)  # optional - gambit
            <Infoset [0] '' for player '2' in game ''>
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node A'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node B'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(6)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(4)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit
            sage: egame = ExtensiveFormGame(g)  # optional - gambit
            sage: egame.info_sets  # optional - gambit
            [[EFG Node "Node A",
              EFG Node "Node B"],
             [EFG Node "Root"]]
        """
        for info_set in gambit_game.infosets:
            info_set_list = [self._gambit_to_sage_node_dict[gambit_node]
                             for gambit_node in info_set.members]
            converted_gambit_game.set_info_set(info_set_list)

    def _gambit_to_sage_setup_attributes(self, converted_root, converted_gambit_game):
        r"""
        A method for ``_gambit_to_sage`` which sets up the attributes of the
        ``ExtensiveFormGame`` created with a Gambit Game. If we wish to test
        the method, we need a gambit game set up::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: move_2.actions[int(2)].label = "E"  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "F"  # optional - gambit
            sage: move_3.actions[int(1)].label = "G"  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node A'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node B'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(6)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(4)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        If we pass the game to an ``ExtensiveFormGame``, the method will be
        used automatically::

            sage: egame = ExtensiveFormGame(g)  # optional - gambit

        We can then see the attributes::

            sage: egame.tree_root  # optional - gambit
            EFG Node "Root"
            sage: egame.tree  # optional - gambit
            Graph on 8 vertices
            sage: egame.nodes  # optional - gambit
            [EFG Node "Node A",
             EFG Node "Node B",
             EFG Node "Root"]
            sage: egame.players  # optional - gambit
            [EFG Player "1", EFG Player "2"]
            sage: egame.info_sets  # optional - gambit
            [[EFG Node "Node A"],
            [EFG Node "Node B"],
            [EFG Node "Root"]]
            sage: egame.leafs  # optional - gambit
            [EFG Leaf with utilities - (Fraction(1, 1), Fraction(6, 1)),
             EFG Leaf with utilities - (Fraction(2, 1), Fraction(1, 1)),
             EFG Leaf with utilities - (Fraction(4, 1), Fraction(0, 1)),
             EFG Leaf with utilities - (Fraction(6, 1), Fraction(2, 1)),
             EFG Leaf with utilities - (Fraction(6, 1), Fraction(2, 1))]
        """
        self.tree_root = converted_root
        self.tree, self.tree_dictionary, self.nodes = converted_gambit_game._grow_tree()
        self.players = converted_gambit_game.players
        self.info_sets = converted_gambit_game.info_sets
        self.leafs = converted_gambit_game.leafs

        self.players.sort(key=attrgetter('name'))
        self.nodes.sort(key=attrgetter('name', 'parent'))
        self.leafs.sort(key=attrgetter('utilities'))
        self.info_sets.sort(key=lambda x: x[0].name)

    def _gambit_to_sage_check_node_names_for_renaming(self):
        r"""
        A method for ``_gambit_to_sage`` which takes any unnamed nodes and sets
        them to False so they get renamed while initalising the
        ``ExtensiveFormGame``

        For this we need a Gambit Game setup::

            sage: import gambit  # optional - gambit
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: move_1 = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: move_1.actions[int(0)].label = "A"  # optional - gambit
            sage: move_1.actions[int(1)].label = "B"  # optional - gambit
            sage: move_2 = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: move_2.actions[int(0)].label = "C"  # optional - gambit
            sage: move_2.actions[int(1)].label = "D"  # optional - gambit
            sage: move_2.actions[int(2)].label = "E"  # optional - gambit
            sage: move_3 = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: move_3.actions[int(0)].label = "F"  # optional - gambit
            sage: move_3.actions[int(1)].label = "G"  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(6)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(4)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(6)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        Then we can setup the game with unnamed nodes and the game should name
        the nodes with the default method::

            sage: egame_without_names = ExtensiveFormGame(g)  # optional - gambit
            sage: egame_without_names.nodes  # optional - gambit
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Tree Root"]

        Then if we label the nodes and set up the game again, we can see that
        the nodes in the game take the new names::

            sage: g.root.label = 'Root'  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node A'  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node B'  # optional - gambit
            sage: egame_with_names = ExtensiveFormGame(g)  # optional - gambit
            sage: egame_with_names.nodes  # optional - gambit
            [EFG Node "Node A",
             EFG Node "Node B",
             EFG Node "Root"]
        """
        for gambit_node in self._gambit_to_sage_node_dict:
            if gambit_node.label is '':
                self._gambit_to_sage_node_dict[gambit_node].name = None


    def set_info_set(self, node_list):
        r"""
        Combines a list of nodes into an information set. It does this by
        removing the nodes to be set from the information sets they are
        currently in (then making the rest of the nodes in those information
        sets in their own individual information sets for each node) and
        then adding the information set to the information sets within the
        game.

        EXAMPLES::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1, 'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3, 'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1, 'D': node_a2}, 'Root 1')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.info_sets
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Root 1"]]
            sage: egame_a.set_info_set([node_a1, node_a2])
            sage: egame_a.info_sets
            [[EFG Node "Node 1",
              EFG Node "Node 2"],
             [EFG Node "Root 1"]]

        Once we've set an information set, we can see it visually on the
        graph::

            sage: egame_a.plot()
            Graphics object consisting of 21 graphics primitives

        Which will be plotted follows:

        .. PLOT::
            :width: 500 px

            player_a1 = EFG_Player('Player 1')
            player_a2 = EFG_Player('Player 2')
            leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            node_a1 = EFG_Node(player_a2, {'A': leaf_a1, 'B': leaf_a2}, 'Node 1')
            node_a2 = EFG_Node(player_a2, {'A': leaf_a3, 'B': leaf_a4}, 'Node 2')
            root_a = EFG_Node(player_a1, {'C': node_a1, 'D': node_a2}, 'Root 1')
            egame_a = ExtensiveFormGame(root_a)
            egame_a.info_sets
            egame_a.set_info_set([node_a1, node_a2])
            p = egame_a.plot()
            sphinx_plot(p)

        On some occasions we might want to plot the tree without showing the
        information sets::

            sage: egame_a.plot(view_info_set_lines = False)
            Graphics object consisting of 20 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_a1 = EFG_Player('Player 1')
            player_a2 = EFG_Player('Player 2')
            leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            node_a1 = EFG_Node(player_a2, {'A': leaf_a1, 'B': leaf_a2}, 'Node 1')
            node_a2 = EFG_Node(player_a2, {'A': leaf_a3, 'B': leaf_a4}, 'Node 2')
            root_a = EFG_Node(player_a1, {'C': node_a1, 'D': node_a2}, 'Root 1')
            egame_a = ExtensiveFormGame(root_a)
            egame_a.info_sets
            egame_a.set_info_set([node_a1, node_a2])
            p = egame_a.plot(view_info_set_lines = False)
            sphinx_plot(p)

        If by setting an information set we leave a node without an information
        set it is automatically put in a set by itself::

            sage: egame_a.set_info_set([node_a1])
            sage: egame_a.info_sets
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Root 1"]]

        If two nodes in an information set don't have the same actions, an
        error is returned::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1: 0, player_b2: 1})
            sage: leaf_b2 = EFG_Leaf({player_b1: 1, player_b2: 0})
            sage: leaf_b3 = EFG_Leaf({player_b1: 2, player_b2: 4})
            sage: leaf_b4 = EFG_Leaf({player_b1: 2, player_b2: 1})
            sage: node_b1 = EFG_Node(player_b2, {'A': leaf_b1, 'B': leaf_b2})
            sage: node_b2 = EFG_Node(player_b2, {'DifferentA': leaf_b3, 'B': leaf_b4})
            sage: root_b = EFG_Node(player_b1, {'C': node_b1, 'D': node_b2})
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: egame_b.set_info_set([node_b1, node_b2])
            Traceback (most recent call last):
            ...
            AttributeError: All nodes in the same information set must have the
            same actions.

        If two nodes in an information set have different players, an error is
        returned::

            sage: player_c1 = EFG_Player('Player 1')
            sage: player_c2 = EFG_Player('Player 2')
            sage: leaf_c1 = EFG_Leaf({player_c1: 0, player_c2: 1})
            sage: leaf_c2 = EFG_Leaf({player_c1: 1, player_c2: 0})
            sage: leaf_c3 = EFG_Leaf({player_c1: 2, player_c2: 4})
            sage: leaf_c4 = EFG_Leaf({player_c1: 2, player_c2: 1})
            sage: node_c1 = EFG_Node(player_c1, {'A': leaf_c1, 'B': leaf_c2})
            sage: node_c2 = EFG_Node(player_c2, {'A': leaf_c3, 'B': leaf_c4})
            sage: root_c = EFG_Node(player_c1, {'C': node_c1, 'D': node_c2})
            sage: egame_c = ExtensiveFormGame(root_c)
            sage: egame_c.set_info_set([node_c1, node_c2])
            Traceback (most recent call last):
            ...
            AttributeError: All nodes in the same information set must have the
            same players.

        If we try to create an information set with a node that isn't in the
        game, an error is returned::

            sage: bad_node = EFG_Node(player_c1, {'A': leaf_c1, 'B': leaf_c2})
            sage: egame_c.set_info_set([node_c1, bad_node])
            Traceback (most recent call last):
            ...
            ValueError: Some of the nodes in this information set are not in
            the game.

        If we try to create an information set where some of the nodes are
        duplicated, the duplicates are removed::

            sage: player_d1 = EFG_Player('Player 1')
            sage: player_d2 = EFG_Player('Player 2')
            sage: leaf_d1 = EFG_Leaf({player_d1: 0, player_d2: 1})
            sage: leaf_d2 = EFG_Leaf({player_d1: 1, player_d2: 0})
            sage: leaf_d3 = EFG_Leaf({player_d1: 2, player_d2: 4})
            sage: leaf_d4 = EFG_Leaf({player_d1: 2, player_d2: 1})
            sage: node_d1 = EFG_Node(player_d2, {'A': leaf_d1,
            ....:                                'B': leaf_d2}, 'Node 1')
            sage: node_d2 = EFG_Node(player_d2, {'A': leaf_d3,
            ....:                                'B': leaf_d4}, 'Node 2')
            sage: root_d = EFG_Node(player_d1, {'C': node_d1,
            ....:                               'D': node_d2}, 'Root')
            sage: egame_d = ExtensiveFormGame(root_d)
            sage: egame_d.set_info_set([node_d1, node_d1, node_d2,
            ....:                       node_d2, node_d2])
            sage: egame_d.info_sets
            [[EFG Node "Node 1", EFG Node "Node 2"], [EFG Node "Root"]]
            sage: egame_d.set_info_set([node_d1, node_d1, node_d1])
            sage: egame_d.info_sets
            [[EFG Node "Node 1"], [EFG Node "Node 2"], [EFG Node "Root"]]
        """

        if len(set([node.player for node in node_list])) != 1:
            raise AttributeError("All nodes in the same information set must "
                                 "have the same players.")
        if len(set([tuple(sorted(node.actions)) for node in node_list])) != 1:
            raise AttributeError("All nodes in the same information set must "
                                 "have the same actions.")

        for node_to_be_set in node_list:
            if node_to_be_set not in self.nodes:
                raise ValueError("Some of the nodes in this information set "
                                 "are not in the game.")
            for info_set in [info_set for info_set in self.info_sets
                             if node_to_be_set in info_set]:
                self.info_sets.remove(info_set)
        self.info_sets.append((sorted(list(set(node_list)), key=lambda x: x.name)))

        for node in self.nodes:
            if not any(node in info_set for info_set in self.info_sets):
                self.info_sets.append([node])

        self.info_sets.sort(key=lambda x: x[0].name)

    def remove_info_set(self, node_list):
        r"""
        Removes an information set and sets all nodes in that information set
        to be in their own information set.

        EXAMPLES::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, 'Node 1')
            sage: node_2 = EFG_Node(player_2, {'A': leaf_3,
            ....:                              'B': leaf_4}, 'Node 2')
            sage: root = EFG_Node(player_1, {'C': node_1,
            ....:                            'D': node_2}, 'Root 1')
            sage: egame = ExtensiveFormGame(root)
            sage: egame.set_info_set([node_1, node_2])
            sage: egame.info_sets
            [[EFG Node "Node 1",
              EFG Node "Node 2"],
             [EFG Node "Root 1"]]
            sage: egame.perfect_info()
            False
            sage: egame.remove_info_set([node_1, node_2])
            sage: egame.info_sets
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Root 1"]]
            sage: egame.perfect_info()
            True

        As this game has perfect information, the Extensive Form tree should
        not show an information sets::

            sage: egame.plot()
            Graphics object consisting of 20 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            leaf_2 = EFG_Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            leaf_3 = EFG_Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            leaf_4 = EFG_Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            egame.set_info_set([node_1, node_2])
            egame.remove_info_set([node_1, node_2])
            p = egame.plot()
            sphinx_plot(p)

        If we try to remove an information set that doesn't exist within the
        game, a specific error is returned::

            sage: egame.remove_info_set([node_1, node_2])
            Traceback (most recent call last):
            ...
            ValueError: Information set to be removed does not exist.

        However if we try to remove a set which has duplicates, it reads it as
        if the duplicates aren't there::

            sage: info_set = [node_1, node_2, node_2, node_2]
            sage: egame.set_info_set(info_set)
            sage: egame.info_sets
            [[EFG Node "Node 1", EFG Node "Node 2"], [EFG Node "Root 1"]]
            sage: egame.remove_info_set(info_set)
            sage: egame.info_sets
            [[EFG Node "Node 1"], [EFG Node "Node 2"], [EFG Node "Root 1"]]
        """
        try:
            duplicates_removed = (sorted(list(set(node_list)),
                                  key=lambda x: x.name))
            self.info_sets.remove(duplicates_removed)
            for node_to_be_readded in duplicates_removed:
                self.info_sets.append([node_to_be_readded])
            self.info_sets.sort(key=lambda x: x[0].name)
        except ValueError:
            raise ValueError("Information set to be removed does not exist.")

    def perfect_info(self):
        r"""
        Returns True or False, depending on whether or not a game has perfect
        information. A game has perfect information if all nodes are contained
        in information sets with no other nodes.

        All games start of with perfect information, adding information sets
        changes this::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, 'Node 1')
            sage: node_2 = EFG_Node(player_2, {'A': leaf_3,
            ....:                              'B': leaf_4}, 'Node 2')
            sage: root = EFG_Node(player_1, {'C': node_1,
            ....:                            'D': node_2}, 'Root 1')
            sage: egame = ExtensiveFormGame(root)
            sage: egame.perfect_info()
            True
            sage: egame.set_info_set([node_1, node_2])
            sage: egame.perfect_info()
            False
            sage: egame.remove_info_set([node_1, node_2])
            sage: egame.perfect_info()
            True
        """
        perfect_info_set = sorted([[node] for node in self.nodes],
                                  key=lambda x: x[0].name)
        return self.info_sets == perfect_info_set

    def plot(self, view_info_set_lines=True, view_info_set_circles=False, circle_size=0.1):
        r"""
        Returns a visual representation of the game.

        INPUT:

        - ``view_info_set_lines`` - If ``True``, shows information sets as
          lines on the graph. Each information set is represented by a
          different color. It is true by default.

        - ``view_info_set_circles`` - If ``True``, shows circles around nodes
          in information sets (unless the information set only contains one
          node). Each information set is represented by a different color. This
          can be useful to view information sets if they overlap. It is
          ``False`` by default.

        - ``circle_size`` - Allows user to change size of the circles around
          nodes. The default value is 0.1, however this may not be suitable for
          all trees.

        Here we create a game and plot it::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1,
            ....:                               'D': node_a2}, 'Root 1')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.set_info_set([node_a1, node_a2])
            sage: egame_a.plot()
            Graphics object consisting of 21 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})
            leaf_2 = EFG_Leaf({player_1: 1, player_2: 0})
            leaf_3 = EFG_Leaf({player_1: 2, player_2: 4})
            leaf_4 = EFG_Leaf({player_1: 2, player_2: 1})
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            egame.set_info_set([node_1, node_2])
            p = egame.plot()
            sphinx_plot(p)

        The plot has the option for whether or not info-sets are visible, the
        default way of viewing information set is by information set lines,
        we can turn this off::

            sage: egame_a.plot(view_info_set_lines = False)
            Graphics object consisting of 20 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})
            leaf_2 = EFG_Leaf({player_1: 1, player_2: 0})
            leaf_3 = EFG_Leaf({player_1: 2, player_2: 4})
            leaf_4 = EFG_Leaf({player_1: 2, player_2: 1})
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            egame.set_info_set([node_1, node_2])
            p = egame.plot(view_info_set_lines = False)
            sphinx_plot(p)

        Another way to view information sets is by information set circles.
        These are circles that are around nodes within information sets.
        This can be useful in situations where the information set lines
        overlap::

            sage: egame_a.plot(view_info_set_circles = True)
            Graphics object consisting of 23 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})
            leaf_2 = EFG_Leaf({player_1: 1, player_2: 0})
            leaf_3 = EFG_Leaf({player_1: 2, player_2: 4})
            leaf_4 = EFG_Leaf({player_1: 2, player_2: 1})
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            egame.set_info_set([node_1, node_2])
            p = egame.plot(view_info_set_circles = True)
            sphinx_plot(p)

        The trees can be plotted with circles but no lines::

            sage: egame_a.plot(view_info_set_lines = False,
            ....: view_info_set_circles = True)
            Graphics object consisting of 22 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})
            leaf_2 = EFG_Leaf({player_1: 1, player_2: 0})
            leaf_3 = EFG_Leaf({player_1: 2, player_2: 4})
            leaf_4 = EFG_Leaf({player_1: 2, player_2: 1})
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            egame.set_info_set([node_1, node_2])
            p = egame.plot(view_info_set_lines = False, view_info_set_circles = True)
            sphinx_plot(p)

        The circles have a default radius of 0.1, which may not be suitable
        for all plots, circle diameter can be changed with ``circle_size``::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            sage: leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            sage: leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            sage: leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            sage: leaf_b5 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            sage: leaf_b6 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            sage: leaf_b7 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            sage: leaf_b8 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            sage: node_b1 = EFG_Node(player_b1, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, 'Node 1')
            sage: node_b2 = EFG_Node(player_b1, {'A': leaf_b3,
            ....:                                'B': leaf_b4}, 'Node 2')
            sage: node_b3 = EFG_Node(player_b1, {'A': leaf_b5,
            ....:                                'B': leaf_b6}, 'Node 3')
            sage: node_b4 = EFG_Node(player_b1, {'A': leaf_b7,
            ....:                                'B': leaf_b8}, 'Node 4')
            sage: node_b5 = EFG_Node(player_b2, {'C': node_b1,
            ....:                                'D': node_b2}, 'Node 5')
            sage: node_b6 = EFG_Node(player_b2, {'C': node_b3,
            ....:                                'D': node_b4}, 'Node 6')
            sage: root_b = EFG_Node(player_b1, {'A': node_b5,
            ....:                               'B': node_b6}, 'Tree Root')
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: egame_b.set_info_set([node_b1, node_b2])
            sage: egame_b.set_info_set([node_b3, node_b4])
            sage: egame_b.set_info_set([node_b5, node_b6])
            sage: egame_b.plot(view_info_set_circles = True, circle_size = 0.2)
            Graphics object consisting of 53 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_b1 = EFG_Player('Player 1')
            player_b2 = EFG_Player('Player 2')
            leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            leaf_b5 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            leaf_b6 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            leaf_b7 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            leaf_b8 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            node_b1 = EFG_Node(player_b1, {'A': leaf_b1,
                                           'B': leaf_b2}, 'Node 1')
            node_b2 = EFG_Node(player_b1, {'A': leaf_b3,
                                           'B': leaf_b4}, 'Node 2')
            node_b3 = EFG_Node(player_b1, {'A': leaf_b5,
                                           'B': leaf_b6}, 'Node 3')
            node_b4 = EFG_Node(player_b1, {'A': leaf_b7,
                                           'B': leaf_b8}, 'Node 4')
            node_b5 = EFG_Node(player_b2, {'C': node_b1,
                                           'D': node_b2}, 'Node 5')
            node_b6 = EFG_Node(player_b2, {'C': node_b3,
                                           'D': node_b4}, 'Node 6')
            root_b = EFG_Node(player_b1, {'A': node_b5,
                                          'B': node_b6}, 'Tree Root')
            egame_b = ExtensiveFormGame(root_b)
            egame_b.set_info_set([node_b1, node_b2])
            egame_b.set_info_set([node_b3, node_b4])
            egame_b.set_info_set([node_b5, node_b6])
            p = egame_b.plot(view_info_set_circles = True, circle_size = 0.3)
            sphinx_plot(p)
        """

        t = copy(self.tree)
        for node in self.nodes:
            actions = node.action_to_node_dict.keys()
            for action in actions:
                t.set_edge_label(node, node.action_to_node_dict[action],
                                 str(node.player.name) + ': ' + action)

        leaf_labels = {leaf: leaf.name + ' - ' + str(leaf.utilities)
                       for leaf in self.leafs}
        # Adding padding so that leaf labels are to the right of the nodes
        leaf_labels = {leaf: (2 * len(leaf_labels[leaf])) * ' ' + leaf_labels[leaf]
                       for leaf in self.leafs}
        t.relabel(leaf_labels)

        node_labels = {node: node.name for node in self.nodes}
        t.relabel(node_labels)

        coloring = [[node.name for node in self.nodes if node.player == player]
                    for player in self.players]
        coloring += [[leaf for leaf in leaf_labels.values()]]


        tree_plot = t.plot(layout='tree', tree_orientation='right',
                           edge_labels=True, tree_root=self.tree_root.name,
                           save_pos=True, axes=False, partition=coloring)

        positions = t.get_pos()

        visible_info_sets = [info_set for info_set in self.info_sets
                             if len(info_set) > 1]
        color_list = rainbow(len(visible_info_sets))
        color_index = 0
        for info_set in visible_info_sets:
            points = sorted([positions[node.name] for node in info_set])
            if view_info_set_lines is True:
                tree_plot += (line2d(points, linestyle="dashed",
                              color=color_list[color_index]))
            if view_info_set_circles is True:
                for point in points:
                    x, y = point
                    info_set_circle = circle((x, y), circle_size,
                                             linestyle="dashed",
                                             color=color_list[color_index])
                    tree_plot += info_set_circle
            color_index += 1
        return tree_plot

    def plot_info_sets(self):
        r"""
        This method returns a DiGraph showing the relationship between
        information sets in the current game. Each vertex on the graph is an
        information set, where vertices are colored according to the players
        for the information sets.

        The following example is a simple 3 node tree with perfect
        information::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1,
            ....:                               'D': node_a2}, 'Root 1')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.plot_info_sets()
            Graphics object consisting of 8 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})
            leaf_2 = EFG_Leaf({player_1: 1, player_2: 0})
            leaf_3 = EFG_Leaf({player_1: 2, player_2: 4})
            leaf_4 = EFG_Leaf({player_1: 2, player_2: 1})
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            p = egame.plot_info_sets()
            sphinx_plot(p)

        The information set graph changes if we change the information sets::

            sage: egame_a.set_info_set([node_a1, node_a2])
            sage: egame_a.plot_info_sets()
            Graphics object consisting of 5 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_1 = EFG_Player('Player 1')
            player_2 = EFG_Player('Player 2')
            leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})
            leaf_2 = EFG_Leaf({player_1: 1, player_2: 0})
            leaf_3 = EFG_Leaf({player_1: 2, player_2: 4})
            leaf_4 = EFG_Leaf({player_1: 2, player_2: 1})
            node_1 = EFG_Node(player_2, {'A': leaf_1, 'B': leaf_2}, 'Node 1')
            node_2 = EFG_Node(player_2, {'A': leaf_3, 'B': leaf_4}, 'Node 2')
            root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            egame = ExtensiveFormGame(root)
            egame.set_info_set([node_1, node_2])
            p = egame.plot_info_sets()
            sphinx_plot(p)

        Here is another example with a bigger tree::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            sage: leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            sage: leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            sage: leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            sage: leaf_b5 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            sage: leaf_b6 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            sage: leaf_b7 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            sage: leaf_b8 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            sage: node_b1 = EFG_Node(player_b1, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, 'Node 1')
            sage: node_b2 = EFG_Node(player_b1, {'A': leaf_b3,
            ....:                                'B': leaf_b4}, 'Node 2')
            sage: node_b3 = EFG_Node(player_b1, {'A': leaf_b5,
            ....:                                'B': leaf_b6}, 'Node 3')
            sage: node_b4 = EFG_Node(player_b1, {'A': leaf_b7,
            ....:                                'B': leaf_b8}, 'Node 4')
            sage: node_b5 = EFG_Node(player_b2, {'C': node_b1,
            ....:                                'D': node_b2}, 'Node 5')
            sage: node_b6 = EFG_Node(player_b2, {'C': node_b3,
            ....:                                'D': node_b4}, 'Node 6')
            sage: root_b = EFG_Node(player_b1, {'A': node_b5,
            ....:                               'B': node_b6}, 'Tree Root')
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: egame_b.set_info_set([node_b5, node_b6])
            sage: egame_b.plot_info_sets()
            Graphics object consisting of 17 graphics primitives

        Which will be plotted as follows:

        .. PLOT::
            :width: 500 px

            player_b1 = EFG_Player('Player 1')
            player_b2 = EFG_Player('Player 2')
            leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            leaf_b5 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            leaf_b6 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            leaf_b7 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            leaf_b8 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            node_b1 = EFG_Node(player_b1, {'A': leaf_b1,
                                           'B': leaf_b2}, 'Node 1')
            node_b2 = EFG_Node(player_b1, {'A': leaf_b3,
                                           'B': leaf_b4}, 'Node 2')
            node_b3 = EFG_Node(player_b1, {'A': leaf_b5,
                                           'B': leaf_b6}, 'Node 3')
            node_b4 = EFG_Node(player_b1, {'A': leaf_b7,
                                           'B': leaf_b8}, 'Node 4')
            node_b5 = EFG_Node(player_b2, {'C': node_b1,
                                           'D': node_b2}, 'Node 5')
            node_b6 = EFG_Node(player_b2, {'C': node_b3,
                                           'D': node_b4}, 'Node 6')
            root_b = EFG_Node(player_b1, {'A': node_b5,
                                          'B': node_b6}, 'Tree Root')
            egame_b = ExtensiveFormGame(root_b)
            egame_b.set_info_set([node_b5, node_b6])
            p = egame_b.plot_info_sets()
            sphinx_plot(p)

        """
        g = copy(self._grow_info_set_graph())

        info_set_dictionary = self._grow_info_set_graph_dictionary()

        for info_set in info_set_dictionary:
            if type(info_set) is list:
                node = info_set[0][0]
            else:
                node = info_set[0]
            actions = node.action_to_node_dict.keys()
            for info_set_2 in info_set_dictionary[info_set]:
                g.set_edge_label(info_set, info_set_2, str(node.player))

        coloring = [[tuple(info_set) for info_set in self.info_sets
                     if info_set[0].player == player]
                    for player in self.players]

        info_set_plot = g.plot(edge_labels=True, axes=False, partition=coloring)

        return info_set_plot

    def _grow_tree(self):
        r"""
        A private method to grow a tree from a given root, returns the tree
        object and a sorted list of all nodes.

        TESTS::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1,
            ....:                               'D': node_a2}, 'Root')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: g, d, nodes = egame_a._grow_tree()
            sage: g
            Graph on 7 vertices
            sage: d == {node_a2: [leaf_a3, leaf_a4],
            ....:       node_a1: [leaf_a1, leaf_a2],
            ....:       root_a: [node_a1, node_a2]}
            True
            sage: nodes
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Root"]

        If the relationship between the nodes does not correspond to a tree, an
        error is returned. As this method is called in the initialisation
        method, the following test is a methodal test and not a true unit
        test::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1: 0, player_b2: 1})
            sage: leaf_b2 = EFG_Leaf({player_b1: 1, player_b2: 0})
            sage: leaf_b3 = EFG_Leaf({player_b1: 2, player_b2: 4})
            sage: leaf_b4 = EFG_Leaf({player_b1: 2, player_b2: 1})
            sage: node_b1 = EFG_Node(player_b2, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, 'Node 1')
            sage: node_b2 = EFG_Node(player_b2, {'A': leaf_b3,
            ....:                                'B': node_b1}, 'Node 2')
            sage: root_b = EFG_Node(player_b1, {'C': node_b1,
            ....:                               'D': node_b2}, 'Root')
            sage: egame_b = ExtensiveFormGame(root_b)
            Traceback (most recent call last):
            ...
            TypeError: Relationship between nodes does not correspond to a tree.
        """
        d, nodes = self._grow_tree_dictionary()
        t = Graph(d)
        if t.is_tree():
            return t, d, nodes
        else:
            raise TypeError("Relationship between nodes does not correspond to a tree.")

    def _grow_tree_dictionary(self):
        r"""
        Returns a dictionary defining the underlying tree as well as a sorted
        list of nodes.

        TESTS::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1,
            ....:                               'D': node_a2}, 'Root')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: d, nodes = egame_a._grow_tree_dictionary()
            sage: nodes
            [EFG Node "Node 1",
             EFG Node "Node 2",
             EFG Node "Root"]
            sage: t = Graph(d)
            sage: t
            Graph on 7 vertices

        If a node if not complete then an error is returned::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1: 0, player_b2: 1})
            sage: leaf_b2 = EFG_Leaf({player_b1: 1, player_b2: 0})
            sage: leaf_b3 = EFG_Leaf({player_b1: 2, player_b2: 4})
            sage: leaf_b4 = EFG_Leaf({player_b1: 2, player_b2: 1})
            sage: node_b1 = EFG_Node(player_b2, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, 'Node 1')
            sage: node_b2 = EFG_Node(args = {'A': leaf_b3,
            ....:                            'B': leaf_b4}, name = 'Node 2')
            sage: root_b = EFG_Node(player_b1, {'C': node_b1,
            ....:                               'D': node_b2}, 'Root 1')
            sage: egame_b = ExtensiveFormGame(root_b)
            Traceback (most recent call last):
            ...
            AttributeError: One or more of the Nodes in the tree are not
            complete.
        """
        to_check = [self.tree_root]
        checked = []
        while to_check:
            checking = to_check.pop()
            for child in checking.children:
                if not isinstance(child, EFG_Leaf):
                    if child._is_complete():
                        child._player_check()
                        to_check.append(child)
                    else:
                        raise AttributeError("One or more of the Nodes in the "
                                              "tree are not complete.")
                checked.append(child)

        d = {node:node.children for node in checked if not isinstance(node, EFG_Leaf)}
        d[self.tree_root] = self.tree_root.children
        return d, sorted(d.keys(), key=attrgetter('name'))

    def _grow_info_set_graph_dictionary(self):
        r"""
        A private method that creates a dictionary needed to plot a graph made
        from the game's information sets. i.e. a dictionary showing which
        information sets are connected to eachother.

        TESTS::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: leaf_a5 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a6 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a7 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a8 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a1, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a1, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: node_a3 = EFG_Node(player_a1, {'A': leaf_a5,
            ....:                                'B': leaf_a6}, 'Node 3')
            sage: node_a4 = EFG_Node(player_a1, {'A': leaf_a7,
            ....:                                'B': leaf_a8}, 'Node 4')
            sage: node_a5 = EFG_Node(player_a2, {'C': node_a1,
            ....:                                'D': node_a2}, 'Node 5')
            sage: node_a6 = EFG_Node(player_a2, {'C': node_a3,
            ....:                                'D': node_a4}, 'Node 6')
            sage: root_a = EFG_Node(player_a1, {'A': node_a5,
            ....:                               'B': node_a6}, 'Tree Root')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: d = egame_a._grow_info_set_graph_dictionary()
            sage: expected_outcome = {(root_a,): [(node_a5,), (node_a6,)],
            ....:                     (node_a5,): [(node_a1,), (node_a2,)],
            ....:                     (node_a6,): [(node_a3,), (node_a4,)]}
            sage: d == expected_outcome
            True

        If we change any information sets, the dictionary changes accordingly::

            sage: egame_a.set_info_set([node_a5, node_a6])
            sage: d = egame_a._grow_info_set_graph_dictionary()
            sage: expected_outcome = {(node_a5, node_a6,): [(node_a1,), (node_a2,),
            ....:                                           (node_a3,), (node_a4,)],
            ....:                     (root_a,): [(node_a5, node_a6,)]}
            sage: d == expected_outcome
            True

            sage: egame_a.set_info_set([node_a3, node_a4])
            sage: egame_a.set_info_set([node_a1, node_a2])
            sage: d = egame_a._grow_info_set_graph_dictionary()
            sage: expected_outcome = {(root_a,): [(node_a5, node_a6,)],
            ....:                     (node_a5, node_a6,): [(node_a1, node_a2,),
            ....:                                           (node_a3, node_a4,)]}
            sage: d == expected_outcome
            True

            sage: egame_a.set_info_set([node_a3, node_a4, node_a1, node_a2])
            sage: d = egame_a._grow_info_set_graph_dictionary()
            sage: expected_outcome =  {(root_a,): [(node_a5, node_a6,)],
            ....:                      (node_a5, node_a6,): [(node_a1, node_a2,
            ....:                                             node_a3, node_a4,)]}
            sage: d == expected_outcome
            True

            sage: egame_a.remove_info_set([node_a5, node_a6])
            sage: d = egame_a._grow_info_set_graph_dictionary()
            sage: expected_outcome =  {(root_a,): [(node_a5,), (node_a6,)],
            ....:                      (node_a5,): [(node_a1, node_a2, node_a3, node_a4,)],
            ....:                      (node_a6,): [(node_a1, node_a2, node_a3, node_a4,)]}
            sage: d == expected_outcome
            True

        """
        d = {}
        for info_set_1 in self.info_sets:
            info_sets = []
            for info_set_2 in self.info_sets:
                if info_set_2 != info_set_1:
                    if any(kid in info_set_2 for nd in info_set_1 for kid in nd.children):
                        info_sets.append(tuple(info_set_2))
                        d[tuple(info_set_1)] = info_sets
        return d

    def _grow_info_set_graph(self):
        """
        A private method to grow information set DiGraph from the information
        sets of the game, returns the DiGraph.

        TESTS::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: leaf_a5 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a6 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a7 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a8 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a1, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a1, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: node_a3 = EFG_Node(player_a1, {'A': leaf_a5,
            ....:                                'B': leaf_a6}, 'Node 3')
            sage: node_a4 = EFG_Node(player_a1, {'A': leaf_a7,
            ....:                                'B': leaf_a8}, 'Node 4')
            sage: node_a5 = EFG_Node(player_a2, {'C': node_a1,
            ....:                                'D': node_a2}, 'Node 5')
            sage: node_a6 = EFG_Node(player_a2, {'C': node_a3,
            ....:                                'D': node_a4}, 'Node 6')
            sage: root_a = EFG_Node(player_a1, {'A': node_a5,
            ....:                               'B': node_a6}, 'Tree Root')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a._grow_info_set_graph()
            Digraph on 7 vertices

        If we change any information sets, the tree changes accordingly::

            sage: egame_a.set_info_set([node_a5, node_a6])
            sage: egame_a._grow_info_set_graph()
            Digraph on 6 vertices

            sage: egame_a.set_info_set([node_a3, node_a4])
            sage: egame_a.set_info_set([node_a1, node_a2])
            sage: egame_a._grow_info_set_graph()
            Digraph on 4 vertices

            sage: egame_a.set_info_set([node_a3, node_a4, node_a1, node_a2])
            sage: egame_a._grow_info_set_graph()
            Digraph on 3 vertices

            sage: egame_a.remove_info_set([node_a5, node_a6])
            sage: egame_a._grow_info_set_graph()
            Digraph on 4 vertices
        """
        d = self._grow_info_set_graph_dictionary()
        return DiGraph(d)

    def _check_node_names_and_find_players(self, generator):
        r"""
        A private method to check the names of the nodes and gives names for
        the ones that do not have names. This also finds all the players.

        This method is embedded in the init method but has been written here to
        improve the readability of the code. The tests for this are methodal
        tests.

        TESTS::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1 : 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1, 'B': leaf_a2})
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3, 'B': leaf_a4})
            sage: root_a = EFG_Node(player_a1, {'C': node_a1, 'D': node_a2})
            sage: node_a1.name is node_a2.name is root_a.name
            True
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: node_a1.name is node_a2.name is root_a.name
            False
            sage: sorted([node_a1.name, node_a2.name])
            ['Node 1', 'Node 2']
            sage: root_a.name
            'Tree Root'
            sage: sorted(egame_a.players, key=lambda x: x.name)
            [EFG Player "Player 1", EFG Player "Player 2"]

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1})
            sage: leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0})
            sage: leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4})
            sage: leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1})
            sage: node_b1 = EFG_Node(player_b2, {'A': leaf_b1, 'B': leaf_b2})
            sage: node_b2 = EFG_Node(player_b2, {'A': leaf_b3,
            ....:                                'B': leaf_b4}, 'Node B')
            sage: root_b = EFG_Node(player_b1, {'C': node_b1, 'D': node_b2})
            sage: node_b1.name is root_b.name
            True
            sage: node_b2.name
            'Node B'
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: node_b1.name is root_b.name
            False
            sage: sorted([node_b1.name, node_b2.name])
            ['Node 1', 'Node B']
            sage: root_b.name
            'Tree Root'
        """
        node_index = 1
        leaf_index = 1
        self.nodes.sort(key=attrgetter('parent', 'actions'))
        for node in self.nodes:
            if node.player not in self.players:
                self.players.append(node.player)
            if node is generator and node.name is None:
                node.name = "Tree Root"
            if node.name is None:
                node.name = "Node %i" % node_index
                node_index += 1
            for child in node.children:
                if isinstance(child, EFG_Leaf) and child.name is None:
                    child.name = "Leaf %i" % leaf_index
                    leaf_index += 1
                    self.leafs.append(child)
                elif isinstance(child, EFG_Leaf):
                    self.leafs.append(child)

    def __repr__(self):
        """
        Representation method for the Extensive Form Game::

        sage: player_a1 = EFG_Player('Player 1')
        sage: player_a2 = EFG_Player('Player 2')
        sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
        sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
        sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
        sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
        sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1, 'B': leaf_a2}, 'Node 1')
        sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3, 'B': leaf_a4}, 'Node 2')
        sage: root_a = EFG_Node(player_a1, {'C': node_a1, 'D': node_a2}, 'Root')
        sage: egame_a = ExtensiveFormGame(root_a)
        sage: egame_a
        Extensive Form Game with the following underlying tree:  {...}

        sage: player_b1 = EFG_Player('Player 1')
        sage: player_b2 = EFG_Player('Player 2')
        sage: leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1}, 'Leaf 1')
        sage: leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0}, 'Leaf 2')
        sage: leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4}, 'Leaf 3')
        sage: leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1}, 'Leaf 4')
        sage: leaf_b5 = EFG_Leaf({player_b1 : 0, player_b2: 1}, 'Leaf 5')
        sage: leaf_b6 = EFG_Leaf({player_b1 : 1, player_b2: 0}, 'Leaf 6')
        sage: leaf_b7 = EFG_Leaf({player_b1 : 2, player_b2: 4}, 'Leaf 7')
        sage: leaf_b8 = EFG_Leaf({player_b1 : 2, player_b2: 1}, 'Leaf 8')
        sage: node_b1 = EFG_Node(player_b1, {'A': leaf_b1, 'B': leaf_b2}, "Node 1")
        sage: node_b2 = EFG_Node(player_b1, {'A': leaf_b3, 'B': leaf_b4}, "Node 2")
        sage: node_b3 = EFG_Node(player_b1, {'A': leaf_b5, 'B': leaf_b6}, "Node 3")
        sage: node_b4 = EFG_Node(player_b1, {'A': leaf_b7, 'B': leaf_b8}, "Node 4")
        sage: node_b5 = EFG_Node(player_b2, {'C': node_b1, 'D': node_b2}, "Node 5")
        sage: node_b6 = EFG_Node(player_b2, {'C': node_b3, 'D': node_b4}, "Node 6")
        sage: root_b = EFG_Node(player_b1, {'A': node_b5, 'B': node_b6})
        sage: egame_b = ExtensiveFormGame(root_b)
        sage: egame_b
        Extensive Form Game with the following underlying tree: {...}
        """

        return "Extensive Form Game with the following underlying tree: " + str(self.tree_dictionary)

    def sage_to_gambit(self):
        r"""
        A method to convert a Sage ``ExtensiveFormGame`` into a Gambit
        ``Game``. Returns the Gambit ``Game`` with all the attributes taken
        from the Sage ``ExtensiveFormGame``. In order to convert a sage
        ``ExtensiveFormGame`` into a Gambit ``Game``, we have to set up the
        sage game as normal, setting up information sets we want, before using
        ``sage_to_gambit``::

            sage: from gambit import Game  # optional - gambit
            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: leaf_a1 = EFG_Leaf({player_a1: 0, player_a2: 1})
            sage: leaf_a2 = EFG_Leaf({player_a1: 1, player_a2: 0})
            sage: leaf_a3 = EFG_Leaf({player_a1: 2, player_a2: 4})
            sage: leaf_a4 = EFG_Leaf({player_a1: 2, player_a2: 1})
            sage: node_a1 = EFG_Node(player_a2, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1,
            ....:                               'D': node_a2}, 'Root')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.info_sets
            [[EFG Node "Node 1"],
             [EFG Node "Node 2"],
             [EFG Node "Root"]]
            sage: egame_a.set_info_set([node_a1, node_a2])
            sage: egame_a.info_sets
            [[EFG Node "Node 1",
              EFG Node "Node 2"],
             [EFG Node "Root"]]
            sage: gambit_egame_a = egame_a.sage_to_gambit()  # optional - gambit

        Attributes can be called from the gambit game using the gambit methods::

            sage: gambit_egame_a.players  # optional - gambit
            [<Player [0] 'Player 1' in game ''>,
             <Player [1] 'Player 2' in game ''>]
            sage: gambit_egame_a.infosets  # optional - gambit
            [<Infoset [0] '(EFG Node "Root")' for player 'Player 1' in game ''>,
             <Infoset [0] '(EFG Node "Node 1",
              EFG Node "Node 2")' for player 'Player 2' in game ''>]
            sage: gambit_egame_a.root  # optional - gambit
            <Node [1] 'Root' in game ''>
            sage: gambit_egame_a.outcomes  # optional - gambit
            [<Outcome [0] 'Leaf 1' in game ''>, <Outcome [1] 'Leaf 2' in game ''>,
             <Outcome [2] 'Leaf 3' in game ''>, <Outcome [3] 'Leaf 4' in game ''>]

        We can also call Gambit's .efg format of the game by simply calling the game::

            sage: gambit_egame_a  # optional - gambit
            EFG 2 R "" { "Player 1" "Player 2" }
            ""
            <BLANKLINE>
            p "Root" 1 1 "(EFG Node \"Root\")" { "C" "D" } 0
            p "Node 1" 2 1 "(EFG Node \"Node 1\",
             EFG Node \"Node 2\")" { "A" "B" } 0
            t "" 1 "Leaf 1" { 0, 1 }
            t "" 2 "Leaf 2" { 1, 0 }
            p "Node 2" 2 1 "(EFG Node \"Node 1\",
             EFG Node \"Node 2\")" { "A" "B" } 0
            t "" 3 "Leaf 3" { 2, 4 }
            t "" 4 "Leaf 4" { 2, 1 }
            <BLANKLINE>

        The following is a test to show that this works for larger trees too::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: leaf_b1 = EFG_Leaf({player_b1 : 0, player_b2: 1}, 'Leaf 1')
            sage: leaf_b2 = EFG_Leaf({player_b1 : 1, player_b2: 0}, 'Leaf 2')
            sage: leaf_b3 = EFG_Leaf({player_b1 : 2, player_b2: 4}, 'Leaf 3')
            sage: leaf_b4 = EFG_Leaf({player_b1 : 2, player_b2: 1}, 'Leaf 4')
            sage: leaf_b5 = EFG_Leaf({player_b1 : 0, player_b2: 1}, 'Leaf 5')
            sage: leaf_b6 = EFG_Leaf({player_b1 : 1, player_b2: 0}, 'Leaf 6')
            sage: leaf_b7 = EFG_Leaf({player_b1 : 2, player_b2: 4}, 'Leaf 7')
            sage: leaf_b8 = EFG_Leaf({player_b1 : 2, player_b2: 1}, 'Leaf 8')
            sage: node_b1 = EFG_Node(player_b1, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, "Node 1")
            sage: node_b2 = EFG_Node(player_b1, {'A': leaf_b3,
            ....:                                'B': leaf_b4}, "Node 2")
            sage: node_b3 = EFG_Node(player_b1, {'A': leaf_b5,
            ....:                                'B': leaf_b6}, "Node 3")
            sage: node_b4 = EFG_Node(player_b1, {'A': leaf_b7,
            ....:                                'B': leaf_b8}, "Node 4")
            sage: node_b5 = EFG_Node(player_b2, {'C': node_b1,
            ....:                                'D': node_b2}, "Node 5")
            sage: node_b6 = EFG_Node(player_b2, {'C': node_b3,
            ....:                                'D': node_b4}, "Node 6")
            sage: root_b = EFG_Node(player_b1, {'A': node_b5, 'B': node_b6})
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: egame_b.set_info_set([node_b5, node_b6])
            sage: egame_b.set_info_set([node_b1, node_b2, node_b3, node_b4])
            sage: gambit_egame_b = egame_b.sage_to_gambit()  # optional - gambit
            sage: gambit_egame_b.players  # optional - gambit
            [<Player [0] 'Player 1' in game ''>, <Player [1] 'Player 2' in game ''>]
            sage: gambit_egame_b.infosets  # optional - gambit
            [<Infoset [0] '(EFG Node "Tree Root")' for player 'Player 1' in game ''>,
             <Infoset [1] '(EFG Node "Node 1", EFG Node "Node 2",
              EFG Node "Node 3", EFG Node "Node 4")' for player 'Player 1' in game ''>,
             <Infoset [0] '(EFG Node "Node 5", EFG Node "Node 6")'
              for player 'Player 2' in game ''>]
            sage: gambit_egame_b.root  # optional - gambit
            <Node [1] 'Tree Root' in game ''>
            sage: gambit_egame_b.outcomes  # optional - gambit
            [<Outcome [0] 'Leaf 1' in game ''>,
             <Outcome [1] 'Leaf 2' in game ''>,
             <Outcome [2] 'Leaf 3' in game ''>,
             <Outcome [3] 'Leaf 4' in game ''>,
             <Outcome [4] 'Leaf 5' in game ''>,
             <Outcome [5] 'Leaf 6' in game ''>,
             <Outcome [6] 'Leaf 7' in game ''>,
             <Outcome [7] 'Leaf 8' in game ''>]
            sage: gambit_egame_b  # optional - gambit
            EFG 2 R "" { "Player 1" "Player 2" }
            ""
            <BLANKLINE>
            p "Tree Root" 1 1 "(EFG Node \"Tree Root\")" { "A" "B" } 0
            p "Node 5" 2 1 "(EFG Node \"Node 5\", EFG Node \"Node 6\")" { "C" "D" } 0
            p "Node 1" 1 2 "(EFG Node \"Node 1\", EFG Node \"Node 2\", EFG Node \"Node 3\",
             EFG Node \"Node 4\")" { "A" "B" } 0
            t "" 1 "Leaf 1" { 0, 1 }
            t "" 2 "Leaf 2" { 1, 0 }
            p "Node 2" 1 2 "(EFG Node \"Node 1\", EFG Node \"Node 2\", EFG Node \"Node 3\",
             EFG Node \"Node 4\")" { "A" "B" } 0
            t "" 3 "Leaf 3" { 2, 4 }
            t "" 4 "Leaf 4" { 2, 1 }
            p "Node 6" 2 1 "(EFG Node \"Node 5\", EFG Node \"Node 6\")" { "C" "D" } 0
            p "Node 3" 1 2 "(EFG Node \"Node 1\", EFG Node \"Node 2\", EFG Node \"Node 3\",
             EFG Node \"Node 4\")" { "A" "B" } 0
            t "" 5 "Leaf 5" { 0, 1 }
            t "" 6 "Leaf 6" { 1, 0 }
            p "Node 4" 1 2 "(EFG Node \"Node 1\", EFG Node \"Node 2\", EFG Node \"Node 3\",
             EFG Node \"Node 4\")" { "A" "B" } 0
            t "" 7 "Leaf 7" { 2, 4 }
            t "" 8 "Leaf 8" { 2, 1 }
            <BLANKLINE>

        This is a test with a tree that isn't symmetrical::

            sage: player_c1 = EFG_Player('Player 1')
            sage: player_c2 = EFG_Player('Player 2')
            sage: leaf_c1 = EFG_Leaf({player_c1 : 0, player_c2: 1}, 'Leaf 1')
            sage: leaf_c2 = EFG_Leaf({player_c1 : 1, player_c2: 0}, 'Leaf 2')
            sage: leaf_c3 = EFG_Leaf({player_c1 : 2, player_c2: 4}, 'Leaf 3')
            sage: leaf_c4 = EFG_Leaf({player_c1 : 2, player_c2: 1}, 'Leaf 4')
            sage: leaf_c5 = EFG_Leaf({player_c1 : 0, player_c2: 1}, 'Leaf 5')
            sage: leaf_c6 = EFG_Leaf({player_c1 : 1, player_c2: 0}, 'Leaf 6')
            sage: node_c1 = EFG_Node(player_c1, {'A': leaf_c1,
            ....:                                'B': leaf_c2}, "Node 1")
            sage: node_c2 = EFG_Node(player_c1, {'A': leaf_c3,
            ....:                                'B': leaf_c4}, "Node 2")
            sage: node_c5 = EFG_Node(player_c2, {'C': node_c1,
            ....:                                'D': node_c2}, "Node 3")
            sage: node_c6 = EFG_Node(player_c2, {'C': leaf_c5,
            ....:                                'D': leaf_c6}, "Node 4")
            sage: root_c = EFG_Node(player_c1, {'A': node_c5,
            ....:                               'B': node_c6}, "Root")
            sage: egame_c = ExtensiveFormGame(root_c)
            sage: egame_c.set_info_set([node_c1, node_c2])
            sage: gambit_egame_c = egame_c.sage_to_gambit()  # optional - gambit
            sage: gambit_egame_c.players  # optional - gambit
            [<Player [0] 'Player 1' in game ''>,
             <Player [1] 'Player 2' in game ''>]
            sage: gambit_egame_c.outcomes  # optional - gambit
            [<Outcome [0] 'Leaf 5' in game ''>,
             <Outcome [1] 'Leaf 6' in game ''>,
             <Outcome [2] 'Leaf 1' in game ''>,
             <Outcome [3] 'Leaf 2' in game ''>,
             <Outcome [4] 'Leaf 3' in game ''>,
             <Outcome [5] 'Leaf 4' in game ''>]
            sage: gambit_egame_c.root  # optional - gambit
            <Node [1] 'Root' in game ''>
            sage: gambit_egame_c  # optional - gambit
            EFG 2 R "" { "Player 1" "Player 2" }
            ""
            <BLANKLINE>
            p "Root" 1 1 "(EFG Node \"Root\")" { "A" "B" } 0
            p "Node 3" 2 1 "(EFG Node \"Node 3\",)" { "C" "D" } 0
            p "Node 1" 1 2 "(EFG Node \"Node 1\", EFG Node \"Node 2\")" { "A" "B" } 0
            t "" 3 "Leaf 1" { 0, 1 }
            t "" 4 "Leaf 2" { 1, 0 }
            p "Node 2" 1 2 "(EFG Node \"Node 1\", EFG Node \"Node 2\")" { "A" "B" } 0
            t "" 5 "Leaf 3" { 2, 4 }
            t "" 6 "Leaf 4" { 2, 1 }
            p "Node 4" 2 2 "(EFG Node \"Node 4\",)" { "C" "D" } 0
            t "" 1 "Leaf 5" { 0, 1 }
            t "" 2 "Leaf 6" { 1, 0 }
            <BLANKLINE>

        This is a test for a game with more than 2 players::

            sage: player_c1 = EFG_Player('Player 1')
            sage: player_c2 = EFG_Player('Player 2')
            sage: player_c3 = EFG_Player('Player 3')
            sage: leaf_c1 = EFG_Leaf({player_c1 : 0, player_c2: 1,
            ....:                     player_c3: -5}, 'Leaf 1')
            sage: leaf_c2 = EFG_Leaf({player_c1 : 1, player_c2: 0,
            ....:                     player_c3: -4}, 'Leaf 2')
            sage: leaf_c3 = EFG_Leaf({player_c1 : 2, player_c2: 4,
            ....:                     player_c3: -3}, 'Leaf 3')
            sage: leaf_c4 = EFG_Leaf({player_c1 : 2, player_c2: 1,
            ....:                     player_c3: -2}, 'Leaf 4')
            sage: node_c1 = EFG_Node(player_c3, {'A': leaf_c1,
            ....:                                'B': leaf_c2}, 'Node 1')
            sage: node_c2 = EFG_Node(player_c2, {'A': leaf_c3,
            ....:                                'B': leaf_c4}, 'Node 2')
            sage: root_c = EFG_Node(player_c1, {'C': node_c1,
            ....:                               'D': node_c2}, 'Root')
            sage: egame_c = ExtensiveFormGame(root_c)
            sage: gambit_egame_c = egame_c.sage_to_gambit()  # optional - gambit
            sage: gambit_egame_c.players  # optional - gambit
            [<Player [0] 'Player 1' in game ''>,
             <Player [1] 'Player 2' in game ''>,
             <Player [2] 'Player 3' in game ''>]
            sage: gambit_egame_c.root  # optional - gambit
            <Node [1] 'Root' in game ''>
            sage: gambit_egame_c.infosets  # optional - gambit
            [<Infoset [0] '(EFG Node "Root")' for player 'Player 1' in game ''>,
             <Infoset [0] '(EFG Node "Node 2",)' for player 'Player 2' in game ''>,
             <Infoset [0] '(EFG Node "Node 1",)' for player 'Player 3' in game ''>]
            sage: gambit_egame_c  # optional - gambit
            EFG 2 R "" { "Player 1" "Player 2" "Player 3" }
            ""
            <BLANKLINE>
            p "Root" 1 1 "(EFG Node \"Root\")" { "C" "D" } 0
            p "Node 1" 3 1 "(EFG Node \"Node 1\",)" { "A" "B" } 0
            t "" 1 "Leaf 1" { 0, 1, -5 }
            t "" 2 "Leaf 2" { 1, 0, -4 }
            p "Node 2" 2 1 "(EFG Node \"Node 2\",)" { "A" "B" } 0
            t "" 3 "Leaf 3" { 2, 4, -3 }
            t "" 4 "Leaf 4" { 2, 1, -2 }
            <BLANKLINE>
        """
        gambit_game = Game.new_tree()

        for player in self.players:
            gambit_game.players.add(player.name)

        info_set_dictionary = self._grow_info_set_graph_dictionary()
        gambit_root = gambit_game.root
        gambit_root.label = self.tree_root.name
        self._sage_to_gambit_set_gambit_move(self.tree_root, self.tree_root,
                                             gambit_root, gambit_game)
        self._sage_to_gambit_add_gambit_outcomes(self.tree_root, gambit_root,
                                                 gambit_game)
        tuple_root = (self.tree_root,)
        current_info_sets = info_set_dictionary[tuple_root]
        parent_dict = {}
        parent_dict[self.tree_root] = gambit_game.root

        while len(current_info_sets)!= 0:
            current_info_sets.sort(key = lambda x:x[0].name)
            self._sage_to_gambit_info_set_list = []
            for info_set in current_info_sets:
                for node in info_set:
                    child_index = self._sage_to_gambit_get_gambit_child_index(node)
                    gambit_node = parent_dict[node.parent].children[child_index]
                    gambit_node.label = node.name
                    if node is info_set[0]:
                        move = self._sage_to_gambit_set_gambit_move(node, info_set,
                                                                    gambit_node, gambit_game)
                    else:
                        gambit_node.append_move(move)
                    self._sage_to_gambit_add_gambit_outcomes(node, gambit_node,
                                                             gambit_game)
                    parent_dict[node] = gambit_node
                self._sage_to_gambit_setup_next_info_sets(info_set, info_set_dictionary)
            current_info_sets = list(set(self._sage_to_gambit_info_set_list))

        return gambit_game

    def _sage_to_gambit_set_gambit_move(self, sage_node, info_set, gambit_node, gambit_game):
        r"""
        A method for ``_sage_to_gambit`` which takes a sage node, and passes
        along its name, actions, and number of children to the corresponding
        gambit branch, and creates a move using this information, and then
        transfers this move to the corresponding Gambit ``Node``.

        The method can be used to create individual branches for a gambit game
        using information from a single sage ``EFG_Node`` per branch, to do
        this we must first set up a Sage ``ExtensiveFormGame``::

            sage: from gambit import Game  # optional - gambit
            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1: 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1: 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, 'Node 1')
            sage: node_2 = EFG_Node(player_2, {'A': leaf_3,
            ....:                              'B': leaf_4}, 'Node 2')
            sage: root = EFG_Node(player_1, {'C': node_1,
            ....:                            'D': node_2}, 'Root 1')
            sage: egame = ExtensiveFormGame(root)

        We then need a Gambit ``Game`` set up with players::

            sage: g = Game.new_tree()  # optional - gambit
            ....:
            sage: g.players.add('Player 1')  # optional - gambit
            <Player [0] 'Player 1' in game ''>
            sage: g.players.add('Player 2')  # optional - gambit
            <Player [1] 'Player 2' in game ''>

        Then we can set up the branches we want within the game, which takes
        information from the sage node, and passes it to the gambit node::

            sage: egame._sage_to_gambit_set_gambit_move(root, root, g.root, g)  # optional - gambit
            <Infoset [0] '(EFG Node "Root 1")' for player 'Player 1' in game ''>
            sage: egame._sage_to_gambit_set_gambit_move(node_1, [node_1, node_2], g.root.children[int(0)], g)  # optional - gambit
            <Infoset [0] '(EFG Node "Node 1", EFG Node "Node 2")' for player 'Player 2' in game ''>


        We can see that the gambit nodes and equivilant sage nodes share the
        same players::

            sage: g.root.player  # optional - gambit
            <Player [0] 'Player 1' in game ''>
            sage: g.root.children[int(0)].player  # optional - gambit
            <Player [1] 'Player 2' in game ''>

        We can also check node 'relationships'::

            sage: g.root.parent  # optional - gambit
            <BLANKLINE>
            sage: g.root.children[int(0)].parent  # optional - gambit
            <Node [1] '' in game ''>
        """

        move = gambit_node.append_move(gambit_game.players[sage_node.player.name],
                                       len(sage_node.children))
        if sage_node is self.tree_root:
            move.label = "(%s)" %self.tree_root
        else:
            move.label = str(tuple(info_set))
        for action_setting_index in range(len(sage_node.actions)):
            move.actions[action_setting_index].label = sage_node.actions[action_setting_index]
        return move


    def _sage_to_gambit_get_gambit_child_index(self, child_of_sage_node):
        r"""
        A method for ``_sage_to_gambit`` which takes takes a node/leaf with a
        parent (which has already been converted into a gambit node) and then
        returns which branch in the gambit game that child should lie on,
        which is determined by which action took it there.

        The method requires an sage extensive form game, and a gambit extensive
        form game set up::

            sage: from gambit import Game  # optional - gambit
            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1: 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1: 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, 'Node 1')
            sage: node_2 = EFG_Node(player_2, {'A': leaf_3,
            ....:                              'B': leaf_4}, 'Node 2')
            sage: root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root')
            sage: egame = ExtensiveFormGame(root)
            sage: g = Game.new_tree()  # optional - gambit
            sage: g.players.add('Player 1')  # optional - gambit
            <Player [0] 'Player 1' in game ''>
            sage: g.players.add('Player 2')  # optional - gambit
            <Player [1] 'Player 2' in game ''>

        We also need to ensure the root node is set up for the gambit game::

            sage: egame._sage_to_gambit_set_gambit_move(root, root, g.root, g)  # optional - gambit
            <Infoset [0] '(EFG Node "Root")' for player 'Player 1' in game ''>

        We can then use the method to check the 'child index'::

            sage: egame._sage_to_gambit_get_gambit_child_index(node_1)
            0
            sage: egame._sage_to_gambit_get_gambit_child_index(node_2)
            1
        """
        sage_node = child_of_sage_node.parent
        for i,action in enumerate(sage_node.actions):
            if sage_node.action_to_node_dict[action] is child_of_sage_node:
                child_index = i
        return child_index

    def _sage_to_gambit_add_gambit_outcomes(self, sage_node, gambit_node, gambit_game):
        r"""
        A method for ``_sage_to_gambit`` which takes an already converted
        gambit node, and searches through its children to determine if any are
        leaves, then if a child is a leaf it adds its payoffs to the game.

        The method requires an sage extensive form game, and a gambit extensive
        form game set up::

            sage: from gambit import Game  # optional - gambit
            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1: 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1: 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, 'Node 1')
            sage: node_2 = EFG_Node(player_2, {'A': leaf_3,
            ....:                              'B': leaf_4}, 'Node 2')
            sage: root = EFG_Node(player_1, {'C': node_1, 'D': node_2}, 'Root 1')
            sage: egame = ExtensiveFormGame(root)
            sage: g = Game.new_tree()  # optional - gambit
            sage: g.players.add('Player 1')  # optional - gambit
            <Player [0] 'Player 1' in game ''>
            sage: g.players.add('Player 2')  # optional - gambit
            <Player [1] 'Player 2' in game ''>

        We also need to ensure the some nodes are set up for the gambit game::

            sage: egame._sage_to_gambit_set_gambit_move(root, root, g.root, g)  # optional - gambit
            <Infoset [0] '(EFG Node "Root 1")' for player 'Player 1' in game ''>
            sage: egame._sage_to_gambit_set_gambit_move(node_1, [node_1, node_2], g.root.children[int(0)], g)  # optional - gambit
            <Infoset [0] '(EFG Node "Node 1", EFG Node "Node 2")' for player 'Player 2' in game ''>

        We can then then use the method to add outcomes to the gambit node::

            sage: g.outcomes  # optional - gambit
            []
            sage: egame._sage_to_gambit_add_gambit_outcomes(node_1, g.root.children[int(0)], g)  # optional - gambit
            sage: g.outcomes  # optional - gambit
            [<Outcome [0] 'Leaf 1' in game ''>,
             <Outcome [1] 'Leaf 2' in game ''>]
        """

        for child in sage_node.children:
            if isinstance(child, EFG_Leaf):
                outcome = gambit_game.outcomes.add(child.name)
                for player_index in range(len(self.players)):
                    player = self.players[player_index]
                    outcome[player_index] = Fraction(str(child[player]))
                leaf_action_index = self._sage_to_gambit_get_gambit_child_index(child)
                gambit_node.children[leaf_action_index].outcome = outcome

    def _sage_to_gambit_setup_next_info_sets(self, info_set, info_set_dictionary):
        r"""
        The method sets up the next list of information sets to be looped
        through in the ``_sage_to_gambit`` method.

        In order to use the method, we must mimic the stages of
        ``_sage_to_gambit`` by first setting up an ``ExtensiveFormGame``::

            sage: from gambit import Game  # optional - gambit
            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1: 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1: 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, 'Node 1')
            sage: node_2 = EFG_Node(player_2, {'A': leaf_3,
            ....:                              'B': leaf_4}, 'Node 2')
            sage: root = EFG_Node(player_1, {'C': node_1,
            ....:                            'D': node_2}, 'Root')
            sage: egame = ExtensiveFormGame(root)

        We then need a Gambit ``Game`` set up with players::

            sage: g = Game.new_tree()  # optional - gambit
            sage: g.players.add('Player 1')  # optional - gambit
            <Player [0] 'Player 1' in game ''>
            sage: g.players.add('Player 2')  # optional - gambit
            <Player [1] 'Player 2' in game ''>

        Then we can set up the branches we want within the game, which takes
        information from the sage node, and passes it to the gambit node::

            sage: egame._sage_to_gambit_set_gambit_move(root, root, g.root, g)  # optional - gambit
            <Infoset [0] '(EFG Node "Root")' for player 'Player 1' in game ''>

        Then we use the method and it'll append the next information
        sets to a list::

            sage: egame._sage_to_gambit_setup_next_info_sets((root,), egame._grow_info_set_graph_dictionary())
            sage: egame._sage_to_gambit_info_set_list
            [(EFG Node "Node 1",), (EFG Node "Node 2",)]
        """
        for key in info_set_dictionary:
            if info_set == key:
                for listed_info_set in info_set_dictionary[key]:
                    self._sage_to_gambit_info_set_list.append(listed_info_set)


    def obtain_nash(self):
        r"""
        A method to obtain Nash Equilibria for an ``ExtensiveFormGame``,
        currently uses Gambit and the LCP algorithm, hence Gambit must be
        installed to run this method. This method returns the Nash Equilibria
        in a format which will show a list of tuples containing Nash Equilibria,
        where in each tuple are tuples with information sets and the
        corresponding probability distributions, mapping actions to their
        probabilities

        To obtain the Nash Equilibria of an ``ExtensiveFormGame``, we firstly
        set up the game as normal::

            sage: from gambit import Game  # optional - gambit
            sage: player_1 = EFG_Player('1')
            sage: player_2 = EFG_Player('2')
            sage: leaf_a1 = EFG_Leaf({player_1: 2, player_2: 0}, 'Leaf 1')
            sage: leaf_a2 = EFG_Leaf({player_1: 3, player_2: 1}, 'Leaf 2')
            sage: leaf_a3 = EFG_Leaf({player_1: 4, player_2: 2}, 'Leaf 3')
            sage: leaf_a4 = EFG_Leaf({player_1: 3, player_2: 5}, 'Leaf 4')
            sage: leaf_a5 = EFG_Leaf({player_1: 4, player_2: 1}, 'Leaf 5')
            sage: node_a1 = EFG_Node(player_1, {'Z': leaf_a2,
            ....:                               'Y': leaf_a3}, 'd')
            sage: node_a2 = EFG_Node(player_2, {'D': leaf_a1,
            ....:                               'C': node_a1}, 'b')
            sage: node_a3 = EFG_Node(player_2, {'B': leaf_a4,
            ....:                               'A': leaf_a5}, 'c')
            sage: node_a4 = EFG_Node(player_1, {'X': node_a2,
            ....:                               'W': node_a3}, 'a')
            sage: example = ExtensiveFormGame(node_a4)

        Then we simply use the obtain_nash method::

            sage: expected_outcome = [((('a',), {'W': 1.0, 'X': 0.0}),
            ....:                      (('d',), {'Y': 0.5, 'Z': 0.5}),
            ....:                      (('b',), {'C': 0.0, 'D': 1.0}),
            ....:                      (('c',), {'A': 0.0, 'B': 1.0})),
            ....:                     ((('a',), {'W': 1.0, 'X': 0.0}),
            ....:                      (('d',), {'Y': 0.5, 'Z': 0.5}),
            ....:                      (('b',), {'C': 0.5, 'D': 0.5}),
            ....:                      (('c',), {'A': 0.0, 'B': 1.0})),
            ....:                     ((('a',), {'W': 0.0, 'X': 1.0}),
            ....:                      (('d',), {'Y': 1.0, 'Z': 0.0}),
            ....:                      (('b',), {'C': 1.0, 'D': 0.0}),
            ....:                      (('c',), {'A': 0.0, 'B': 1.0}))]
            sage: expected_outcome == example.obtain_nash()   # optional - gambit
            True

        Here is an example with a different tree::

            sage: leaf_b1 = EFG_Leaf({player_1: 1, player_2: 5})
            sage: leaf_b2 = EFG_Leaf({player_1: 5, player_2: 2})
            sage: leaf_b3 = EFG_Leaf({player_1: 9, player_2: 1})
            sage: leaf_b4 = EFG_Leaf({player_1: 3, player_2: 0})
            sage: leaf_b5 = EFG_Leaf({player_1: 2, player_2: 7})
            sage: leaf_b6 = EFG_Leaf({player_1: 1, player_2: 5})
            sage: node_b3 = EFG_Node(player_1, {'f': leaf_b2,
            ....:                               'g': leaf_b3}, 'Node A')
            sage: node_b2 = EFG_Node(player_2, {'d': leaf_b1,
            ....:                               'e': node_b3}, 'Node B')
            sage: node_b4 = EFG_Node(player_2, {'h': leaf_b5,
            ....:                               'i': leaf_b6}, 'Node C')
            sage: node_b1 = EFG_Node(player_1, {'a': node_b2, 'b': leaf_b4,
            ....:                               'c': node_b4}, 'Tree Root')
            sage: example_b = ExtensiveFormGame(node_b1)
            sage: expected_outcome = [((('Tree Root',), {'a': 0.0, 'b': 1.0, 'c': 0.0}),
            ....:                      (('Node A',), {'f': 0.5, 'g': 0.5}),
            ....:                      (('Node B',), {'d': 1.0, 'e': 0.0}),
            ....:                      (('Node C',), {'h': 1.0, 'i': 0.0}))]
            sage: example_b.obtain_nash() == expected_outcome  # optional - gambit
            True

        The following is a test to show that this works for larger trees too::

            sage: player_c1 = EFG_Player('Player 1')
            sage: player_c2 = EFG_Player('Player 2')
            sage: leaf_c1 = EFG_Leaf({player_c1 : 0, player_c2: 2}, 'Leaf 1')
            sage: leaf_c2 = EFG_Leaf({player_c1 : 2, player_c2: 0}, 'Leaf 2')
            sage: leaf_c3 = EFG_Leaf({player_c1 : 2, player_c2: 8}, 'Leaf 3')
            sage: leaf_c4 = EFG_Leaf({player_c1 : 4, player_c2: 1}, 'Leaf 4')
            sage: leaf_c5 = EFG_Leaf({player_c1 : 0, player_c2: 1}, 'Leaf 5')
            sage: leaf_c6 = EFG_Leaf({player_c1 : 1, player_c2: 0}, 'Leaf 6')
            sage: leaf_c7 = EFG_Leaf({player_c1 : 2, player_c2: 4}, 'Leaf 7')
            sage: leaf_c8 = EFG_Leaf({player_c1 : 2, player_c2: 1}, 'Leaf 8')
            sage: node_c1 = EFG_Node(player_c1, {'A': leaf_c1,
            ....:                                'B': leaf_c2}, "Node 1")
            sage: node_c2 = EFG_Node(player_c1, {'A': leaf_c3,
            ....:                                'B': leaf_c4}, "Node 2")
            sage: node_c3 = EFG_Node(player_c1, {'A': leaf_c5,
            ....:                                'B': leaf_c6}, "Node 3")
            sage: node_c4 = EFG_Node(player_c1, {'A': leaf_c7,
            ....:                                'B': leaf_c8}, "Node 4")
            sage: node_c5 = EFG_Node(player_c2, {'C': node_c1,
            ....:                                'D': node_c2}, "Node 5")
            sage: node_c6 = EFG_Node(player_c2, {'C': node_c3,
            ....:                                'D': node_c4}, "Node 6")
            sage: root_c = EFG_Node(player_c1, {'A': node_c5, 'B': node_c6})
            sage: example_c = ExtensiveFormGame(root_c)
            sage: example_c.set_info_set([node_c5, node_c6])
            sage: example_c.set_info_set([node_c1, node_c2])
            sage: example_c.set_info_set([node_c3, node_c4])
            sage: expected_outcome = [((('Tree Root',), {'A': 1.0, 'B': 0.0}),
            ....:                      (('Node 1', 'Node 2'), {'A': 0.0, 'B': 1.0}),
            ....:                      (('Node 3', 'Node 4'), {'A': 0.5, 'B': 0.5}),
            ....:                      (('Node 5', 'Node 6'), {'C': 0.0, 'D': 1.0}))]
            sage: expected_outcome == example_c.obtain_nash()  # optional - gambit
            True

        Curently we cannot obtain Nash Equilibria for games with more than two
        players, so an error will be returned in these instances::

            sage: player_d1 = EFG_Player('Player 1')
            sage: player_d2 = EFG_Player('Player 2')
            sage: player_d3 = EFG_Player('Player 3')
            sage: leaf_d1 = EFG_Leaf({player_d1 : 0, player_d2: 1,
            ....:                     player_d3: -5}, 'Leaf 1')
            sage: leaf_d2 = EFG_Leaf({player_d1 : 1, player_d2: 0,
            ....:                     player_d3: -4}, 'Leaf 2')
            sage: leaf_d3 = EFG_Leaf({player_d1 : 2, player_d2: 4,
            ....:                     player_d3: -3}, 'Leaf 3')
            sage: leaf_d4 = EFG_Leaf({player_d1 : 2, player_d2: 1,
            ....:                     player_d3: -2}, 'Leaf 4')
            sage: node_d1 = EFG_Node(player_d3, {'A': leaf_d1,
            ....:                                'B': leaf_d2}, 'Node 1')
            sage: node_d2 = EFG_Node(player_d2, {'A': leaf_d3,
            ....:                                'B': leaf_d4}, 'Node 2')
            sage: root_d = EFG_Node(player_d1, {'C': node_d1,
            ....:                               'D': node_d2}, 'Root')
            sage: egame_d = ExtensiveFormGame(root_d)
            sage: egame_d.obtain_nash()  # optional - gambit
            Traceback (most recent call last):
            ...
            NotImplementedError: Nash equilibrium for games with more than 2
            players have not been implemented yet.
        """
        from gambit.nash import ExternalLCPSolver
        if Game is None:
            raise NotImplementedError("gambit is not installed")
        if len(self.players) > 2:
            raise NotImplementedError("Nash equilibrium for games with more "
                                      "than 2 players have not been "
                                      "implemented yet.")

        gambit_efg = self.sage_to_gambit()
        solver = ExternalLCPSolver()
        lcp_output = solver.solve(gambit_efg)
        nasheq = Parser(lcp_output).format_gambit_efg_tree(gambit_efg)
        # I am not 100% convinced this is done in the right way.
        return nasheq

    def is_constant_sum(self):
        r"""
        This is a method which returns whether a game is constant sum or not
        (i.e. the sum of all the utilites for each leaf is the same for every
        leaf in the game).
        The following game is a constant sum game::

            sage: player_a1 = EFG_Player('Player 1')
            sage: player_a2 = EFG_Player('Player 2')
            sage: player_a3 = EFG_Player('Player 3')
            sage: leaf_a1 = EFG_Leaf({player_a1 : 0, player_a2: 5,
            ....:                     player_a3: -5}, 'Leaf 1')
            sage: leaf_a2 = EFG_Leaf({player_a1 : 1, player_a2: 3,
            ....:                     player_a3: -4}, 'Leaf 2')
            sage: leaf_a3 = EFG_Leaf({player_a1 : 2, player_a2: 1,
            ....:                     player_a3: -3}, 'Leaf 3')
            sage: leaf_a4 = EFG_Leaf({player_a1 : 1, player_a2: 1,
            ....:                     player_a3: -2}, 'Leaf 4')
            sage: node_a1 = EFG_Node(player_a3, {'A': leaf_a1,
            ....:                                'B': leaf_a2}, 'Node 1')
            sage: node_a2 = EFG_Node(player_a2, {'A': leaf_a3,
            ....:                                'B': leaf_a4}, 'Node 2')
            sage: root_a = EFG_Node(player_a1, {'C': node_a1,
            ....:                               'D': node_a2}, 'Root')
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.is_constant_sum()
            True

        The following game is not a constant sum game::

            sage: player_b1 = EFG_Player('Player 1')
            sage: player_b2 = EFG_Player('Player 2')
            sage: player_b3 = EFG_Player('Player 3')
            sage: leaf_b1 = EFG_Leaf({player_a1 : 0, player_b2: 5,
            ....:                                    player_b3: -5}, 'Leaf 1')
            sage: leaf_b2 = EFG_Leaf({player_a1 : 1, player_b2: 3,
            ....:                                    player_b3: -4}, 'Leaf 2')
            sage: leaf_b3 = EFG_Leaf({player_a1 : 2, player_b2: 1,
            ....:                                    player_b3: -3}, 'Leaf 3')
            sage: leaf_b4 = EFG_Leaf({player_a1 : 1, player_b2: 0,
            ....:                                    player_b3: -2}, 'Leaf 4')
            sage: node_b1 = EFG_Node(player_a3, {'A': leaf_b1,
            ....:                                'B': leaf_b2}, 'Node 1')
            sage: node_b2 = EFG_Node(player_a2, {'A': leaf_b3,
            ....:                                'B': leaf_b4}, 'Node 2')
            sage: root_b = EFG_Node(player_a1, {'C': node_b1,
            ....:                               'D': node_b2}, 'Root')
            sage: egame_b = ExtensiveFormGame(root_b)
            sage: egame_b.is_constant_sum()
            False
        """
        return len(set([sum(l.utilities) for l in self.leafs])) == 1


class EFG_Node():
    def __init__(self, player=False, args={}, name=None):
        r"""
        An ``EFG_Node`` is specifically for creating an ``ExtensiveFormGame``.
        Node ``args`` must be in a dictionary format, consisting of the actions
        and the children of that node::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: child_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = EFG_Node(player_1, {'Action1': child_1,
            ....:                                   'Action2': child_2}, 'Mother')
            sage: mother_node.actions
            ['Action1', 'Action2']
            sage: mother_node.children
            [EFG Leaf with utilities - (0, 1),
             EFG Leaf with utilities - (1, 0)]

        If we then create a second node, who has :code:`mother_node` as one of
        its children, then the parent of :code:`mother_node` will be set to
        that node::

            sage: sis_node = EFG_Node()
            sage: mother_node.parent
            False
            sage: grandmother_node = EFG_Node(False, {'ActionA': mother_node,
            ....:                                     'ActionB': sis_node}, 'Node A')
            sage: mother_node.parent
            EFG Node "Node A"

        Nodes automatically have player set to ``False``, if we do not assign a
        player upon initalisation, we can simply set it up afterwards::

            sage: grandmother_node = EFG_Node()
            sage: grandmother_node.player
            False
            sage: grandmother_node.player = player_1
            sage: grandmother_node.player
            EFG Player "Player 1"

        If we try to pass an ``args`` that isn't a dictionary, an error is
        returned::

            sage: grandmother_node = EFG_Node(args = 5)
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an args in the form of a dictionary.

            sage: grandmother_node = EFG_Node(args = 'This is a string')
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an args in the form of a dictionary.

            sage: grandmother_node = EFG_Node(args = matrix([[1, 1], [1, 1]]))
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an args in the form of a dictionary.

            sage: sis_node = EFG_Node()
            sage: grandmother_node = EFG_Node(args = sis_node)
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an args in the form of a dictionary.
        """

        self.action_to_node_dict = {}
        self.player = player
        self.name = name
        self.actions = []
        self.children = False
        self.parent = False
        self.is_root = False

        if type(args) is dict:
            self.action_to_node_dict = args
            self.actions = args.keys()
            self.children = args.values()
            for child in self.children:
                child.parent = self

        else:
            raise TypeError("Node must be passed an args in the form of a "
                            "dictionary.")


    def __repr__(self):
        """
        Representation method for the Node::

            sage: repr_node = EFG_Node()
            sage: repr_node
            Unnamed EFG Node
            sage: repr_node.name = "A named Node"
            sage: repr_node
            EFG Node "A named Node"
            sage: repr_node = EFG_Node(name = "A different name")
            sage: repr_node
            EFG Node "A different name"
        """
        if self.name is None:
            return 'Unnamed EFG Node'
        else:
            return 'EFG Node "' + self.name + '"'

    def _is_complete(self):
        r"""
        This is a private method to ensure all nodes have the following
        If we create a node where their children aren't specified and no parent
        is set, the node is considered incomplete::

            sage: b = EFG_Node()
            sage: b._is_complete == True
            False

        However, when we do specify all those attributes, the node is then
        considered complete::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: child_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = EFG_Node(player_1, {'Action1': child_1,
            ....:                                   'Action2': child_2}, 'Node B')
            sage: sis_node = EFG_Node()
            sage: grandmother_node = EFG_Node(False, {'ActionA': mother_node,
            ....:                                     'ActionB': sis_node}, 'Node A')
            sage: mother_node._is_complete()
            True

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: child_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = EFG_Node(False, {'Action1': child_1,
            ....:                                'Action2': child_2}, 'Node B')
            sage: sis_node = EFG_Node()
            sage: grandmother_node = EFG_Node(False, {'ActionA': mother_node,
            ....:                                     'ActionB': sis_node}, 'Node A')
            sage: mother_node._is_complete()
            False

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: child_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = EFG_Node(player_1, {'Action1': child_1,
            ....:                                   'Action2': child_2}, 'Node B')
            sage: mother_node._is_complete()
            False

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: child_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = EFG_Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = EFG_Node(player_1, {'Action1': child_1,
            ....:                                   'Action2': child_2}, 'Node B')
            sage: mother_node.children = False
            sage: sis_node = EFG_Node()
            sage: grandmother_node = EFG_Node(False, {'ActionA': mother_node,
            ....:                                     'ActionB': sis_node}, 'Node A')
            sage: mother_node._is_complete()
            False
        """
        return all([self.parent , self.actions, self.children, self.player])

    def is_terminal(self):
        """
        Checks to see if a node is terminal (i.e. only has leafs as children)
        and if it is, the method return true::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = EFG_Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = EFG_Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = EFG_Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: leaf_5 = EFG_Leaf({player_1 : 0, player_2: 1}, 'Leaf 5')
            sage: node_3 = EFG_Node(player_1, {'C': leaf_3,
            ....:                              'D': leaf_4}, "Node 3")
            sage: node_2 = EFG_Node(player_2, {'A': leaf_1,
            ....:                              'B': leaf_2}, "Node 2")
            sage: node_1 = EFG_Node(player_2, {'A': leaf_5,
            ....:                              'B': node_3}, "Node 1")
            sage: root = EFG_Node(player_1, {'E': node_1, 'F': node_2}, "Root")

        Once we have build nodes we can check if each one is terminal::

            sage: root.is_terminal()
            False
            sage: node_1.is_terminal()
            False
            sage: node_2.is_terminal()
            True
            sage: node_3.is_terminal()
            True
        """
        leaf_count = 0
        for child in self.children:
            if isinstance(child, EFG_Leaf):
                leaf_count += 1
        if leaf_count is len(self.children):
            return True
        else:
            return False

    def _player_check(self):
        r"""
        A private method primarily used later for creating of Extensive Form Games::

            sage: grandmother_node = EFG_Node()
            sage: mother_node = EFG_Node()
            sage: grandmother_node.player
            False
            sage: grandmother_node.player = mother_node
            sage: grandmother_node._player_check()
            Traceback (most recent call last):
            ...
            TypeError: Cannot assign a non Player as a player to a Node.
        """
        if self.player is not False:
            if not isinstance(self.player, EFG_Player):
                raise TypeError("Cannot assign a non Player as a player to a Node.")


class EFG_Leaf():
    def __init__(self, payoffs, name=None):
        """
        An ``EFG_Leaf`` is specifically for creating an ``ExtensiveFormGame``.
        Each Leaf within the ``ExtensiveFormGame`` must be unique instances of
        this class and contain a ``payoffs`` dictionary which maps players to
        outcomes.

        In order to correctly create a Leaf, we must first create players in
        which to use::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')

        The payoffs must be in dictionary form such that the keys are players,
        and the values are of the following types: ``Rational``, ``RealLiteral``,
        ``Integer``, ``float``, ``int``, ``Fraction`` or ``Decimal`` ::

            sage: leaf_1 = EFG_Leaf({player_1: 0, player_2: 1})

        We can see which players the Leaf has utilites for::

            sage: leaf_1.players
            [EFG Player "Player 1", EFG Player "Player 2"]

        We can also simply check the utilites for that Leaf::

            sage: leaf_1.utilities
            (0, 1)

        We can check the utility for a specific player in two ways::

            sage: leaf_1.payoffs[player_1]
            0
            sage: leaf_1[player_1]
            0
            sage: leaf_1[player_2]
            1

        If we try to create leaves where the keys in the payoff dictionary
        aren't players, we get an error::

            sage: node_1 = EFG_Node(); node_2 = EFG_Node()
            sage: leaf_1 = EFG_Leaf({node_1: 0, node_2: 1})
            Traceback (most recent call last):
            ...
            TypeError: EFG_Leaf input must be in dictionary form with each key
            as a EFG_Player.

        Simarly, if we try to create a leaf with utilities that aren't the
        correct type, we get an error::

            sage: leaf_1 = EFG_Leaf({player_1: 'Win', player_2: 'Lose'})
            Traceback (most recent call last):
            ...
            TypeError: EFG_Leaf utilities must be of the following types:
            Rational, RealLiteral, Integer, float, int,
            Fraction or Decimal.

        The payoffs can only be inputted in dictionary form::

            sage: leaf_1 = EFG_Leaf([0, 1])
            Traceback (most recent call last):
            ...
            TypeError: EFG_Leaf input must be in dictionary form.
        """

        if type(payoffs) is not dict:
            raise TypeError("EFG_Leaf input must be in dictionary form.")

        self.payoffs =  payoffs
        self.name = name
        self.players = sorted(payoffs.keys(), key=lambda x:x.name)
        self.parent = False
        self.utilities = tuple([self[player] for player in sorted(self.players, key=lambda x:x.name)])

        for utility in self.utilities:
            if (type(utility) is float) or (type(utility) is int) \
            or (isinstance(utility, Rational)) \
            or (isinstance(utility, RealLiteral)) \
            or (isinstance(utility, Integer)) \
            or isinstance(utility, Fraction) \
            or isinstance(utility, Decimal):
                utility = float(utility)
            else:
                raise TypeError("EFG_Leaf utilities must be of the following types: "
                                "Rational, RealLiteral, Integer, float, int, "
                                "Fraction or Decimal.")

        for player in self.players:
            if not isinstance(player, EFG_Player):
                raise TypeError("EFG_Leaf input must be in dictionary form "
                                "with each key as a EFG_Player.")

    def __repr__(self):
        """
        Representation method for the leaf::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1: 0, player_2: 1}, 'end_leaf')
            sage: leaf_1
            EFG Leaf with utilities - (0, 1)

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_1 = EFG_Leaf({player_1: 5, player_2: 1}, 'end_leaf')
            sage: leaf_1
            EFG Leaf with utilities - (5, 1)

            sage: player_1 = EFG_Player('Vince')
            sage: player_2 = EFG_Player('Hannah')
            sage: player_3 = EFG_Player('James')
            sage: leaf_1 = EFG_Leaf({player_1: 5, player_2: 1,
            ....:                    player_3: 10}, 'end_leaf')
            sage: leaf_1
            EFG Leaf with utilities - (1, 10, 5)
        """
        return 'EFG Leaf with utilities - ' + str(self.utilities)


    def __delitem__(self, key):
        """
        This method is one of a collection that aims to make a leaf
        instance behave like a dictionary.

        Here we set up deleting an element of the payoffs dictionary::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_dict = EFG_Leaf({player_1: 0, player_2: 1})
            sage: del(leaf_dict[player_2])
            sage: leaf_dict.payoffs
            {EFG Player "Player 1": 0}
        """
        self.payoffs.pop(key, None)

    def __getitem__(self, key):
        """
        This method is one of a collection that aims to make a leaf
        instance behave like a dictionary.

        Here we allow for querying a key::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_dict = EFG_Leaf({player_1: 0, player_2: 1})
            sage: leaf_dict[player_1]
            0
            sage: del(leaf_dict[player_1])
            sage: leaf_dict[player_1]
            Traceback (most recent call last):
            ...
            KeyError: EFG Player "Player 1"

        """
        return self.payoffs[key]

    def __iter__(self):
        """
        This method is one of a collection that aims to make a game
        instance behave like a dictionary.

        Here we allow for iteration over the leaf to correspond to
        iteration over keys of the utility dictionary::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_dict = EFG_Leaf({player_1: 0, player_2: 1})
            sage: for key in leaf_dict:
            ....:     print "The player: {}, has payoff {}".format(key, leaf_dict[key])
            The player: EFG Player "Player 2", has payoff 1
            The player: EFG Player "Player 1", has payoff 0
        """
        return iter(self.payoffs)

    def __setitem__(self, key, value):
        """
        This method is one of a collection that aims to make a game
        instance behave like a dictionary.

        Here we set up setting the value of a key::

            sage: player_1 = EFG_Player('Player 1')
            sage: player_2 = EFG_Player('Player 2')
            sage: leaf_dict = EFG_Leaf({player_1: 0, player_2: 1})
            sage: leaf_dict[player_1]
            0
            sage: leaf_dict[player_1] = 2
            sage: leaf_dict[player_1]
            2
        """
        self.payoffs[key] = value


class EFG_Player():
    def __init__(self, name = None):
        """
        An ``EFG_Player`` is specifically for creating an ``ExtensiveFormGame``.

        We can create players and assign them names, either during
        initalisation or afterwards::

            sage: ben_player = EFG_Player('Benjamin')
            sage: ben_player.name
            'Benjamin'
            sage: tim_player = EFG_Player()
            sage: tim_player
            Unnamed EFG Player
            sage: tim_player.name = 'Tim'
            sage: tim_player.name
            'Tim'

        Players can be assigned to an ``EFG_Node`` either during initalisation
        or afterwards::

            sage: jack_1 = EFG_Node()
            sage: jack_1.player = EFG_Player('Jack')
            sage: jack_1.player
            EFG Player "Jack"

        If a node is not specificed a player, then this should return false::

            sage: sam_1 = EFG_Node()
            sage: sam_1.player
            False
            sage: sam_player = EFG_Player('Sam')
            sage: sam_1.player = sam_player
            sage: sam_1.player
            EFG Player "Sam"
            sage: sam_2 = EFG_Node(sam_player)
            sage: sam_2.player
            EFG Player "Sam"

        A player can be reassigned for Nodes::

            sage: andy_player = EFG_Player('Andy')
            sage: sam_2.player = andy_player
            sage: sam_2.player
            EFG Player "Andy"
        """

        self.name = name

    def __repr__(self):
        """
        Representation method for the player::

            sage: apple_1 = EFG_Player('Apple')
            sage: apple_1
            EFG Player "Apple"

            sage: apple_2 = EFG_Player()
            sage: apple_2
            Unnamed EFG Player
        """

        if self.name is None:
            return "Unnamed EFG Player"
        else:
            return 'EFG Player "' + str(self.name) + '"'

    def __hash__(self):
        """
        Makes the player class hashable::

            sage: apple_1 = EFG_Player('Apple')
            sage: banana_1 = EFG_Leaf({apple_1 : 0})
            sage: banana_1.payoffs
            {EFG Player "Apple": 0}
        """
        return hash(self.name)
