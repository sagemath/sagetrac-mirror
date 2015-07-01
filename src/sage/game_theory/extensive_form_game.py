"""
Extensive Form Games with N players.

This module implements a class for extensive form games
[NN2007]_. Graphical representations of the games are implemented and solution
algorithms are being developed (with an interface to gambit).

A well known example that can be implemented as an extensive form game is the
battle of the sexes. Consider two players, Celine and Bob. The two are deciding
on how to spend their evening, they can either watch Sports or go and see a
Comedy. Bob would prefer to see a Comedy, Celine would prefer to watch a Sports
movie. Depending on who choses what, there are different payoffs, this can be
demonstrated in tree Form.

.. PLOT::
    :width: 500 px

    player_1 = Player('Bob')
    player_2 = Player('Celine')
    leaf_1 = Leaf({player_1: 2, player_2: 3})
    leaf_2 = Leaf({player_1: 0, player_2: 0})
    leaf_3 = Leaf({player_1: 1, player_2: 1})
    leaf_4 = Leaf({player_1: 3, player_2: 2})
    node_3 = Node({'Sports': leaf_1, 'Comedy': leaf_2}, 'c', player_2)
    node_2 = Node({'Sports': leaf_3, 'Comedy': leaf_4}, 'b', player_2)
    node_1 = Node({'Sports': node_3, 'Comedy': node_2}, 'a', player_1)
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    battle_of_the_sexes.set_info_set([node_2, node_3])
    p = battle_of_the_sexes.plot(view_info_sets=True)
    sphinx_plot(p)

We can see there are three nodes, one for Bob, two for Celine. Connecting those
nodes are actions. These actions represent choices made by one player, and the
actions then lead on to a node of another player.  So Bob either chooses Sports
or Comedy, and Celine chooses Sports or Comedy. The location of the payoffs
correspond to the leaf of the underlying tree and they show the outcome for each
player.  So if Bob chooses Sports, and Celine chooses Sports, we see the payoff
is (2, 3), which represents Bob getting a payoff of 2, and Celine getting a
payoff of 3. Note that the green line between Celine's two nodes indicate that
they are in the same 'information set', in other words, Celine does not know
what Bob has picked. Thus the following game corresponds to a different
situation which is easier for both players to coordinate:

.. PLOT::
    :width: 500 px

    player_1 = Player('Bob')
    player_2 = Player('Celine')
    leaf_1 = Leaf({player_1: 2, player_2: 3})
    leaf_2 = Leaf({player_1: 0, player_2: 0})
    leaf_3 = Leaf({player_1: 1, player_2: 1})
    leaf_4 = Leaf({player_1: 3, player_2: 2})
    node_3 = Node({'Sports': leaf_1, 'Comedy': leaf_2}, 'c', player_2)
    node_2 = Node({'Sports': leaf_3, 'Comedy': leaf_4}, 'b', player_2)
    node_1 = Node({'Sports': node_3, 'Comedy': node_2}, 'a', player_1)
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    p = battle_of_the_sexes.plot(view_info_sets=True)
    sphinx_plot(p)

The first game (with the information set) corresponds to the following normal
form game (which are also implemented in Sage)::

    sage: A = matrix([[3, 1], [0, 2]])
    sage: B = matrix([[2, 1], [0, 3]])
    sage: battle_of_the_sexes = NormalFormGame([A, B])
    sage: battle_of_the_sexes
    Normal Form Game with the following utilities: {(0, 1): [1, 1], (1, 0): [0, 0], (0, 0): [3, 2], (1, 1): [2, 3]}

To generate an extensive form game we need to generate the nodes and assign them
to players as well as describing the actions they have and to which node each
action goes. As such it makes sense to start with the terminal nodes of the
tree, but the initial step is to create players as each node will map players to
utilities::

    sage: player_1, player_2 = Player('Bob'), Player('Celine')
    sage: player_1, player_2
    (Bob, Celine)

Once we have done this, we are ready to create our leafs::

    sage: leaf_1 = Leaf({player_1: 2, player_2: 3})
    sage: leaf_1
    (2, 3)
    sage: leaf_2 = Leaf({player_1: 0, player_2: 0})
    sage: leaf_2
    (0, 0)
    sage: leaf_3 = Leaf({player_1: 1, player_2: 1})
    sage: leaf_3
    (1, 1)
    sage: leaf_4 = Leaf({player_1: 3, player_2: 2})
    sage: leaf_4
    (3, 2)

We can then create the parents of these leafs, the general :code:`Node` class
takes 3 arguments: a dictionary mapping actions to other nodes, a name for the
node and finally the player who makes the decision at this node::

    sage: node_1 = Node({'Sports': leaf_3, 'Comedy': leaf_4}, 'b', player_2)
    sage: node_1
    b
    sage: node_1.player
    Celine
    sage: node_1.name
    'b'
    sage: node_2 = Node({'Sports': leaf_1, 'Comedy': leaf_2}, 'c', player_2)
    sage: node_2
    c
    sage: node_2.player
    Celine
    sage: node_2.name
    'c'

Finally, we create the root of the tree::

    sage: root = Node({'Sports': node_2, 'Comedy': node_1}, 'a', player_1)
    sage: root
    a
    sage: root.player
    Bob
    sage: root.name
    'a'

The extensive form game can then be created by passing this root (which
recursively has all required information)::

    sage: battle_of_the_sexes = ExtensiveFormGame(root)
    sage: battle_of_the_sexes
    <sage.game_theory.extensive_form_game.ExtensiveFormGame instance ...

By default all nodes are in their own information set. If we plot the tree we
see this::

    sage: battle_of_the_sexes.plot(view_info_sets=True)
    Graphics object consisting of 23 graphics primitives

Here is the output (this is the same tree as above):

.. PLOT::
    :width: 500 px

    player_1 = Player('Bob')
    player_2 = Player('Celine')
    leaf_1 = Leaf({player_1: 2, player_2: 3})
    leaf_2 = Leaf({player_1: 0, player_2: 0})
    leaf_3 = Leaf({player_1: 1, player_2: 1})
    leaf_4 = Leaf({player_1: 3, player_2: 2})
    node_3 = Node({'Sports': leaf_1, 'Comedy': leaf_2}, 'c', player_2)
    node_2 = Node({'Sports': leaf_3, 'Comedy': leaf_4}, 'b', player_2)
    node_1 = Node({'Sports': node_3, 'Comedy': node_2}, 'a', player_1)
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    p = battle_of_the_sexes.plot(view_info_sets=True)
    sphinx_plot(p)

An extensive form game where all nodes are in their own information set is said
to have 'perfect information'. The game we have so far still has perfect
information::

    sage: battle_of_the_sexes.info_sets
    [[a], [b], [c]]
    sage: battle_of_the_sexes.perfect_info()
    True

To set the information sets as described above we use the :code:`set_info_set`
method::

    sage: battle_of_the_sexes.set_info_set([node_1, node_2])
    sage: battle_of_the_sexes.info_sets
    [[a], [b, c]]

Now the game does not have perfect information::

    sage: battle_of_the_sexes.perfect_info()
    False

Information sets are demonstrated visually on the graph we plot by setting
```view_info_sets``` to be ```True``` while plotting::

    sage: battle_of_the_sexes.plot(view_info_sets = True)
    Graphics object consisting of 23 graphics primitives

Which will be plotted as follows:

.. PLOT::
    :width: 500 px

    player_1 = Player('Bob')
    player_2 = Player('Celine')
    leaf_1 = Leaf({player_1: 2, player_2: 3})
    leaf_2 = Leaf({player_1: 0, player_2: 0})
    leaf_3 = Leaf({player_1: 1, player_2: 1})
    leaf_4 = Leaf({player_1: 3, player_2: 2})
    node_3 = Node({'Sports': leaf_1, 'Comedy': leaf_2}, 'c', player_1)
    node_2 = Node({'Sports': leaf_3, 'Comedy': leaf_4}, 'b', player_1)
    node_1 = Node({'Sports': node_3, 'Comedy': node_2}, 'a', player_2)
    battle_of_the_sexes = ExtensiveFormGame(node_1)
    battle_of_the_sexes.set_info_set([node_2, node_3])
    p = battle_of_the_sexes.plot(view_info_sets = True)
    sphinx_plot(p)
"""
from sage.graphs.all import Graph
from sage.plot.line import line2d
from sage.graphs.generic_graph import GenericGraph
from operator import attrgetter

class ExtensiveFormGame():
    r"""
    An object representing an Extensive Form Game. Primarily used to compute the
    Nash Equilibria.

    INPUT:

    - ``generator`` - Can be an instance of the Node class which serves as the
      root of the tree.
    """
    def __init__(self, generator, name=False, extensive_root=False):
        r"""
        Initializes an Extensive Form game and checks the inputs.

        EXAMPLES:

        A game with 2 players and 8 nodes::

            sage: player_a1 = Player('Player 1')
            sage: player_a2 = Player('Player 2')
            sage: leaf_a1 = Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a2 = Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a3 = Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a4 = Leaf({player_a1 : 2, player_a2: 1})
            sage: leaf_a5 = Leaf({player_a1 : 0, player_a2: 1})
            sage: leaf_a6 = Leaf({player_a1 : 1, player_a2: 0})
            sage: leaf_a7 = Leaf({player_a1 : 2, player_a2: 4})
            sage: leaf_a8 = Leaf({player_a1 : 2, player_a2: 1})
            sage: node_a1 = Node({'A': leaf_a1, 'B': leaf_a2}, player = player_a1)
            sage: node_a2 = Node({'A': leaf_a3, 'B': leaf_a4}, player = player_a1)
            sage: node_a3 = Node({'A': leaf_a5, 'B': leaf_a6}, player = player_a1)
            sage: node_a4 = Node({'A': leaf_a7, 'B': leaf_a8}, player = player_a1)
            sage: node_a5 = Node({'C': node_a1, 'D': node_a2}, player = player_a2)
            sage: node_a6 = Node({'C': node_a3, 'D': node_a4}, player = player_a2)
            sage: root_a = Node({'A': node_a5, 'B': node_a6}, player = player_a1)
            sage: egame_a1 = ExtensiveFormGame(root_a)
            sage: egame_a1.tree
            Graph on 15 vertices

        The generated tree has a variety of attributes::

            sage: egame_a1.players
            [Player 1, Player 2]
            sage: egame_a1.nodes
            [Node 1, Node 2, Node 3, Node 4, Node 5, Node 6, Tree Root]
            sage: egame_a1.info_sets
            [[Node 1], [Node 2], [Node 3], [Node 4], [Node 5], [Node 6], [Tree Root]]

        It is possible to have games with more than two players::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: player_3 = Player('Player 3')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1, player_3: -5}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0, player_3: -4}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4, player_3: -3}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1, player_3: -2}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_3)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.players
            [Player 1, Player 2, Player 3]
            sage: egame_1.nodes
            [Node 1, Node 2, Root 1]
            sage: egame_1.leafs
            [(0, 1, -5), (1, 0, -4), (2, 4, -3), (2, 1, -2)]
            sage: egame_1.tree
            Graph on 7 vertices

        If we do not name our nodes, unique names are
        automatically set during the initialisation::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_a = Leaf({player_1 : 0, player_2: 1})
            sage: leaf_b = Leaf({player_1 : 1, player_2: 0})
            sage: leaf_c = Leaf({player_1 : 2, player_2: 4})
            sage: leaf_d = Leaf({player_1 : 2, player_2: 1})
            sage: node_a = Node({'A': leaf_a, 'B': leaf_b}, player = player_2)
            sage: node_b = Node({'A': leaf_c, 'B': leaf_d}, player = player_2)
            sage: root_a = Node({'C': node_a, 'D': node_b}, player = player_1)
            sage: node_a.name is node_b.name is root_a.name
            True
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: node_a.name is node_b.name is root_a.name
            False
            sage: sorted([node_a.name, node_b.name])
            ['Node 1', 'Node 2']
            sage: root_a.name
            'Tree Root'

        If the input is a root, it needs children, actions and players::

            sage: false_root = Node({'C': node_1, 'D': node_2}, 'False Root')
            sage: egame_2 = ExtensiveFormGame(false_root)
            Traceback (most recent call last):
            ...
            AttributeError: Root node has no player.

            sage: false_root = Node(['Action1', 'Action2'])
            sage: false_root.player = player_1
            sage: egame_2 = ExtensiveFormGame(false_root)
            Traceback (most recent call last):
            ...
            AttributeError: Root node has no children.

            sage: false_root = Node([])
            sage: egame_2 = ExtensiveFormGame(false_root)
            Traceback (most recent call last):
            ...
            AttributeError: Root node has no actions.


        If we try to put an object that isn't a graph or a node in the game,
        we'll also return an error::

            sage: egame_2 = ExtensiveFormGame(player_1)
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an input in the form of a Node or a Graph object.

            sage: egame_2 = ExtensiveFormGame([node_1, node_2])
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an input in the form of a Node or a Graph object.

            sage: egame_2 = ExtensiveFormGame(leaf_1)
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an input in the form of a Node or a Graph object.


        Similarly, we cannot create a tree with a Node with a player attribute
        which isn't an instance of the ```Player``` class::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', leaf_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            Traceback (most recent call last):
            ...
            TypeError: Cannot assign a non Player as a player to a Node.

        """
        self.nodes = []

        if isinstance(generator, Node):
            if generator.actions is False:
                raise AttributeError("Root node has no actions.")
            elif generator.children is False:
                raise AttributeError("Root node has no children.")
            elif generator.player is False:
                raise AttributeError("Root node has no player.")
            else:
                self.tree_root = generator
                self.tree = self.grow_tree()
                self.nodes = self.grow_tree_dictionary().keys()
                self.nodes.sort(key=lambda x: x.actions[0])
                self.players = []
                self.info_sets = [[node] for node in self.nodes]
                self.leafs = []

                self._check_node_names_and_find_players(generator)

                self.players.sort(key=lambda x: x.name)
                self.info_sets.sort(key=lambda x: x[0].name)
                self.nodes.sort(key=lambda x: x.name)
                #self.leafs.sort(key=attrgetter('name', 'payoffs'))

        else:
            raise TypeError("Extensive form game must be passed an input in the form of a Node or a Graph object.")

    def set_info_set(self, node_list):
        """
        We can assign information set to  a set of nodes::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.info_sets
            [[Node 1], [Node 2], [Root 1]]
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.info_sets
            [[Node 1, Node 2], [Root 1]]

        Once we've set an info_set, we can see it visually on the graph::

            sage: egame_1.plot()
            Graphics object consisting of 20 graphics primitives
            sage: egame_1.plot(view_info_sets = True)
            Graphics object consisting of 23 graphics primitives

        If two nodes don't have the same actions, an error is returned::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, player = player_2)
            sage: node_2 = Node({'DifferentA': leaf_3, 'B': leaf_4}, player = player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, player = player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            Traceback (most recent call last):
            ...
            AttributeError: All nodes in the same information set must have the same actions.

        If two nodes have different players, an error is returned::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, player = player_1)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, player = player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, player = player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            Traceback (most recent call last):
            ...
            AttributeError: All nodes in the same information set must have the same players.
        """
        num_of_same_players = 1
        num_of_same_actions = 1
        previous_actions = []
        previous_player = True
        for node in node_list:
            if node.player == previous_player:
                num_of_same_players += 1
            if node.actions == previous_actions:
                num_of_same_actions += 1
            previous_player = node.player
            previous_actions = node.actions
        if num_of_same_players is not len(node_list):
            raise AttributeError("All nodes in the same information set must have the same players.")
        if num_of_same_actions is not len(node_list):
            raise AttributeError("All nodes in the same information set must have the same actions.")

        for node_to_be_set in node_list:
            for info_set in self.info_sets:
                for node_in_a_set in info_set:
                        if node_in_a_set is node_to_be_set and len(info_set) is not 1:
                            raise ValueError("Cannot assign information sets to nodes already in information sets")
            self.info_sets.remove([node_to_be_set])
        self.info_sets.append(sorted(node_list, key=lambda x: x.name))
        self.info_sets.sort(key=lambda x: x[0].name)

    def remove_info_set(self, node_list):
        """
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.info_sets
            [[Node 1, Node 2], [Root 1]]
            sage: egame_1.perfect_info()
            False
            sage: egame_1.remove_info_set([node_1, node_2])
            sage: egame_1.info_sets
            [[Node 1], [Node 2], [Root 1]]
            sage: egame_1.perfect_info()
            True
        """
        self.info_sets.remove(node_list)
        for node_to_be_readded in node_list:
            self.info_sets.append([node_to_be_readded])
        self.info_sets.sort(key=lambda x: x[0].name)

    def perfect_info(self):
        """
        All games start of with perfect information, adding information sets
        changes this::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.perfect_info()
            True
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.perfect_info()
            False
            sage: egame_1.remove_info_set([node_1, node_2])
            sage: egame_1.perfect_info()
            True
        """
        perfect_info_set = [[node] for node in self.nodes]
        return len(self.info_sets) == len(self.nodes)

    def grow_tree(self):
        """
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.grow_tree()
            Graph on 7 vertices
        """
        d = self.grow_tree_dictionary()
        t = Graph(d)
        if t.is_tree():
            return t
        else:
            raise TypeError("Graph isn't tree")

    def plot(self, view_info_sets=False):
        """
        Returns a visual representation of the game::
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.plot()
            Graphics object consisting of 20 graphics primitives

        The plot has the option for whether or not info-sets are visible::

            sage: egame_1.plot(view_info_sets = True)
            Graphics object consisting of 23 graphics primitives
        """

        keylist = []
        t = self.grow_tree()
        for node in self.nodes:
            keylist = node.node_input.keys()
            for key in keylist:
                t.set_edge_label(node, node.node_input[key], key)

        tree_plot = t.plot(layout='tree', tree_orientation='right',
                        edge_labels=True, tree_root=self.tree_root,
                        save_pos=True, axes=False)
        positions = t.get_pos()

        past_info_node = self.info_sets[0][0]
        if view_info_sets is True:
            for info_set in self.info_sets:
                past_info_node = info_set[0]
                for node in info_set:
                    for key in positions.keys():
                        if node is key:
                            tree_plot += (line2d([positions[past_info_node],
                                          positions[node]], linestyle="dashed",
                                          color='green'))
                            past_info_node = node
        return tree_plot

    def grow_tree_dictionary(self):
        """
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1', player_2)
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2', player_2)
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1', player_1)
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: t = Graph(egame_1.grow_tree_dictionary())
            sage: t
            Graph on 7 vertices
        """
        to_check = [self.tree_root]
        checked = []
        while to_check:
            checking = to_check.pop()
            for child in checking.children:
                if not isinstance(child, Leaf):
                    if child._is_complete():
                        child._player_check()
                        to_check.append(child)
                    else:
                        raise AttributeError("One or more of the Nodes in tree are not complete.")
                checked.append(child)

        d = {node:node.children for node in checked if not isinstance(node, Leaf)}
        d[self.tree_root] = self.tree_root.children
        return d

    def _check_node_names_and_find_players(self, generator):
        """
        A method to check the names of the nodes and gives names for the ones
        that do not have names. This also finds all the players.

        This method is embedded in the init method but has been written here to
        improve the readability of the code. The tests for this are functional
        tests.

        TESTS::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_a = Leaf({player_1 : 0, player_2: 1})
            sage: leaf_b = Leaf({player_1 : 1, player_2: 0})
            sage: leaf_c = Leaf({player_1 : 2, player_2: 4})
            sage: leaf_d = Leaf({player_1 : 2, player_2: 1})
            sage: node_a = Node({'A': leaf_a, 'B': leaf_b}, player = player_2)
            sage: node_b = Node({'A': leaf_c, 'B': leaf_d}, player = player_2)
            sage: root_a = Node({'C': node_a, 'D': node_b}, player = player_1)
            sage: node_a.name is node_b.name is root_a.name
            True
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: node_a.name is node_b.name is root_a.name
            False
            sage: sorted([node_a.name, node_b.name])
            ['Node 1', 'Node 2']
            sage: root_a.name
            'Tree Root'
            sage: sorted(egame_a.players, key=lambda x: x.name)
            [Player 1, Player 2]

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_a = Leaf({player_1 : 0, player_2: 1})
            sage: leaf_b = Leaf({player_1 : 1, player_2: 0})
            sage: leaf_c = Leaf({player_1 : 2, player_2: 4})
            sage: leaf_d = Leaf({player_1 : 2, player_2: 1})
            sage: node_a = Node({'A': leaf_a, 'B': leaf_b}, player = player_2)
            sage: node_b = Node({'A': leaf_c, 'B': leaf_d}, 'Node B', player = player_2)
            sage: root_a = Node({'C': node_a, 'D': node_b}, player = player_1)
            sage: node_a.name is root_a.name
            True
            sage: node_b.name
            'Node B'
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: node_a.name is node_b.name is root_a.name
            False
            sage: sorted([node_a.name, node_b.name])
            ['Node 1', 'Node B']
            sage: root_a.name
            'Tree Root'
        """
        node_index = 1
        leaf_index = 1
        for node in self.nodes:
            if node.player not in self.players:
                self.players.append(node.player)
            if node is generator and node.name is False:
                node.name = "Tree Root"
            if node.name is False:
                node.name = "Node %i" % node_index
                node_index += 1
            for child in node.children:
                if isinstance(child, Leaf) and child.name is False:
                    child.name = "Leaf %i" % leaf_index
                    leaf_index += 1
                    self.leafs.append(child)
                elif isinstance(child, Leaf):
                    self.leafs.append(child)


class Node():
    def __init__(self, node_input, name=False, player=False):
        """
        Node input will be in a dictionary format, consisting of the actions and
        the children of that node::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Mother')
            sage: mother_node.actions
            ['Action1', 'Action2']
            sage: mother_node.children
            [(0, 1), (1, 0)]

        If we then create a second node, who has :code:`mother_node` as one of
        its children, then the parent of :code:`mother_node` will be set to that
        node::

            sage: sisternode = Node(['inputhere'])
            sage: mother_node.parent
            False
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node.parent
            Node A

        Nodes can also be created without specifying children or parents by just
        passing the list of actions.  This so that nodes can be passed via a
        tree Sage Graph object to the extensive form game class::

            sage: grandmother_node = Node(['ActionA', 'ActionB'])
            sage: grandmother_node.children
            False
            sage: grandmother_node.parent
            False


        Nodes automatically have player set to false, we can then assign a
        player to that node::

            sage: grandmother_node.player
            False
            sage: grandmother_node.player = player_1
            sage: grandmother_node.player
            Player 1

        If we try to pass an node_input that isn't a dictionary or a list, an
        error is returned::

            sage: grandmother_node = Node(5)
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an node_input in the form of a dictionary or a list.

            sage: grandmother_node = Node('This is a string')
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an node_input in the form of a dictionary or a list.

            sage: grandmother_node = Node(matrix([[1, 1], [1, 1]]))
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an node_input in the form of a dictionary or a list.

            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node(sisternode)
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an node_input in the form of a dictionary or a list.
        """
        self.node_input = node_input
        self.player = player
        self.name = name
        self.actions = False
        self.children = False
        self.parent = False
        self.is_root = False

        if type(node_input) is dict:
            self.actions = node_input.keys()
            self.children = node_input.values()
            for child in self.children:
                child.parent = self

        elif type(node_input) is list:
            if node_input == []:
                self.actions = False
            else:
                self.actions = node_input

        else:
            raise TypeError("Node must be passed an node_input in the form of a dictionary or a list.")

    def __repr__(self):
        """
        Representation method for the Node::

            sage: repr_node = Node(['inputhere'])
            sage: repr_node
            False
            sage: repr_node.name = "A named Node"
            sage: repr_node
            A named Node
            sage: repr_node = Node(['inputhere'], "A different name")
            sage: repr_node
            A different name
        """
        if self.name is False:
            return "False"
        else:
            return self.name

    def _is_complete(self):
        """
        If we create a node where their children aren't specified and no parent
        is set, the node is considered incomplete::

            sage: b = Node(['Action1', 'Action2'])
            sage: b._is_complete == True
            False

        However, when we do specify all those attributes, the node is then
        considered complete::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B', player_1)
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode},
            ....:                         'Node A')
            sage: mother_node._is_complete()
            True

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode},
            ....:                         'Node A')
            sage: mother_node._is_complete()
            False

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2},
            ....:                    'Node B', player_1)
            sage: sisternode = Node(['inputhere'])
            sage: mother_node._is_complete()
            False

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2},
            ....:                    'Node B', player_1)
            sage: mother_node.children = False
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode},
            ....:                         'Node A')
            sage: mother_node._is_complete()
            False
        """
        return all([self.parent , self.actions, self.children, self.player])

    def _player_check(self):
        """
        A check primarily used later for creating of Extensive Form Games::
            sage: grandmother_node = Node([])
            sage: mother_node = Node([])
            sage: grandmother_node.player
            False
            sage: grandmother_node.player = mother_node
            sage: grandmother_node._player_check()
            Traceback (most recent call last):
            ...
            TypeError: Cannot assign a non Player as a player to a Node.
        """
        if self.player is not False:
                    if not isinstance(self.player, Player):
                        raise TypeError("Cannot assign a non Player as a player to a Node.")


class Leaf():
    def __init__(self, payoffs, name = False):
        """
        We can check payoffs of any leaf::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_1.payoffs[player_1]
            0
            sage: leaf_1[player_1]
            0
            sage: leaf_1[player_2]
            1
            sage: leaf_1.players
            [Player 1, Player 2]

        The payoffs must be in dictionary form such that the keys are players,
        and the values are either float or integers::

            sage: node_1 = Node(['input']); node_2 = Node(['input'])
            sage: leaf_1 = Leaf({node_1: 0, node_2: 1})
            Traceback (most recent call last):
            ...
            TypeError: The payoffs within Leaf must be in dictionary form with players as keys, and numbers as payoffss.

            sage: leaf_1 = Leaf([0, 1])
            Traceback (most recent call last):
            ...
            TypeError: The payoffs within Leaf must be in dictionary form with players as keys, and numbers as payoffss.
        """
        if type(payoffs) is not dict:
            raise TypeError("The payoffs within Leaf must be in dictionary form with players as keys, and numbers as payoffss.")

        self.payoffs =  payoffs
        self.name = name
        self.players = sorted(payoffs.keys(), key=lambda x:x.name)
        self.parent = False

        for player in self.players:
            if not isinstance(player, Player):
                raise TypeError("The payoffs within Leaf must be in dictionary form with players as keys, and numbers as payoffss.")

    def __repr__(self):
        """
        Representation method for the leaf::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1}, 'end_leaf')
            sage: leaf_1
            (0, 1)

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 5, player_2: 1}, 'end_leaf')
            sage: leaf_1
            (5, 1)

            sage: player_1 = Player('Vince')
            sage: player_2 = Player('Hannah')
            sage: player_3 = Player('James')
            sage: leaf_1 = Leaf({player_1: 5, player_2: 1, player_3:10}, 'end_leaf')
            sage: leaf_1
            (1, 10, 5)
        """
        return str(tuple([self[player] for player in
                        sorted(self.players, key=lambda x:x.name)]))


    def __delitem__(self, key):
        """
        This method is one of a collection that aims to make a leaf
        instance behave like a dictionary.

        Here we set up deleting an element of the payoffs dictionary::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_dict = Leaf({player_1: 0, player_2: 1})
            sage: del(leaf_dict[player_2])
            sage: leaf_dict.payoffs
            {Player 1: 0}
        """
        self.payoffs.pop(key, None)

    def __getitem__(self, key):
        """
        This method is one of a collection that aims to make a leaf
        instance behave like a dictionary.

        Here we allow for querying a key::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_dict = Leaf({player_1: 0, player_2: 1})
            sage: leaf_dict[player_1]
            0
            sage: del(leaf_dict[player_1])
            sage: leaf_dict[player_1]
            Traceback (most recent call last):
            ...
            KeyError: Player 1

        """
        return self.payoffs[key]

    def __iter__(self):
        """
        This method is one of a collection that aims to make a game
        instance behave like a dictionary.

        Here we allow for iteration over the leaf to correspond to
        iteration over keys of the utility dictionary::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_dict = Leaf({player_1: 0, player_2: 1})
            sage: for key in leaf_dict:
            ....:     print "The player: {}, has payoff {}".format(key, leaf_dict[key])
            The player: Player 2, has payoff 1
            The player: Player 1, has payoff 0
        """
        return iter(self.payoffs)

    def __setitem__(self, key, value):
        """
        This method is one of a collection that aims to make a game
        instance behave like a dictionary.

        Here we set up setting the value of a key::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_dict = Leaf({player_1: 0, player_2: 1})
            sage: leaf_dict[player_1]
            0
            sage: leaf_dict[player_1] = 2
            sage: leaf_dict[player_1]
            2

        """
        self.payoffs[key] = value


class Player():
    def __init__(self, name):
        """
        We can use Player() to assign players to nodes::

            sage: jack_1 = Node([0, 1])
            sage: jack_1.player = Player('Jack')
            sage: jack_1.player
            Jack

        If a node is not specificed a player, then this should return false::

            sage: sam_1 = Node([0, 1])
            sage: sam_1.player
            False
            sage: sam_player = Player('Sam')
            sage: sam_1.player = sam_player
            sage: sam_1.player
            Sam
            sage: sam_2 = Node([0, 1], player = sam_player)
            sage: sam_2.player
            Sam

        A player can be reassigned for Nodes::
            sage: andy_player = Player('Andy')
            sage: sam_2.player = andy_player
            sage: sam_2.player
            Andy

        We can create players and assign them names::

            sage: ben_player = Player('Benjamin')
            sage: ben_player.name
            'Benjamin'
        """
        self.name = name

    def __repr__(self):
        """
        Representation method for the player::

            sage: apple_1 = Player('Apple')
            sage: apple_1
            Apple
        """
        if self.name is False:
            return "False"
        else:
            return self.name

    def __hash__(self):
        """
        Makes the player class hashable::

            sage: apple_1 = Player('Apple')
            sage: banana_1 = Leaf({apple_1 : 0})
            sage: banana_1.payoffs
            {Apple: 0}
        """
        return hash(self.name)
