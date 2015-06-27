"""
This module implements a class for Extensive Form Games.

A well known example that can be implememnted as an Extensive Form Game is the battle of  the sexes.
Consider two players, Celine and Bob. The two are deciding on a movie, they can either see a Sports or a Comedy.
Bob would prefer to see a Comedy, Celine would prefer to watch a Sports movie.
Depending on who choses what, there are different payoffs, this can be demonstrated in tree Form.

[some form of diagram here]

We can see there are three nodes, one for Bob, two for Celine. Connecting those nodes are actions.
These actions represent choices made by one player, and the actions then lead on to a node of another player.
So Bob either choses Sports or Comedy, and then Celine choses Sports or Comedy.
However we can see that the payoffs differ depending on what Bob chose first.
The location where a payoff is located is called a Leaf, and  they show the outcome for each player.
So if Bob choses Sports, and Celine choses Sports, we see the payoff is (2, 3), which represents Bob getting a payoff of 2,
and Celine getting a payoff of 3.

We can create the game in sage two ways, one by passing it a root, and the other by passing it a tree::

    sage: player_1 = Player('Bob')
    sage: player_2 = Player('Celine')
    sage: leaf_1 = Leaf({player_1: 2, player_2: 3})
    sage: leaf_2 = Leaf({player_1: 0, player_2: 0})
    sage: leaf_3 = Leaf({player_1: 1, player_2: 1})
    sage: leaf_4 = Leaf({player_1: 3, player_2: 2})
    sage: node_3 = Node({'Sports': leaf_1, 'Comedy': leaf_2}, 'c', player_1)
    sage: node_2 = Node({'Sports': leaf_3, 'Comedy': leaf_4}, 'b', player_1)
    sage: node_1 = Node({'Sports': node_3, 'Comedy': node_2}, 'a', player_2)
    sage: battle_of_the_sexes = ExtensiveFormGame(node_1)

This can then be shown visually in graph form, as is shown here::

    sage: battle_of_the_sexes.plot()
    [tree here]

In Extensive Form Games we can create information sets, when creating a game all the nodes are in their own individual information sets.
This is called perfect information, and we can check this using sage::

    sage: battle_of_the_sexes.info_sets
    [a, b, c]
    sage: battle_of_the_sexes.perfect_info()
    True

If we wanted to make it so Celine does not know the actions Bob has taken, we can put Celine's nodes in one information set::

    sage: battle_of_the_sexes.set_info_set[node_2, node_3]
    sage: battle_of_the_sexes.info_sets
    [a, [b, c]]

Now the game does not have perfect information::

    sage: battle_of_the_sexes.perfect_info()
    False

Information sets are demonstrated visually on the graph we plot (hopefully)::

    sage: battle_of_the_sexes.plot()
    [tree here]
"""
from sage.graphs.all     import *
from sage.structure.sage_object import SageObject
from sage.plot.graphics      import show_default, Graphics
from sage.rings.all     import *
import networkx
from sage.graphs.generic_graph import GenericGraph


class ExtensiveFormGame():
    def __init__(self, game_input , name = False, extensive_root = False):
        """
        Just a really big test to see if anything breaks::

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
            sage: node_a1 = Node({'A': leaf_a1, 'B': leaf_a2})
            sage: node_a2 = Node({'A': leaf_a3, 'B': leaf_a4})
            sage: node_a3 = Node({'A': leaf_a5, 'B': leaf_a6})
            sage: node_a4 = Node({'A': leaf_a7, 'B': leaf_a8})
            sage: node_a5 = Node({'C': node_a1, 'D': node_a2})
            sage: node_a6 = Node({'C': node_a3, 'D': node_a4})
            sage: node_a1.player = player_a1
            sage: node_a2.player = player_a1
            sage: node_a3.player = player_a1
            sage: node_a4.player = player_a1
            sage: node_a5.player = player_a2
            sage: node_a6.player = player_a2
            sage: root_a = Node({'A': node_a5, 'B': node_a6})
            sage: root_a.player = player_a1
            sage: egame_a1 = ExtensiveFormGame(root_a)
            sage: egame_a1.players
            sage: egame_a1.tree
            Graph on 15 vertices
            sage: egame_a1.plot()
            Graphics object consisting of 44 graphics primitives
            sage: egame_a1.info_sets
            [[Node 1], [Node 2], [Node 3], [Node 4], [Node 5], [Node 6], [Tree Root]]
            sage: egame_a1.set_info_set([node_a3, node_a5])

        ExtensiveFormGame can take input in the form of a Node::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.players
            [Player 2, Player 2, Player 1]
            sage: egame_1.tree
            Graph on 7 vertices
            sage: egame_1.plot()
            Graphics object consisting of 20 graphics primitives

        Temporary test for Node naming::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_a = Leaf({player_1 : 0, player_2: 1})
            sage: leaf_b = Leaf({player_1 : 1, player_2: 0})
            sage: leaf_c = Leaf({player_1 : 2, player_2: 4})
            sage: leaf_d = Leaf({player_1 : 2, player_2: 1})
            sage: node_a = Node({'A': leaf_a, 'B': leaf_b})
            sage: node_b = Node({'A': leaf_c, 'B': leaf_d})
            sage: node_a.player = player_2
            sage: node_b.player = player_2
            sage: root_a = Node({'C': node_a, 'D': node_b})
            sage: root_a.player = player_1
            sage: node_b.name
            sage: egame_a = ExtensiveFormGame(root_a)
            sage: egame_a.plot()
            Graphics object consisting of 20 graphics primitives


        If the game_input is a root, it needs children, actions and players::

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


        If we try to put an object that isn't a graph or a node in the game, we'll also return an error::

            sage: egame_2 = ExtensiveFormGame(player_1)
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an game_input in the form of a Node or a Graph object.

            sage: egame_2 = ExtensiveFormGame([node_1, node_2])
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an game_input in the form of a Node or a Graph object.

            sage: egame_2 = ExtensiveFormGame(leaf_1)
            Traceback (most recent call last):
            ...
            TypeError: Extensive form game must be passed an game_input in the form of a Node or a Graph object.



        Similarly, we cannot create a tree with a Node with a player attribute which isn't
        an instance of the ```Player``` class::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = leaf_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            Traceback (most recent call last):
            ...
            TypeError: Cannot assign a non Player as a player to a Node.

        """
        self.nodes = []
        self.check = []

        if isinstance(game_input, Node):
            if game_input.actions is False:
                raise AttributeError("Root node has no actions.")
            elif game_input.children is False:
                raise AttributeError("Root node has no children.")
            elif game_input.player is False:
                raise AttributeError("Root node has no player.")
            else:
                self.tree_root = game_input
                self.tree = self.grow_tree()
                self.nodes = sorted((self.grow_tree_dictionary()).keys())
                self.players = []
                self.info_sets = sorted((self.grow_tree_dictionary()).keys())
                for i in self.nodes:
                    self.players.append(i.player)
        else:
            raise TypeError("Extensive form game must be passed an game_input in the form of a Node or a Graph object.")

    def set_info_set(self, nodelist):
        """
        We can assign information set to  a set of nodes::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.info_sets
            [[Node 1], [Node 2], [Root 1]]
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.info_sets

        If two nodes don't have the same actions, an error is returned::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2})
            sage: node_2 = Node({'DifferentA': leaf_3, 'B': leaf_4})
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2})
            sage: root_1.player = player_1
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
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2})
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4})
            sage: node_1.player = player_1
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2})
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            Traceback (most recent call last):
            ...
            AttributeError: All nodes in the same information set must have the same players.
        """
        j = 1
        previousplayer = True
        for i in nodelist:
            if i.player == previousplayer:
                j += 1
            previousplayer = i.player
        if j is not len(nodelist):
            raise AttributeError("All nodes in the same information set must have the same players.")

        j = 1
        previousactions = []
        for i in nodelist:
            if i.actions == previousactions:
                j += 1
            previousactions = i.actions

        if j is not len(nodelist):
            raise AttributeError("All nodes in the same information set must have the same actions.")


        for i in self.info_sets:
            for j in nodelist:
                if type(i) is list:
                    for k in i:
                        if k is j:
                            raise ValueError("Cannot assign information sets to nodes already in information sets")
                self.info_sets.remove(j)
        self.info_sets.append(nodelist)

    def perfect_info(self):
        """
        All games start of with perfect information, adding information sets changes this::
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.perfect_info()
            True
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.perfect_info()
            False
        """
        if self.info_sets == self.nodes:
            return True
        else:
            return False

    def remove_info_set(self, nodelist):
        """
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1 : 0, player_2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player_1 : 1, player_2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player_1 : 2, player_2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player_1 : 2, player_2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.info_sets
            sage: egame_1.perfect_info()
            False
            sage: egame_1.remove_info_set([node_1, node_2])
            sage: egame_1.info_sets
            sage: egame_1.perfect_info()
            True
        """
        self.info_sets.remove(nodelist)
        for j in nodelist:
            self.info_sets.append(j)



    def grow_tree(self):
        """
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
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

    def plot(self):
        keylist = []
        t = self.grow_tree()
        for i in self.nodes:
            keylist = i.node_input.keys()
            for j in keylist:
                t.set_edge_label(i, i.node_input[j], j)
        return t.plot(layout='tree',tree_orientation='right', edge_labels=True, tree_root = self.tree_root)

    def grow_tree_dictionary(self):
        """
            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_2 = Leaf({player_1: 1, player_2: 0})
            sage: leaf_3 = Leaf({player_1: 2, player_2: 4})
            sage: leaf_4 = Leaf({player_1: 2, player_2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: node_1.player = player_2
            sage: node_2.player = player_2
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: root_1.player = player_1
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: t = Graph(egame_1.grow_tree_dictionary())
            sage: t
            Graph on 7 vertices
        """
        to_check = [self.tree_root]  # Put the one node we have in a list of things we need to check
        checked = []  # A list of nodes that we have checked
        while to_check:  # A while loop to keep going until there is nothing left in the `to_check` list
            checking = to_check.pop()  # The pop command returns the node AND removes it from the list
            for child in checking.children:  # Loop through that node's children
                if not isinstance(child, Leaf):  # If it's not a Leaf ...
                    if child._is_complete():
                        child._player_check()
                        to_check.append(child)  # ... append it the the list of nodes we need to check
                    else:
                        raise AttributeError("One or more of the Nodes in tree are not complete.")
                checked.append(child)  # Put the child in the list of checked nodes

        nodevalue = 1
        for i in self.nodes:
            self.check.append(i.name)
            if i is self.tree_root and i.name is False:
                i.name = "Tree Root"
            if i.name is False:
                i.name = "Node %i" %nodevalue
                nodevalue += 1

        leafvalue = 1
        for i in checked:
            if isinstance(i, Leaf) and i.name is False:
                i.name = "Leaf %i" %leafvalue
                leafvalue += 1

        # Create the dictionary
        d = {node:node.children for node in checked if not isinstance(node, Leaf)}  # Build the dictionary mapping the leafs to their children
        # The above does not include the root
        d[self.tree_root] = self.tree_root.children  # The above does not actually include the original root so we need to include it
        return d


class Node():
    def __init__(self, node_input, name=False, player=False):
        """
        Node input will be in a dictionary format, consisting of the actions and the children of that node::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Mother')
            sage: mother_node.actions
            ['Action1', 'Action2']
            sage: mother_node.children
            [Player 1: 0 Player 2: 1 , Player 1: 1 Player 2: 0 ]

        If we then create a second node, who has :code:`mother_node` as one of its children,
        then the parent of :code:`mother_node` will be set to that node::

            sage: sisternode = Node(['inputhere'])
            sage: mother_node.parent
            False
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node.parent
            Node A

        Nodes can also be created without specifying children or parents by just passing the list of actions.
        This so that nodes can be passed via a tree Sage Graph object to the extensive form game class::

            sage: grandmother_node = Node(['ActionA', 'ActionB'])
            sage: grandmother_node.children
            False
            sage: grandmother_node.parent
            False


        Nodes automatically have player set to false, we can then assign a player to that node::

            sage: grandmother_node.player
            False
            sage: grandmother_node.player = player_1
            sage: grandmother_node.player
            Player 1

        If we try to pass an node_input that isn't a dictionary or a list, an error is returned::

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
        if self.name is False:
            return "False"
        else:
            return self.name

    def attributes(self):
        """
        We can use this function to check the attributes of each singular node, the following is what would happen if no attibutes are assigned::

            sage: laura_1 = Node(['inputhere'])
            sage: laura_1.attributes()
            "The Node has the following attributes. Actions: ['inputhere']. Children: False. Parent: False. Player: False."
        """
        return "The Node has the following attributes. Actions: %s. Children: %s. Parent: %s. Player: %s." %(self.actions, self.children, self.parent, self.player)

    def _is_complete(self):
        """
        If we create a node where their children aren't specified and no parent is set, the node is considered incomplete::

            sage: b = Node(['Action1', 'Action2'])
            sage: b._is_complete == True
            False

        However, when we do specify all those attributes, the node is then considered complete::

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: mother_node.player = player_1
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node._is_complete()
            True

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node._is_complete()
            False

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: mother_node.player = player_1
            sage: sisternode = Node(['inputhere'])
            sage: mother_node._is_complete()
            False

            sage: player_1 = Player('Player 1')
            sage: player_2 = Player('Player 2')
            sage: child_1 = Leaf({player_1: 0, player_2: 1}, 'Child 1')
            sage: child_2 = Leaf({player_1: 1, player_2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: mother_node.player = player_1
            sage: mother_node.children = False
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node._is_complete()
            False
        """
        return all([self.parent , self.actions, self.children, self.player])

    def _player_check(self):
        """
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
            sage: leaf_1.payoffs[player_2]
            1
            sage: leaf_1.players
            [Player 2, Player 1]

        The payoffs must be in dictionary form such that the keys are players, and the values are either float or intergers::

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
        self.players = sorted(payoffs.keys())
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

            sage: d = {'Missty':2, 'Auray':5, 'Toby':4}
            sage: d == {'Auray':5, 'Missty':2, 'Toby':4}
            True
        """
        return str(tuple([self[plry] for plry in sorted(self.players, key=lambda x:x.name)]))


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

        We can create players and assign them names::

            sage: ben_player = Player('Benjamin')
            sage: ben_player.name
            'Benjamin'
        """
        self.name = name

    def __repr__(self):
        """
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
        TESTS::

            sage: apple_1 = Player('Apple')
            sage: banana_1 = Leaf({apple_1 : 0})
            sage: banana_1.payoffs
            {Apple: 0}
        """
        return hash(self.name)
