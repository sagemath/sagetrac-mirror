"""
This file will contain a class for extensive form games.
"""
from sage.graphs.all     import *
from sage.structure.sage_object import SageObject
from sage.plot.graphics      import show_default, Graphics
import networkx
from sage.graphs.generic_graph import GenericGraph
class ExtensiveFormGame():
    def __init__(self, argument , name = False):
        """
        ExtensiveFormGame can take input in the form of a Node::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: leaf_1 = Leaf({player1 : 0, player2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player1 : 1, player2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player1 : 2, player2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player1 : 2, player2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.tree
            sage: type(egame_1.tree)
            Some tree here
            sage: egame_1.nodes
            [An extensive form game node - Root 1, An extensive form game node - Node 2, An extensive form game node - Node 1]
            sage: egame_1.players
            [False, False, False]

        Or it can take a Graph object, so long as it is a tree::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: leaf_1 = Leaf({player1 : 0, player2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player1 : 1, player2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player1 : 2, player2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player1 : 2, player2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: tree_1 = Graph({root_1:[node_1, node_2], node_1:[leaf_1, leaf_2], node_2:[leaf_3, leaf_4]})
            sage: egame_2 = ExtensiveFormGame(tree_1)
            sage: egame_2.nodes
            [An extensive form game node - Root 1, An extensive form game node - Node 2, An extensive form game node - Node 1]
            sage: egame_2.players
            [False, False, False]


        In the above examples, the two games created should be equal::

            sage: egame_1 == egame_2
            True

        If we input either a Graph that is not a tree, an error is returned::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: leaf_1 = Leaf({player1 : 0, player2: 1}, 'Leaf 1')
            sage: leaf_2 = Leaf({player1 : 1, player2: 0}, 'Leaf 2')
            sage: leaf_3 = Leaf({player1 : 2, player2: 4}, 'Leaf 3')
            sage: leaf_4 = Leaf({player1 : 2, player2: 1}, 'Leaf 4')
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: root_1 = Node({'C': node_1, 'D': node_2}, 'Root 1')
            sage: tree_2 = Graph({root_1: [node_1, node_2], node_1:[leaf_1, node_2], node_2:[leaf_3, leaf_4]})
            sage: egame_2 = ExtensiveFormGame(tree_2)
            Traceback (most recent call last):
            ...
            TypeError: Graph inputted is not a tree.

        If we try to put an empty tree in as a game, we'll also return an error::

            sage: tree_3 = Graph()
            sage: egame_2 = ExtensiveFormGame(tree_3)
            Traceback (most recent call last):
            ...
            ValueError: Graph inputted is empty.
        """
        self.nodes = []
        self.info_sets = []
        if isinstance(argument, Node):
            self.tree_root = argument
            self.tree = self.grow_tree()
            self.nodes = (self.grow_tree_dictionary()).keys()
            self.players = []
            for i in self.nodes:
                self.players.append(i.player)
        self.keycheck = {}
        if isinstance(argument, Graph):
            if argument.is_tree is False:
                raise TypeError("Graph inputted is not a tree.")
            dictionary = argument.to_dictionary()
            if dictionary == {}:
                raise ValueError("Graph inputted is empty.")
            nodelist = dictionary.keys()
            for i in nodelist:
                if isinstance(i, Node):
                    for child in i.children:
                        if isinstance(child, Node) or isinstance(child, Leaf):
                            i.children = dictionary[i]
                        else:
                            raise TypeError("Nodes in tree must be connected to either leaves or other nodes.")
                else:
                    pass
                    #raise TypeError("Tree keys must be Nodes.")
            self.keycheck = dictionary

    def set_info_set(self, nodelist):
        """
        We can assign information set to  a set of nodes::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: leaf_1 = Leaf({player1: 0, player2: 1})
            sage: leaf_2 = Leaf({player1: 1, player2: 0})
            sage: leaf_3 = Leaf({player1: 2, player2: 4})
            sage: leaf_4 = Leaf({player1: 2, player2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2}, 'Node 1')
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4}, 'Node 2')
            sage: root_1 = Node({'C': node_1, 'D': node_2})
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            sage: egame_1.info_sets
            [[An extensive form game node - Node 1, An extensive form game node - Node 2]]

        If two nodes don't have the same actions, an error is returned::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: leaf_1 = Leaf({player1: 0, player2: 1})
            sage: leaf_2 = Leaf({player1: 1, player2: 0})
            sage: leaf_3 = Leaf({player1: 2, player2: 4})
            sage: leaf_4 = Leaf({player1: 2, player2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2})
            sage: node_2 = Node({'DifferentA': leaf_3, 'B': leaf_4})
            sage: node_1.player = player1
            sage: node_2.player = player1
            sage: root_1 = Node({'C': node_1, 'D': node_2})
            sage: egame_1 = ExtensiveFormGame(root_1)
            sage: egame_1.set_info_set([node_1, node_2])
            Traceback (most recent call last):
            ...
            AttributeError: All nodes in the same information set must have the same actions.

        If two nodes have different players, an error is returned::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: leaf_1 = Leaf({player1: 0, player2: 1})
            sage: leaf_2 = Leaf({player1: 1, player2: 0})
            sage: leaf_3 = Leaf({player1: 2, player2: 4})
            sage: leaf_4 = Leaf({player1: 2, player2: 1})
            sage: node_1 = Node({'A': leaf_1, 'B': leaf_2})
            sage: node_2 = Node({'A': leaf_3, 'B': leaf_4})
            sage: node_1.player = player1
            sage: node_2.player = player2
            sage: root_1 = Node({'C': node_1, 'D': node_2})
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

        self.info_sets.append(nodelist)

    def grow_tree(self):
        d = self.grow_tree_dictionary()
        t = Graph(d)
        if t.is_tree():
            return t
        else:
           raise TypeError("Graph isn't tree")

    def plot_tree(self):
        keylist = []
        t = self.grow_tree()
        for i in self.nodes:
            keylist = i.argument.keys()
            for j in keylist:
                t.set_edge_label(i, i.argument[j], j)
        t.show(layout='tree',tree_orientation='right', edge_labels=True, tree_root = self.tree_root)

    def grow_tree_dictionary(self):
        to_check = [self.tree_root]  # Put the one node we have in a list of things we need to check
        checked = []  # A list of nodes that we have checked
        while to_check:  # A while loop to keep going until there is nothing left in the `to_check` list
            checking = to_check.pop()  # The pop command returns the node AND removes it from the list
            for child in checking.children:  # Loop through that node's children
                if not isinstance(child, Leaf):  # If it's not a Leaf ...
                    to_check.append(child)  # ... append it the the list of nodes we need to check
                checked.append(child)  # Put the child in the list of checked nodes

        # Create the dictionary
        d = {node:node.children for node in checked if not isinstance(node, Leaf)}  # Build the dictionary mapping the leafs to their children
        # The above does not include the root
        d[self.tree_root] = self.tree_root.children  # The above does not actually include the original root so we need to include it
        return d


class Node():
    def __init__(self, argument, name = False, player = False, is_root = False):
        """
        Node input will be in a dictionary format, consisting of the actions and the children of that node::

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: child_1 = Leaf({player1: 0, player2: 1}, 'Child 1')
            sage: child_2 = Leaf({player1: 1, player2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Mother')
            sage: mother_node.actions
            ['Action1', 'Action2']
            sage: mother_node.children
            [An extensive form game leaf - Child 1,
             An extensive form game leaf - Child 2]

        If we then create a second node, who has :code:`mother_node` as one of its children,
        then the parent of :code:`mother_node` will be set to that node::

            sage: sisternode = Node(['inputhere'])
            sage: mother_node.parent
            False
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node.parent
            An extensive form game node - Node A

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
            sage: grandmother_node.player = player1
            sage: grandmother_node.player
            Player 1

        If we try to pass an argument that isn't a dictionary or a list, an error is returned::

            sage: grandmother_node = Node(5)
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an argument in the form of a dictionary or a list.

            sage: grandmother_node = Node('This is a string')
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an argument in the form of a dictionary or a list.

            sage: grandmother_node = Node(matrix([[1, 1], [1, 1]]))
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an argument in the form of a dictionary or a list.

            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node(sisternode)
            Traceback (most recent call last):
            ...
            TypeError: Node must be passed an argument in the form of a dictionary or a list.
        """
        self.argument = argument
        self.player = player
        self.name = name
        self.actions = False
        self.children = False
        self.parent = False
        self.is_root = False

        if type(argument) is dict:
            self.actions = argument.keys()
            self.children = argument.values()
            for child in self.children:
                child.parent = self

        elif type(argument) is list:
            self.actions = argument

        else:
            raise TypeError("Node must be passed an argument in the form of a dictionary or a list.")

    def __repr__(self):

        s = 'An extensive form game node'
        if self.name:
            s += ' - ' + str(self.name)
        return s

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

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: child_1 = Leaf({player1: 0, player2: 1}, 'Child 1')
            sage: child_2 = Leaf({player1: 1, player2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: mother_node.player = player1
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node._is_complete()
            True

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: child_1 = Leaf({player1: 0, player2: 1}, 'Child 1')
            sage: child_2 = Leaf({player1: 1, player2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node._is_complete()
            False

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: child_1 = Leaf({player1: 0, player2: 1}, 'Child 1')
            sage: child_2 = Leaf({player1: 1, player2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: mother_node.player = player1
            sage: sisternode = Node(['inputhere'])
            sage: mother_node._is_complete()
            False

            sage: player1 = Player('Player 1')
            sage: player2 = Player('Player 2')
            sage: child_1 = Leaf({player1: 0, player2: 1}, 'Child 1')
            sage: child_2 = Leaf({player1: 1, player2: 0}, 'Child 2')
            sage: mother_node = Node({'Action1': child_1, 'Action2': child_2}, 'Node B')
            sage: mother_node.player = player1
            sage: mother_node.children = False
            sage: sisternode = Node(['inputhere'])
            sage: grandmother_node = Node({'ActionA':mother_node, 'ActionB':sisternode}, 'Node A')
            sage: mother_node._is_complete()
            False
        """
        return all([self.parent , self.actions, self.children, self.player])


class Leaf():
    def __init__(self, argument, name = False):
        """
        We can check payoffs of any leaf::

            sage: player_1 = Player('player 1')
            sage: player_2 = Player('player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1})
            sage: leaf_1.payoffs[player_1]
            0
            sage: leaf_1.payoffs[player_2]
            1
            sage: leaf_1.players
            [player 2, player 1]

        The payoffs must be in dictionary form such that the keys are players, and the values are either float or intergers::

            sage: node_1 = Node(['input']); node_2 = Node(['input'])
            sage: leaf_1 = Leaf({node_1: 0, node_2: 1})
            Traceback (most recent call last):
            ...
            TypeError: The payoffs within Leaf must be in dictionary form with players as keys, and numbers as arguments.

            sage: leaf_1 = Leaf([0, 1])
            Traceback (most recent call last):
            ...
            TypeError: The payoffs within Leaf must be in dictionary form with players as keys, and numbers as arguments.
        """

        if type(argument) is not dict:
            raise TypeError("The payoffs within Leaf must be in dictionary form with players as keys, and numbers as arguments.")

        self.argument = argument
        self.payoffs = argument
        self.name = name
        self.players = argument.keys()
        self.parent = False

        for player in self.players:
            if not isinstance(player, Player):
                raise TypeError("The payoffs within Leaf must be in dictionary form with players as keys, and numbers as arguments.")

    def __repr__(self):
        """
        Representation method for the leaf::

            sage: player_1 = Player('player 1')
            sage: player_2 = Player('player 2')
            sage: leaf_1 = Leaf({player_1: 0, player_2: 1}, 'end_leaf')
            sage: leaf_1
            An extensive form game leaf - end_leaf
        """

        s = 'An extensive form game leaf'
        if self.name:
            s += ' - ' + str(self.name)
        return s


class Player():
    def __init__(self, name):
        """
        We can use Player() to assign players to nodes::

            sage: jack_1 = Node([0, 1])
            sage: jack_1.player = (Player('Jack'))
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