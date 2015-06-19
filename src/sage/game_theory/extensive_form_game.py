"""
This file will contain a class for extensive form games.
"""

class ExtensiveFormGame():
    """

    Class to be able to build a tree. To do this we first need to build nodes::

        sage: leaf_1 = Leaf((2, 3))
        sage: leaf_2 = Leaf((0, 0))
        sage: leaf_3 = Leaf((1, 1))
        sage: leaf_4 = Leaf((3, 2))
        sage: celine_1 = Node({'Sports': leaf_1, 'Comedy': leaf_2})
        sage: celine_2 = Node({'Sports': leaf_3, 'Comedy': leaf_4})
        sage: root = Root({'Sports': celine_1, 'Comedy': celine_2})


    We then want to assign the nodes to players::

        sage: root.player = Player('Vince')
        sage: celine_1.player = Player('Celine')
        sage: celine_2.player = Player('Celine')

    We see that :code:`celine_1` has the following attributes::

        sage: celine_1.actions
        ['Sports, 'Comedy']
        sage: celine_1.children
        [leaf_1, leaf_2]
        sage: celine_1.player
        Celine

    We can then create the game::

        sage: g1 = ExtensiveFormGame(root)
        sage: g1
        An extensive form game...

    This creates an underlying tree object as well as a dictionary (or something) that keeps track of all the information::

        sage: g.tree
        The game tree...  # This is just a sage tree object

    After creating the game, :code:`celine_1` has the following extra attributes:

        sage: celine_1.actions
        ['Sports, 'Comedy']
        sage: celine_1.children
        [leaf_1, leaf_2]
        sage: celine_1.parent
        Vince
        sage: celine_1.player
        True
        sage: celine_1.is_complete
        True

    If we try to create a tree from a single node that is not a root we get an error::

        sage: g1 = ExtensiveFormGame(celine_1)
        ERROR

    We can also create a game without specifying the children/parents of each node but by passing a tree with all the structural information.


        sage: leaf_1 = Leaf((2, 3))
        sage: leaf_2 = Leaf((0, 0))
        sage: leaf_3 = Leaf((1, 1))
        sage: leaf_4 = Leaf((3, 2))
        sage: celine_1 = Node(['Sports', 'Comedy'])
        sage: celine_2 = Node(['Sports', 'Comedy'])
        sage: root = Root(['Sports', 'Comedy'])
        sage: tree = Graph({root:[celine_1, celine_2], celine_1:[leaf_1, leaf_2], celine_2:[leaf_3, leaf_4]})

    We need to specify the players::

        sage: root.player = Player('Vince')
        sage: celine_1.player = Player('Celine')
        sage: celine_2.player = Player('Celine')

    Our nodes now have the following information::

        sage: celine_1.actions
        ['Sports, 'Comedy']
        sage: celine_1.player
        True

    The tree now contains all the structural information and the nodes contain the other information::

        sage: g2 = ExtensiveFormGame(tree)
        sage: g2
        An extensive form game...

    After creating the game the nodes have all the information::

        sage: celine_1.actions
        ['Sports, 'Comedy']
        sage: celine_1.children
        [leaf_1, leaf_2]
        sage: celine_1.parent
        Vince
        sage: celine_1.player
        True
        sage: celine_1.is_complete
        True

    We see that these two trees although created differently are equal::

        sage: g1 == g2
        True

    If we try to create a tree from a root but the root was generated with no specified children, we also get an error:

        sage: root = Root(['Sports', 'Comedy'])
        sage: g1 = ExtensiveFormGame(root)
        ERROR

    If we create


    """


class Node(something):
    def __init__(self, input):
        """

        """
        self.player= False
        
    def attributes():
        """
        We can use this function to check the attributes of each singular node, the following is what would happen if no attibutes are assigned::
            sage: laura_1 = Node()
            sage: laura_1.attributes 
            The node has the following attributes. Actions: False. Children: False. Parent: False. Player: False. 
        """

class Root(Node):
    """
    A root is just another type of node, so we can get attributes, however Parent will always be false. Attempting to add a parent will return an error::
        sage: jess_1 = Node([Green, Yellow]); jess_2 = Node([Green, Yellow]); jess_3 = Node([Green, Yellow])
        sage: bethan_1 = Root({'Red': jess_1, 'Blue': jess_2}))
        sage: bonnie_1 = Root({'Black': jess_1, 'White': jess_3})
        ERROR

    We cannot have more than one Root in a game, so if we try to connect a second Root to a connected set of nodes that already have a root, an error will be displayed::
        sage: 
    """

class Player():
    def __init__(self, name):
        """
        We can use Player() to assign players to nodes::
        sage: jack_1 = Node()
        sage: jack_1.player = Player('Jack')
        sage: jack_1.player
        Jack

        If a node is not specificed a player, then this should return false::
        sage: sam_1 = Node()
        sage: sam_1.player
        False
        """

