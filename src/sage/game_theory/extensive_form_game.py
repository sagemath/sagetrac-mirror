"""
This file will contain a class for extensive form games.
"""

class ExtensiveFormGame():
    """
    """


class Node():
    def __init__(self, input, name = 'False'):
        """
        Node input can be read to determine actions of node and children of node.
        
            sage: child_1 = Leaf((0, 1))
            sage: child_2 = Leaf((1, 0))
            sage: mothernode = Node({'Action1': child_1, 'Action2'L child_2})
            sage: mothernode.actions
            ['Action1', 'Action2']
            sage: mothernode.children
            [child_1, child_2]

        If we then create a second node, who has :code:`mothernode` as one of it's children, then the parent of :code:`mothernode` will be set to that node::
        
            sage: sisternode = Node()
            sage: grandmothernode = Node({'ActionA':mothernode, 'ActionB':sisternode})
            sage: mothernode.parent
            grandmothernode


        We can also set names for the Node, which will then be returned instead::
        
            sage: grandmothernode = Node({'ActionA':mothernode, 'ActionB':sisternode}, 'Node A')
            sage: mothernode.parent
            Node A

        Node input can only be in dictionary form.

            sage: falsenode = Node([Fight, Flight])
            Traceback (most recent call last):
            ...
            TypeError: Node can only have its input in dictionary form.
        """
        #self.player = False
        
    
    def attributes():
        """
        We can use this function to check the attributes of each singular node, the following is what would happen if no attibutes are assigned::
           
            sage: laura_1 = Node()
            sage: laura_1.attributes 
            The node has the following attributes. Actions: False. Children: False. Parent: False. Player: False. 
        """


    def assigninfoset(type):
        """
        If a node has no parents, and a root node hasn't been set, a node can become a root node::
        
            sage: andy_1 = Node()
            sage: andy_1.change_type(Root)


        If the node has parents, an error message will be returned::
        
            sage: andy_1 = Node()
            sage: andy_2 = Node()
            sage: dave_1 = Node({'A': andy_1, 'B': andy_2})
            sage: andy_1.change_type(Root)
            Traceback (most recent call last):
            ...
            AttributeError: Node with parents cannot be a root.


        If the node is connected to a set of nodes that have a root associated with them, an error is returned::
           
            sage: jill_1 = node({'A'; helen_1 'B'; helen_2})
            sage: helen_1 = Root()
            sage: jill_1 = change_type(Root)
            Traceback (most recent call last):
            ...
            AttributeError: Extensive Form Game cannot have two roots


        A node can also be changed into a leaf if it has no parents, children, or actions. (i.e it is a blank node)::

            sage: jones_1 = Node()
            sage: jones_1.change_type(Leaf)


        If a node has any attribues, an error is returned::

            sage: williams_1 = Node({'A', 'B'})
            sage: williams_1.change_type(Leaf)
            Traceback (most recent call last):
            ...
            AttributeError: Node has attributes, cannot be leaf.
        """


class Leaf():
    def __init__(self, payoffs):
        """
        We can check payoffs of any leaf.
            sage: leaf_1 = Leaf([0,1])
            sage: leaf_1.payoffs
            [0,1]
        """

    def changepayoffs(something):
        """
        We can change the payoff of a leaf after setting it
            sage: leaf_1.changepayoffs([2,4])
            sage: leaf_1.payoffs
            [2,4]
        """

class Root(Node):
    """
    A root is just another type of node, so we can get attributes, however Parent will always be false. Attempting to add a parent will return an error::
        sage: bethan_1 = Root({'Red': jess_1, 'Blue': jess_2}))
        sage: bethan_1.attribues
        some output of attributes
        

    We cannot have more than one Root in a game, so if we try to connect a second Root to a connected set of nodes that already have a root, an error will be displayed::
        sage: jess_1 = Node([Green, Yellow]); jess_2 = Node([Green, Yellow]); jess_3 = Node([Green, Yellow])
        sage: bonnie_1 = Root({'Black': jess_1, 'White': jess_3})
        Traceback (most recent call last):
        ...
        AttributeError: Extensive Form Game cannot have two roots
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

