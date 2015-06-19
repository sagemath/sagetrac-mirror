"""
This file will contain a class for extensive form games.
"""

class ExtensiveFormGame():
    """
    """
    def checkplayers(self):
        """
        I may not include this
        """

    def assignplayers(self, players):
        """
        Another example of what's to come::
            sage: h = ExtensiveTree({0:[1,2], 1:[3,4], 2:[5,6]})
            sage: f = ExtensiveFormGame()
            sage: g = Players({"Player One":[0,3,4,5,6], "Player 2":[1, 2]})
            sage: f.assignplayers(g)
            sage: f.checkplayers()
            The game has the following players ...               
    
    
        If any of the nodes being assigned do not exist, assignplayers should return an error::
            sage: g = Players({"Player One":[0,3,4,5,6], "Player 2":[1, 2]})
            sage: h = ExtensiveTree({0:[1,2], 1:[3,4], 2:[5,6]})
            sage: f = ExtensiveFormGame()
            sage: f.assignplayers(g)
            Traceback (most recent call last):
            ...            
            SomeError: Tried to assign players to non existant nodes
        """
    def assigninfoset(self, informationset):
        """
        Here we want to assign an info set::
            sage: g = ExtensiveFormGame()
            sage: f = InfoSet(2)
            sage: g.assigninfoset(f)
            The information set '' has been assigned to the nodes ''.

        It should be able to reject any information set given to nodes with different numbers of successors::
            sage: g = ExtensiveFormGame()
            sage: a = 1; b = 2
            sage: f = InfoSet(a)
            sage: g.assigninfoset(f)
            Traceback (most recent call last):
            ...
            SomeError: Tried to assign information set to nodes with unequal successors
        """
    def assigntree(self, extensivetree):
        """
        Here we want to allow a tree stucture to be assigned to the game::
            sage: g = ExtensiveFormGame()
            sage: t = ExtensiveTree:({0:[1,2], 1:[3,4], 2:[5,6,7]})
            sage: g.assigntree(t)
            sage: g.gameinfo()
            An extensive form game with tree format ({0:[1,2], 1:[3,4], 2:[5,6,7]})...
        """
    def gameinfo(self):
        """
        Here is the sort of thing we'll want to return:
            sage: g = ExtensiveFormGame()
            sage: g.gameinfo()
            An extensive form game with tree format ({}), players: (Player One:[nodes], Player Two:[nodes], ...), information sets: (one [0,1], two [nodes]), root node [this is it].
        """

class ExtensiveTree():
    """
    
    Potential input and output using Networkx dictionary format::
        sage: g = ExtensiveTree({0:[1,2], 1:[3,4], 2:[5,6]})
        sage: print g
        Tree with (output text to be determined)


    Should be able to read other formats, potentially adjacency matrix::
        sage: m = matrix([0,1])
        sage: g = ExtensiveTree(m)
        sage: print g
        Tree with ...


    If we do end up accepting adjacency matrix, it should be able to tell when it's not in a tree format::
        sage: m = matrix([0,1])
        sage: g = ExtensiveTree(m)
        sage: print g
        Traceback (most recent call last):
        ...
        SomeError: Input not tree.

    """
    def __init__(self, etinput = None):
        self.etinput = etinput
        
        
class Players():
    """
    """
    def __init__(self, pinput = 'none'):
        self.pinput = pinput
        
    

class ChanceNode():
    """
    """

class InfoSet():
    """
    Information set class, which will simply be grouping nodes together
    """
    def __init__(self, iset = None):
        self.iset = iset
        
