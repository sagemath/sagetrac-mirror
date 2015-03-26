import six
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.categories.sets_cat import Sets
from sage.combinat.tableaux.abstract_tableau import (AbstractTableau,
     BadShapeTableau, SkewTableau, StraightTableau)
from sage.misc.classcall_metaclass import ClasscallMetaclass

class AbstractTableaux(UniqueRepresentation, Parent):
    r"""
    Class of all abstract tableaux.
    """
    Element = AbstractTableau

    def __init__(self):
        r"""
        Initialize the parent.
        """
        Parent.__init__(self, category=Sets())

    def _repr_(self):
        r"""
        Return the representation string.
        
        OUTPUT:
        
        A string.
        """
        return "Abstract Tableaux"

    def _element_constructor_(self, t):
        r"""
        Constructs an object from ``t`` as an element of ``self``,
        if possible.
        
        INPUT:
        
        - ``t`` -- Data which can be interpreted as a tableau
        
        OUTPUT:
        
        - The corresponding tableau object
        """
        if not t in self:
            raise ValueError("%s is not an element of %s."%(t, self))

        return self.element_class(self, t)

class BadShapeTableaux(AbstractTableaux):
    Element = BadShapeTableau

    @staticmethod
    def _element_constructor_(cls, dct, check=True):
        r"""
        Construct a new BadShapeTableau, optionally validating input.
        
        INPUT:
        
        - ``dct`` -- a dictionary (or more generally something
          passable to ``dict``) whose keys are pairs of integers
        - ``check`` -- (default: ``True``) if ``True``, then check that
          the keys of ``dct`` are in fact pairs of integers
        """
        dct = deepcopy(dict(dct))

        if check:
            if not all(x in ZZ and y in ZZ for x, y
                       in six.iterkeys(dct)):
                raise ValueError('keys must be pairs of integers')

        return cls.element_class(cls, dct)

    def _repr_(self):
        r"""
        Return the representation string.
        
        OUTPUT:
        
        A string.
        """
        return "Bad Shape Tableaux"

class SkewTableaux(BadShapeTableaux):
    Element = SkewTableau
    
    def _repr_(self):
        r"""
            Return the representation string.
            
            OUTPUT:
            
            A string.
            """
        return "Skew Tableaux"


class StraightTableaux(SkewTableaux):
    Element = StraightTableau

    def _repr_(self):
        r"""
            Return the representation string.
            
            OUTPUT:
            
            A string.
            """
        return "Straight Tableaux"