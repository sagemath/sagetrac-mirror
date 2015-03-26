
class AbstractTableaux(UniqueRepresentation, Parent):
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
                        

                        
class BadShapeTableaux(AbstractTableaux):
    Element = BadShapeTableau
    
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
    Element = SkewTableau
    
    def _repr_(self):
        r"""
            Return the representation string.
            
            OUTPUT:
            
            A string.
            """
        return "Straight Tableaux"