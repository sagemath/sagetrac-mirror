"""
Commutativity Axiom
"""

from axiom import Axiom


class Commutative(Axiom):

    class ParentMethods:
        
        def is_commutative(self):
            """
            Return True, since commutative magmas are commutative.
        
            EXAMPLES::
            
            sage: Parent(QQ,category=CommutativeRings()).is_commutative()
            True
            """
            return True



class AdditiveCommutative(Axiom):
    pass
    
