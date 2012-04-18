

class Map(object):
    """
    Map from the set of right cosets of Gamma0(N) in SL_2(Z) to a
    coefficient module.
    """
    def __init__(self, manin_relations, defining_data):
        self._manin_relations = manin_relations
        raise NotImplementedError
        self._codomain = ?

    def __add__(self, right):
        """
        Return difference self + right.
        """
        raise NotImplementedError

    def __sub__(self, right):
        """
        Return difference self - right.
        """
        raise NotImplementedError

    def __mul__(self, left):
        """
        Return left*self, where left acts on the codomain.
        """
        raise NotImplementedError

    def __repr__(self):
        return "Map from the set of right cosets of Gamma0(%s) in SL_2(Z) to %s"%(
            self._manin_relations.level(), self._codomain)

    
    

    
