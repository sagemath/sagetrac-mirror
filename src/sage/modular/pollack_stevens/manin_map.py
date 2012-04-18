"""
Manin Maps.

This is a class that represents maps from a set of right coset
representatives to a coefficient module.  This is a basic building
block for implementing modular symbols, and provides basic arithmetic
and right action of matrices.


"""



class ManinMap(object):
    """
    Map from a set of right coset representatives of Gamma0(N) in
    SL_2(Z) to a coefficient module that satisfies the Manin
    relations.
    """
    def __init__(self, manin_relations, defining_data):
        self._manin_relations = manin_relations
        raise NotImplementedError
        #self._codomain = ?

    def __add__(self, right):
        """
        Return difference self + right, where self and right are
        assumed to have identical codomains and Manin relations.
        """
        raise NotImplementedError

    def __sub__(self, right):
        """
        Return difference self - right, where self and right are
        assumed to have identical codomains and Manin relations.
        """
        raise NotImplementedError

    def __mul__(self, right):
        """
        Return scalar multiplication self*right, where right is in the
        base ring of the codomain.
        """
        if isinstance(right, sage.matrix.matrix_integer_2x2.Matrix_integer_2x2):
            return self._right_action(right)
        raise NotImplementedError

    def __repr__(self):
        return "Map from the set of right cosets of Gamma0(%s) in SL_2(Z) to %s"%(
            self._manin_relations.level(), self._codomain)

    def __call__(self, A):
        # NOTE: we should use the high performance specialized type
        #   sage.matrix.matrix_integer_2x2.Matrix_integer_2x2
        # for our matrices.
        raise NotImplementedError

    def apply(self, f):
        """
        Returns Manin map given by x |--> f(self(x)), where f is
        anything that can be called with elements of the coefficient
        module.

        This might be used to normalize, reduce modulo a prime, change
        base ring, etc.
        """
        raise NotImplementedError

    def __iter__(self):
        """
        Returns iterator over the values of this map on the reduced
        representatives.

        This might be used to compute the valuation.
        """
        raise NotImplementedError

    def _right_action(self, gamma):
        """
        Returns self | gamma, where gamma is a 2x2 integer matrix.
        
        The action is defined by (self | gamma)(D) = self(gamma D)|gamma

        For the action by a single element gamma to be well defined,
        gamma must normalize Gamma_0(N).  However, this right action
        can also be used to define Hecke operators, in which case each
        individual self | gamma is not a modular symbol on Gamma_0(N),
        but the sum over acting by the appropriate double coset
        representatives is.

        INPUT:

        - ``gamma`` - 2 x 2 matrix which acts on the values of self

        OUTPUT:

        - ManinMap
        """
        raise NotImplementedError
    

    
