# The parent
############

class pAdicFieldLattice(pAdicRingBaseGeneric):
    # Internal variables:
    #  . self._working_precision_cap
    #    a cap for the working precision
    #    meaning that the precision lattice always contains p^(self._working_precision_cap)
    #  . self._elements
    #    list of weak references of elements in this parent
    #  . self._precision_lattice
    #    a matrix over ZZ representing the lattice of precision
    #    (its columns are indexed by self._elements)

    def __init__(self, p, working_precision_cap, print_mode):
        # Initialize variables here

    def _echelonize(self):
        # Echelonize the matrix giving the precision lattice

    def _add_element(self, elt, prec):
        # A new element in the parent has just been created
        # We should:
        #  . add a weak reference to it to the list self._elements
        #  . add a column to self._precision_lattice
        #    and update this matrix according to prec 
        #    (and possibly the working precision cap)
        # NOTE: prec be either an integer or a formal linear combinaison
        # of the precision on the other elements of this parent

    def _del_element(self, elt):
        # The element elt has just been garbage collected
        # We should:
        #  . remove it to the list self._elements
        #  . remove the corresponding column of self._precision_lattice

    def precision_absolute(self, elt):
        # Return the (optimal) absolute precision of the element elt
        # This precision can be read off on self._precision_lattice:
        # it is the smallest valuation of an entry of the column of
        # self._precision_lattice corresponding to the element elt

    def working_precision(self, elt):
        # Return the working precision of the element elt
        # This precision can be read off on the precision lattice:
        # it is the smallest integer n for which the precision
        # lattice contains p^n*[elt] where [elt] denotes the 
        # basis vector corresponding to elt

    def _element_constructor_(self, x, prec):
        # We ask for the creation of an element in this parent
        # We should:
        #  . create this element (called elt hereafter) from x
        #    Note: elt is an instance of the class pAdicLatticeElement below
        #  . call the method _new_element(elt, prec)
        #  . install a callback so that when elt will be collected
        #    by the garbage collector, the method _del_element will
        #    be called

QpLP = pAdicFieldLattice


# The elements
##############

class pAdicLatticeElement(Element):
    # Internal variable:
    #  . self._approximation
    #    an approximation of this p-adic number
    #    it is defined modulo p^(working_precision)

    def working_precision(self):
        return self.parent().working_precision(self)

    def approximation(self):
        # We should:
        #  . reduce self._approximation modulo p^(self.working_precision())
        #    (and update self._approximation accordingly)
        #  . return the result

    def precision_absolute(self):
        return self.parent().precision_absolute(self)

    def valuation(self, secure=False):
        # We should:
        #  . compute the valuation of self.approximation()
        #  . compare it to the self.precision_absolute()
        #  . if the former is less than the latter, we return the former
        #  . otherwise, if secure is False, we return the latter
        #               if secure is True, we raise an error

    def precision_relative(self, secure=False):
        return self.precision_absolute() - self.valuation(secure=secure)

    def _repr_(self):
        # Print like this:
        #   self.approximation() + O(p^self.precision_absolute())


    def _add_(self, other):
        # We should:
        #  . compute app = self.approximation() + other.approximation()
        #  . create a new element from app and the precision given by the differential:
        #    dapp = dself + dother

    def _mul_(self, other):
        # We should:
        #  . compute app = self.approximation() * other.approximation()
        #  . create a new element from app and the precision given by the differential:
        #    dapp = self*dother + other*dself
