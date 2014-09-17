from sage.combinat.permutation import Permutation, Permutations
from sage.misc.abstract_method import abstract_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.sage_object import SageObject

def bijection(func):
    """
    TESTS::

        sage: from sage.misc.bijection import bijection
        sage: class A:
                @bijection
                def foo(self, sigma):
                    return sigma
        sage: s = A().foo(Permutation([1,2,3])); s
        <class 'sage.misc.bijection.BijectionFromPermutation'>
        sage: ~s
        <class 'sage.misc.bijection.BijectionFromPermutation'>
    """
    def wrapper(self, bij):
        return func(self, Bijection(bij))
    return wrapper


class Bijection(SageObject):

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, param):
        """
        TESTS::

            sage: from sage.misc.bijection import Bijection
            sage: b = Bijection({1:2,2:1})
            sage: b(1)
            2
            sage: b((~b)(1))
            1

            sage: c = Bijection(Permutation([3,1,2]))
            sage: (~c)(1)
            2

            sage: d = Bijection([1,3,2])
            sage: d(3)
            2
            sage: (~d)(3)
            2
        """
        if isinstance(param, dict):
            return BijectionFromDict(param)
        elif param in Permutations():
            return BijectionFromPermutation(Permutation(param))
        else:
            raise NotImplemented

    @abstract_method
    def __call__(self, obj):
        """
        Method to encode `f(x)`.
        """

    @abstract_method
    def __invert__(self):
        """
        Method to encode `f^{-1}`
        """

class BijectionFromDict(Bijection):

    def __init__(self, dct):
        self._dct_ = dct

    def __call__(self, obj):
        return self._dct_[obj]

    def __invert__(self):
        dct = self._dct_
        dctInv = {dct[k]:k for k in dct.keys()}
        return BijectionFromDict(dctInv)

class BijectionFromPermutation(Bijection):

    def __init__(self, sigma):
        self._sigma_ = sigma

    def __call__(self, obj):
        return self._sigma_[obj-1]

    def __invert__(self):
        return BijectionFromPermutation(self._sigma_.inverse())

