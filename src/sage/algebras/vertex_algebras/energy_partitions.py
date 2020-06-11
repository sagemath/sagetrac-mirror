"""
Energy Partitions

AUTHORS:

- Reimundo Heluani (06-09-2020): Initial implementation.

.. linkall

"""


#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.combinat.partition import Partitions, RegularPartitions, Partition
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.all import QQ, ZZ, NN
from sage.combinat.partitions import ZS1_iterator_nk
from sage.functions.other import ceil

class EnergyPartition(Partition):
    """
    EnergyPartitions Element class
    """
    def energy(self):
        """
        Returns the energy of this element.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: v = EnergyPartitions(1/2,4,length=2).an_element(); v
            [4, 1]
            sage: v.energy()
            4

        TESTS::

            sage: EnergyPartitions(1/2)([]).energy()
            0
        """
        return self.length()*(self.parent().w-1) + self.size()

class EnergyPartitions(Partitions):
    r"""
    Base class for Energy Partitions.

    INPUT:

    - ``w`` -- a positive rational number; the energy of the partition
      `[1]`
    - ``n`` -- a non-negative rational number or ``None`` (default :
      ``None``); if ``None`` it returns the class of all Energy
      Partitions, if a positive rational number it returns the class
      of Energy Partitions with total energy ``n``

    in addition the following keyword arguments are recognized:

    - ``regular`` -- a positive integer number; the maximum number of
      times that a part can appear plus one
    - ``length`` -- a non-negative integer number; the number of parts

    EXAMPLES:

    We lists some partitions::

        sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
        sage: EnergyPartitions(1/2,3).list()
        [[1, 1, 1, 1, 1, 1], [2, 1, 1, 1], [3, 1], [2, 2]]
        sage: EnergyPartitions(1/2,3,regular=2).list()
        [[3, 1]]
        sage: EnergyPartitions(1/2,3,length=3).list()
        []
        sage: EnergyPartitions(1/2,5,regular=2,length=2).list()
        [[5, 1], [4, 2]]
        sage: EnergyPartitions(1/3,3).list()
        [[1, 1, 1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1], [3, 1, 1], [2, 2, 1]]
        sage: EnergyPartitions(1,3).list()
        [[3], [2, 1], [1, 1, 1]]

    Check the energy of some partitions::

        sage: V = EnergyPartitions(1/2,regular=2); V
        2-Regular Energy Partitions with weight 1/2
        sage: V.an_element()
        [1]
        sage: V([1]).energy()
        1/2
        sage: V([5,2]).energy()
        6
        sage: V([]).energy()
        0

    Notice that the energy depends on the weight::

        sage: v = EnergyPartitions(1).an_element(); v
        [3, 2, 1]
        sage: v.energy()
        6

        sage: v = EnergyPartitions(1/2).an_element(); v
        [3, 2, 1]
        sage: v.energy()
        9/2
        sage: v.parent()
        Energy Partitions with weight 1/2

    TESTS::

        sage: V([5,5,2])
        Traceback (most recent call last):
        ...
        ValueError: [5, 5, 2] is not an element of 2-Regular Energy Partitions with weight 1/2
        sage: EnergyPartitions(1,0).list()
        [[]]
        sage: EnergyPartitions(1/2,-1)
        Traceback (most recent call last):
        ...
        ValueError: n must be a non-negative rational or be equal to `None`
    """


    @staticmethod
    def __classcall_private__(cls, w = None,n = None, **kwargs):
        regular = kwargs.get('regular',None)
        if regular is not None:
            if regular not in NN or regular <= 0:
                raise ValueError("regular must be a positive integer.")
        length = kwargs.get('length',None)
        if length is not None:
            if length not in NN or length < 0:
                raise ValueError("length must be a non-negative integer.")

        if w not in QQ or w <= 0:
            raise ValueError("w must be a positive rational number, got {}"\
                                .format(w))
        if n is None:
            if len(kwargs) > 0:
                if len(kwargs) == 1:
                    if 'regular' in kwargs:
                        return RegularEnergyPartitions_all(w,regular)
                raise ValueError("n must be a non-negative rational with these"\
                                 " keyword arguments")
            return EnergyPartitions_all(w)
        elif n in QQ and n >= 0:
            if len(kwargs) == 0:
                return EnergyPartitions_n(w,n)

            if len(kwargs) == 1:
                if 'length' in kwargs:
                    return EnergyPartitions_nk(w,n,length)

                elif 'regular' in kwargs:
                    return RegularEnergyPartitions_n(w, n, regular)
            elif len(kwargs) == 2:
                if 'length' in kwargs and 'regular' in kwargs:
                    return RegularEnergyPartitions_nk(w,n,length,regular)
            raise ValueError("Unrecognized combination of keyword arguments")

        raise ValueError("n must be a non-negative rational or be equal to "\
                         "`None`")

    def __init__(self, w, is_infinite=False):
        super(EnergyPartitions,self).__init__(is_infinite)
        self.w = QQ(w)

    Element = EnergyPartition

class EnergyPartitions_all(EnergyPartitions):
    def __init__(self,w):
        """
        Base class for all Energy Partitions.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2)[0:10]
            [[],
             [1],
             [1, 1],
             [1, 1, 1],
             [2],
             [1, 1, 1, 1],
             [2, 1],
             [1, 1, 1, 1, 1],
             [2, 1, 1],
             [3]]

        TESTS::

            sage: EnergyPartitions(0)
            Traceback (most recent call last):
            ...
            ValueError: w must be a positive rational number, got 0
        """
        EnergyPartitions.__init__(self,w,is_infinite=True)

    def _repr_(self):
        """
        A visual representation of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2)
            Energy Partitions with weight 1/2
        """
        return "Energy Partitions with weight {}".format(self.w)

    def an_element(self):
        """
        An element of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1).an_element()
            [3, 2, 1]
        """
        return self.element_class(self, [3,2,1])

    def __iter__(self):
        """
        Iterate over all Energy Partitions.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1)[0:5]
            [[], [1], [2], [1, 1], [3]]
            sage: EnergyPartitions(1/2)[0:5]
            [[], [1], [1, 1], [1, 1, 1], [2]]
        """
        yield self.element_class(self,[])
        n = self.w
        while True:
            for p in EnergyPartitions_n(self.w,n):
                yield self.element_class(self,p)
            n += 1/self.w.denominator()

class EnergyPartitions_n(EnergyPartitions):
    def __init__(self, w, n):
        """
        The class of Energy Partitions with a fixed energy.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: V = EnergyPartitions(1/2,3); V
            Energy Partitions of 3 with weight 1/2
            sage: V.list()
            [[1, 1, 1, 1, 1, 1], [2, 1, 1, 1], [3, 1], [2, 2]]
            sage: V = EnergyPartitions(1/2,3, regular=2); V
            2-Regular Energy Partitions of energy 3 with weight 1/2
            sage: V.list()
            [[3, 1]]
        """
        EnergyPartitions.__init__(self,w)
        self.n = QQ(n)

    def __contains__(self,x):
        """
        Whether this element is an Energy Partition in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: [3,1] in EnergyPartitions(1/2,3)
            True
            sage: [3,1] in EnergyPartitions(1,3)
            False

        TESTS::

            sage: [] in EnergyPartitions(1,3)
            False
            sage: [] in EnergyPartitions(1,0)
            True
        """
        return x in EnergyPartitions(self.w) and \
                sum(x) == self.n - (self.w-1)*len(x)

    def _repr_(self):
        """
        A visual representation of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1,2)
            Energy Partitions of 2 with weight 1
        """
        return "Energy Partitions of {} with weight {}".format(self.n, self.w)

    def an_element(self):
        """
        A Partition in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,4).an_element()
            [4, 1]
            sage: EnergyPartitions(1,3).an_element()
            [3]
        """
        if self.n == 0:
            lst = []
        elif self.n < self.w or self.w.denominator()*self.n not in NN:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("There are no Energy Partitions of energy {}"\
                                " and weight {}".format(self.n,self.w))
        else:
            k = 1
            while True:
                t = self.n - k*self.w
                if t in ZZ:
                    lst = [t+1]+[1]*(k-1)
                    break
                k += 1
        return self.element_class(self,lst)

    def cardinality(self):
        """
        The number of partitions in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/3,5).cardinality()
            11
            sage: EnergyPartitions(1,8).cardinality() == Partitions(8).cardinality()
            True

        TESTS::

            sage: EnergyPartitions(1/2,0).cardinality()
            1
        """
        ieqs = [[0,1,0],[0,0,1]]
        nn = self.n.numerator()
        nd = self.n.denominator()
        wn = self.w.numerator()
        wd = self.w.denominator()
        eqns = [[-nn*wd, wn*nd,wd*nd]]
        P = Polyhedron(ieqs=ieqs, eqns=eqns)
        return sum(Partitions(y,max_length=x).cardinality()\
                    for x,y in P.integral_points())

    def __iter__(self):
        """
        Iterate over the set of these partitions.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,5).list()
            [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
             [2, 1, 1, 1, 1, 1, 1, 1],
             [3, 1, 1, 1, 1, 1],
             [2, 2, 1, 1, 1, 1],
             [4, 1, 1, 1],
             [3, 2, 1, 1],
             [2, 2, 2, 1],
             [5, 1],
             [4, 2],
             [3, 3]]

        TESTS::

            sage: EnergyPartitions(1/2,0).list()
            [[]]
        """
        ieqs = [[0,1,0],[0,0,1]]
        wn = self.w.numerator()
        wd = self.w.denominator()
        eqns = [[-self.n*wd, wn,wd]]
        P = Polyhedron(ieqs=ieqs, eqns=eqns)
        for x,y in P.integral_points():
            for p in ZS1_iterator_nk(y,x):
                v = [i+1 for i in p]
                adds = [1]*(x - len(v))
                yield self.element_class(self,v+adds)

class EnergyPartitions_nk(EnergyPartitions):
    def __init__(self, w, n, k):
        """
        The class of Energy Partitions with prescribed energy and
        length.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,5,length=4).list()
            [[4, 1, 1, 1], [3, 2, 1, 1], [2, 2, 2, 1]]
            sage: EnergyPartitions(1,8,length=2).list()
            [[7, 1], [6, 2], [5, 3], [4, 4]]
        """
        EnergyPartitions.__init__(self,w)
        self.n = QQ(n)
        if k not in ZZ or k < 0:
            raise ValueError("length must be a non-negative integer")
        self.k = k

    def __contains__(self,x):
        """
        Whether ``x`` is an Energy Partition in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: [3,3] in EnergyPartitions(1/2,5,length=2)
            True
            sage: [3,3] in EnergyPartitions(1/2,5,length=3)
            False

        TESTS::

            sage: [] in EnergyPartitions(1/2,0,length=0)
            True
            sage: [] in EnergyPartitions(1/2,0,length=3)
            False
            sage: [0,0,0] in EnergyPartitions(1/2,0,length=3)
            False
        """
        return x in EnergyPartitions(self.w) and \
                sum(x) == self.n - (self.w-1)*self.k and len(x) == self.k

    def _repr_(self):
        """
        A visual representation of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1,8,length=2)
            Energy Partitions of 8 of length 2 with weight 1
        """
        return "Energy Partitions of {} of length {} with weight {}"\
                .format(self.n, self.k, self.w)

    def an_element(self):
        """
        An Energy Partition in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1,8,length=2).an_element()
            [7, 1]
            sage: EnergyPartitions(1/3,4,length=3).an_element()
            [4, 1, 1]
        """
        rest = self.n - self.k*self.w
        if rest not in NN:
                from sage.categories.sets_cat import EmptySetError
                raise EmptySetError
        if rest == 0:
            lst = [1]*self.k
        else:
            lst = [rest + 1] + [1]*(self.k -1)
        return self.element_class(self,lst)

    def __iter__(self):
        """
        Iterate over all Energy Partitions in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/3,4,length=3).list()
            [[4, 1, 1], [3, 2, 1], [2, 2, 2]]
        """
        rest = self.n - self.k*self.w
        if rest not in NN:
            return
        for p in ZS1_iterator_nk(rest, self.k):
            v = [i + 1 for i in p]
            adds = [1] * (self.k - len(v))
            yield self.element_class(self, v + adds)

    def cardinality(self):
        """
        The number of Energy Partitions in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,28,length=4).cardinality()
            206

        TESTS::

            sage: EnergyPartitions(1/2,0,length=0).cardinality()
            1
            sage: EnergyPartitions(1/2,0,length=1).cardinality()
            0
            sage: EnergyPartitions(1/2,1,length=0).cardinality()
            0
        """
        rest = self.n -self.k*self.w
        if rest not in NN:
            return 0
        return Partitions(ZZ(rest),max_length=self.k).cardinality()

class RegularEnergyPartitions(EnergyPartitions):
    def __init__(self, w, ell, is_infinite=False):
        """
        Base class for Regular Energy Partitions.

        INPUT:

        - ``w`` -- a positive rational; the weight.
        - ``is_infinite`` a boolean; wether this class of partitions
          contains infinitely many elements.
        """
        EnergyPartitions.__init__(self, w, is_infinite)
        self._ell = ell

    def __contains__(self,x):
        """
        Whether this element is a Regular Energy Partition.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: [] in EnergyPartitions(1/2,regular=2)
            True
            sage: [2,1] in EnergyPartitions(1/2,regular=2)
            True
            sage: [1,1] in EnergyPartitions(1/2,regular=2)
            False
        """
        return EnergyPartitions.__contains__(self,x) and \
                RegularPartitions.__contains__(self,x)

    def _fast_iterator(self, n, max_part):
        """
        A recursive iterator to go over regular partitions.
        """
        if n == 0:
            yield []
            return

        if n < max_part:
            max_part = n
        bdry = self._ell - 1

        for i in reversed(range(1, max_part+1)):
            for p in self._fast_iterator(n-i, i):
                if p.count(i) < bdry:
                    yield [i] + p


class RegularEnergyPartitions_all(RegularEnergyPartitions):
    def __init__(self, w, ell):
        """
        The class of all Regular Energy Partitions of this weight.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,regular=3)
            3-Regular Energy Partitions with weight 1/2
        """
        RegularEnergyPartitions.__init__(self, w, ell, bool(ell > 1))

    def _repr_(self):
        """
        The name of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,regular=3)
            3-Regular Energy Partitions with weight 1/2
        """
        return "{}-Regular Energy Partitions with weight {}".format(self._ell,
                                                                        self.w)
    def __iter__(self):
        """
        Iterate over Regular Energy Partitions of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,regular=3)[0:10]
            [[1],
             [1, 1],
             [2],
             [2, 1],
             [2, 1, 1],
             [3],
             [3, 1],
             [2, 2],
             [3, 1, 1],
             [2, 2, 1]]

        TESTS::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,regular=1).list()
            [[]]
        """
        if self._ell == 1:
            yield self.element_class(self,[])
            return

        n = self.w
        while True:
            for p in RegularEnergyPartitions_n(self.w,n,self._ell):
                yield(self.element_class(self,p))
            n += 1/self.w.denominator()



class RegularEnergyPartitions_n(RegularEnergyPartitions, EnergyPartitions_n):
    def __init__(self, w, n, ell):
        """
        Regular Energy Partitions of a fixed energy.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,4,regular=2)
            2-Regular Energy Partitions of energy 4 with weight 1/2
            sage: EnergyPartitions(1/2,4,regular=2).list()
            [[4, 1], [3, 2]]
        """
        RegularEnergyPartitions.__init__(self, w, ell)
        EnergyPartitions_n.__init__(self,w,n)

    def _repr_(self):
        """
        The name of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,4,regular=2)
            2-Regular Energy Partitions of energy 4 with weight 1/2
        """
        return "{}-Regular Energy Partitions of energy {} with weight {}".\
                    format(self._ell,self.n,self.w)

    def __contains__(self,x):
        """
        Whether this class contains ``x``.

        EXAMPLES::

        sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
        sage: [3,2,1] in EnergyPartitions(1/2,9/2,regular=2)
        True
        """
        return RegularPartitions.__contains__(self,x) and \
               EnergyPartitions_n.__contains__(self,x)

    def __iter__(self):
        """
        Iterate over all Energy Partitions in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,17/2,regular=2).list()
            [[7, 2, 1], [6, 3, 1], [5, 4, 1], [5, 3, 2], [9]]
        """
        ieqs = [[0,1,0],[0,0,1]]
        wn = self.w.numerator()
        wd = self.w.denominator()
        eqns = [[-self.n*wd, wn,wd]]
        P = Polyhedron(ieqs=ieqs, eqns=eqns)
        for x,y in P.integral_points():
            for p in self._fast_iterator(y,y):
                if 0 <= x - len(p) < self._ell:
                    v = [i+1 for i in p]
                    adds = [1]*(x - len(v))
                    yield self.element_class(self,v+adds)

    def cardinality(self):
        """
        The number of Energy Partitions in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,17/2,regular=2).cardinality()
            5
        """
        if self._ell > self.n:
            return EnergyPartitions_n.cardinality(self)
        return ZZ.sum(1 for x in self)

class RegularEnergyPartitions_nk(RegularEnergyPartitions, EnergyPartitions_nk):
    def __init__(self,w,n,k,ell):
        """
        The class of Regular Energy Partitions of prescribed energy and
        length.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,17/2,length=3,regular=2)
            2-Regular Energy Partitions of energy 17/2 and length 3 with weight 1/2
        """
        RegularEnergyPartitions.__init__(self,w,ell)
        EnergyPartitions_nk.__init__(self,w,n,k)

    def _repr_(self):
        """
        The name of this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,17/2,length=3,regular=2)
            2-Regular Energy Partitions of energy 17/2 and length 3 with weight 1/2
        """
        return "{}-Regular Energy Partitions of energy {} and length {} with "\
               "weight {}".format(self._ell, self.n, self.k, self.w)

    def __contains__(self,x):
        """
        Whether this class contains ``x``.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: [5,4,1] in EnergyPartitions(1/2,17/2,length=3,regular=2)
            True

        TESTS::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: [] in EnergyPartitions(1/2,0,length=0,regular=1)
            True
        """
        return RegularPartitions.__contains__(self,x) and \
               EnergyPartitions_nk.__contains__(self,x)

    def __iter__(self):
        """
        Iterate over all Energy Partitions in this class.

        EXAMPLES::

        sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
        sage: EnergyPartitions(1/2,17/2,length=3,regular=2).list()
        [[7, 2, 1], [6, 3, 1], [5, 4, 1], [5, 3, 2]]
        """
        if self.k == 0:
            if self.n == 0:
                yield self.element_class(self,[])
            return None

        rest = self.n - self.k*self.w
        if rest not in NN:
            return
        for p in self._fast_iterator(rest,rest):
            if self.k - len(p) < self._ell:
                v = [i + 1 for i in p]
                adds = [1] * (self.k - len(v))
                yield self.element_class(self, v + adds)

    def _fast_iterator(self, n, max_part, depth=0):
        """
        A recursive iterator which returns a list (from
        RegularPartitions_truncated)

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: P = EnergyPartitions(1/2,7,length=2,regular=2)
            sage: list(P._fast_iterator(5, 5))
            [[5], [4, 1], [3, 2]]
            sage: list(P._fast_iterator(5, 3))
            [[3, 2]]
            sage: list(P._fast_iterator(5, 6))
            [[5], [4, 1], [3, 2]]
        """
        if n == 0 or depth >= self.k:
            yield []
            return

        # Special case
        if depth + 1 == self.k:
            if max_part >= n:
                yield [n]
            return

        if n < max_part:
            max_part = n
        bdry = self._ell - 1

        for i in reversed(range(1, max_part+1)):
            for p in self._fast_iterator(n-i, i, depth+1):
                if p.count(i) < bdry:
                    yield [i] + p

    def cardinality(self):
        """
        The number of Energy Partitions in this class.

        EXAMPLES::

            sage: from sage.algebras.vertex_algebras.energy_partitions import EnergyPartitions
            sage: EnergyPartitions(1/2,17/2,length=3,regular=2).cardinality()
            4
        """
        if self._ell > self.n:
            return EnergyPartitions_nk.cardinality(self)
        return ZZ.sum(1 for x in self)

