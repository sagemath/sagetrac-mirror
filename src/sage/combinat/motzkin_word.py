from sage.structure.list_clone import ClonableArray
from .backtrack import GenericBacktracker

class MotzkinWord(ClonableArray):
    def __init__(self, parent, word):
        for i in range(1,len(word)+1):
            if sum(word[1:i])<0:
                raise ValueError("A Motzkinword is not allowed to go beneath the starting level, error at position {}".format(i))
            if word[i] notin {1,-1,0}:
                raise ValueError("Only steps in {1,-1,0} are allowed, however you used {}.".format(word[i]))
        ClonableArray.__init__(self, parent,word)
    def __repr__(self):
        return "MotzkinWord {}".format(self.word)


class MotzkinWords(UniqueRepresentation,Parent):
    def __init__(self,n):
        self._n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())
    def __repr__(self):
        return "Motzkinwords of lenght [{}]".format(self._n)
    def __iter__(self):
        n = 0
        yield self.element_class(self, [])
        while True:
            for k in range(1, n+1):
                for x in DyckWords_size(k, n-k):
                    yield self.element_class(self, list(x))
            n += 1


class MotzkinWordBacktracker(GenericBacktracker):

    def __init__(self, k1, k2, k3):
        GenericBacktracker.__init__(self, [], (0, 0))
        # note that the comments in this class think of our objects as
        # Dyck paths, not words; having k1 opening parens and k2 closing
        # parens corresponds to paths of length k1 + k2 ending at height
        # k1 - k2.
        k1 = ZZ(k1)
        k2 = ZZ(k2)
        k3 = ZZ(k3)
        self.n = k1 + k2 + k3
        self.endht = k1 - k2

    def _rec(self, path, state):
        len, ht = state

        if len < self.n - 1:
            # if length is less than n-1, new path won't have length n, so
            # don't yield it, and keep building paths

            # if the path isn't too low and is not touching the x-axis, we can
            # yield a path with a downstep at the end
            if ht > (self.endht - (self.n - len)) and ht > 0:
                yield path + [-1], (len + 1, ht - 1), False

            # if the path isn't too high, it can also take an upstep
            if ht < (self.endht + (self.n - len)):
                yield path + [1], (len + 1, ht + 1), False

            # there can always be a horizontal step
            yield path + [0], (len + 1, ht), False

        else:
            # length is n - 1, so add a single step (up, down or horizontal,
            # according to current height and endht), don't try to
            # construct more paths, and yield the path
            if ht < self.endht:
                yield path + [1], None, True
            elif ht > self.endht:
                yield path + [-1], None, True
            else:
                yield path + [0], None, True


class MotzkinWords_size(MotzkinWords):
    def __init__(self, k1, k2, k3):
        self.k1 = ZZ(k1)
        self.k2 = ZZ(k2)
        self.k3 = ZZ(k3)
        MotzkinWords.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        return "Motzkin words with %s upsteps, %s downsteps and %s horziontalsteps" % (self.k1, self.k2, self.k3)

#    def __contains__(self, x):
#        return is_a(x, self.k1, self.k2, self.k3)

    def __iter__(self):
        if self.k1 == 0:
            yield self.element_class(self, [])
        elif (self.k2 == 0 and self.k3==0):
            yield self.element_class(self, [1]*self.k1)
        elif (self.k1 == 0 and self.k2==0):
            yield self.element_class(self, [0]*self.k1)
        else:
            for w in DyckWordBacktracker(self.k1, self.k2, self.k3):
                yield self.element_class(self, w)

#    def cardinality(self):
#        from sage.arith.all import binomial
#        return (self.k1 - self.k2 + 1) * binomial(self.k1 + self.k2, self.k2) // (self.k1 + 1)
