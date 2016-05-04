"""
Gosper iterator

A class which serves as a stateful iterable for computing the terms of the continued fraction of `(a*x+b)/(c*x+d)`,
where `a, b, c, d` are integers, and `x` is a continued fraction.

EXAMPLES:
::
    sage: from sage.rings.continued_fraction_gosper import gosper_iterator
    sage: x = continued_fraction(pi)
    sage: it = iter(gosper_iterator(3,2,3,1,x))
    sage: Word(it, length='infinite')
    word: 1,10,2,2,1,4,1,1,1,97,4,1,2,1,2,45,6,4,9,1,27,2,6,1,4,2,3,1,3,1,15,2,1,1,2,1,1,2,32,1,...
"""

from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR
from sage.rings.continued_fraction import ContinuedFraction_periodic

class gosper_iterator:


    def __init__(self, a, b, c, d, x):
        """
        Constructs the class.
        INPUT:

        - ``a, b, c, d`` -- Integer coefficients of the transformation.
        - ``x`` -- An instance of a continued fraction.

        OUTPUT:

        - The instance of gosper_iterator class.

        TESTS:
        ::
            sage: a = Integer(randint(-10,10)); b = Integer(randint(-10,10));
            sage: c = Integer(randint(-10,10)); d = Integer(randint(-10,10));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: x = continued_fraction(([1,2],[3,4])); i = iter(gosper_iterator(a,b,c,d,x))
            sage: l = list(i)
            sage: preperiod_length = i.output_preperiod_length
            sage: preperiod = l[:preperiod_length]
            sage: period = l[preperiod_length:]
            sage: continued_fraction((preperiod, period), x.value()) == continued_fraction((a*x.value()+b)/(c*x.value()+d))
            True
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d

        self.cf = x
        self.x = iter(x)

        self.states = []

        self.currently_emitted = 0
        self.currently_read = 0

        # Rational or quadratic case
        if isinstance(self.cf, ContinuedFraction_periodic):
            self.input_preperiod_length = self.cf.preperiod_length()
        # Infinite case
        else:
            self.input_preperiod_length = +Infinity

        self.output_preperiod_length = 0

    def __iter__(self):
        """
        Returns the iterable instance of the class. Is called upon `iter(gosper_iterator(a,b,c,d,x))`.

        TESTS:
        ::
            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: ig = iter(gosper_iterator(a,b,c,d,continued_fraction(pi))); icf = iter(continued_fraction((a*pi+b)/(c*pi+d)));
            sage: lg = [next(ig) for _ in xrange(10)]; lcf = [next(icf) for _ in xrange(10)];
            sage: lg == lcf
            True
        """
        return self

    def next(self):
        """
        Returns the next term of the transformation.

        TESTS:
        ::
            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: ig = iter(gosper_iterator(a,b,c,d,continued_fraction(pi))); icf = iter(continued_fraction((a*pi+b)/(c*pi+d)));
            sage: for i in range(10):
            ....:     assert next(ig) == next(icf)
        """
        limit = 100
        while True:
            if self.currently_read >= self.input_preperiod_length:
                current_state = {
                'a': self.a,
                'b': self.b,
                'c': self.c,
                'd': self.d,
                'next': self.cf.quotient(self.currently_read),
                'currently_emitted': self.currently_emitted,
                }
                for state in self.states:
                    if self.compare_dicts(state, current_state, ['currently_emitted']):
                        self.output_preperiod_length = state['currently_emitted']
                        raise StopIteration
                self.states.append(current_state)
                if len(self.states) > 100:
                    print "ERROR: Stopping iteration, danger of memory overflow."
                    raise StopIteration

            if (self.c == 0 and self.d == 0):
                raise StopIteration

            ub = self.bound(self.a, self.c)
            lb = self.bound(self.a + self.b, self.c + self.d)
            s = -self.bound(self.c, self.d)

            if ub == lb and s<1:
                self.emit(ub)
                return Integer(ub)
            else:
                self.ingest()

            limit -= 1
            if limit < 1:
                print "ERROR: Next loop iteration ran too many times."
                raise StopIteration


    def emit(self, q):
        """
        Changes the state of the iterator, correspondingly to emitting the term `q`.

        TESTS:
        ::
            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gi = gosper_iterator(a,b,c,d,continued_fraction(pi))
            sage: for i in range(10):
            ....:     gi.emit(i)
            sage: gi.currently_emitted
            10
        """
        self.currently_emitted += 1
        # This is being computed for the case when no states are being saved (still reading preperiod).
        if self.currently_read <= self.input_preperiod_length:
            self.output_preperiod_length = self.currently_emitted
        a = self.a
        b = self.b
        self.a = self.c
        self.b = self.d
        self.c = a - q*self.c
        self.d = b - q*self.d


    def ingest(self):
        """
        Changes the state of the iterator, correspondingly to ingesting another term from the input continued fraction.

        TESTS:
        ::
            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gi = gosper_iterator(a,b,c,d,continued_fraction(pi))
            sage: for i in range(10):
            ....:     gi.ingest()
            sage: gi.currently_read
            10
        """
        try:
            p = next(self.x)
            self.currently_read += 1
            a = self.a
            c = self.c
            self.a = a*p + self.b
            self.b = a
            self.c = c*p + self.d
            self.d = c
        except StopIteration:
            self.b = self.a
            self.d = self.c

    @staticmethod
    def bound(n,d):
        """
        Helper function for division. Returns infinity if denominator is zero.

        TESTS:
        ::
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gosper_iterator.bound(1,0)
            +Infinity
        """
        if d == 0:
            return +Infinity
        else:
            return (RR(n)/RR(d)).floor()

    @staticmethod
    def compare_dicts(d1, d2, ignore_keys):
        """
        Helper function, used to compare two dictionaries, ignoring the keys in `ignore_keys`.

        TESTS:
        ::
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gosper_iterator.compare_dicts({"compared": 1, "ignored": 2},{"compared": 1, "ignored": 3}, ["ignored"])
            True
        """
        d1_filtered = dict((k, v) for k,v in d1.iteritems() if k not in ignore_keys)
        d2_filtered = dict((k, v) for k,v in d2.iteritems() if k not in ignore_keys)
        return d1_filtered == d2_filtered