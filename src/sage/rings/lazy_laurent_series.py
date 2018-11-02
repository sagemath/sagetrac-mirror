from .infinity import infinity

from sage.structure.element import AlgebraElement
from sage.structure.richcmp import richcmp

class LazyLaurentSeries(AlgebraElement):

    def __init__(self, parent, coefficient=None, approximate_valuation=0, constant=None):
        AlgebraElement.__init__(self, parent)

        self._parent = parent
        self._ring = parent.base_ring()

        self._coefficient_function = coefficient
        self._approximate_valuation = approximate_valuation
        self._constant = constant

        self._cache = dict()

    def __nonzero__(self):
        """
        Return True if this power series is not equal to 0.

        EXAMPLES::

            sage: R.<q> = ZZ[[ ]]; R
            Laurent Series Ring in q over Integer Ring
            sage: f = 1 + 3*q + O(q^10)
            sage: f.is_zero()
            False
            sage: (0 + O(q^2)).is_zero()
            True
            sage: R(0).is_zero()
            True
            sage: (0 + O(q^1000)).is_zero()
            True
        """
        return self.valuation() == infinity

    def _repr_(self):
        """
        Return the string representation of this power series.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3))._repr_()
            't^2 + O(t^3)'

        ::

            sage: R.<t> = QQ[[]]
            sage: 1 / (1+2*t +O(t^5))
            1 - 2*t + 4*t^2 - 8*t^3 + 16*t^4 + O(t^5)

        ::

            sage: R.<t> = LaurentSeriesRing(QQ, sparse=True)
            sage: 1 / (1+2*t +O(t^5))
            1 - 2*t + 4*t^2 - 8*t^3 + 16*t^4 + O(t^5)
            sage: -13/2 * t^3  + 5*t^5 + O(t^10)
            -13/2*t^3 + 5*t^5 + O(t^10)
        """
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        X = self._parent.variable_name()

        try:
            n = self.valuation()
        except ValueError:
            n = self._approximate_valuation

        if self._constant is None:
            m = n + 7 # long enough
        elif self._constant[0] != 0:
            m = self._constant[1] + 3
        else:
            m = self._constant[1]

        s = ' '
        first = True
        while n < m:
            x = repr(self.coefficient(n))
            if x != '0':
                if not first:
                    s += ' + '
                if not atomic_repr and n > 0 and (x[1:].find('+') != -1 or x[1:].find('-') != -1):
                    x = '({})'.format(x)
                if n > 1 or n < 0:
                    var = '*%s^%s'%(X,n)
                elif n == 1:
                    var = '*%s'%X
                else:  # n == 0
                    var = ''
                s += '{}{}'.format(x,var)
                first = False
            n += 1

        s = s.replace(" + -", " - ").replace(" 1*"," ").replace(" -1*", " -")[1:]

        if len(s) == 0:
            s = '0'

        if self._constant is None or self._constant[1] > m or self._constant[0] != 0:
            s += ' + {}'.format('...')

        return s

    def coefficient(self, n):
        if self._approximate_valuation == infinity:
            return self._ring.zero()
        elif n < self._approximate_valuation:
            return self._ring.zero()
        elif self._constant is not None and n >= self._constant[1]:
            return self._constant[0]

        cache = self._cache
        try:
            c = cache[n]
        except:
            c = self._coefficient_function(self, n)
            cache[n] = c

        return c

    def _mul_(self, other):

        def mul(s, n):
            c = self._ring.zero()
            for k in range(self._approximate_valuation, n - other._approximate_valuation + 1):
                c += self.coefficient(k) * other.coefficient(n-k)
            return c

        a = self._approximate_valuation + other._approximate_valuation

        c = None
        if self._constant is not None and other._constant is not None:
            if self._constant[0] == 0 and other._constant[0] == 0:
                c = (self._constant[0], self._constant[1] + other._constant[1] - 1)

        return LazyLaurentSeries(self.parent(), coefficient=mul,
                               approximate_valuation=a, constant=c)

    def _add_(self, other):

        def add(s, n):
            return self.coefficient(n) + other.coefficient(n)

        a = min(self._approximate_valuation, other._approximate_valuation)

        c = None
        if self._constant is not None and other._constant is not None:
            c = (self._constant[0] + other._constant[0],
                 max(self._constant[1], other._constant[1]))

        return LazyLaurentSeries(self.parent(), coefficient=add,
                               approximate_valuation=a, constant=c)

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*y + 1/x); a                           # indirect doctest
            (-x^2/(8*x^5 + x^2 + 1/2))*y + (2*x^3 + x)/(16*x^5 + 2*x^2 + 1)
            sage: a*(2*y + 1/x)
            1
        """
        try:
            v = self.valuation()
        except ValueError:
            raise ZeroDivisionError('Might be zero')

        a0 = self.coefficient(v)
        zero = self._ring.zero()

        def inv(s, n):
            if n == -v:
                return ~a0
            c = zero
            for k in range(-v, n):
                c += s.coefficient(k) * self.coefficient(n + v - k)
            return -c / a0

        return LazyLaurentSeries(self.parent(), coefficient=inv,
                               approximate_valuation=-v, constant=None)

    def valuation(self, verify=False):
        """
        Return the valuation of the series.

        If not decidable from the known coefficients, raise ``ValueError`` unless
        ``verify`` is ``True``
        """
        if self._constant is None:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n in cache:
                    if cache[n] != 0:
                        return n
                    n += 1
                elif verify:
                    if self.coefficient(n) != 0:
                        return n
                    n += 1
                else:
                    raise ValueError('Unknown')

        n = self._approximate_valuation
        m = self._constant[1]
        while n <= m:
            if self.coefficient(n) != 0:
                return n
            n += 1
        return infinity






