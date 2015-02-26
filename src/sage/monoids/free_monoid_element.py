"""
Monoid Elements

AUTHORS:

- David Kohel (2005-09-29)

Elements of free monoids are represented internally as lists of
pairs of integers.
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.structure.element import MonoidElement

def is_FreeMonoidElement(x):
    return isinstance(x, FreeMonoidElement)

class FreeMonoidElement(MonoidElement):
    """
    Element of a free monoid.

    EXAMPLES::

        sage: a = FreeMonoid(5, 'a').gens()
        sage: x = a[0]*a[1]*a[4]**3
        sage: x**3
        a0*a1*a4^3*a0*a1*a4^3*a0*a1*a4^3
        sage: x**0
        1
        sage: x**(-1)
        Traceback (most recent call last):
        ...
        TypeError: bad operand type for unary ~: 'FreeMonoid_class_with_category.element_class'
    """
    def __init__(self, F, x, check=True):
        """
        Create the element `x` of the FreeMonoid `F`.

        This should typically be called by a FreeMonoid.
        """
        MonoidElement.__init__(self, F)
        if isinstance(x, (int, long, Integer)):
            if x == 1:
                self._element_list = []
            else:
                raise TypeError("Argument x (= %s) is of the wrong type."%x)
        elif isinstance(x, list):
            if check:
                x2 = []
                for v in x:
                    if not isinstance(v, tuple) and len(v) == 2:
                        raise TypeError("x (= %s) must be a list of 2-tuples or 1."%x)
                    if not (isinstance(v[0], (int,long,Integer)) and \
                            isinstance(v[1], (int,long,Integer))):
                        raise TypeError("x (= %s) must be a list of 2-tuples of integers or 1."%x)
                    if len(x2) > 0 and v[0] == x2[len(x2)-1][0]:
                        x2[len(x2)-1] = (v[0], v[1]+x2[len(x2)-1][1])
                    else:
                        x2.append(v)
                self._element_list = x2
            else:
                self._element_list = list(x)  # make copy, so user can't accidentally change monoid.

        else:
            # TODO: should have some other checks here...
            raise TypeError("Argument x (= %s) is of the wrong type."%x)

    def __iter__(self):
        """
        Returns an iterator which yields tuples of variable and exponent.

        EXAMPLES::

            sage: a = FreeMonoid(5, 'a').gens()
            sage: list(a[0]*a[1]*a[4]**3*a[0])
            [(a0, 1), (a1, 1), (a4, 3), (a0, 1)]
        """
        gens=self.parent().gens()
        return ((gens[index], exponent) \
                for (index, exponent) in self._element_list)

##     def __cmp__(left, right):
##         """
##         Compare two free monoid elements with the same parents.

##         The ordering is the one on the underlying sorted list of
##         (monomial,coefficients) pairs.

##         EXAMPLES:
##             sage: R.<x,y> = FreeMonoid(2)
##             sage: x < y
##             True
##             sage: x * y < y * x
##             True
##             sage: x * y * x^2 < x * y * x^3
##             True
##         """
##         return cmp(left._element_list, right._element_list)

    def _repr_(self):
        s = ""
        v = self._element_list
        x = self.parent().variable_names()
        for i in range(len(v)):
            if len(s) > 0: s += "*"
            g = x[int(v[i][0])]
            e = v[i][1]
            if e == 1:
                s += "%s"%g
            else:
                s += "%s^%s"%(g,e)
        if len(s) == 0: s = "1"
        return s

    def _latex_(self):
        r"""
        Return latex representation of self.

        EXAMPLES::

            sage: F = FreeMonoid(3, 'a')
            sage: z = F([(0,5),(1,2),(0,10),(0,2),(1,2)])
            sage: z._latex_()
            'a_{0}^{5}a_{1}^{2}a_{0}^{12}a_{1}^{2}'
            sage: F, (alpha,beta,gamma) = FreeMonoid(3, 'alpha,beta,gamma').objgens()
            sage: latex(alpha*beta*gamma)
            \alpha\beta\gamma
        """
        s = ""
        v = self._element_list
        x = self.parent().latex_variable_names()
        for i in range(len(v)):
            g = x[int(v[i][0])]
            e = v[i][1]
            if e == 1:
                s += "%s"%(g,)
            else:
                s += "%s^{%s}"%(g,e)
        if len(s) == 0: s = "1"
        return s

    def __call__(self, *x, **kwds):
        """
        EXAMPLES::

            sage: M.<x,y,z>=FreeMonoid(3)
            sage: (x*y).subs(x=1,y=2,z=14)
            2
            sage: (x*y).subs({x:z,y:z})
            z^2
            sage: M1=MatrixSpace(ZZ,1,2)
            sage: M2=MatrixSpace(ZZ,2,1)
            sage: (x*y).subs({x:M1([1,2]),y:M2([3,4])})
            [11]

            sage: M.<x,y> = FreeMonoid(2)
            sage: (x*y).substitute(x=1)
            y

            sage: M.<a> = FreeMonoid(1)
            sage: a.substitute(a=5)
            5

        AUTHORS:

        - Joel B. Mohler (2007-10-27)
        """
        if kwds and x:
            raise ValueError("must not specify both a keyword and positional argument")

        P = self.parent()

        if kwds:
            x = self.gens()
            gens_dict = dict([(name, i) for i, name in enumerate(P.variable_names())])
            for key, value in kwds.iteritems():
                if key in gens_dict:
                    x[gens_dict[key]] = value

        if isinstance(x[0],tuple):
            x = x[0]

        if len(x) != self.parent().ngens():
            raise ValueError("must specify as many values as generators in parent")

        # I don't start with 0, because I don't want to preclude evaluation with
        #arbitrary objects (e.g. matrices) because of funny coercion.
        one = P.one()
        result = None
        for var_index, exponent in self._element_list:
            # Take further pains to ensure that non-square matrices are not exponentiated.
            replacement = x[var_index]
            if exponent > 1:
                c = replacement ** exponent
            elif exponent == 1:
                c = replacement
            else:
                c = one

            if result is None:
                result = c
            else:
                result *= c

        if result is None:
            return one

        return result

    def _mul_(self, y):
        """
        Multiply two elements ``self`` and ``y`` of the
        free monoid.

        EXAMPLES::

            sage: a = FreeMonoid(5, 'a').gens()
            sage: x = a[0] * a[1] * a[4]**3
            sage: y = a[4] * a[0] * a[1]
            sage: x*y
            a0*a1*a4^4*a0*a1
        """
        M = self.parent()
        z = M(1)
        x_elt = self._element_list
        y_elt = y._element_list
        if not x_elt:
            z._element_list = y_elt
        elif not y_elt:
            z._element_list = x_elt
        else:
            k = len(x_elt)-1
            if x_elt[k][0] != y_elt[0][0]:
                z._element_list = x_elt + y_elt
            else:
                m = (y_elt[0][0], x_elt[k][1]+y_elt[0][1])
                z._element_list = x_elt[:k] + [ m ] + y_elt[1:]
        return z

    def __len__(self):
        """
        Return the degree of the monoid element ``self``, where each
        generator of the free monoid is given degree `1`.

        For example, the length of the identity is `0`, and the
        length of `x_0^2x_1` is `3`.

        EXAMPLES::

            sage: F = FreeMonoid(3, 'a')
            sage: z = F(1)
            sage: len(z)
            0
            sage: a = F.gens()
            sage: len(a[0]**2 * a[1])
            3
        """
        s = 0
        for x in self._element_list:
            s += x[1]
        return s

    def __cmp__(self,y):
##         """
##         The comparison operator, defined via x = self:
##             x < y <=> x.__cmp__(y) == -1
##             x == y <=> x.__cmp__(y) == 0
##             x > y <=> x.__cmp__(y) == 1
##         It is not possible to use __cmp__ to define a
##         non-totally ordered poset.
##         Question: How can the operators <, >, ==, !=,
##         <=, and >= be defined for a general poset?
##         N.B. An equal operator __equal__ may or may not
##         have been introduced to define == and != but can
##         not be used in conjunction with __cmp__.
##        """
        if not isinstance(y,FreeMonoidElement) or y.parent() != self.parent():
            #raise TypeError, "Argument y (= %s) is of the wrong type."%y
            return 1
        n = len(self)
        m = len(y)
        if n < m:
            return -1
        elif m < n:
            return 1
        elif n == 0:
            return 0 # n = m = 0 hence x = y = 1
        x_elt = self._element_list
        y_elt = y._element_list
        for i in range(len(x_elt)):
            k = x_elt[i][0]
            l = y_elt[i][0]
            if k < l:
                return -1
            elif k > l:
                return 1
            e = x_elt[i][1]
            f = y_elt[i][1]
            if e < f:
                # x_elt is longer so compare next index
                if x_elt[i+1][0] < l:
                    return -1
                else:
                    return 1
            elif f < e:
                # y_elt is longer so compare next index
                if k < y_elt[i+1][0]:
                    return -1
                else:
                    return 1
        return 0 # x = self and y are equal


    def _acted_upon_(self, x, self_on_left):
        """
        Currently, returns the action of the integer 1 on this
        element.

        EXAMPLES::

            sage: M.<x,y,z>=FreeMonoid(3)
            sage: 1*x
            x
        """
        if x == 1:
            return self
        return None

    def to_word(self, alph=None):
        """
        Return ``self`` as a word.

        INPUT:

        - ``alph`` -- (optional) the alphabet which the result should
          be specified in

        EXAMPLES::

            sage: M.<x,y,z> = FreeMonoid(3)
            sage: a = x * x * y * x
            sage: w = a.to_word(); w
            word: xxyx
            sage: w.to_monoid_element() == a
            True

        .. SEEALSO::

            :meth:`to_list`
        """
        from sage.combinat.words.finite_word import Words
        gens = self.parent().gens()
        if alph is None:
            alph = gens
        alph = map(str, alph)
        W = Words(alph)
        return W(sum([ [alph[gens.index(i[0])]] * i[1] for i in list(self) ], []))

    def to_list(self, indices=False):
        """
        Return ``self`` as a list of generators.

        If ``self`` equals `x_{i_1} x_{i_2} \cdots x_{i_n}`, with
        `x_{i_1}, x_{i_2}, \ldots, x_{i_n}` being some of the
        generators of the free monoid, then this method returns
        the list `[x_{i_1}, x_{i_2}, \ldots, x_{i_n}]`.

        If the optional argument ``indices`` is set to ``True``,
        then the list `[i_1, i_2, \ldots, i_n]` is returned instead.

        EXAMPLES::

            sage: M.<x,y,z> = FreeMonoid(3)
            sage: a = x * x * y * x
            sage: w = a.to_list(); w
            [x, x, y, x]
            sage: M.prod(w) == a
            True
            sage: w = a.to_list(indices=True); w
            [0, 0, 1, 0]
            sage: a = M.one()
            sage: a.to_list()
            []

        .. SEEALSO::

            :meth:`to_word`
        """
        if not indices:
            return sum( ([i[0]] * i[1] for i in list(self)), [])
        gens = self.parent().gens()
        return sum( ([gens.index(i[0])] * i[1] for i in list(self)), [])

