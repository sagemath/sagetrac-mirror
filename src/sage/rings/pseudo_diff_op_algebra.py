r"""
Pseudo-differential operators

An implementation of pseudo-differential operators. These are
formal expansions of the form $\sum_n=-N^\infty u_{-n}(x) \partial_x^{-n}$.
Multiplication is such that $[\partial_x, u(x)] =\frac{\partial u}{\partial x}$.

AUTHORS:

    - Timo Kluck (2010-11-15): initial version

EXAMPLES::

We first create a pseudo-differential operator algebra with coefficients in the symbolic ring::
    sage: x=var('x')
    sage: from sage.rings.pseudo_diff_op_algebra import *
    sage: A.<d> = PseudoDiffOpAlgebra(SR,x); A
    Algebra of pseudo-differential operators over Symbolic Ring acting on x
    
The formal inverse of d behaves as it should::
    sage: d^-1 * d
    1
    sage: d * d^-1
    1

And we have the usual commutation relation::
    sage: u=function('u',x)
    sage: bracket(d,u)     # long time
    D[0](u)(x)

This illustrates associativity of the product::
    sage: d*(d^-1 * (u * d^-3)) # long time
    (u(x))d^(-3) + O(d^-9)

Note that [d^3, u] is not the same as [d, [d, [d, u]]]::
    sage: bracket(d,bracket(d,bracket(d,u))) # long time
    D[0, 0, 0](u)(x)
    sage: bracket(d^3,u)  # long time
    (3*D[0](u)(x))d^(2) + (3*D[0, 0](u)(x))d + D[0, 0, 0](u)(x)    

.. WARNING::

    This algebra is only interesting when used over the symbolic ring, because
    multiplication involves differentation. However, in the current
    implementation, this can be very slow. The reason for this is that we use
    LaurentSeriesRing, which checks whether any coefficient vanishes after
    every operation. For symbolic expressions, checking for vanishing spawns a
    Maxima process and this can be very slow.
"""

#*****************************************************************************
#       Copyright (C) 2010-2011 William Stein <wstein@gmail.com>
#                     2010-2011 Timo Kluck <tkluck@infty.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex
from sage.rings.integer import Integer
from sage.rings.ring import Algebra
from sage.structure.element import AlgebraElement
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.infinity import Infinity
from sage.rings.arith import binomial
from sage.calculus.functional import diff

def bracket(a,b):
    r"""
    Returns the commutator bracket between two elements.
    Equivalent to a*b - b*a.
    """
    return a*b - b*a

class PseudoDiffOpAlgebra(Algebra):
    r"""
    An algebra of pseudo-differential operators.
    
    Pseudo-differential operators  are formal expansions of the form
    $\sum_n=-N^\infty u_{-n}(x) \partial_x^{-n}$. Multiplication is such that
    $[\partial_x, u(x)] =\frac{\partial u}{\partial x}$.
    
    EXAMPLES::
        sage: x=var('x')
        sage: from sage.rings.pseudo_diff_op_algebra import *
        sage: A.<d> = PseudoDiffOpAlgebra(SR,x)
        sage: d^-1 * d
        1
        sage: d * d^-1
        1
        sage: u=function('u',x)
        sage: bracket(d,u) # long time
        D[0](u)(x)
        sage: d*(d^-1 * (u * d^-3)) # long time
        (u(x))d^(-3) + O(d^-9)

    AUTHORS:
        Timo Kluck (2010-11-15)    
    """
    
    def __init__(self,base,variable,prec=10,*kwds, **args):
        r"""
        Initializes an algebra of pseudo-differential operators.
        
        INPUT:
        
        -   ``base`` - the ring of coefficients
        
        -   ``variable`` - the variable with respect to which the operators
            should differentiate
            
        -   ``prec`` - the default order up to which to evaluate expansions, default is 10
        
        -   ``names`` - the name of the generator variable
        
        EXAMPLES::
        
            sage: x=var('x')
            sage: from sage.rings.pseudo_diff_op_algebra import *
            sage: A.<d> = PseudoDiffOpAlgebra(SR,x); A
            Algebra of pseudo-differential operators over Symbolic Ring acting on x

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        Algebra.__init__(self,base,*kwds,**args)
        self._var = variable
        self._prec = prec
        self.__symbolAlgebra = LaurentSeriesRing(base,'z_inverse')
        
    def prec(self):
        r"""
        Returns the default order up to which to evaluate expansions
        
        EXAMPLES::
        We first create a pseudo-differential operator algebra with coefficients in the symbolic ring::
            sage: x=var('x')
            sage: from sage.rings.pseudo_diff_op_algebra import *
            sage: A.<d> = PseudoDiffOpAlgebra(SR,x,20); A
            Algebra of pseudo-differential operators over Symbolic Ring acting on x

        This illustrates associativity of the product, but now with precision equal to 20::
            sage: d*(d^-1 * (u * d^-3)) # long time
            (u(x))d^(-3) + O(d^-19)

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        return self._prec
        
    def var(self):
        r"""
        Returns the variable upon which the differential operators act.

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        
    def _repr_(self):
        r"""
        Prints out an algebra of pseudo-differential operators.
        
        EXAMPLES::
        
            sage: x=var('x')
            sage: from sage.rings.pseudo_diff_op_algebra import *
            sage: A.<d> = PseudoDiffOpAlgebra(SR,x); A
            Algebra of pseudo-differential operators over Symbolic Ring acting on x

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        return "Algebra of pseudo-differential operators over %s acting on %s" % (str(self.__symbolAlgebra.base_ring()), str(self._var))
        
    def _latex_(self):
        r"""
        Returns a latex representation.
        
        EXAMPLES::
        
            sage: x=var('x')
            sage: from sage.rings.pseudo_diff_op_algebra import *
            sage: A.<d> = PseudoDiffOpAlgebra(SR,x); latex(A)
            \texttt{Algebra of pseudo-differential operators over } \text{SR} \texttt{ acting on } x

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        return r"\texttt{Algebra of pseudo-differential operators over } %s \texttt{ acting on } %s" % (latex(self.__symbolAlgebra.base_ring()), latex(self._var))
    
           
    def _coerce_map_from_(self, S):
        r"""
        A coerce map exists iff a coerce map to the symbol algebra exists.
        
        Note that the symbol algebra is not an accessible property, because the
        implementation with z_inverse as a generator is a bit messy. This means 
        that if the user wants to construct a pseudo-differential operator, s/he
        cannot do that by constructing a symbol and coerceing it. However, the
        generator is accessible via gen(0), so the user will only need coercion
        from the base ring to PseudoDiffOpAlgebra. Coercion via the symbol
        ring provides this.

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        return self.__symbolAlgebra.has_coerce_map_from(S)

    def _element_constructor_(self, x):
        r"""
        Returns a pseudo-differential operator with symbol x
        
        INPUT:
        
        -    ``x`` - the symbol of x as an element of the symbol algebra
        
        OUTPUT:
        
        An instance of PseudoDiffOpAlgebra_element with symbol equal to x

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """       
        
        x = self.__symbolAlgebra(x)
        return PseudoDiffOpAlgebraElement(self, x)

    def is_commutative(self):
        r"""
        Returns False, since for example $[\partial, u] = u^\prime \neq 0$.

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        return False
    
    def gen(self, n):
        r"""
        Returns a generator for the algebra. There is only one generator.
        
        EXAMPLES::
            sage: from sage.rings.pseudo_diff_op_algebra import *
            sage: A.<d> = PseudoDiffOpAlgebra(SR,x)
            sage: A.gen(0)
            d

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        
        if(n!=0):
            raise IndexError, "Generator n not defined."
        z_inverse = self.__symbolAlgebra.gen(0)
        return PseudoDiffOpAlgebraElement(parent=self,symbol=1/z_inverse,is_gen=True)
        
    def ngens(self):
        r"""
        Returns 1; the only generator is $\partial_{\mathrm{var}}$.

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """
        return 1


class PseudoDiffOpAlgebraElement(AlgebraElement):
    r"""
    A pseudo-differential operator.
    
    Pseudo-differential operators  are formal expansions of the form
    $\sum_n=-N^\infty u_{-n}(x) \partial_x^{-n}$. Multiplication is such that
    $[\partial_x, u(x)] =\frac{\partial u}{\partial x}$.
    
    EXAMPLES::
        sage: x=var('x')
        sage: from sage.rings.pseudo_diff_op_algebra import *
        sage: A.<d> = PseudoDiffOpAlgebra(SR,x)
        sage: type(d)
        <class 'sage.rings.pseudo_diff_op_algebra.PseudoDiffOpAlgebraElement'>
        sage: d^-1 * d
        1
        sage: d * d^-1
        1
        sage: u=function('u',x)
        sage: bracket(d,u)  # long time
        D[0](u)(x)
        sage: d*(d^-1 * (u * d^-3))  # long time
        (u(x))d^(-3) + O(d^-9)
        
    AUTHORS:
    
        Timo Kluck (2010-11-15)    
    """
    def __init__(self, parent, symbol,is_gen=False):
        r"""
        Initializes a pseudo-differential operator with symbol ``symbol``
        
        INPUT:
        
        -    ``symbol`` - the symbol of x as an element of the symbol algebra
        
        -   ``is_gen`` - whether this symbol is equal to z_inverse^(-1) (for
            internal use only!)
            
        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """

        AlgebraElement.__init__(self,parent)
        self._symbol = symbol
        self._is_gen = is_gen
        
    def _add_(left, right):
        r"""
        Adds two pseudo-differential operators. This is just addition of the
        symbols.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        
        new_symbol = left._symbol + right._symbol
        return PseudoDiffOpAlgebraElement(left.parent(), new_symbol)
        
    def _sub_(left, right):
        r"""
        Substracts two pseudo-differential operators. This is just addition of the
        symbols.

        AUTHORS:
            Timo Kluck (2010-11-15)    
        """

        new_symbol = left._symbol - right._symbol
        return PseudoDiffOpAlgebraElement(left.parent(), new_symbol)
    
    
    def _mul_(left, right):
        r"""
        Multiply two pseudo-differential operators.
        
        Multiplication is such that $[\partial_x, u(x)] =\frac{\partial u}{\partial x}$,
        and such that $\partial_x^{-1} \cdot \partial_x = 1$. More precisely,
        we define multiplication by the formula 
        \[
        a \partial_x^n * b \partial_x^m = \sum_{k=0}^{\infty} \choose{n}{k} a \frac{\partial b}{\partial x} \partial_x^{m+m-k}
        \]
        
        EXAMPLES::
            sage: x=var('x')
            sage: from sage.rings.pseudo_diff_op_algebra import *
            sage: A.<d> = PseudoDiffOpAlgebra(SR,x)
            sage: type(d)
            <class 'sage.rings.pseudo_diff_op_algebra.PseudoDiffOpAlgebraElement'>
            sage: d^-1 * d
            1
            sage: d * d^-1
            1
            sage: u=function('u',x)
            sage: bracket(d,u)   # long time
            D[0](u)(x)
            sage: d*(d^-1 * (u * d^-3))   # long time
            (u(x))d^(-3) + O(d^-9)

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        # some shortcuts
        #if(left.is_zero() or right.is_zero()):
        #    return PseudoDiffOpAlgebraElement(left.parent(), left._symbol.parent()(0))
        #if(right.is_gen()):
        #    return PseudoDiffOpAlgebraElement(left.parent(), left._symbol>>1)
        #if(right.is_monomial() and left.exponents()==[0]):
        #    return PseudoDiffOpAlgebraElement(left.parent(), left._symbol*right._symbol)
            
        
        z_inverse = left._symbol.parent().gen(0)
        
        new_prec = (left._symbol * right._symbol).prec()
        if(not new_prec < Infinity):
            new_prec=left.parent().prec()
        new_symbol=left._symbol.parent()(0)
        for a,n in zip(left.coefficients(), left.exponents()):
            for b,m in zip(right.coefficients(), right.exponents()):
                k=0
                while(k-m-n < new_prec):
                    new_symbol += binomial(n,k) * a * diff(b, left.parent().var(), k) * z_inverse**(k-m-n)
                    k+=1
                if binomial(n,k)!=0 and diff(b,left.parent().var(),k)!=0:
                    new_symbol = new_symbol.add_bigoh(new_prec)
                         
            
        return PseudoDiffOpAlgebraElement(left.parent(), new_symbol)

    def __pow__(self, exponent):
        r"""
        Raises a pseudo-differential operator to an integer power.
        
        .. NOTE::
        
            Raising any operator but $\partial_x$ to a negative power raises a
            value error. I am not sure if that is mathematically right.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        
        if not isinstance(exponent, Integer) or \
                isinstance(exponent, int):
                    try:
                        exponent = Integer(exponent)
                    except TypeError:
                        raise TypeError, "Pseudo-differential operators can only be raised to integer powers."
                        
        if exponent == 0:
            return self.parent()(1)
        
        if self._is_gen:
            z_inverse = self._symbol.parent().gen(0)
            return self.parent()(z_inverse**-exponent)    
        else:
            if exponent < 0:
                raise ValueError, "Pseudo-differential operators can not be raised to negative powers (except for the generator)."
            else:
                #inductively
                return self * self ** (exponent-1)
    
    def __call__(self, fn):
        """
        Apply this operator to a function.
        
        This can only be done for operators with only a positive part. Otherwise,
        a ValueError is raised.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        if self._symbol.degree() > 0:
            raise ValueError, "Can only apply the positive part of a pseudo-differential operator to a function."
        
        ret = 0
        var = self.parent().var()
        for c, e in zip(self.coefficients(), self.exponents()):
            ret += c * diff(fn, var, e)
        return ret
        
    def derivative(self, *args):
        r"""
        Returns the pseudo-differential operator obtained by differentiation
        of the coefficients.
        
        This is most useful when the base ring is the symbolic ring, or a
        power series ring, or similar.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        # differentiate a pseudo diff op's coefficients
        d = self.parent().gen(0)
        new_c = [diff(c,*args) for c in self.coefficients()]
        return sum(c * d ** e for c,e in zip(new_c,self.exponents()))
    
    def symbol(self):
        return self._symbol
        
    def pos_part(self):
        r"""
        Returns the positive part of the pseudo-differential operator.
        
        The positive part is the summands with non-negative exponents of $\partial_x$.
        These are the ones that can be applied to a function of x.
        
        .. NOTE::
            The current implementation requires a patch of laurent series,
            trac #10272

        AUTHORS:
        
            Timo Kluck (2010-11-15)
        """
        new_symbol = self._symbol.truncate(1)
        return PseudoDiffOpAlgebraElement(self.parent(), new_symbol)
        
    def neg_part(self):
        r"""
        Returns the negative part of the pseudo-differential operator.
        
        The negative part is the summands with negative exponents of $\partial_x$.
        It is given by self - self.pos_part()
        
        .. NOTE::
            The current implementation requires a patch of laurent series,
            trac #10272

        AUTHORS:
        
            Timo Kluck (2010-11-15)
        """
        return self - self.pos_part()
    
    def coefficients(self):
        r"""
        Returns a list of coefficients of self.
        
        To obtain the corresponding exponents, use self.exponents()
        
        See the documentation for LaurentSeries.coefficients() for details.

        AUTHORS:
        
            Timo Kluck (2010-11-15)
        """
        return self._symbol.coefficients()
    
    def exponents(self):
        r"""
        Returns a list of exponents of self.
        
        To obtain the corresponding coefficients, use self.coefficients()

        See the documentation for LaurentSeries.exponents() for details.
        
        AUTHORS:
        
            Timo Kluck (2010-11-15)
        """
        return [-e for e in self._symbol.exponents()]
    
    def dict(self):
        r"""
        Returns a dictionary with keys=>value equal to exponent=>coefficient.
        
        See the documentation for LaurentSeries.dict() for details.
        
        AUTHORS:
        
            Timo Kluck (2010-11-15)
        """
        return dict(zip(self.exponents(),self.coefficients()))
            
    def simplify(self):
        r"""
        Returns the pseudo-differential operator obtained by simplifying
        the coefficients.
        
        This is most useful when the base ring is the symbolic ring.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        from sage.calculus.functional import simplify
        new_coeff = [simplify(c) for c in self.coefficients()]
        z_inverse = self._symbol.parent().gen(0)
        new_symbol = self._symbol.parent()(0)
        new_symbol = new_symbol.add_bigoh(self._symbol.prec())
        for c,e in zip(new_coeff,self.exponents()):
            new_symbol += c * z_inverse**-e
        return PseudoDiffOpAlgebraElement(self.parent(), new_symbol)
    
    def expand(self): 
        r"""
        Returns the pseudo-differential operator obtained by expanding
        the coefficients.
        
        This is most useful when the base ring is the symbolic ring.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        from sage.calculus.functional import expand
        new_coeff = [expand(c) for c in self.coefficients()]
        z_inverse = self._symbol.parent().gen(0)
        new_symbol = self._symbol.parent()(0)
        new_symbol = new_symbol.add_bigoh(self._symbol.prec())
        for c,e in zip(new_coeff,self.exponents()):
            new_symbol += c * z_inverse**-e
        return PseudoDiffOpAlgebraElement(self.parent(), new_symbol)
    
    def __getitem__(self,n):
        r"""
        returns the nth coefficient of this pseudo-differential operator.
        """
        
        return self._symbol[-n]
    
    def is_gen(self):
        return self._is_gen
        
    def is_monomial(self):
        r"""
        Returns True when this pseudo-differential operator consists of a single
        term in the expansion in powers of $\partial_x$.
        
        NOTE:
            self.is_monomial() returns False if its precision is not Infinity.

        AUTHORS:
        
            Timo Kluck (2010-11-15)    
        """
        return len(self._symbol.coefficients())==1 and self._symbol.prec()==Infinity
        
    def is_zero(self):
        return len(self._symbol.coefficients())==0
    
    def degree(self):
        r"""
        Returns the maximum integer $n$ for which the coefficient of $\partial_x^n$ 
        is non-zero.
        """
        return -self._symbol.valuation()    

    def add_bigoh(self, prec):
        return PseudoDiffOpAlgebraElement(self.parent(), self._symbol.add_bigoh(-prec))

    def O(self,prec):
        return self.add_bigoh(self, prec)
        
    def prec(self):
        return -self._symbol.prec()
    
    def common_prec(self,other):
        return -self._symbol.common_prec(other._symbol)
    
    def _repr_(self):
        list = []
        for (c,e) in zip(self._symbol.coefficients(),self._symbol.exponents()):
            if(c==0):
                continue
            if(c==1):
                coeff=""
            elif(e!=0):
                coeff = "(%s)" % str(c)
            else:
                coeff = str(c)
            if(-e==1 or e==0):
                ex=""
            else:
                ex="^(%s)" % str(-e)
            if(e==0):
                v=""
            else:
                v="d"
            if(c==1 and e==0):
                coeff="1"
                
            list.append("%(coeff)s%(var)s%(ex)s" % {"coeff": coeff, "var": v,"ex": ex})
        
        prec=self._symbol.prec()    
        if(prec<Infinity):
            if(prec==-1):
                list.append("O(d)")
            elif(prec==0):
                list.append("O(1)")
            else:
                list.append("O(d^%d)" % -prec)
                
                
        if(len(list)==0):
            return "0"
        else:
            return " + ".join(list)
            
    def _latex_(self):
        list = []
        v=self.parent().var()._latex_()
        for (c,e) in zip(self._symbol.coefficients(),self._symbol.exponents()):
            if(c==0):
                continue
            if(c==1):
                coeff=""
            elif(e!=0):
                coeff = "\\left(%s\\right)" % c._latex_()
            else:
                coeff = c._latex_()
            if(-e==1 or e==0):
                ex=""
            else:
                ex="^{%s}" % str(-e)
            if(e==0):
                d=""
            else:
                d="\\partial_{%s}" % v
            if(c==1 and e==0):
                coeff="1"

            list.append("%(coeff)s%(d)s%(ex)s" % {"coeff": coeff, "d": d, "ex": ex})

        prec=self._symbol.prec()    
        if(prec<Infinity):
            if(prec==-1):
                list.append("\\mathcal{O}(%(var)s)" % {"var":v})
            elif(prec==0):
                list.append("\\mathcal{O}(1)")
            else:
                list.append("\\mathcal{O}(\\partial_{%(var)s}^{%(ex)s})" % {"var":v, "ex": -prec})               


        if(len(list)==0):
            return "0"
        else:
            return " + ".join(list)
