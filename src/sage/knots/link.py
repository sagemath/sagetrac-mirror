r"""
Link class
"""
#*****************************************************************************
#  Copyright (C) 2014
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.free_group import FreeGroupElement
from sage.groups.braid import Braid
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.groups.braid import BraidGroup
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.integer_mod import Mod

class Link:
    r"""
    The base class for Link, taking input in three formats namely Briadword, gauss_code, dt_code
    """
    def __init__(self, input = None, gauss_code = None, dt_code = None):
        if type(input) == Braid:
            self._braid = input
            self._gauss_code = None
            self._dt_code = None

        elif gauss_code != None:
            self._braid = None
            self._gauss_code = gauss_code
            self._dt_code = None

        elif dt_code != None:
            self._braid = None
            self._gauss_code = None
            self._dt_code = dt_code

        else:
            raise Exception("Invalid input")

    def braidword(self):
        r"""
        Returns the braidword

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Braidword representation of the INPUT

        EXAMPLES::
        """

        if self._braid != None:
            return list(self._braid.Tietze())

        elif self._gauss_code != None:
            return "Not implemented Error"

        elif self._dt_code != None:
            return "Not Implemented Error"

    def braid(self):
        r"""
        Returns the braid

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Braid representation of the INPUT

        EXAMPLES::
        """

        if self._braid != None:
            return self._braid

        elif self._gauss_code != None:
            return "Not implemented Error"

        elif self._dt_code != None:
            return "Not Implemented Error"

    def gauss_code(self):
        r"""
        Returns the gauss_code

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Gauss code representation of the INPUT

        EXAMPLES::
        """
        if self._gauss_code != None:
            return self._gauss_code

        elif self._braid != None:
            L = self.braid()
            B = L.parent()
            l = self.braidword()
            L = Link(B(l)).dt_code()
            gc = Link(dt_code = L).gauss_code()
            return gc

        elif self._dt_code != None:
            dt = self._dt_code
            gauss = []
            y = [None for i in range(2*len(dt))]
            x = [0 for i in range(2*len(dt))]
            for i in range(len(dt)):
                x[2*i] = 2*i + 1
                x[2*i + 1] = dt[i]
            for i in range(len(dt)):
                if x[2*i+1] > 0:
                    y[2*i+1] = 'under'
                    y[2*i] = 'over'
                elif x[2*i+1] < 0:
                    y[2*i+1] = 'over'
                    y[2*i] = 'under'
            for i in range(1,len(x)+1):
                for j in range(0,len(x)):
                    if abs(x[j]) == i:
                        if y[j] == 'under':
                            gauss.append(-(j//2 + 1))
                        elif y[j] == 'over':
                            gauss.append(j//2 + 1)
            return gauss

    def dt_code(self):
        r"""
        Returns the dt_code

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - DT Code representation of the INPUT

        EXAMPLES::
        """
        if self._dt_code != None:
            return self._dt_code

        elif self._braid != None:
            b = list(self._braid.Tietze())
            N = len(b)
            label = [0 for i in range(2*N)]
            string = 1
            next_label = 1
            type1 = 0
            crossing = 0
            while(next_label <= 2*N):
                string_found = 0
                for i in range(crossing, N):
                    if(abs(b[i]) == string or abs(b[i]) == string - 1):
                        string_found = 1
                        crossing = i
                        break
                if(string_found == 0):
                    for i in range(0,crossing):
                        if(abs(b[i]) == string or abs(b[i]) == string - 1):
                            string_found = 1
                            crossing = i
                            break
                if(label[2*crossing + next_label%2] == 1):
                    return "error"
                else:
                    label[2*crossing + next_label%2] =  next_label
                    next_label = next_label + 1
                if(type1 == 0):
                    if(b[crossing] < 0):
                        type1 = 1
                    else:
                        type1 = -1
                else:
                    type1 = -1 * type1
                    if((abs(b[crossing]) == string and b[crossing] * type1 > 0) or (abs(b[crossing]) != string and b[crossing] * type1 < 0)):
                        if(next_label%2 == 1):
                            label[2*crossing] = label[2*crossing] * -1
                if(abs(b[crossing]) == string):
                    string = string + 1
                else:
                    string = string - 1
                crossing = crossing + 1

            code = [0 for i in range(N)]
            for i in range(N):
                for j in range(N):
                    if label[2*j+1] == 2*i+1:
                        code[i] = label[2*j]
                        break
            return code

         #this needs to further worked upon, the signs are not taken in consideration yet
        elif self._gauss_code != None:
            gc = self._gauss_code
            l = [0 for i in range(len(gc))]
            for i in range(len(gc)):
                k = abs(gc[i])
                if l[2*(k-1)] == 0:
                    l[2*(k-1)] = i + 1
                else:
                    l[2*k-1] = i + 1
            l = [l[i] for i in range(len(gc)) if l[i]%2 == 0]
            return l

    def _braidwordcomponents(self):
        r"""
        Returns the braid components in an array

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Array containing the components is returned

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponents()
            [[-1, 1], [3, 3]]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponents()
            [[-1, 1, 1, 1], [3], [5, 7, 6]]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponents()
            [[-2, 1, 1], [4, 4], [6]]
        """
        b = self.braid()
        ml = list(b.Tietze())
        if ml == []:
            raise Exception("The braid remains the same with no components")
        else:
            l = list(set([abs(k) for k in ml]))
            missing1 = list(set(range(min(l),max(l)+1)) - set(l))
            if len(missing1) == 0:
                return [ml]
            else:
                missing = sorted(missing1)
                x = [[] for i in range(len(missing) + 1)]
                for i in range(len(missing)):
                    for j in range(len(ml)):
                        if(ml[j] != 0 and abs(ml[j]) < missing[i]):
                            x[i].append(ml[j])
                            ml[j] = 0
                        elif(ml[j] != 0 and abs(ml[j]) > missing[-1]):
                            x[-1].append(ml[j])
                            ml[j] = 0
                y2 = [x for x in x if x != []]
                return y2

    def _braidwordcomponentsvector(self):
        r"""
        From braidwordcomponents it is converted to a vector

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Vector containing non-zero values

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponentsvector()
            [-1, 1, 3, 3]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponentsvector()
            [-1, 1, 1, 1, 3, 5, 7, 6]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponentsvector()
            [-2, 1, 1, 4, 4, 6]
        """
        y2 = self._braidwordcomponents()
        if len(y2) == 1:
            return y2[0]
        else:
            y3 = []
            for i in range(len(y2)):
                y = y2[i]
                for j in range(len(y)):
                    y3.append(y[j])
            return y3

    def homology_generators(self):
        r"""
        Returns the homology generators

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - The homology generators relating to the braid word representation

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.homology_generators()
            [1, 0, 3]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.homology_generators()
            [1, 2, 3, 0, 0, 0, 0]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.homology_generators()
            [0, 2, 0, 4, 0]
        """
        x4 = self._braidwordcomponentsvector()
        hom_gen = []
        for j in range(len(x4)-1):
            a = abs(x4[j])
            for i in range(j+1, len(x4)):
                    if(a == abs(x4[i])):
                        hom_gen.append(i)
                        break
            else:
                hom_gen.append(0)
        return hom_gen

    def Seifert_Matrix(self):
        r"""
        Returns the Seifert Matrix

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Returns the Seifert Matrix of the link.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.Seifert_Matrix()
            [ 0  0]
            [ 0 -1]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.Seifert_Matrix()
            [ 0  0  0]
            [ 1 -1  0]
            [ 0  1 -1]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.Seifert_Matrix()
            [-1  0]
            [ 0 -1]
        """
        x5 = self._braidwordcomponentsvector()
        h = self.homology_generators()
        hl = len(h)
        A = matrix(ZZ, hl, hl)
        for i in range(hl):
            if h[i] != 0:
                for j in range(i,hl):
                        if i == j:
                            A[i,j] = -cmp((x5[i] + x5[h[i]]),0)
                        elif (h[i] > h[j]):
                            A[i,j] = 0
                            A[j,i] = 0
                        elif (h[i] <  j):
                            A[i,j] = 0
                            A[j,i] = 0
                        elif (h[i] == j):
                            if(x5[j] > 0):
                                A[i,j] = 0
                                A[j,i] = 1
                            else:
                                A[i,j] = -1
                                A[j,i] = 0
                        elif abs(abs(x5[i]) - abs(x5[j])) > 1:
                            A[i,j] =  0
                        elif (abs(x5[i]) - abs(x5[j]) == 1):
                            A[i,j] = 0
                            A[j,i] = -1
                        elif (abs(x5[j])- abs(x5[i]) == 1):
                            A[i,j] = 1
                            A[j,i] = 0
                        else: # for debugging
                            A[i,j] = 2
                            A[j,i] = 2
            else:
                for k in range(hl):
                    A[k,i] = 0
                    A[i,k] = 0
        k = []
        for i in range(hl):
                if h[i] == 0:
                    k.append(i)
        for i in reversed(k):
                A = A.delete_rows([i])
                A = A.delete_columns([i])
        return A


    def link_number(self):
        r"""
        Returns the link number

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Link number of the link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.link_number()
            4
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.link_number()
            5
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.link_number()
            1
        """
        p = self.braid().permutation()
        return len(p.to_cycles())

    def is_knot(self):
        r"""
        Returns true if the link is knot

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - True or False

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([1,3,1,-3]))
            sage: L.is_knot()
            False
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([1, 2, 3, 4, 5, 6]))
            sage: L.is_knot()
            True
        """
        if self.link_number() == 1:
            return True
        else:
            return False

    def genus(self):
        r"""
        Returns the genus of the link

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Genus of the Link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.genus()
            0
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.genus()
            0
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.genus()
            1
        """
        b = self.braidword()
        if b == []:
            return 0
        else:
            B = self.braid().parent()
            x = self._braidwordcomponents()
            q = []
            genus = 0
            s = [Link(B(x[i])).smallest_equivalent() for i in range(len(x))]
            t = [Link(B(s[i])).link_number() for i in range(len(s))]
            for i in range(len(s)):
                if s[i] == []:
                    s[i].append(-2)
            for i in range(len(s)):
                q1 = (abs(k)+1 for k in s[i])
                q2 = max(q1)
                q.append(q2)
            g = [((2 - t[i]) + len(x[i]) - q[i])/2 for i in range(len(x))]
            for i in range(len(g)):
                genus = genus + g[i]
            return genus

    def smallest_equivalent(self):
        r"""
        Returns the braidword

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Smallest equivalent of the given braid word representation.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(5)
            sage: L = link.Link(B([-2, 4, 2, 4]))
            sage: L.smallest_equivalent()
            [-1, 3, 1, 3]
            sage: L = link.Link(B([-1, 1]))
            sage: L.smallest_equivalent()
            []
        """
        b = list(self.braid().Tietze())
        if not b:
            return list(b)
        else:
            b1 = min([abs(k) for k in b])
            for i in range(len(b)):
                if b[i] > 0:
                    b[i] = b[i] - b1 + 1
                else:
                    b[i] = b[i] + b1 - 1
            return b

    def signature(self):
        r"""
        Returns the signature of the link

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Signature of the Link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.signature()
            -1
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.signature()
            -2
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.signature()
            -2
        """
        m = 2*(self.Seifert_Matrix() + self.Seifert_Matrix().transpose())
        e = m.eigenvalues()
        sum = 0
        s = []
        for i in range(len(e)):
            s.append(cmp(e[i],0))
            sum = sum + s[i]
        return sum

    def alexander_polynomial(self, var ='t'):
        r"""
        Returns the alexander polynomial of the link

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Alexander Polynomial of the Link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.alexander_polynomial()
            0
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.alexander_polynomial()
            t^2 - 2*t + 1
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.alexander_polynomial()
            t^2 - t + 1
        """
        R = PolynomialRing(ZZ, var)
        t = R.gen()
        m2 = self.Seifert_Matrix() - t* (self.Seifert_Matrix().transpose())
        return m2.determinant()

    def knot_determinant(self):
        r"""
        Returns the determinant of the knot

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Determinant of the Knot

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 2, 1, 2]))
            sage: L.knot_determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.knot_determinant()
            3
            sage: L = link.Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.knot_determinant()
            65
        """
        if self.is_knot() == True:
            a = self.alexander_polynomial()
            return abs(a(-1))
        else:
            return "is defined for a knot"

    def arf_invariant(self):
        r"""
        Returns the arf invariant only if the link is knot

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Arf invariant of knot

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 2, 1, 2]))
            sage: L.arf_invariant()
            0
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 3, 1, 2, 1, 4]))
            sage: L.arf_invariant()
            0
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.arf_invariant()
            1
        """

        if self.is_knot() == True:
            a = self.alexander_polynomial()
            if ((Mod(a(-1),8) == 1) or (Mod(a(-1),8) == 7)):
                return 0
            else:
                return 1
        else:
            raise Exception("Arf invariant is defined only for knots")
