r"""
The Pre-Lie Operad

AUTHORS:

- Florent Hivert, Frederic Chapoton (2011)
"""

#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@lri.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.functions.other import factorial
from sage.categories.all import (
    OperadsWithBasis, GradedHopfAlgebrasWithBasis, tensor )
from sage.combinat.free_module import (
    CombinatorialFreeModule,
    CombinatorialFreeModuleElement )
from sage.combinat.rooted_tree import LabelledRootedTrees
from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.split_nk import SplitNK
from sage.operads.free_operad import FreeOperad

# PreLie operad : rooted trees

# bases combinatoires : les etiquettes d'un arbre
# x.labels()

# bases combinatoires : la racine d'un arbre
# x.label()

# bases combinatoires : la largeur d'un arbre
# len(x)

# bases combinatoires : les sous-arbres t=B(t0,..,tk)
# x[0],x[1],...

class PreLieOperad(CombinatorialFreeModule):
    r"""
    An example of an operad with basis: the PreLie operad

    This class illustrates a minimal implementation of an operad with basis.
    """

    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = PreLieOperad(QQ); A
            The Pre-Lie operad over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, LabelledRootedTrees(),
            latex_prefix = "",
            monomial_cmp = LabelledRootedTrees().graded_cmp(),
            category = (OperadsWithBasis(R), GradedHopfAlgebrasWithBasis(R)))

    def _repr_(self):
        """
        EXAMPLES::

            sage: PreLieOperad(QQ) #indirect doctest
            The Pre-Lie operad over Rational Field
        """
        return "The Pre-Lie operad over %s"%(self.base_ring())

    def species(self):
        """
        The species of rooted trees

        EXAMPLES::

            sage: f=PreLieOperad(QQ).species()
            sage: f.generating_series().coefficients(5)
            [0, 1, 1, 3/2, 8/3]
        """
        from sage.combinat.species.library import (
            SingletonSpecies, SetSpecies, CombinatorialSpecies)
        X=SingletonSpecies()
        E=SetSpecies()
        R=CombinatorialSpecies()
        R.define(X*E(R))
        return R

    @cached_method
    def one_basis(self, letter = "@"):
        """
        Returns the tree with one vertex, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.one_basis("a")
            a[]
            sage: A.one("a")
            B[a[]]
        """
        return self.basis().keys()([],label=letter)

    def degree_on_basis(self,t):
        """
        Returns the degree of a rooted tree in the Pre-Lie operad.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: A.degree_on_basis(LT([LT([],'b')], label='a'))
            2
        """
        return t.node_number()

    def labelling_on_basis(self,t):
        """
        Put canonical labels on a tree in the Pre-Lie operad.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: A.labelling_on_basis(LT([LT([],'b')], label='a'))
            B[1[2[]]]
        """
        return self.basis()[t.canonical_labelling()]

    def unlabelling_on_basis(self,t):
        """
        Removes the labels of a tree in the Pre-Lie operad.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: A.unlabelling_on_basis(LT([LT([],'b')], label='a'))
            B[[[]]]
        """
        return self.basis()[t.shape()]

    def one(self, letter = "@"):
        """
        I overload the one of the operad so that it can also serve as the one
        of the Hopf algebra

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.one("a")
            B[a[]]
            sage: A.one()
            B[@[]]
        """
        # super(PreLieOperad, self).one(letter) bug of cached_method
        return self.monomial(self.one_basis(letter))

    def some_elements(self):
        """
        Returns some elements of the operad

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.some_elements()
            [B[@[]], B[@[b[]]], B[@[b[@[b[]]]]] + B[@[b[], @[b[]]]], B[@[b[@[]]]] + B[@[@[], b[]]]]
        """
        x = self.one() < self.one("b")
        return [self.one(), x, x < x, x < self.one()]

    # procedures de composition
    # nota bene : these algorithms only work if all labels of y are distinct !

    def singleGraft(self,y,x,graftingFunction, path_prefix = ()):
        """
        This returns the rooted tree obtained from a rooted tree y, a
        rooted tree x and a grafting function graftingFunction from
        range(len(x)) to the set of paths in y.

        .. WARNING:

        these algorithms only work if all labels of y are distinct !

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.singleGraft(y,x,dict([[0,(0,)]]))
            a[b[d[]]]
            sage: t=LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s=LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: A.singleGraft(t,s,dict([[0,(0,)],[1,(1,)]]))
            a[b[d[]], c[e[]]]

        """
        y1=self.basis().keys()([self.singleGraft(suby,x,graftingFunction,
                                                 path_prefix+(i,))
                                for i,suby in enumerate(y)],label=y.label())
        with y1.clone() as y2:
            for k in range(len(x)):
                if graftingFunction[k]==path_prefix:
                    y2.append(x[k])
        return y2

    def composition_on_basis_in_root(self,x,y):
        """
        This returns a list of rooted trees obtained from a rooted tree x
        and a rooted tree y by the composition x o_i y where i is the root
        of x.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.composition_on_basis_in_root(x,y)
            [a[b[], d[]], a[b[d[]]]]
        """
        return [self.singleGraft(y,x,graftingFunction)
                for graftingFunction in
                CartesianProduct(*([list(y.paths())]*len(x)))]

    def composition_on_basis_list(self,x,y,i):
        """
        This returns a list of rooted trees obtained from a rooted tree x
        and a rooted tree y by the composition x o_i y. i should be a label of x.
        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.composition_on_basis_list(x,y,'c')
            [a[b[], d[]], a[b[d[]]]]
        """
        if not(i in x.labels()):
            raise ValueError, "The composition index is not present in the first argument."
        elif x.label()==i:
            return self.composition_on_basis_in_root(x,y)
        else:
            resu=[]
            for j in range(len(x)):
                if i in x[j].labels():
                    for sx in self.composition_on_basis_list(x[j],y,i):
                        with x.clone() as x1:
                            x1[j]=sx;
                        resu.append(x1);
            return(resu)

    def composition_on_basis(self,x, y, i):
        """
        This computes the composition x o_i y as a sum of rooted trees,
        for rooted trees x and y. i should be a label of x.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.composition_on_basis(x,y,'a')
            Traceback (most recent call last):
            ...
            ValueError: The composition index is not present.
            sage: A.composition(A(x), A(y), 'd')
            B[c[a[b[]]]]
            sage: A.composition(A(x), A(y), 'c')
            B[a[b[], d[]]] + B[a[b[d[]]]]
        """
        if not(i in x.labels()):
            raise ValueError, "The composition index is not present."
        else:
            return sum(self.basis()[t] for t in self.composition_on_basis_list(x,y,i))

    def product_on_basis(self, x, y):
        """
        The product of the Hopf algebra is the composition at the root.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = A(LT([LT([],'b')], label='@'))
            sage: x = A(LT([LT([],'d')], label='@'))
            sage: x*y          # indirect doctest
            B[@[b[], d[]]] + B[@[b[d[]]]]
        """
        return sum(self.monomial(t) for t in
                   self.composition_on_basis_in_root(x, y))


    def coproduct_on_basis(self, x):
        """
        The coproduct of the Hopf algebra is un-shuffle.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = A(LT([LT([],'a'), LT([],'b'), LT([],'c')], label='@'))
            sage: y.coproduct()    # random # indirect doctest
            B[@[]] # B[@[a[], b[], c[]]] + B[@[a[]]] # B[@[b[], c[]]] + B[@[b[], c[]]] # B[@[a[]]] + B[@[a[], b[]]] # B[@[c[]]] + B[@[a[], b[], c[]]] # B[@[]] + B[@[b[]]] # B[@[a[], c[]]] + B[@[a[], c[]]] # B[@[b[]]] + B[@[c[]]] # B[@[a[], b[]]]
        """
        LT = self.basis().keys()
        def list_to_tree(lst):
            return self.monomial(LT([x[i] for i in lst], label = '@'))
        return sum(tensor([list_to_tree(u), list_to_tree(v)])
                   for k in range(len(x)+1)
                   for (u,v) in SplitNK(len(x), k))

    def operad_generator_basis(self, fstlabel = 0, sndlabel = 1):
        """
        Returns the generator of self

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.operad_generator_basis()
            0[1[]]
        """
        LT = self.basis().keys()
        return LT([LT([],label=sndlabel)], label=fstlabel)


    def operad_generator(self, fstlabel = 0, sndlabel = 1):
        """
        Returns the generator of self

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.operad_generator()
            B[0[1[]]]
        """
        return self.monomial(self.operad_generator_basis(fstlabel, sndlabel))

    def operad_generators(self):
        """
        EXAMPLES::

        sage: PreLieOperad(QQ).operad_generators()
        Finite family {'pre_Lie_product': B[1[2[]]]}
        """
        from sage.sets.family import Family
        return Family(dict([("pre_Lie_product",
                             self.operad_generator(1,2))]))

    def pre_Lie_product(self, x, y):
        """
        This computes the pre-Lie product.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b')], label='a'))
            sage: y = A(LT([LT([],'d')], label='c'))
            sage: A.pre_Lie_product(x, y)
            B[a[b[c[d[]]]]] + B[a[b[], c[d[]]]]
        """
        return self.operad_generator().compose(x, 0).compose(y, 1)

    # def normalize_on_basis(self, tr):
    #     """
    #     EXAMPLES::

    #         sage: A = PreLieOperad(QQ)
    #         sage: LT = A.basis().keys()
    #         sage: x = LT([LT([],'b'), LT([],'c')], label='a')
    #         sage: A.normalize_on_basis(x)
    #         B[a[b[], c[]]]
    #         sage: x = LT([LT([],'c'), LT([],'b')], label='a')
    #         sage: A.normalize_on_basis(x)
    #         B[a[b[], c[]]]
    #     """
    #     return self.basis()[tr.cayley_normalize()]

    # @lazy_attribute
    # def normalize(self):
    #     """
    #     EXAMPLES::

    #         sage: A = PreLieOperad(QQ)
    #         sage: LT = A.basis().keys()
    #         sage: x = A(LT([LT([],'b'), LT([],'c')], label='a'))
    #         sage: A.normalize(x)
    #         B[a[b[], c[]]]
    #         sage: x = A(LT([LT([],'c'), LT([],'b')], label='a'))
    #         sage: A.normalize(x)
    #         B[a[b[], c[]]]
    #     """
    #     return self.module_morphism(self.normalize_on_basis, codomain = self)

    def operad_morphism_on_basis(self,t,cod):
        """
        Defines a morphism from the PreLie operad to the target operad

        the target operad has to possess a method called pre_Lie_product

        the argument should not have repeated labels

        EXAMPLES::

            sage: PL=PreLieOperad(QQ)
            sage: RT=PL.basis().keys()
            sage: tr=RT([],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL)==PL(tr)
            True
            sage: tr=RT([RT([],label='b')],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL)==PL(tr)
            True
            sage: tr=RT([RT([],label='b'),RT([],label='c')],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL)==PL(tr)
            True
            sage: tr=RT([RT([RT([],label='c')],label='b'),RT([],label='e'),RT([],label='d')],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL)==PL(tr)
            True
        """
        targetPreLieProduct=cod.pre_Lie_product
        width = len(t)
        if width==0 :
            return cod.one(t.label())
        elif width==1 :
            return targetPreLieProduct(cod.one(t.label()),self.operad_morphism_on_basis(t[0],cod))
        else:
            t_red = self.basis().keys()([t[i] for i in range(width-1)],label = t.label())
            somme1=targetPreLieProduct(self.operad_morphism_on_basis(t_red,cod),self.operad_morphism_on_basis(t[width-1],cod))
            somme2=cod(0)
            for j in range(width-1):
                with t_red.clone() as tj:
                    tj[j]=self.basis().keys()([],label=t_red[j].label())
                    somme2 += cod.composition(self.operad_morphism_on_basis(tj,cod),targetPreLieProduct(self.operad_morphism_on_basis(t_red[j],cod),self.operad_morphism_on_basis(t[width-1],cod)),t_red[j].label())
            return somme1-somme2

    def corolla(self,n,x,y,N=10):
        r"""
        evaluate the corolla with `n` leaves, with `x` in the root and `y`
        in the leaves

        The result is computed up to order `N` (included).

        INPUT:

        - `n`: an integer

        - `x`, `y`: two tree-indexed series

        - `N` (optional): an integer (default: 10)

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL  = PreLieOperad(ZZ)
            sage: PLT = PL.basis().keys()

            sage: a = PL(PLT([PLT([],label=None)],label=None))
            sage: PL.corolla(1,a,a)
            B[None[None[None[None[]]]]] + B[None[None[], None[None[]]]]

            sage: b = PL.one()
            sage: PL.corolla(3,b,b,4)
            B[@[@[], @[], @[]]]

            sage: PL.corolla(2,a,b,4)
            B[None[None[@[], @[]]]] + B[None[@[], @[], None[]]] + 2*B[None[@[], None[@[]]]]
        """
        PL = x.parent()
        if n+1 > N:
            return PL(0)
        else:
            xx = x.truncate(N+1-n)
            yy = y.truncate(N+1-n)
            PLT = PL.basis().keys()
            crln = PL(PLT([PLT([],i+1) for i in range(n)],0))
            resu = crln.compose(xx,0)
            for i in range(1,n+1):
                resu = PL.composition_truncated(resu,yy,i,N)
            return resu

    def sum_corolla_prelie(self,n,x,y,N=10):
        r"""
        evaluate the sum of all corollas with up to `n` leaves, with
        `x` in the root and `y` in the leaves

        The result is computed up to order `N` (included).

        INPUT:

        - `n`: an integer

        - `a`, `b`: two tree-indexed series

        - `N` (optional): an integer (default: 10)

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL = PreLieOperad(QQ)
            sage: a = PL.one("@")
            sage: b = PL.one("O")
            sage: PL.sum_corolla_prelie(3,a,b,3)
            B[@[]] + B[@[O[]]] + 1/2*B[@[O[], O[]]]
            sage: PL.sum_corolla_prelie(3,a,b,10)
            B[@[]] + B[@[O[]]] + 1/2*B[@[O[], O[]]] + 1/6*B[@[O[], O[], O[]]]
        """
        PL = x.parent()
        return x + sum(PL.corolla(i,x,y,N)/factorial(i) for i in range(1,n+1))

    def diese_product(self,x,y,N=10):
        r"""
        evaluate the `\#` product of `x` and `y` up to order `N`

        The `\#` product is an associative product on tree-indexed `\#`
            series. This is similar to the Baker-Campbell-Hausdorff
            formula.

        The tree-indexed series `x` is in the root. The tree-indexed
        series `y` is in the leaves.

        The result is computed up to order `N` (included).

        INPUT:

        - `x`, `y`: two tree-indexed series

        - `N` (optional): an integer (default: 10)

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL  = PreLieOperad(QQ)
            sage: PLT = PL.basis().keys()

            sage: a = PL(PLT([PLT([],label='@')],label='@'))
            sage: PL.diese_product(a,a,2)
            2*B[@[@[]]]
            sage: PL.diese_product(a,a,4)
            2*B[@[@[]]] + B[@[@[@[@[]]]]] + B[@[@[], @[@[]]]]
            sage: b = PL.one("O")
            sage: PL.diese_product(b,a,5)
            B[O[]] + B[@[@[]]] + B[O[@[@[]]]] + 1/2*B[O[@[@[]], @[@[]]]]
        """
        PL = x.parent()
        yy = y.truncate(N+1)
        xx = x.truncate(N+1)
        return yy + xx + sum(PL.corolla(i,xx,yy,N)/factorial(i) for i in range(1,N))

    def inverse_diese(self,x,N=10):
        r"""
        computes the inverse of `x` for the # product up to order `N`

        The result is computed up to order `N` (included).

        INPUT:

        - `x`: a tree-indexed series

        - `N` (optional): an integer (default: 10)

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL = PreLieOperad(QQ)
            sage: a = PL.one()
            sage: PL.inverse_diese(a,2)
            -B[@[]] + B[@[@[]]]
        """
        PL = x.parent()
        inverse = PL.zero()
        for i in range(1,N+1):
            inverse = inverse - PL.diese_product(x,inverse,i)
        return inverse

    class Element(CombinatorialFreeModuleElement):
        def __lt__(self, other):
            r"""
            Shortcut for the prelie product

            EXAMPLES::

                sage: A = PreLieOperad(QQ)
                sage: A.one("x") < A.one("x")
                B[x[x[]]]

            .. warning:: Due to priority rules for operators, term must be put
                within parentheses inside sum, product... For example you must
                write::

                    sage: a = A.one('a'); b = A.one('b'); c = A.one('c')
                    sage: (a<b) + c
                    B[c[]] + B[a[b[]]]

                Indeed ``a<b + c`` is understood as ``a< (b + c)``

                    sage: (a<b + c) - (a < (b + c))
                    0
            """
            parent = self.parent()
            assert(parent == other.parent())
            return parent.pre_Lie_product(self, other)

