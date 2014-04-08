from sage.misc.cachefunc import cached_method
from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ordered_tree import LabelledOrderedTrees
# NAP (Non Assocative Permutative) operad : based on rooted trees

class NapOperad(CombinatorialFreeModule):
    r"""
    An example of an operad with basis: the Nap operad

    This class illustrates a minimal implementation of an operad with basis.
    """

    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = NapOperad(QQ); A
            An example of an operad with basis: the Nap operad over Rational Field
            sage: TestSuite(A).run()

        """
        CombinatorialFreeModule.__init__(self, R, LabelledOrderedTrees(), category = OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: NapOperad(QQ)         # indirect doctest
            An example of an operad with basis: the Nap operad over Rational Field
        """
        return "An example of an operad with basis: the Nap operad over %s"%(self.base_ring())

    def species(self):
        """
        The species of rooted trees

        EXAMPLES::

            sage: f=NapOperad(QQ).species()
            sage: f.generating_series().coefficients(5)
            [0, 1, 1, 3/2, 8/3]
        """
        from sage.combinat.species.library import *
        X=SingletonSpecies()
        E=SetSpecies()
        R=CombinatorialSpecies()
        R.define(X*E(R))
        return R

    @cached_method
    def one_basis(self,letter):
        """
        Returns the tree with one vertex, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = NapOperad(QQ)
            sage: A.one_basis("a")
            a[]
        """
        return self.basis().keys()([],label=letter)

    def composition_on_basis_in_root(self,x,y):
        """
        This returns a rooted tree
        obtained from a rooted tree x and a rooted tree y by the
        composition x o_i y where i is the root of x.

        EXAMPLES::

            sage: NAP=NapOperad(QQ)
            sage: LT=NAP.basis().keys()
            sage: t=LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s=LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: NAP.composition_on_basis_in_root(t,s)
            f[b[], c[], d[], e[]]

        TESTS::

            sage: toto = LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: titi = LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: test1=NAP.composition_on_basis_in_root(toto,titi)

        """
        return self.basis().keys()(list(x)+list(y),label=y.label())

    def composition_on_basis(self,x,y,i):
        """
        This returns a rooted tree obtained from a rooted tree x
        and a rooted tree y by the composition x i y. i should be a label of x.

        EXAMPLES::

            sage: NAP=NapOperad(QQ)
            sage: LT=NAP.basis().keys()
            sage: t=LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s=LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: NAP.composition_on_basis(t,s,"b")
            a[f[d[], e[]], c[]]
        """
        if not(i in x.labels()):
            return "The composition index is not present in the first argument."
        elif x.label()==i:
            return self.composition_on_basis_in_root(x,y)
        else:
            j=[k for k in range(len(x)) if i in x[k].labels()][0]
            with x.clone() as x1:
                x1[j]=self.composition_on_basis(x[j],y,i);
            return x1
