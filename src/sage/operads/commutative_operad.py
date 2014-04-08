from sage.misc.cachefunc import cached_method
from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
class CommutativeOperad(CombinatorialFreeModule):
    r"""
    The Commutative operad
    """

    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = CommutativeOperad(QQ); A
            The Commutative operad over Rational Field
            sage: TestSuite(A).run()

        """
        CombinatorialFreeModule.__init__(self, R, Words(), category = OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: CommutativeOperad(QQ)       # indirect doctest
            The Commutative operad over Rational Field
        """
        return "The Commutative operad over %s"%(self.base_ring())

    def species(self):
        """
        The species of non-empty sets

        EXAMPLES::

            sage: f=CommutativeOperad(QQ).species()
            sage: f.generating_series().coefficients(5)
            [0, 1, 1/2, 1/6, 1/24]
        """
        from sage.combinat.species.library import *
        return SetSpecies().restricted(min=1)

    @cached_method
    def one_basis(self,letter):
        """
        Returns the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def degree_on_basis(self,t):
        """
        Returns the degree of a word `t` in the Commutative operad.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.degree_on_basis(m)
            4
        """
        return t.length()

    def map_labels(self,t,f):
        """
        Maps the function `f` on the word `t`.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.map_labels(m,lambda u:u)
            word: 1234
        """
        return self.basis().keys()(sorted([f(u) for u in t]))

    def labelling_on_basis(self,t):
        """
        Put canonical labels on a word in the Commutative operad.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.labelling_on_basis(m)
            B[word: 1234]
        """
        return self.basis()[self.basis().keys()(sorted([1+i for i in range(t.length())]))]

    def unlabelling_on_basis(self,t):
        """
        Removes the labels of a tree in the Commutative operad.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.unlabelling_on_basis(m)
            B[word: 1111]
        """
        return self.basis()[self.basis().keys()([1 for i in range(t.length())])]

    def grafts(self,x,y,i):
        """
        Auxiliary procedure: inserts a word y at position i in a word x
        and returns a word

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.grafts(Words("acb"), Words("de"),"c")
            word: abde

        """
        if x[0]==i:
            return self.basis().keys()(sorted(y+x[1:]))
        else:
            return self.basis().keys()(sorted(x[:1]+self.grafts(x[1:],y,i)))

    def composition_on_basis(self,x,y,i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis(Words("acb"), Words("de"),"c")
            B[word: abde]
        """
        if not(i in x):
            return "The composition index is not present."
        else:
            return self.basis()[self.grafts(x,y,i)]

    def commutative_product(self, x, y):
        """
        This computes the commutative product.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: W = A.basis().keys()
            sage: x = A(W('ab'))
            sage: y = A(W('dc'))
            sage: A.commutative_product(x, y)
            B[word: abcd]
        """
        gen=self.basis()[self.basis().keys()([0,1])]
        return gen.compose(x,0).compose(y,1)

    def operad_generators(self):
        """
        EXAMPLES::

        sage: CommutativeOperad(QQ).operad_generators()
        Finite family {'commutative_product': B[word: 12]}
        """
        from sage.sets.family import Family
        return Family(dict([("commutative_product",
                             self.basis()[self.basis().keys()([1,2])])]))

    def operad_morphism_on_basis(self,t,codomain):
        """
        defines a morphism from the Commutative operad to the target operad

        the target operad has to possess a method called commutative_product

        the argument should not have repeated labels

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
        """
        targetProduct=codomain.commutative_product
        n = len(t)
        if n==1 :
            return codomain.one(t[0])
        else:
            return targetProduct(self.operad_morphism_on_basis(t[0],codomain),self.operad_morphism_on_basis(t[1:],codomain))
