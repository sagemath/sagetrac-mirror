from math import*
from cartan_type import CartanType
class SuperBracket(x,y)

        """
        Returns the Lie super bracket [x, y] = x y -(-1)^{|x||y|} y x of x and y.

        INPUT:

    x, y -- elements of self 

EXAMPLES

                sage: F = SuperAlgebrasWithBasis?(QQ).example()
                sage: F
                An example of an superalgebra with basis: the free superalgebra on the generators ('a', 'b', 'c') over Rational Field
                sage: a,b,c = F.super_algebra_generators()
                sage: F.super_bracket(a,b)
                B[word: ab] - (-1)^{degree(a)*degree(b)}B[word: ba]

        This measures the default of commutation between x and y.
        F endowed with the super_bracket operation is a Lie superalgebra;

in particular, it satisfies super Jacobi's identity

            sage: F.super_bracket( F.super_bracket(a,b), c) + F.super_bracket(F.super_bracket(b,c),a) + F.super_bracket(F.super_bracket(c,a),b)
            0

        """
    def even_odd(self)
	if self.x != cartan_type:
           return even

	   else:
		  odd

    def degree(self,x):

         if self.x is even:

               0
         else:

               1

return x*y -(-1)^{degree(x)*degree(y)}* y*x
