r"""
Convergent bounded series on ultrametric disks


Throughout this tutorial, we work over :math:`\mathbb Q_5` but 
everything presented below extends to any field equipped with a discrete 
valuation for which it is complete.

    sage: K = Qp(5); K
    5-adic Field with capped relative precision 20

We denote by :math:`v_K` the valuation on K and we assume that it takes 
integral values.


Mathematical background
-----------------------

Let nu be a real number. The disk of log radius nu, denoted in this 
worksheet by :math:`D_\nu`, is defined as the set of all points in 
an algebraic closure of K whose valuation is bigger than \nu. When
:math:`\nu = 0`, the disk :math:`D_0` is the open unit disk.

A series  

.. MATH::

    f(x) = \sum_{i=0}^\infty a_i x^i

converges and is bounded on :math:`D_\nu` if and only if the quantity:

.. MATH::

    v_\nu(f) = \inf_i \,\, v(a_i) + \nu i

is finite. When it occurs, :math:`v_\nu(f)` is called the *Gauss 
valuation* of f. It satisfies the relations 

.. MATH::

    v_\nu(f+g) \geq \min (v_\nu(f), v_\nu(g))

    v_\nu(fg) = v_\nu(f) + v_\nu(g)

We shall denote by :math:`S_\nu` the ring of convergent bounded series 
on the disk :math:`D_\nu`.

From now on, we assume that :math:`\nu` is a rational number. With this 
extra assumption, it is easy to see that the infimum that defines the 
Gauss valuation is reached for some index i. The smallest such 
:math:`i` is called the Weierstrass degree of :math:`f` and will be 
denoted below by :math:`deg_\nu(f)`. By convention, we agree that the 
Weierstrass degree of 0 is :math:`-\infty`. We then have the usual 
relation:

.. MATH::

    deg_\nu(fg) =  deg_\nu(f) + deg_\nu(g).

We have a Weierstrass preparation theorem: any :math:`f \in D_\nu` can 
be written as a produit :math:`f = f_0 \times u` where :math:`f_0` is a 
*polynomial* of degree :math:`deg_\nu(f)` and :math:`u` is a *unit* is 
:math:`S_\nu`.

We also have an Euclidean division on :math:`S_\nu`: if :math:`a` and 
:math:`b` and two elements of :math:`D_\nu` with :math:`b \neq 0`, then 
there exist :math:`q` and :math:`r` such that :math:`a = bq + r` and 
:math:`r` is a *polynomial* of degree :math:`< \text{deg}_\nu(b)`. 
Furthermore the above properties determine uniquely :math:`q` and 
:math:`r`.

We emphasize that the Weierstrass decomposition and the Euclidean 
division both depend strongly on the log radius :math:`\nu`.


Construction of rings and elements
----------------------------------

We can construct the ring :math:`S_\nu` by using the function 
*BoundedSeriesRing* . The syntax is the following:

    sage: S1.<t> = BoundedSeriesRing(K,1); S1
    Bounded Convergent Series Ring in t on {val > 1} over 5-adic Field with capped relative precision 20

If the log radius is omitted, the computer uses the default value 0
(corresponding to :math:`D_0` which is the unit open disk).

    sage: S.<x> = BoundedSeriesRing(K); S
    Bounded Convergent Series Ring in x on {val > 0} over 5-adic Field with capped relative precision 20

    sage: S.log_radius()
    0

The rings we have just constructed are also equipped with a default 
x-adic precision: it's the precision at which series in this parent will 
be truncated when non exact operations (like inversion or euclidean 
divisions) are performed.

The default value is 20, but we can change it by setting the attribute 
*default_prec* :

    sage: S.default_prec()
    20

    sage: S.<x> = BoundedSeriesRing(K, default_prec=10)
    sage: S.default_prec()
    10

As usual the above syntax has defined the generators t and x 

    sage: x
    (1 + O(5^20))*x
    sage: x.parent()
    Bounded Convergent Series Ring in x on {val > 0} over 5-adic Field with capped relative precision 20

    sage: t.parent()
    Bounded Convergent Series Ring in t on {val > 1} over 5-adic Field with capped relative precision 20

The simplest way to define an element is just to write it:

    sage: f = 25 + 5*x + x^3 + 10*x^8; f
    5^2 + O(5^22) + (5 + O(5^21))*x + (1 + O(5^20))*x^3 + (2*5 + O(5^21))*x^8

The above f is in fact an actual polynomial without any O(x^n). We can 
check it using the method *prec* , which returns the x-adic precision.

    sage: f.prec()
    +Infinity

Note that, even if f belongs to the ring :math:`S_0`, its log radius of 
convergence is infinite since it is a polynomial. Hence, you have to be 
very careful with the behaviour of the different methods *log_radius* .

    sage: f.log_radius()
    +Infinity

    sage: f.parent().log_radius()
    0

We now can add a O(x^n) to f using the method *change_prec* :

    sage: g = f.change_prec(15); g
    5^2 + O(5^22) + (5 + O(5^21))*x + (1 + O(5^20))*x^3 + (2*5 + O(5^21))*x^8 + O(x^15)

The meaning of this new term O(x^15) needs extra explanations: it 
replaces an unknown sum of the form

.. MATH:

    \sum_{i=15}^\infty a_i x^i

where all :math:`a_i` are *integral* (*i.e.* have a nonnegative 
valuation).

Using the attribute *valuation_final_terms* we can require that they all 
have a valuation not less than some given value, which may be different 
from 0:

    sage: h = f.change_prec(15, valuation_final_terms=-2); h
    5^2 + O(5^22) + (5 + O(5^21))*x + (1 + O(5^20))*x^3 + (2*5 + O(5^21))*x^8 + O(5^(-2)*x^15)

Warning: the value of valuation_final_terms corresponds to the valuation 
of the final terms in the ring :math:`S_\nu` and not to the valuation of 
individual coefficients. Both notions coincide when the log radius 
vanishes, but they do not in general. Here is an example.

    sage: t.change_prec(15, valuation_final_terms=-2)
    (1 + O(5^20))*t + O(5^(-17)*t^15)

In the example above, the term O(5^(-17)*t^15) hides a sum 

.. MATH:

    \sum_{i=15}^\infty a_i x^i

with :math:`v_K(a_i) + \nu i \geq -2` that is :math:`v_K(a_i) \geq 
-2-i`. In particulier the valuation of :math:`a_{15}` has to be at least 
-17. This explains the magic exponent -17 that has appeared above.


Operations on elements
----------------------

Sums, oppsites, products are implemented. Here is a small demo 
(look also at precisions):

    sage: g + h
    2*5^2 + O(5^22) + (2*5 + O(5^21))*x + (2 + O(5^20))*x^3 + (4*5 + O(5^21))*x^8 + O(5^(-2)*x^15)
    sage: g * h
    5^4 + O(5^24) + (2*5^3 + O(5^23))*x + (5^2 + O(5^22))*x^2 + (2*5^2 + O(5^22))*x^3 + (2*5 + O(5^21))*x^4 + (1 + O(5^20))*x^6 + (4*5^3 + O(5^23))*x^8 + (4*5^2 + O(5^22))*x^9 + (4*5 + O(5^21))*x^11 + O(5^(-2)*x^15)
    sage: h^3
    5^6 + O(5^26) + (3*5^5 + O(5^25))*x + (3*5^4 + O(5^24))*x^2 + (5^3 + 3*5^4 + O(5^23))*x^3 + (5^3 + 5^4 + O(5^23))*x^4 + (3*5^2 + O(5^22))*x^5 + (3*5^2 + O(5^22))*x^6 + (3*5 + O(5^21))*x^7 + (5^5 + 5^6 + O(5^25))*x^8 + (1 + 2*5^4 + 2*5^5 + O(5^20))*x^9 + (5^3 + 5^4 + O(5^23))*x^10 + (2*5^3 + 2*5^4 + O(5^23))*x^11 + (2*5^2 + 2*5^3 + O(5^22))*x^12 + (5 + 5^2 + O(5^21))*x^14 + O(5^(-6)*x^15)

If the two summands (or the two factors) do not converge on the same 
disk, coercions (*i.e.* restrictions in that case) are made 
automatically:

    sage: g1 = g.change_log_radius(1); g1
    5^2 + O(5^22) + (5 + O(5^21))*x + (1 + O(5^20))*x^3 + (2*5 + O(5^21))*x^8 + O(x^15)
    sage: g.parent()
    Bounded Convergent Series Ring in x on {val > 0} over 5-adic Field with capped relative precision 20
    sage: g1.parent()
    Bounded Convergent Series Ring in x on {val > 1} over 5-adic Field with capped relative precision 20
    sage: ans = g1 * g; ans
    5^4 + O(5^24) + (2*5^3 + O(5^23))*x + (5^2 + O(5^22))*x^2 + (2*5^2 + O(5^22))*x^3 + (2*5 + O(5^21))*x^4 + (1 + O(5^20))*x^6 + (4*5^3 + O(5^23))*x^8 + (4*5^2 + O(5^22))*x^9 + (4*5 + O(5^21))*x^11 + O(5^2*x^15)
    sage: ans.parent()
    Bounded Convergent Series Ring in x on {val > 1} over 5-adic Field with capped relative precision 20

The units in :math:`S_\nu` are exactly those series whose Weierstrass 
degree vanishes. For such elements, one can compute their inverse as 
follows:

    sage: unit = 1 + x + 5*x^2; unit
    1 + O(5^20) + (1 + O(5^20))*x + (5 + O(5^21))*x^2
    sage: inv = unit.inverse(); inv
    1 + O(5^20) + (4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x + (1 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^2 + (4 + 5 + O(5^20))*x^3 + (1 + 2*5 + O(5^20))*x^4 + (4 + 3*5 + 2*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^5 + (1 + O(5^20))*x^6 + (4 + 5^2 + 2*5^3 + O(5^20))*x^7 + (1 + 3*5 + 3*5^2 + 2*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^8 + (4 + 2*5 + 5^3 + 3*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^9 + O(x^10)
    sage: inv * unit
    1 + O(5^20) + O(x^10)

We remark that the answer is truncated to x-adic precision given by the 
default precision of the parent. Nevertheless, if the series carries 
itself it own precision, the answer will be truncated to the latter one 
(even if it is larger than the default precision of the parent):

    sage: unit.change_prec(15).inverse()
    1 + O(5^20) + (4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x + (1 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^2 + (4 + 5 + O(5^20))*x^3 + (1 + 2*5 + O(5^20))*x^4 + (4 + 3*5 + 2*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^5 + (1 + O(5^20))*x^6 + (4 + 5^2 + 2*5^3 + O(5^20))*x^7 + (1 + 3*5 + 3*5^2 + 2*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^8 + (4 + 2*5 + 5^3 + 3*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^9 + (1 + 5 + 5^2 + 4*5^4 + O(5^20))*x^10 + (4 + 4*5 + 4*5^3 + 4*5^4 + O(5^20))*x^11 + (1 + 4*5 + 2*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 3*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + O(5^20))*x^12 + (4 + 5 + 2*5^2 + 4*5^3 + O(5^20))*x^13 + (1 + 2*5 + 3*5^2 + 2*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + O(5^20))*x^14 + O(x^15)

When the series in not invertible, an error is raised

    sage: f.inverse()
    Traceback (most recent call last):
    ...
    ZeroDivisionError: This series in not invertible in the ring of bounded convergent series on {val >= 0}

The methods *gauss_valuation* (or simply *valuation*) and 
*weierstrass_degree* (or simply *degree*) compute respectively the Gauss 
valuation and the Weierstrass degree of a bounded series.

    sage: g.gauss_valuation()
    0
    sage: g.weierstrass_degree()
    3

It may append in some cases that the known informations are not enough 
to determine the Gauss valuation or the Weierstrass degree: for 
instance, the element h defined above may have Gauss valuation 0, -1
or -2.

By default, in that case, the highest possibility is returned (i.e. one 
does as if there was no big-oh).

    sage: h.gauss_valuation()
    0

Nevertheless if we set the attribute *secure=True* , an error is raised 
if the Gauss valuation can't be determined.

    sage: g.gauss_valuation(secure=True)
    0
    sage: h.gauss_valuation(secure=True)
    Traceback (most recent call last):
    ...
    PrecisionError: Unable to determine for sure the Gauss valuation

The method *is_secure* tells whether a given series has well determined 
Gauss valuation and Weierstrass degree or not.

    sage: g.is_secure()
    True
    sage: h.is_secure()
    False

The newton polygon of an element 

.. MATH:

    \sum_{i=0}^\infty a_i x^i

is defined as the convex hull in the plane of points of coordinates
:math:`(i, v_K(a_i))` together with the two following extra points at 
infinity: :math:`(0, +\infty)` and :math:`(\infty, -\nu\:\infty)`.

We can construct these newton polygons as follows:

    sage: f.newton_polygon()
    Infinite Newton polygon with 3 vertices: (0, 2), (1, 1), (3, 0) ending by an infinite line of slope 0

When the Newton polygon is not determined due to a lack of precision 
(*e.g.* is the series is not secure), an error is raised.

    sage: h.newton_polygon()
    Traceback (most recent call last):
    ...
    PrecisionError: The Newton polygon is not determined

The method *weierstrass_preparation* computes the decomposition of a 
series as a product of a polynomial and a unit.

    sage: g0, unit = g.weierstrass_preparation()
    sage: g0
    5^2 + 5^6 + O(5^7) + (5 + 4*5^5 + 3*5^6 + O(5^7))*x + (3*5^4 + 5^5 + O(5^6))*x^2 + (1 + O(5^20))*x^3
    sage: unit
    1 + 4*5^4 + O(5^6) + (2*5^3 + O(5^5))*x + (3*5^3 + 4*5^4 + O(5^5))*x^2 + (3*5^2 + 4*5^3 + O(5^4))*x^3 + (2*5 + O(5^3))*x^5 + O(x^12)
    sage: g0 * unit
    5^2 + O(5^7) + (5 + O(5^7))*x + (1 + O(5^6))*x^3 + (2*5 + O(5^3))*x^8 + O(x^12)
    sage: g0 * unit == g
    True

Remark that the product g0*unit has much less precision than g itself. 
In other terms, computing the decomposition of the Weierstrass 
preparation theorem may generate important losses of precision... 
something that should be kept in mind.

Let's remark furthermore that the O(x^n) may decrease the p-adic 
precision of the coefficients of the polynomial part (and as well the 
first coefficients of the unit part). For instance, compare the two 
following results:

    sage: g0, _ = g.weierstrass_preparation(); g0
    5^2 + 5^6 + O(5^7) + (5 + 4*5^5 + 3*5^6 + O(5^7))*x + (3*5^4 + 5^5 + O(5^6))*x^2 + (1 + O(5^20))*x^3
    sage: f0, _ = f.weierstrass_preparation(); f0
    5^2 + 5^6 + 4*5^7 + 4*5^8 + 2*5^9 + 3*5^11 + 2*5^12 + 5^15 + 3*5^16 + 3*5^17 + 3*5^18 + 3*5^20 + 4*5^21 + O(5^22) + (5 + 4*5^5 + 3*5^6 + 4*5^7 + 2*5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 3*5^13 + 4*5^14 + 4*5^15 + 3*5^17 + 4*5^18 + 2*5^19 + O(5^21))*x + (3*5^4 + 5^5 + 2*5^9 + 4*5^10 + 3*5^11 + 2*5^12 + 5^13 + 4*5^14 + 3*5^15 + 2*5^16 + 2*5^17 + 5^18 + 2*5^19 + 2*5^20 + 4*5^21 + 3*5^22 + 2*5^23 + O(5^24))*x^2 + (1 + O(5^20))*x^3

If the series we started with is not secure, the computes returns an 
answer like this:

    sage: h0, unit = h.weierstrass_preparation()
    sage: h0
    5^2 + O(5^5) + (5 + O(5^5))*x + (1 + O(5^20))*x^3 + O(unknown)
    sage: unit
    1 + O(5^4) + O(unknown)

The O(unknown) means that there is no guarantee on any coefficient... 
but nonetheless, the returned answer is somehow possible.

If we set the attribute *secure=True* , the computer raises an error 
instead of trying to guess what is the best answer.

    sage: h.weierstrass_preparation(secure=True)
    Traceback (most recent call last):
    ...
    PrecisionError: This series is not secure


Euclidean division is available *via* the usual method *quo_rem* :

    sage: a = S.random_element(prec=20)
    sage: q, r = a.quo_rem(f)

We check the required properties:

    sage: r.prec() is Infinity
    True
    sage: r.polynomial().degree() < f.degree()
    True
    sage: a == q*f + r
    True

If the divisor is not secure, there is no guarantee on the result. As 
above, it appears on the writing *via* a term O(unknown).

    sage: q, r = (h+1).quo_rem(h)
    sage: q
    1 + O(5^4) + O(unknown)
    sage: r
    1 + O(5^5) + O(unknown)

As above, if the attribute *secure* is set to True, an error is raised 
instead of this:

    sage: (h+1).quo_rem(h, secure=True)
    Traceback (most recent call last):
    ...
    PrecisionError: Unable to determine for sure the Gauss valuation

If a is an element of the disk :math:`D_\nu`, we can evaluate a series 
belonging to :math:`S_\nu` a as follows:

    sage: PolRing.<var> = PolynomialRing(K)
    sage: L.<pi> = K.extension(var^2 + 5)
    sage: gpi = g(pi); gpi
    pi^4 + 3*pi^10 + pi^12 + O(pi^15)

If a does not belong to the domain of convergence, an error is raised:

    sage: g(1+pi)
    Traceback (most recent call last):
    ...
    ValueError: Image does not define a valid morphism

Composition of series are also implemented. The syntax is analoguous:

    sage: g1 = g(x^2+5); g1
    2*5^2 + 5^3 + 2*5^9 + O(5^15) + (5 + 3*5^2 + 5^8 + 3*5^9 + O(5^14))*x^2 + (3*5 + 5^7 + 5^8 + 2*5^9 + O(5^13))*x^4 + (1 + 2*5^6 + 2*5^7 + 4*5^8 + O(5^12))*x^6 + (3*5^6 + 5^8 + O(5^11))*x^8 + (2*5^4 + 2*5^5 + 4*5^6 + O(5^10))*x^10 + (5^3 + 5^4 + 2*5^5 + O(5^9))*x^12 + (5^2 + 3*5^3 + O(5^8))*x^14 + (2*5 + O(5^7))*x^16 + O(x^30)

The (log) radius of convergence of the composite is automatically computed:

    sage: g2 = g(x/25); g2
    5^2 + O(5^22) + (5^-1 + O(5^19))*x + (5^-6 + O(5^14))*x^3 + (2*5^-15 + O(5^5))*x^8 + O(5^(-30)*x^15)
    sage: g2.log_radius()
    2
    sage: g2.parent()
    Bounded Convergent Series Ring in x on {val > 2} over 5-adic Field with capped relative precision 20

If the composite does not converge anywhere, an error is raised:

    sage: g(x+1)   # 1 does not belong to the radius of convergence
    ...            # thus, we can't translate the Taylor expansion in 1
    Traceback (most recent call last):
    ...
    ValueError: The composite does not converge


Morphisms from Bounded Series Ring
----------------------------------

As usual in Sage, we can use the method *hom* (of the parent) to 
construct a morphism from a bounded series ring to either a finite 
extension of the base ring (case of evaluation) or another bounded 
series ring (case of composition) by specifying the image of the 
variable. Here is a simple example:

    sage: phi = S.hom(pi); phi
    Ring morphism:
      From: Bounded Convergent Series Ring in x on {val > 0} over 5-adic Field with capped relative precision 20
      To:   Eisenstein Extension of 5-adic Field with capped relative precision 20 in pi defined by (1 + O(5^20))*var^2 + (O(5^21))*var + (5 + O(5^21))
      Defn: x |--> pi + O(pi^41)
    sage: phi(x)
    pi + O(pi^41)
    sage: phi(g)
    pi^4 + 3*pi^10 + pi^12 + O(pi^15)
    sage: phi(g) == gpi
    True

    sage: psi = S.hom(x^2 + 5); psi
    Ring endomorphism of Bounded Convergent Series Ring in x on {val > 0} over 5-adic Field with capped relative precision 20
      Defn: x |--> 5 + O(5^21) + (1 + O(5^20))*x^2
    sage: psi(g)
    2*5^2 + 5^3 + 2*5^9 + O(5^15) + (5 + 3*5^2 + 5^8 + 3*5^9 + O(5^14))*x^2 + (3*5 + 5^7 + 5^8 + 2*5^9 + O(5^13))*x^4 + (1 + 2*5^6 + 2*5^7 + 4*5^8 + O(5^12))*x^6 + (3*5^6 + 5^8 + O(5^11))*x^8 + (2*5^4 + 2*5^5 + 4*5^6 + O(5^10))*x^10 + (5^3 + 5^4 + 2*5^5 + O(5^9))*x^12 + (5^2 + 3*5^3 + O(5^8))*x^14 + (2*5 + O(5^7))*x^16 + O(x^30)
    sage: psi(g) == g1
    True

By default, morphisms constructed as before act trivially on the base 
ring. We can nevertheless specify another action using the attribute 
*morphism_on_coefficients* . As an example, we can construct the 
"Frobenius" on :math:`\mathbb Z_{5^3}[x]` as follows:

    sage: F.<a> = Qq(5^3)
    sage: sigma = F.frobenius_endomorphism()
    sage: S.<x> = BoundedSeriesRing(F)
    sage: Frob = S.hom(x^5, morphism_on_coefficients=sigma)
    sage: Frob
    Ring endomorphism of Bounded Convergent Series Ring in x on {val > 0} over Unramified Extension of 5-adic Field with capped relative precision 20 in a defined by (1 + O(5^20))*x^3 + (O(5^20))*x^2 + (3 + O(5^20))*x + (3 + O(5^20))
      Defn: x |--> (1 + O(5^20))*x^5
            Action on coefficients: Frob

    sage: Frob(a*x)
    ((2*a^2 + 4*a + 4) + (2*a^2 + 3*a + 4)*5 + a*5^2 + 2*a*5^3 + (a^2 + 3*a + 2)*5^4 + (a^2 + 2)*5^5 + 2*a*5^6 + 2*a*5^7 + (2*a^2 + 4*a + 4)*5^8 + (3*a^2 + 4*a + 1)*5^9 + 2*a^2*5^10 + (2*a^2 + a)*5^11 + (2*a^2 + a)*5^12 + (4*a^2 + 3*a + 4)*5^13 + (3*a^2 + a + 2)*5^14 + (3*a^2 + 4*a + 2)*5^15 + (3*a^2 + a + 2)*5^16 + (a^2 + 2*a + 3)*5^17 + (4*a^2 + 3*a + 3)*5^18 + (4*a^2 + 2*a + 4)*5^19 + O(5^20))*x^5
    sage: Frob(a*x) == sigma(a)*x^5
    True

When the image of the variable is a power of itself, instead of giving 
this image to the method hom, we can just give the exponent.

    sage: Frob2 = S.hom(5, morphism_on_coefficients=sigma); Frob2
    Ring endomorphism of Bounded Convergent Series Ring in x on {val > 0} over Unramified Extension of 5-adic Field with capped relative precision 20 in a defined by (1 + O(5^20))*x^3 + (O(5^20))*x^2 + (3 + O(5^20))*x + (3 + O(5^20))
      Defn: x |--> (1 + O(5^20))*x^5
            Action on coefficients: Frob
    sage: Frob2 == Frob
    True

There is actually a direct (and simpler) way to construct the Frobenius 
endomorphism on a ring of bounded series: it is given by the method 
*frobenius_endomorphism* .

    sage: S.frobenius_endomorphism()
    Ring endomorphism of Bounded Convergent Series Ring in x on {val > 0} over Unramified Extension of 5-adic Field with capped relative precision 20 in a defined by (1 + O(5^20))*x^3 + (O(5^20))*x^2 + (3 + O(5^20))*x + (3 + O(5^20))
      Defn: x |--> (1 + O(5^20))*x^5
            Action on coefficients: Frob
    sage: S.frobenius_endomorphism() == Frob
    True


AUTHOR::

- Xavier Caruso
"""

#############################################################################
#    Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation 
from sage.structure.category_object import normalize_names

import sage.categories.basic as categories

from integer import Integer
from rational_field import QQ 
from sage.structure.element import Element 
from ring import Field, CommutativeAlgebra
from infinity import Infinity 
import sage.misc.latex as latex

from polynomial.polynomial_ring_constructor import PolynomialRing 
from polynomial.polynomial_element import Polynomial 
from sage.structure.parent_gens import ParentWithGens

from bounded_series_ring_element import BoundedSeries 
from bounded_series_ring_morphism import BoundedSeriesBaseringInjection
from bounded_series_ring_morphism import BoundedSeriesRestriction


def BoundedSeriesRing(base_ring, log_radius=0, name=None, names=None,
                      sparse=False, default_prec=None):
    if default_prec is None:
        default_prec = 20

    if name is None:
        name = names
    if name is None:
        raise TypeError, "You must specify the name of the indeterminate of the Bounded Power series ring."
    try:
        name = normalize_names(1, name)[0]
    except IndexError:
        raise NotImplementedError("Multivariate bounded power series rings are not implemented yet.")
    except TypeError:
        raise TypeError, "illegal variable name"
    
    if not (isinstance(base_ring, Field)): # and base_ring.is_cdvf()):  <-- TODO: fix this
        raise TypeError("base_ring must be a complete discrete valuation field")
    return BoundedSeriesRing_generic(base_ring, log_radius, name, default_prec, sparse=sparse)


def is_BoundedSeriesRing(R):
    return isinstance(R, BoundedSeriesRing_generic)


class BoundedSeriesRing_generic(CommutativeAlgebra, UniqueRepresentation):
    @staticmethod
    def __classcall__(cls, base_ring, log_radius=0, name=None, default_prec=20, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse bounded series are not implemented")
            else:
                from bounded_series_poly import BoundedSeries_poly
                element_class = BoundedSeries_poly
        return super(BoundedSeriesRing_generic,cls).__classcall__(cls, base_ring, log_radius, name, default_prec, sparse, element_class)
        
    def __init__(self, base_ring, log_radius=0, name=None, default_prec=20, sparse=False, element_class=None, category=None):
        self.__is_sparse = sparse
        self._series_class = element_class
        self._by_one = False
        self._log_radius = log_radius
        self._default_prec = default_prec
        self._polynomial_ring = PolynomialRing(base_ring, name=name)

        # Algebra.__init__ also calls __init_extra__ of Algebras(...).parent_class, which tries to provide a conversion from the base ring, if it does not exist. This is for algebras that only do the 
        # generic stuff in their initialisation. But here, we want to use PolynomialBaseringInjection. Hence, we need to wipe the memory and construct the conversion from scratch.
        CommutativeAlgebra.__init__(self, base_ring, names=name, normalize=True, category=category)
        self.__generator = self._series_class(self, self._polynomial_ring.gen(), valuation_final_terms=Infinity, prec=Infinity, is_gen=True)
        self._base_inject = BoundedSeriesBaseringInjection(base_ring,self)
        self._coercions_log_radius = { }

        if log_radius in QQ:
            self._refine_category_(categories.EuclideanDomains())

    #def __reduce__(self):

    def _coerce_map_from_(self, S):
        base = self.base_ring()
        if base is S:
            return self._base_inject
        elif base.has_coerce_map_from(S):
            return self._base_inject * base.coerce_map_from(S) 
        elif isinstance(S, BoundedSeriesRing_generic):
            log_radius = S.log_radius()
            if self.base_ring() is S.base_ring() and self.variable_name() is S.variable_name() and self._log_radius >= log_radius:
                key = (self.default_prec(), log_radius)
                try:
                    return self._coercions_log_radius[key]
                except KeyError:
                    map = BoundedSeriesRestriction(S, self, check=False)
                    self._coercions_log_radius[key] = map
                    return map

    def _element_constructor_(self, x=None, valuation_final_terms=Infinity, prec=Infinity, check=True, is_gen=False):
        if is_gen:
            x = [ 0,1 ]
            prec = Infinity
            valuation_final_terms = Infinity
        elif isinstance(x, Element):
            if isinstance(x, BoundedSeries):
                if x.log_radius() > self.log_radius():
                    raise TypeError("Impossible to decrease the logarithmic radius of convergence")
                prec = x.prec()
                valuation_final_terms = x.valuation_final_terms()
                if valuation_final_terms is not None:
                    valuation_final_terms += prec * (self.log_radius() - x.log_radius())
                x = x.polynomial()
            #elif isinstance(x, PowerSeries):
            if isinstance(x, Polynomial):
                pass
            else:
                x = [ x ]
        elif isinstance(x, int):
            x = [ x ]
        return self._series_class(self, x, valuation_final_terms, prec)

    def _repr_(self):
        s = "Bounded Convergent Series Ring in %s on {val > %s} over %s"%(self.variable_name(), self.log_radius(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def is_sparse(self):
        return self.__is_sparse

    def is_dense(self):
        return not self.__is_sparse

    def log_radius(self):
        return self._log_radius

    def default_prec(self):
        return self._default_prec

    def _latex_(self):
        return "%s\left<%s\right>"%(latex.latex(self.base_ring()), self.latex_variable_names()[0])

    #def _is_valid_homomorphism_(self, codomain, im_gens):

    def base_extend(self, R):
        raise NotImplementedError
    
    def change_ring(self, R):
        raise NotImplementedError

    def change_log_radius(self, log_radius):
        return self.__class__(self.base_ring(), log_radius, self.variable_name(), sparse=self.is_sparse(), default_prec=self.default_prec())

    def change_default_prec(self, default_prec):
        return self.__class__(self.base_ring(), self.log_radius(), self.variable_name(), sparse=self.is_sparse(), default_prec=default_prec)

    def is_exact(self):
        return False

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "generator n>0 not defined"
        return self.__generator

    def ngens(self):
        return 1

    def random_element(self, prec=None, degree=None, valuation=0, *args, **kwds):
        from sage.functions.other import ceil
        if prec is None:
            prec = self.default_prec()
        if degree is None or degree >= prec:
            degree = prec - 1
        if degree is Infinity:
            degree = self.default_prec() - 1
        val = valuation
        log_radius = self._log_radius
        coeffs = [ ]
        base = self.base_ring()
        integer_base = base.integer_ring()
        for i in range(degree+1):
            coeffs.append(base(integer_base.random_element()) << ceil(val))
            val -= log_radius
        if prec is Infinity:
            valuation_final_terms = Infinity
        else:
            valuation_final_terms = valuation
        return self._series_class(self, coeffs, valuation_final_terms=valuation_final_terms, prec=prec)

    def __cmp__(self, other):
        return self is other

    def uniformizer(self):
        if self._log_radius not in QQ:
            raise ValueError("This ring does not have a uniformizer")
        base = self.base_ring()
        unif = base.uniformizer()
        image_val = unif.valuation()
        fraction = self._log_radius / image_val
        num = fraction.numerator()
        denom = fraction.denominator()
        _, n, v = num.xgcd(denom)
        if n < 0:
            n += denom
            v -= num
        coeffs = n * [ base(0) ] + [ unif ** v ]
        return self._series_class(self, coeffs)

    def is_atomic_repr(self):
        return False

    def is_commutative(self):
        return True

    def is_field(self, proof = True):
        return False

    def is_finite(self):
        return False

    def characteristic(self):
        return self.base_ring().characteristic()

    #def integers(self):
    #    return BoundedByOneSeriesRing_generic(self.base_ring(), self._log_radius(), self.variable_name(), self._default_prec, self.__is_sparse, self._series_class)

    def hom(self, im_gen, morphism_on_coefficients=None):
        from sage.categories.homset import Hom
        from bounded_series_ring_morphism import BoundedSeriesHomomorphism_im_gens
        if isinstance(im_gen, (int, Integer)):
            parent = self
        elif isinstance(im_gen, list):
            parent = im_gen[0].parent()
        else:
            parent = im_gen.parent()
        homset = Hom(self, parent)
        return BoundedSeriesHomomorphism_im_gens(homset, im_gen, morphism_on_coefficients)

    def frobenius_endomorphism(self, n=1, b=None):
        from sage.categories.homset import Hom
        from bounded_series_ring_morphism import BoundedSeriesHomomorphism_im_gens
        base = self.base_ring()
        if b is None:
            b = base.characteristic()
            if b == 0: b = base.prime()
        morphism_on_coefficients = base.frobenius_endomorphism(n=n)
        homset = Hom(self, self)
        return BoundedSeriesHomomorphism_im_gens(homset, b**n, morphism_on_coefficients)
