r"""
Certain functions used in the sieve step in solving an S-unit equation.

AUTHORS:

- Alejandra Alvarado, Angelos Koutsianas, Beth Malmskog, Christopher Rasmussen, Christelle Vincent, Mckenzie West (2017-01-10): original version

EXAMPLES:

"""

#*************************************************************************************************
#       Copyright (C) 2017 Alvarada, Koutsianas, Malmskog, Rasmussen, Vincent, West
#
# This program is free software: you can redistricute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*************************************************************************************************

# We introduce some terminology.
#
# If tau is an S-unit, and rho_0, ..., rho_t is a (normalized) set of generators for the S-units,
# then tau can almost be uniquely described by its **exponent vector**. This is a (t+1)-tuple of integers
# (a_0, a_1, ..., a_t), such that
#
#         t
#       -----
#        | |         a_j
# tau =  | |   rho_j
#        | |
#       j = 0
#
# If K has exactly w roots of unity, and we insist that 0 <= a_0 < w, then the exponent vector is unique.
#
# Let m be an integer. An exponent vector modulo m is a tuple (a_0, a_1, ..., a_t) where each entry, except possibly the first,
# is given by a residue class in ZZ / m ZZ rather than an integer. Sadly, we may refer to 'exponent vectors' when we really mean
# 'exponent vectors modulo m'
#
# More terminology. A **residue field vector,** or rfv for short.
#
# Suppose q is a prime, not appearing in S, which splits completely, say into OK-ideals Q[0], ..., Q[d-1].
#
# There is a map Phi_q from the S-units of K to the d-fold product of copies of F_q, namely:
#
# tau |-----> ( tau mod Q[0], tau mod Q[1], ..., tau mod Q[d-1] )
#
# The image of this map is the **residue field vector** for tau.
#
# If we fix a residue field vector in characteristic q, we can try to determine which exponent vectors modulo (q-1)
# are associated to this residue field vector.
#
# Since the exponent vector for tau is unique to tau, we will also talk about the residue field vector associated
# to an exponent vector (when we really mean 'associated to tau', the S-unit determined by the exponent vector.)
#
# NOTE: If tau0 + tau1 = 1, then this equation also holds modulo Q[j] for any j.
# Thus, the rfv for tau0 and tau1 must sum to (1,1,...,1).
#
# Thus, given a residue field vector v, we define the **complement vector** to v as the vector w in F_q^d
# such that v + w = (1, 1, ..., 1).
#
# In the upcoming functions, we often refer to a dictionary which uses the exponent vectors as keys
# and the residue field vectors as values. In our original code, we often referred to this as the
# 'rfv dictionary.'
#
# The main purpose of such dictionaries is towards finding all possible solutions to the S-unit equation.
# Thus, when we conclude a dictionary entry cannot possibly be used to create a solution, it will be dropped.
#
# More terminology: Let a be a exponent vector modulo (q-1), and let tau be the associated S-unit.
#
# We say an exponent vector b is a **complementary exponent vector modulo (q-1)** to a
# if the residue field vectors associated to a and b (modulo q) are complementary.

# We want to consider exponent vectors with respect to two different moduli,
# and determine whether they are compatible, in the sense that they could
# both descend from the same characteristic zero exponent vector.
#
# Suppose a0 is an exponent vector modulo m0 (likely m0 = q0-1 for some prime q0)
# and a1 is an exponent vector modulo m1. Let g = gcd(m0, m1).
#
# We say a0 and a1 are **compatible** if:
#    - a0[0] == a1[0] (this is the torsion exponent, so it has to match exactly)
#    - for each i > 0, a0[i] == a1[i] modulo g.
#
# Let S0 = [q0, q1, ..., qk] be a set of primes, and suppose that for each qj in S0,
# we have the dictionary comp_exp_vec[qj], which has as its keys possible exponent vectors modulo qj-1
# and which has as its values a list of complementary exponent vectors modulo qj-1.
#
# A **compatible system** associated to S0 is a pair of lists:
# [v0, .., vk]
# [w0, .., wk]
# satisfying all the following conditions:
#       * vj is an exponent vector modulo qj-1 for all j
#       * wj is an exponent vector modulo qj-1 for all j
#       * vi and vj are compatible for all i, j
#       * wi and wj are compatible for all i, j
#       * wj lies in comp_exp_vec[ vj ] for all j

from sage.rings.fast_arith import prime_range
from sage.arith.all import gcd, factor, lcm, CRT
from sage.misc.misc_c import prod
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.number_field.number_field import NumberField_absolute
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ

def completely_split_primes(K, Bound = 200):
    r"""
    Returns a list of rational primes which split completely in the number field `K`.

    INPUT:

    - ``K`` -- a number field
    - Bound -- a positive integer bound (default 200)

    OUTPUT:

    A list of all primes `p < Bound` which split completely in `K`.

    NOTE:

    This function does not screen against primes in S, specified for determining S-Unit equation solutions.

    EXAMPLES:

        sage: K.<xi> = NumberField(x^3 - 3*x + 1)
        sage: completely_split_primes(K, 100)
        [17, 19, 37, 53, 71, 73, 89]

    """

    f = K.defining_polynomial()
    split_primes = []
    for p in prime_range(Bound):
        Fp = GF(p)
        FpT = PolynomialRing(Fp,'T')
        g = FpT(K.defining_polynomial())
        if len(factor(g)) == K.degree():
            split_primes.append(p)
    return split_primes

def split_primes_large_lcm(K, Bound):
    r"""
    Returns a list L of completely split primes in K such that lcm {q - 1 : q in L } exceeds 2*Bound + 1

    INPUT:

    - ``K`` -- a number field
    - ``Bound`` -- a positive integer

    OUTPUT:

    A list L of rational primes q, each splitting completely in K, and such that
    lcm {q - 1 : q in L } is greater than or equal to 2 * Bound + 1.

    NOTES:

    - A series of compatible exponent vectors for the primes in L will lift to **at most**
      one integer exponent vector whose entries a_i satisfy |a_i| <= Bound.

    - The ordering of this set is not very intelligent for the purposes of the later
      sieving processes.

    EXAMPLES::

        sage: K.<xi> = NumberField(x^3 - 3*x + 1)
        sage: split_primes_large_lcm(K, 200)
        [17, 19, 37, 53]

    With a tiny bound, SAGE may ask you to increase the bound.

    ::

        sage: K.<xi> = NumberField(x^2 + 163)
        sage: split_primes_large_lcm(K, 8)
        Traceback (most recent call last):
        ...
        ValueError: Not enough split primes found. Increase bound.

    """

    split_prime_list = completely_split_primes(K, 4*Bound + 4)
    lcm_list = []
    L = 1
    while L < 2*Bound + 1:
        if split_prime_list == []:
            # Need More Primes!
            raise ValueError('Not enough split primes found. Increase bound.')
        q = split_prime_list.pop(0)
        L = lcm(L, q-1)
        lcm_list.append(q)
    return lcm_list

def sieve_ordering(SUK, q):
    r"""
    Returns ordered data for running sieve on the primes in `SUK` over the rational prime `q`.

    INPUT:

    - ``SUK`` -- the S-unit group of a number field K
    - ``q``   -- a rational prime number which splits completely in `K`

    OUTPUT:

    A list of lists, [ideals_over_q, residue_fields, rho_images, product_rho_orders], where
        - ideals_over_q is a list of the `d = [K:\mathbb{Q}]` ideals in `K` over `q`
        - residue_fields[i] is the residue field of ideals_over_q[i]
        - rho_images[i] is a list of the generators in `rho` modulo ideals_over_q[i]
        - product_rho_orders[i] is the product of the multiplicative orders of the elements in rho_images[i]

    The list Q_ideal is sorted so that the product of orders is smallest for Q_ideal[0], as this will
    make the later sieving steps more efficient.

    NOTES:

    - The primes of S may not lie over over q.

    EXAMPLES:

        sage: K.<xi> = NumberField(x^3 - 3*x + 1)
        sage: SUK = K.S_unit_group(S=3)
        sage: sieve_data = sieve_ordering(SUK, 19)
        sage: sieve_data[0]
        [Fractional ideal (-2*xi^2 + 3),
        Fractional ideal (xi - 3),
        Fractional ideal (2*xi + 1)]

        sage: sieve_data[1]
        [Residue field of Fractional ideal (-2*xi^2 + 3), 
        Residue field of Fractional ideal (xi - 3), 
        Residue field of Fractional ideal (2*xi + 1)]

        sage: sieve_data[2]
        [[18, 17, 6, 8], [18, 11, 2, 4], [18, 13, 8, 10]]

        sage: sieve_data[3]
        [972, 972, 3888]

    """

    K = SUK.number_field()
    rho = SUK.gens_values()
    d = K.absolute_degree()
    primes_over_q = K.primes_above(q)
    # q must split completely.
    if len(primes_over_q) != d:
        raise ValueError('The prime q is not completely split.')

    for P in SUK.primes():
        if P in K.primes_above(q):
            raise ValueError('There is a prime in S over q.')

    # We initialize some unsorted lists:

    ideals_over_q_unsorted = []
    # will contain the d ideals over q. In the comments, we use Q[i] as shorthand for ideals_over_q_unsorted[i]

    residue_fields_unsorted = []
    # ith entry is the residue field of Q[i].
    # Note that these are all isomorphic to GF(q), but each has a different reduction map from OK.

    rho_images_unsorted = []
    # The ith entry is a list of the images of the generators in rho under the natural reduction map
    #     OK ---> residue_fields_unsorted[i] == OK / Q[i]

    product_rho_orders_unsorted = []
    # Viewed as elements in the multiplicative group GF(p)^\times, each generator has some multiplicative order.
    # product_rho_orders_Q_unsorted[i] records the product of these multiplicative orders.

    # We will sort the ideals so that the 0th ideal has the minimum value for this product.

    for Qi in primes_over_q:
        ideals_over_q_unsorted.append( Qi )
        residue_fields_unsorted.append( Qi.residue_field() )
        rho_mod_Qi = [ Qi.residue_field()( rho_j) for rho_j in rho ]
        rho_images_unsorted.append( rho_mod_Qi )
        product_of_orders_of_images = prod( rho_ij.multiplicative_order() for rho_ij in rho_mod_Qi )
        product_rho_orders_unsorted.append( product_of_orders_of_images )

    # zip() will turn the four separate lists into one list of 4-tuples.
    # the sort command will sort in the first entry -- i.e., the product of rho orders.

    q_data_zip = zip( product_rho_orders_unsorted, ideals_over_q_unsorted, residue_fields_unsorted, rho_images_unsorted )
    q_data_zip.sort()

    # The remaining code simply regurgitates the data into four separate lists.

    ideals_over_q = []
    residue_fields = []
    rho_images = []
    product_rho_orders = []

    for qdata_item in q_data_zip:
        product_rho_orders.append( qdata_item[0] )
        ideals_over_q.append( qdata_item[1] )
        residue_fields.append( qdata_item[2] )
        rho_images.append( qdata_item[3] )

    return ideals_over_q, residue_fields, rho_images, product_rho_orders

def bounded_integer_lifts(r, m, bound):
    r"""
    Returns all integers up to a given bound which are equivalent to r modulo m.

    INPUT:

    - ``r`` -- an integer, representing a residue class modulo `m`
    - ``m`` -- an integer; the modulus under consideration
    - ``bound`` -- a positive integer

    OUTPUT:

    A list of integers `x` satisfying both `x \equiv r \mod m` and `abs(x) \leq bound`.

    EXAMPLES:

        sage: bounded_integer_lifts(2, 7, 13)
        [2, 9, -5, -12]

        sage: bounded_integer_lifts(4, 12, 1)
        []

    """

    r0 = r % m
    lifts = []
    r_plus = r0
    # first find the lifts in the range [0, bound]
    while r_plus <= bound:
        lifts.append(r_plus)
        r_plus += m

    r_minus = r0 - m
    # now, find the lifts in the range [-bound, 0)
    while r_minus >= -bound:
        lifts.append(r_minus)
        r_minus -= m

    return lifts

def bounded_vector_lifts( exponent_vector, m, bound ):
    r"""
    Given an exponent vector modulo m, construct the possible lifts which agree modulo m and with no entry exceeding     in absolute value.

    INPUT:

    - exponent_vector -- an exponent vector (to be viewed as an exponent vector modulo `m`)
    - m               -- a positive integer > 1
    - bound           -- a positive integer, bounding the absolute value of entries in lifts

    OUTPUT:

    A list of all exponent vectors with integer entries which satisfy the following criteria:
        - the 0th entry matches the 0th entry of exponent_vector
        - for each $j > 0$, the $j$th entry is congruent to the $j$th entry of exponent_vector modulo $m$
        - all entries, except possibly the first, are bounded in absolute value by `bound`.

    EXAMPLES:

        sage: bounded_vector_lifts((2,7), 16, 44)
        [(2, 7), (2, 23), (2, 39), (2, -9), (2, -25), (2, -41)]

    """

    length = len(exponent_vector)
    if length == 0:
        return [ () ]
    elif length == 1:
        return [(exponent_vector[0],)]
    elif length > 1:
        # We work recursively. First, lift the final entry of the vector.
        final_entry_lifts = bounded_integer_lifts( exponent_vector[-1], m, bound )
        # Second, get the lifts of the shorter vector:
        start_of_exponent_vector = exponent_vector[:-1]
        lifted_vectors = []
        for start_of_lift in bounded_vector_lifts( start_of_exponent_vector, m, bound ):
            for final_entry in final_entry_lifts:
                new_vector = start_of_lift + (final_entry,)
                lifted_vectors.append( new_vector )
        return lifted_vectors

def clean_rfv_dict( rfv_dictionary ):
    r"""
    Given a residue field vector dictionary, removes some impossible keys and entries.

    INPUT:

    - rfv_dictionary -- a dictionary whose keys are exponent vectors and whose values are
                        residue field vectors

    OUTPUT:

    None. But it removes some keys from the input dictionary.

    NOTES:

    - The keys of a residue field vector dictionary are exponent vectors modulo (q-1) for some prime q.
    - The values are residue field vectors. It is known that the entries of a residue field vector
      which comes from a solution to the S-unit equation cannot have 1 in any entry.

    EXAMPLES:
    """

    garbage = []
    for a in rfv_dictionary:
        if 1 in rfv_dictionary[a]:
            garbage.append(a)
    for a in garbage:
        trash = rfv_dictionary.pop(a,0)

def construct_rfv_to_ev( rfv_dictionary, q, d, verbose_flag = False ):
    r"""
    Returns a reverse lookup dictionary, to find the exponent vectors associated to a given residue field vector.

    INPUTS:

    - rfv_dictionary -- a dictionary whose keys are exponent vectors and whose values are the associated residue field vectors
    - ``q`` -- a prime (assumed to split completely in the relevant number field)
    - ``d`` -- the number of primes in `K` above the rational prime `q`
    - verbose_flag -- a boolean flag to indicate more detailed output is desired

    OUTPUT:

    A dictionary `P` whose keys are residue field vectors and whose values are lists of all exponent vectors
    which correspond to the given residue field vector.

    NOTES:

    - For example, if rfv_dictionary[ e0 ] = r0, then P[ r0 ] is a list which contains e0.
    - During construction, some residue field vectors can be eliminated as coming from
      solutions to the S-unit equation. Such vectors are dropped from the keys of the dictionary P.

    EXAMPLES:

    """

    # Step 0. Initialize P:
    # The keys in P are just the possible first entries of a residue field vector.
    # The values (all empty lists now) will be added in the next step.

    P = {}
    P = { (v,) : [] for v in xrange(2, q) }

    # Step 1. Populate the empty lists in P[ (v,) ].
    # Loop through the keys in rfv_dictionary. For each, look at the output rf_vector.
    # Find the key in P which matches the first entry of the rf_vector.
    # Dump the **rest** of the rf_vector into a pair [exp_vec, rf_vec[1:]],
    # and append this pair into the dictionary P at the key (rf_vec[0], ).

    # Now, P[ (v,) ] = [ [a_0, e_0], [a_1, e_1], ...]
    #
    # The relationship between v, a_i, and e_i is as follows:
    #
    # a_i is an exponent vector, whose associated residue field vector is the
    # concatenation of v with e_i.

    for exponent_vector in rfv_dictionary:
        residue_field_vector = rfv_dictionary[exponent_vector]
        rf_vector_start = (residue_field_vector[0], )
        rf_vector_end = residue_field_vector[1:]
        P[rf_vector_start].append( [exponent_vector, rf_vector_end] )

    if verbose_flag:
        print "Populated P. Currently it has ", len(P), "keys."

    # Step 2: We build a new dictionary, P_new, from P.
    #
    # This is a step that will be repeated, once for each of the d primes over q.
    #
    # P is a dictionary whose keys are tuples of length m, representing the beginning of known residue field vectors.
    #
    # For any such beginning `s`,
    #
    # P[s] = [ [a_0, e_0], [a_1, e_1], ...]
    #
    # where for any exponent vector a_i, the associated residue field vector is the concatenation s + e_i.
    #
    # The dictionary P_new is constructed from the dictionary P. The new keys will be tuples of length m + 1.
    #
    # During the construction, we look for impossible entries for S-unit solutions, and drop them from the dictionary as needed.

    for j in xrange( d-1 ):
        if verbose_flag:
            print "Constructing ", j, " th place of the residue field vectors, out of ", d-1, " total."
        P_new = {}
        garbage = {}

        # we loop over each key of P.
        for rf_vector_start in P:

            # each key of P provides q-2 possible keys for P_new, which we introduce and assign an empty list.
            for w in xrange(2, q):
                new_rf_vector_start = tuple( list( rf_vector_start ) + [w] )
                P_new[ new_rf_vector_start ] = []

            # we populate P_new[ new_rf_vector_start ] using P[rf_vector_start]
            for exponent_vector, rf_vector_end in P[ rf_vector_start ]:
                new_rf_vector_end = rf_vector_end[1:]
                w = rf_vector_end[0]
                new_rf_vector_start = tuple( list( rf_vector_start ) + [w] )
                P_new[ new_rf_vector_start ].append( [exponent_vector, new_rf_vector_end] )

        if verbose_flag:
            print "P_new is populated with ", len(P_new), " keys."

        # we now loop over the keys of P_new, looking for incompatible entries.

        for rf_vector_start in P_new:
            # the final entry of rf_vector_start or rf_vector_complement_start must be < (q+3)/2.
            # No loss to insist that it is rf_vector_start.
            if rf_vector_start[-1] < (q+3)/2:
                # we find the complement to rf_vector_start:
                rf_vector_complement_start = tuple( [ q+1-j for j in rf_vector_start] )
                if P_new[ rf_vector_start ] == [] or P_new[ rf_vector_complement_start ] == []:
                    # these can't be solutions. Mark them for deletion.
                    garbage[ rf_vector_start ] = True
                    garbage[ rf_vector_complement_start ] = True

        # garbage removal
        for rf_vector_start in garbage:
            trash = P_new.pop(rf_vector_start, 0)

        if verbose_flag:
            print "After removing incompatible entries, P_new is down to ", len(P_new), " keys."

        # Time to move on to the next dictionary.
        P = P_new.copy()

    # Now, we just clean up P.
    for residue_field_vector in P:
        # at this instant, P[ residue_field_vector ] is a list of pairs: [ [a0,e0], ... ]
        # We only care about the exponent vectors a0,...
        P[residue_field_vector] = [ a[0] for a in P[residue_field_vector] ]

    if verbose_flag:
        print "Returning dictionary P with ", len(P), " keys."

    return P.copy()

def construct_comp_exp_vec( rfv_to_ev_dict, q ):
    r"""
    Constructs a dictionary associating complement vectors to residue field vectors.

    INPUT:

    - rfv_to_ev_dict -- a dictionary whose keys are residue field vectors and whose values are lists of exponent vectors with the associated residue field vector.
    - ``q`` -- the characteristic of the residue field

    OUTPUT:

    - comp_exp_vec_dict  -- a dictionary whose typical key is an exponent vector a, and whose associated value is a list of complementary exponent vectors to a.

    """

    comp_exp_vec_dict = {}
    for residue_field_vector in rfv_to_ev_dict:
        rf_vector_complement = tuple( [q + 1 - j for j in residue_field_vector] )
        exponent_vector_list = rfv_to_ev_dict[ residue_field_vector ][:]
        exponent_vector_complement_list = rfv_to_ev_dict[ rf_vector_complement ][:]
        for exponent_vector in exponent_vector_list:
            comp_exp_vec_dict[ exponent_vector ] = exponent_vector_complement_list
    return comp_exp_vec_dict

def drop_vector(ev, p, q, complement_ev_dict):
    r"""
    Determines if ev may be removed from the complement dictionary during construction.

    INPUTS:

        - ev -- an exponent vector modulo `p - 1`
        - ``p`` -- the prime such that ev is an exponent vector modulo `p-1`
        - ``q`` -- a prime, distinct from `p`, that is a key in the complement_ev_dict
        - complement_ev_dict -- a dictionary of dictionaries, whose keys are primes
          complement_ev_dict[q] is a dictionary whose keys are exponent vectors modulo `q-1`
          and whose values are lists of complementary exponent vectors modulo `q-1`

    OUTPUT:

    Returns True if ev may be dropped from the complement exponent vector dictionary, and False if not.

    NOTES:

    - If ev is not compatible with any of the vectors modulo q-1, then it can no longer correspond to a solution
      of the S-unit equation. It returns True to indicate that it should be removed.

    EXAMPLES:


    """
    # returns True if it is OK to drop exp_vec given the current comp_exp_vec dictionary associated to some q.
    # returns False otherwise
    # loop over the possible compatible vectors in the other modulus
    for compatible_exp_vec in compatible_vectors(ev, p-1, q-1):
        # do they appear in the other dictionary?
        if compatible_exp_vec in complement_ev_dict[q]:
            # OK, but the complements need to be compatible, too!
            ev_complement_list = complement_ev_dict[p][ ev ]
            for ev_comp in ev_complement_list:
                for compatible_cv in compatible_vectors( ev_comp, p-1, q-1 ):
                    if compatible_cv in complement_ev_dict[q][compatible_exp_vec]:
                        return False
    return True

def construct_complement_dictionaries(split_primes_list, SUK, verbose_flag = False):
    r"""
    A function to construct the complement exponent vector dictionaries.

    INPUT:
        - split_primes_list -- a list of rational primes which split completely in the number field `K`
        - SUK -- the S-unit group for a number field K
        - verbose_flag (optional) -- a boolean to provide additional feedback

    OUTPUT:
        - comp_exp_vec, a dictionary of dictionaries. The keys coincide with the primes in split_primes_list
          For each q, comp_exp_vec[q] is a dictionary whose keys are exponent vectors modulo q-1,
          and whose values are lists of exponent vectors modulo q-1

          If w is an exponent vector in comp_exp_vec[q][v], then the residue field vectors modulo q for
          v and w sum to [1,1,...,1]

    NOTES:
        - The data of comp_exp_vec will later be lifted to ZZ to look for true S-Unit equation solutions.
        - During construction, the various dictionaries are compared to each other several times to
          eliminate as many mod q solutions as possible.
        - The authors acknowledge a helpful discussion with Norman Danner which helped formulate this code.

    EXAMPLES:

        sage: f = x^2 + 5
        sage: H = 10
        sage: K.<xi> = NumberField(f)
        sage: SUK = K.S_unit_group( S = K.primes_above(10) )
        sage: split_primes_list = [3, 7]
        sage: construct_complement_dictionaries( split_primes_list, SUK )
        {3: {(0, 1, 0): [(1, 0, 0), (0, 1, 0)], (1, 0, 0): [(1, 0, 0), (0, 1, 0)]},
        7: {(0, 1, 0): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        (0, 1, 2): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        (0, 3, 2): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        (0, 3, 4): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        (0, 5, 0): [(0, 1, 2), (0, 3, 4), (0, 5, 0)],
        (0, 5, 4): [(1, 0, 0), (1, 4, 4), (1, 2, 2)],
        (1, 0, 0): [(0, 5, 4), (0, 3, 2), (0, 1, 0)],
        (1, 0, 2): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        (1, 0, 4): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        (1, 2, 0): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        (1, 2, 2): [(0, 5, 4), (0, 3, 2), (0, 1, 0)],
        (1, 2, 4): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        (1, 4, 0): [(1, 0, 4), (1, 4, 2), (1, 2, 0)],
        (1, 4, 2): [(1, 2, 4), (1, 4, 0), (1, 0, 2)],
        (1, 4, 4): [(0, 5, 4), (0, 3, 2), (0, 1, 0)]}}


    """

    # we define a custom function to flatten tuples for use in a later step.
    # see the definition of ev_iterator, below.

    def ev_flatten(vec):
        # turns (a, (b1,...,bn)) to (a, b1, ..., bn)
        return tuple([vec[0]] + list(vec[1]))

    # We initialize some data.

    rho = SUK.gens_values()
    rho_length = len(rho)
    rho_images_dict = {}
    rho_orders_dict = {}

    K = SUK.number_field()
    for q in split_primes_list:
        ideals_over_q, residue_fields, rho_images, product_rho_orders = sieve_ordering(SUK, q)
        rho_images_dict[q] = rho_images
        rho_orders_dict[q] = product_rho_orders

    nK = K.absolute_degree()
    w0 = rho[0].multiplicative_order()

    # We build a dictionary of dictionaries.
    # rfv_to_ev[q] is the 'mod q' residue field vector to exponent vector dictionary.

    rfv_to_ev = {}

    # We build a second dictionary of dictiories.
    # comp_exp_vec[q] is the dictionary mod q which assigns to each exponent vector
    # a list of 'complementary' exponent vectors.

    comp_exp_vec = {}

    q0 = split_primes_list[0]

    if verbose_flag:
        print "Using the following primes: ", split_primes_list
    import itertools
    for q in split_primes_list:
        rho_images = rho_images_dict[q]
        if verbose_flag:
            print "q = ", q

        def epsilon_q(a, i):
            # a is an exponent vector
            # i is an index for one of the primes over q
            # returns the value of rho_j^a_j inside the
            # residue field of Qi. (Necessarily isomorphic to F_q.)
            # rho_images[i][j] == rho[j] modulo Q[i]
            eps_value = rho_images[i][0]**a[0] % q
            for j in xrange(1, rho_length):
                eps_value = eps_value * rho_images[i][j]**a[j] % q
            return eps_value

        if verbose_flag:
            print "The evaluation function epsilon has beend defined using rho_images = ", rho_images

        # Now, we run through the vectors in the iterator, but only keep the ones
        # which are compatible with the previously constructed dictionaries. That is,
        # in order to keep an exp_vec mod q, there must exist a compatible exp_vec mod p
        # in the keys of the rfv_to_ev[p] dictionary for each completely split prime
        # p appearing prior to q in split_primes_list.

        if q == q0:
            # for the first prime, there is no filtering possible, and we just build the exponent vector
            # iterator.

            # This should consist of all vectors (a0,...,a_{t-1}), where
            # a0 is in the range 0 .. w_0 - 1 and
            # aj is in the range 0 .. q - 2   (for j > 0)

            lumpy_ev_iterator = itertools.product( xrange(w0), itertools.product( xrange(q-1), repeat = rho_length - 1))
            ev_iterator = itertools.imap(ev_flatten, lumpy_ev_iterator)

            # With the iterator built, we construct the exponent vector to residue field dictionary.

            ev_to_rfv_dict = { ev : [epsilon_q(ev, i) for i in xrange(nK) ] for ev in ev_iterator }

            if verbose_flag:
                print "The residue field dictionary currently has ", len(ev_to_rfv_dict), " exponent vector keys."

        else:
            ev_to_rfv_dict = {}
            # We use compatibility requirements to keep the size of the dictionary down.
            # Later on, we'll compare all dictionaries pairwise. But for now, we just
            # check against the first.

            # That is, rather than loop over every possible exponent vector mod q-1,
            # we only consider those evs which are compatible with the mod q0 - 1 vectors.

            # Loop over exponent vectors modulo q0 - 1
            for exp_vec_mod_q0 in comp_exp_vec[q0]:
                # Loop only over exponent vectors modulo q-1 which are compatible with exp_vec_mod_q0
                for exp_vec in compatible_vectors(exp_vec_mod_q0, q0-1, q-1):
                    # fill the dictionary with the residue field vectors using the evaluation function.
                    ev_to_rfv_dict[exp_vec] = [epsilon_q(exp_vec, i) for i in xrange(nK) ]

        if verbose_flag:
            print "The residue field dictionary currently has ", len(ev_to_rfv_dict), " exponent vector keys."

        # At this point, we now have a dictionary ev_to_rfv_dict, which attaches
        # to each exponent vector a 'residue field vector,' which is a tuple of the
        # nK values epsilon_q(a,0),...,epsilon_q(a,nK-1).

        clean_rfv_dict( ev_to_rfv_dict )

        if verbose_flag:
            print "clean_rfv_dict executed."
            print "The residue field dictionary currently has ", len(ev_to_rfv_dict), " exponent vector keys."

        # We essentially construct an inverse dictionary: one whose keys are residue field vectors,
        # and whose values are the exponent vectors that yield each key

        rfv_to_ev[q] = construct_rfv_to_ev( ev_to_rfv_dict, q, nK, verbose_flag )

        if verbose_flag:
            print "construct_rfv_to_ev executed."
            print "The rfv_to_ev dictionary currently has ", len(rfv_to_ev[q]), "rfv keys."

        comp_exp_vec[q] = construct_comp_exp_vec( rfv_to_ev[q], q )

        if verbose_flag:
            print "construct_comp_exp_vec executed."

        if verbose_flag:
            print "Size of comp_exp_vec[q]: ", len(comp_exp_vec[q])

        # Now that we have a new dictionary, we compare all the dictionaries pairwise,
        # looking for opportunities to remove 'impossible' solutions.

        for p in [qi for qi in comp_exp_vec.keys() if qi != q]:

            if verbose_flag:
                print "Comparing dictionaries for p = ", p, "and q = ", q, "."

            old_size_p = len(comp_exp_vec[p])

            if verbose_flag:
                print "Size of comp_exp_vec[p] is: ", old_size_p, "."
                cv_size = ( (q-1)/py_gcd(p-1, q-1) )**( rho_length - 1 )
                print "Length of compatible_vectors: ", cv_size, "."
                print "Product: ", old_size_p * cv_size

            for exp_vec in comp_exp_vec[p].copy():
                if drop_vector(exp_vec, p, q, comp_exp_vec):
                    trash = comp_exp_vec[p].pop(exp_vec)

            if verbose_flag:
                print "Shrunk dictionary p from ", old_size_p, " to ", len(comp_exp_vec[p])

            # Now, repeat, but swap p and q.

            old_size_q = len(comp_exp_vec[q])

            if verbose_flag:
                print "Size of comp_exp_vec[q] is: ", old_size_q, "."
                cv_size = ( (p-1)/py_gcd(p-1, q-1) )**( rho_length - 1 )
                print "Length of compatible_vectors: ", cv_size, "."
                print "Product: ", old_size_q * cv_size

            for exp_vec in comp_exp_vec[q].copy():
                if drop_vector(exp_vec, q, p, comp_exp_vec):
                    trash = comp_exp_vec[q].pop(exp_vec)

            if verbose_flag:
                print "Shrunk dictionary q from ", old_size_q, " to ", len(comp_exp_vec[q])

    return comp_exp_vec

def compatible_classes(a, m0, m1):
    r"""
    Given a congruence class a modulo m0, returns those b modulo m1 such that x = a mod m0, x = b mod m1 has a solution

    INPUTS:
        - ``a`` -- an integer
        - ``m0`` -- a positive integer
        - ``m1`` -- a positive integer

    OUTPUT:

    A list of integers b in the range 0..(m1-1) such that the Chinese Remainder Theorem problem

        x = a mod m0
        x = b mod m1

    has a solution.

    NOTES:

    - For efficiency, the solutions are not computed.
    - A necessary and sufficient condition is that a and b are congruent modulo g = gcd(m0, m1)

    EXAMPLES:

        sage: compatible_classes(2, 18, 27)
        [2, 11, 20]

    Use CRT to check the output

        sage: CRT(2, 2, 18, 27)
        2
        sage: CRT(2, 11, 18, 27)
        38
        sage: CRT(2, 20, 18, 27)
        20

    """

    g = gcd(m0, m1)
    a0 = a % g
    return [a0 + b0*g for b0 in xrange(m1/g) ]

def compatible_vectors_check( a0, a1, m0, m1):
    r"""
    Given exponent vectors with respect to two moduli, determines if they are compatible.

    INPUTS:

        - ``a0`` -- an exponent vector modulo `m0`
        - ``a1`` -- an exponent vector modulo `m1`
        - ``m0`` -- a positive integer giving the modulus of `a0`
        - ``m1`` -- a positive integer giving the modulus of `a1`

    OUTPUT:

    True if there is an integer exponent vector a satisfying

        a[0] == a0[0] == a1[0]
        a[1:] == a0[1:] modulo m0
        a[1:] == a1[1:] modulo m1

    and False otherwise.

    NOTES:

    - Exponent vectors must agree exactly in the first coordinate.
    - If exponent vectors are different lengths, an error is raised.

    EXAMPLES:

        sage: a0 = (3, 1, 8, 11)
        sage: a1 = (3, 5, 6, 13)
        sage: a2 = (5, 5, 6, 13)
        sage: a3 = (3, 1, 8)
        sage: compatible_vectors_check(a0, a1, 12, 22)
        True
        sage: compatible_vectors_check(a0, a2, 12, 22)
        False
        sage: compatible_vectors_check(a3, a0, 12, 22)
        Traceback (most recent call last):
        ...
        ValueError: Exponent vectors a0 and a1 are not the same length.
    """

    g = gcd(m0, m1)
    a0_mod_g = (x % g for x in a0)

    length = len( a0 )
    if length != len( a1 ):
        raise ValueError("Exponent vectors a0 and a1 are not the same length.")

    if a0[0] != a1[0]:
        # exponent vectors must agree exactly in the 0th coordinate.
        return False
    else:
        for j in xrange(1, length):
            if ( a0[j]-a1[j] ) % g != 0:
                return False

    # all conditions hold
    return True

def compatible_vectors(a, m0, m1):
    r"""
    Given an exponent vector a modulo m0, returns a list of exponent vectors for the modulus m1, such that a lift to the lcm modulus exists.

    INPUT:

    - a  -- an exponent vector for the modulus m0
    - m0 -- a positive integer (specifying the modulus for a)
    - m1 -- a positive integer (specifying the alternate modulus)

    OUTPUT:

    - compatible_list -- a list of exponent vectors modulo m1 which are compatible with a.

    NOTES:

    - Exponent vectors must agree exactly in the 0th position in order to be compatible.

    EXAMPLES:

        sage: a = (3, 1, 8, 1)
        sage: compatible_vectors(a, 18, 12)
        [(3, 1, 2, 1), 
        (3, 1, 2, 7), 
        (3, 1, 8, 1), 
        (3, 1, 8, 7),
        (3, 7, 2, 1), 
        (3, 7, 2, 7), 
        (3, 7, 8, 1), 
        (3, 7, 8, 7)]

    The order of the moduli matters.

        sage: len(compatible_vectors(a, 18, 12))
        8
        sage: len(compatible_vectors(a, 12, 18))
        27

    """

    # to start, recall that the 0th entry must be an exact match.
    compatible_list = [ ( a[0], ) ]

    # we now build a new list, extending the length of the compatible vectors.
    compatible_list_new = []

    for entry in a[1:]:
        compatible_entries = compatible_classes( entry, m0, m1 )
        for compatible_vector in compatible_list:
            for new_entry in compatible_entries:
                compatible_list_new.append( tuple(list(compatible_vector) + [new_entry]) )
        compatible_list = compatible_list_new
        compatible_list_new = []

    return compatible_list

def compatible_systems( split_prime_list, complement_exp_vec_dict ):
    r"""
    Given dictionaries of complement exponent vectors for various primes that split in K, compute all possible compatible systems.

    INPUT:

    - split_prime_list -- a list of rational primes that split completely in `K`
    - complement_exp_vec_dict -- a dictionary of dictionaries. The keys are primes from split_prime_list.

    OUTPUT:

    - system_list -- a list of compatible systems of exponent vectors.

    NOTES:

    - For any `q` in split_prime_list, complement_exp_vec_dict[q] is a dictionary whose keys are exponent vectors modulo `q-1`
      and whose values are lists of exponent vectors modulo `q-1` which are complementary to the key.

    - an item in system_list has the form [ [v0, w0], [v1, w1], ..., [vk, wk] ], where
        - qj = split_prime_list[j]
        - vj and wj are complementary exponent vectors modulo qj - 1
        - the pairs are all simultaneously compatible.

    - Let H = lcm( qj - 1 : qj in split_primes_list ). Then for any compatible system, there is at most one pair of integer
      exponent vectors [v, w] such that
        - every entry of v and w is bounded in absolute value by H
        - for any qj, v and vj agree modulo (qj - 1)
        - for any qj, w and wj agree modulo (qj - 1)

    EXAMPLES:

        sage: f = x^2 + 5
        sage: H = 10
        sage: K.<xi> = NumberField(f)
        sage: SUK = K.S_unit_group( S=K.primes_above(10) )
        sage: split_primes_list = [3, 7]
        sage: comp_exp_vec = construct_complement_dictionaries( split_primes_list, SUK, False )
        sage: compatible_systems( split_primes_list, comp_exp_vec )
        [[[(0, 1, 0), (0, 1, 0)], [(0, 3, 4), (0, 1, 2)]], [[(0, 1, 0), (0, 1, 0)], [(0, 3, 4), (0, 3, 4)]], [[(0, 1, 0), (0, 1, 0)], [(0, 3, 4), (0, 5, 0)]],
         [[(1, 0, 0), (1, 0, 0)], [(1, 2, 4), (1, 0, 4)]], [[(1, 0, 0), (1, 0, 0)], [(1, 2, 4), (1, 4, 2)]], [[(1, 0, 0), (1, 0, 0)], [(1, 2, 4), (1, 2, 0)]],
         [[(0, 1, 0), (0, 1, 0)], [(0, 1, 2), (0, 1, 2)]], [[(0, 1, 0), (0, 1, 0)], [(0, 1, 2), (0, 3, 4)]], [[(0, 1, 0), (0, 1, 0)], [(0, 1, 2), (0, 5, 0)]],
         [[(0, 1, 0), (1, 0, 0)], [(0, 5, 4), (1, 0, 0)]], [[(0, 1, 0), (1, 0, 0)], [(0, 5, 4), (1, 4, 4)]], [[(0, 1, 0), (1, 0, 0)], [(0, 5, 4), (1, 2, 2)]],
         [[(1, 0, 0), (1, 0, 0)], [(1, 4, 2), (1, 2, 4)]], [[(1, 0, 0), (1, 0, 0)], [(1, 4, 2), (1, 4, 0)]], [[(1, 0, 0), (1, 0, 0)], [(1, 4, 2), (1, 0, 2)]],
         [[(1, 0, 0), (1, 0, 0)], [(1, 0, 4), (1, 2, 4)]], [[(1, 0, 0), (1, 0, 0)], [(1, 0, 4), (1, 4, 0)]], [[(1, 0, 0), (1, 0, 0)], [(1, 0, 4), (1, 0, 2)]],
         [[(0, 1, 0), (1, 0, 0)], [(0, 3, 2), (1, 0, 0)]], [[(0, 1, 0), (1, 0, 0)], [(0, 3, 2), (1, 4, 4)]], [[(0, 1, 0), (1, 0, 0)], [(0, 3, 2), (1, 2, 2)]],
         [[(1, 0, 0), (0, 1, 0)], [(1, 0, 0), (0, 5, 4)]], [[(1, 0, 0), (0, 1, 0)], [(1, 0, 0), (0, 3, 2)]], [[(1, 0, 0), (0, 1, 0)], [(1, 0, 0), (0, 1, 0)]],
         [[(1, 0, 0), (1, 0, 0)], [(1, 2, 0), (1, 2, 4)]], [[(1, 0, 0), (1, 0, 0)], [(1, 2, 0), (1, 4, 0)]], [[(1, 0, 0), (1, 0, 0)], [(1, 2, 0), (1, 0, 2)]],
         [[(0, 1, 0), (1, 0, 0)], [(0, 1, 0), (1, 0, 0)]], [[(0, 1, 0), (1, 0, 0)], [(0, 1, 0), (1, 4, 4)]], [[(0, 1, 0), (1, 0, 0)], [(0, 1, 0), (1, 2, 2)]],
         [[(0, 1, 0), (0, 1, 0)], [(0, 5, 0), (0, 1, 2)]], [[(0, 1, 0), (0, 1, 0)], [(0, 5, 0), (0, 3, 4)]], [[(0, 1, 0), (0, 1, 0)], [(0, 5, 0), (0, 5, 0)]],
         [[(1, 0, 0), (0, 1, 0)], [(1, 2, 2), (0, 5, 4)]], [[(1, 0, 0), (0, 1, 0)], [(1, 2, 2), (0, 3, 2)]], [[(1, 0, 0), (0, 1, 0)], [(1, 2, 2), (0, 1, 0)]],
         [[(1, 0, 0), (1, 0, 0)], [(1, 4, 0), (1, 0, 4)]], [[(1, 0, 0), (1, 0, 0)], [(1, 4, 0), (1, 4, 2)]], [[(1, 0, 0), (1, 0, 0)], [(1, 4, 0), (1, 2, 0)]],
         [[(1, 0, 0), (1, 0, 0)], [(1, 0, 2), (1, 0, 4)]], [[(1, 0, 0), (1, 0, 0)], [(1, 0, 2), (1, 4, 2)]], [[(1, 0, 0), (1, 0, 0)], [(1, 0, 2), (1, 2, 0)]],
         [[(1, 0, 0), (0, 1, 0)], [(1, 4, 4), (0, 5, 4)]], [[(1, 0, 0), (0, 1, 0)], [(1, 4, 4), (0, 3, 2)]], [[(1, 0, 0), (0, 1, 0)], [(1, 4, 4), (0, 1, 0)]]]

    """

    S0 = split_prime_list[:]

    if len(S0) == 0:
        return []
    if len(S0) == 1:
        system_list = []
        q = S0[0]
        for exponent_vector in complement_exp_vec_dict[q]:
            for complementary_vector in complement_exp_vec_dict[q][ exponent_vector ]:
                pair = [ [exponent_vector, complementary_vector] ]
                system_list.append( pair )
        return system_list
    if len(S0) > 1:
        system_list = []
        S1 = S0[:-1]
        num_primes = len(S1)
        old_systems = compatible_systems( S1, complement_exp_vec_dict )
        q = S0[-1]
        for exp_vec in complement_exp_vec_dict[q]:
            for comp_vec in complement_exp_vec_dict[q][ exp_vec ]:
                CompatibleSystem = True
                for old_system in old_systems:
                    for j in xrange(num_primes):
                        qj = S1[j]
                        exp_vec_qj = old_system[j][0]
                        comp_vec_qj = old_system[j][1]
                        CompatibleSystem = compatible_vectors_check( exp_vec, exp_vec_qj, q-1, qj-1 )
                        if CompatibleSystem:
                            CompatibleSystem = compatible_vectors_check( comp_vec, comp_vec_qj, q-1, qj-1 )
                        if not CompatibleSystem:
                            # no reason to finish the j loop.
                            break
                    if CompatibleSystem:
                        # build the new system and append it to the list.
                        new_system = old_system + [ [exp_vec, comp_vec] ]
                        system_list.append( new_system )
        return system_list

def compatible_system_lift( compatible_system, split_primes_list ):
    r"""
    Given a compatible system of exponent vectors and complementary exponent vectors, return a lift to the integers.

    INPUTS:

        - compatible_system -- a list of pairs [ [v0, w0], [v1, w1], .., [vk, wk] ]
          where [vi, wi] is a pair of complementary exponent vectors modulo qi - 1, and all pairs are compatible.
        - split_primes_list -- a list of primes [ q0, q1, .., qk ]

    OUTPUT:

    A pair of vectors [v, w] satisfying:
        - v[0] == vi[0] for all i
        - w[0] == wi[0] for all i
        - v[j] == vi[j] modulo qi - 1 for all i and all j > 0
        - w[j] == wi[j] modulo qi - 1 for all i and all j > 0
        - every entry of v and w is bounded by L/2 in absolute value, where L is the
          least common multiple of {qi - 1 : qi in split_primes_list }

    EXAMPLES:

        sage: f = x^2 + 5
        sage: H = 10
        sage: K.<xi> = NumberField(f)
        sage: S = K.primes_above(10)
        sage: SUK = K.S_unit_group(S=S)
        sage: split_primes_list = [3, 7]
        sage: comp_exp_vec = construct_complement_dictionaries( split_primes_list, SUK, False )
        sage: comp_sys = compatible_systems( split_primes_list, comp_exp_vec )[0]

    compatible_systems returns a long list of compatible systems; we just grab the first

        sage: compatible_system_lift( comp_sys, split_primes_list )
        [(0, 3, -2), (0, 1, 2)]

    """

    if len(split_primes_list) != len( compatible_system ):
        raise ValueError("The number of primes does not match the length of the given exponent vectors.")

    m = len( split_primes_list )
    t = len( compatible_system[0][0] )

    # the first entries are already determined.
    exponent_vector_lift = ( ZZ(compatible_system[0][0][0]), )
    complement_vector_lift = ( ZZ(compatible_system[0][1][0]), )

    # fill in exponent_vector_lift
    moduli_list = [q-1 for q in split_primes_list]
    L = lcm( moduli_list )

    for i in xrange(1,t):
        exp_coord_residues = [ compatible_system[j][0][i] for j in xrange(m) ]
        comp_coord_residues = [ compatible_system[j][1][i] for j in xrange(m) ]

        ev_lift_coordinate = CRT( exp_coord_residues, moduli_list)
        cv_lift_coordinate = CRT( comp_coord_residues, moduli_list)

        # these values lie in the range [0, L-1], so we must shift them if they are bigger than L/2.

        if ev_lift_coordinate > L/2:
            ev_lift_coordinate -= L
        if cv_lift_coordinate > L/2:
            cv_lift_coordinate -= L

        exponent_vector_lift = exponent_vector_lift + ( ev_lift_coordinate, )
        complement_vector_lift = complement_vector_lift + (cv_lift_coordinate, )

    return [ exponent_vector_lift, complement_vector_lift ]

def solutions_from_systems( SUK, bound, cs_list, split_primes_list ):
    r"""
    Lifts compatible systems to the integers and returns the S-unit equation solutions the lifts yield.

    INPUT:

        - SUK -- the group of S-units where we search for solutions
        - bound -- a bound for the entries of all entries of all lifts
        - cs_list -- a list of compatible systems of exponent vectors modulo `q-1` for
                     various primes `q`
        - split_primes_list -- a list of primes giving the moduli of the exponent vectors in cs_list

    OUTPUT:

        - solutions -- a list of solutions to the S-unit equation. Entries of solutions have the form:

            [ ev, cv, iota_exp, iota_comp ]

          where:

        - ev is an exponent vector over the integers
        - cv is an exponent vector over the integers
        - every entry of ev is less than or equal to bound in absolute value
        - every entry of cv is less than or equal to bound in absolute value
        - iota_exp is the S-unit corresponding to ev
        - iota_comp is the S-unit corresponding to cv
        - iota_exp + iota_comp == 1

    EXAMPLES:

        sage: f = x^2 - 15
        sage: K.<xi> = NumberField(f)
        sage: S = K.primes_above(2)
        sage: SUK = K.S_unit_group(S=S)
        sage: split_primes_list = [7, 17]
        sage: comp_exp_vec = construct_complement_dictionaries( split_primes_list, SUK )
        sage: compatible_systems_list = compatible_systems( split_primes_list, comp_exp_vec )
        sage: solutions_from_systems( SUK, 20, compatible_systems_list, split_primes_list )
        [[(1, 0, 0), (0, 0, 1), -1, 2], [(1, 2, 0), (0, 1, 3), -8*xi - 31, 8*xi + 32],
         [(1, -2, 0), (0, -1, 3), 8*xi - 31, -8*xi + 32], [(0, -1, 3), (1, -2, 0), -8*xi + 32, 8*xi - 31],
         [(0, 0, 1), (1, 0, 0), 2, -1], [(0, -1, -3), (0, 1, -3), -1/8*xi + 1/2, 1/8*xi + 1/2],
         [(0, 0, -1), (0, 0, -1), 1/2, 1/2], [(0, 1, 3), (1, 2, 0), 8*xi + 32, -8*xi - 31],
         [(0, 1, -3), (0, -1, -3), 1/8*xi + 1/2, -1/8*xi + 1/2]]

    """

    solutions = []

    for system in cs_list:
        lift = compatible_system_lift( system, split_primes_list )
        ev = lift[0]
        cv = lift[1]
        t = len( ev )
        ValidLift = True
        for x in ev[1:]:
        # coordinates must be less than or equal to H in absolute value
            if abs(x) > bound:
                ValidLift = False
                break
        for x in cv[1:]:
            if abs(x) > bound:
                ValidLift = False
                break
        if ValidLift:
            # the entries are all below the bound, so there is nothing left to do
            # except construct the elements and see if they are solutions to
            # the S-unit equation
            iota_exp = SUK.exp( ev )
            iota_comp = SUK.exp( cv )
            if iota_exp + iota_comp == 1:
                sol = [ ev, cv, iota_exp, iota_comp ]
                solutions.append( sol )

    return solutions

def clean_sfs( sfs_list ):
    r"""
    Given a list of S-unit equation solutions, remove trivial redundancies.

    INPUT:

        - sfs_list -- a list of solutions to the S-unit equation

    OUTPUT:

        - new_sfs_list -- a list of solutions to the S-unit equation

    NOTES:

        - The function looks for cases where x + y = 1 and y + x = 1 appear
          as separate solutions, and removes one.

    EXAMPLES:

        sage: f = x^2 - 15
        sage: K.<xi> = NumberField(f)
        sage: S = K.primes_above(2)
        sage: SUK = K.S_unit_group(S=S)
        sage: split_primes_list = [7, 17]
        sage: comp_exp_vec = construct_complement_dictionaries( split_primes_list, SUK )
        sage: compatible_systems_list = compatible_systems( split_primes_list, comp_exp_vec )
        sage: suniteq_sols = solutions_from_systems( SUK, 20, compatible_systems_list, split_primes_list )
        sage: clean_sfs( suniteq_sols )
        [[(1, 0, 0), (0, 0, 1), -1, 2], [(1, 2, 0), (0, 1, 3), -8*xi - 31, 8*xi + 32],
         [(1, -2, 0), (0, -1, 3), 8*xi - 31, -8*xi + 32],
         [(0, -1, -3), (0, 1, -3), -1/8*xi + 1/2, 1/8*xi + 1/2],
         [(0, 0, -1), (0, 0, -1), 1/2, 1/2]]

    """
    # given the output from solutions_from_systems,
    # look for trivial redundancies: swapping exp_vec, comp_vec, particularly.
    new_sfs = []
    for entry in sfs_list:
        swapped_entry = [entry[1], entry[0], entry[3], entry[2]]
        repeat = False
        if entry in new_sfs or swapped_entry in new_sfs:
            repeat = True
        if not repeat:
            new_sfs.append( entry )
    return new_sfs

