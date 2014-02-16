r"""
The Matroid Union takes matroids as arguments and returns their union.

Let Mi = (Si,Ii)
Matroid Union M = M1 v M2 v M3 v ... v Mk is defined as (S1 U S2 U S3 U ... U Sk,I1 v I2 v I3 v ... v Ik )
where I1 v I2 v I3 v ... v Ik = {i1 U i2 U i3 U ... U ik | ij \in Ij }

EXAMPLES::

	sage: M1 = matroids.Uniform(3,6)
	sage: M2 = matroids.Uniform(2,4)
	sage: M = MatroidUnion(M1,M2)
	sage: M.is_valid()
	True

	sage: N = MatroidUnion(M2,M1)
	sage: M == N
	True

	sage: P = MatroidUnion([M1,M2])
	sage: Q = MatroidUnion([M1,M2],M1)


"""
from matroid import Matroid
from itertools import combinations
import sage.matrix.matrix_space as matrix_space
from sage.matrix.constructor import Matrix
from sage.graphs.all import Graph, graphs
import sage.matrix.matrix
from sage.rings.all import ZZ, QQ, FiniteField, GF
import sage.matroids.matroid
import sage.matroids.basis_exchange_matroid
from minor_matroid import MinorMatroid
from dual_matroid import DualMatroid
from rank_matroid import RankMatroid
from circuit_closures_matroid import CircuitClosuresMatroid
from basis_matroid import BasisMatroid
from linear_matroid import LinearMatroid, RegularMatroid, BinaryMatroid, TernaryMatroid, QuaternaryMatroid
import sage.matroids.utilities
from networkx import NetworkXError


def MatroidUnion(matroids=[],*args): #Only the first argument can be a iterable

	#Checking if the arguments given are valid

	if isinstance(matroids,set):
		matroids = list(matroids)

	if not isinstance(matroids,list):
		matroids = [matroids]

	matroids.extend(list(args))

	if len(matroids) < 2:
		raise ValueError("Atleast two arguments expected")

	for i in range(len(matroids)):
		if not (isinstance(matroids[i],Matroid) and matroids[i].is_valid()):
			raise ValueError("Argument "+str(i+1)+" doesn't seem to be a matroid")
		# if not (matroids[i].base_ring == ZZ or matroids[i].base_ring().is_field()):
		# 	raise NotImplementedError("MatroidUnion is only implemented for fields or integer ring")
		# TODO: Should add more checks,as some matroids don't support independent_r_sets function

	#Creating a matroid using the bases of each Matroid 

	prev_ind_sets = set([])

	prev_matroid = matroids[0]
	rank = prev_matroid.rank()
	max_rank = rank; #Max rank is used to calculate the rank of union of matroids
	ground_set = matroids[0].groundset() 

	for it in prev_matroid.independent_r_sets(rank):
		prev_ind_sets.add(it) #Adding bases of first matroid


	temp_ind_sets = set([]) #Using set to eliminate duplicates,maintains the to be added ind. sets
	

	index=1
	
	while index < len(matroids):
		curr_matroid = matroids[index]
		rank = curr_matroid.rank()
		ground_set = ground_set.union(curr_matroid.groundset()) #Extending groundset

		for it in curr_matroid.independent_r_sets(rank):
			for s in prev_ind_sets:
				n = frozenset(set(s).union(it)) 
				max_rank = max(len(n),max_rank)
				temp_ind_sets.add(n)
			#prev_ind_sets.discard(s)

		prev_ind_sets=prev_ind_sets.union(temp_ind_sets)


		index = index + 1

	bases_of_union = set([]) #Bases for the union of matroid

	for it in prev_ind_sets:
		if(len(it)==max_rank):
			bases_of_union.add(it)

	M = BasisMatroid(groundset=ground_set,bases=bases_of_union)
	return M

