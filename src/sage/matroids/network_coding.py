r"""
A scalar linear network coding achievability prover 

<Paragraph description>

AUTHORS:

- Jayant Apte (2015-07-08): initial version

EXAMPLES::

<Lots and lots of examples>

REFERENCES
==========


..  [Yeung] Raymond W. Yeung. 2008. Information Theory and Network Coding (1 ed.). Springer Publishing Company, Incorporated.
..  [Oxley] James Oxley, "Matroid Theory, Second Edition". Oxford University Press, 2011.

"""

#*****************************************************************************
#       Copyright (C) 2015 Jayant Apte <jayant91089@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject
from sage.matroids.network_coding_helpers import *

class NCinstance(SageObject):
     r"""
     
     The network coding instance class
     
     Network Coding is a paradigm for communication over networks 
     (assumed to be directed acyclic hypergraphs) with error-free links.
     A network code for a network coding instance of size `n` is a 
     collection of `n` discrete random variables, satisfying constraints 
     imposed on entropies of subsets by the network coding instance. 
     
     .. TODO::


     """
     def __init__(self, size, nsrc, constraints):
         assert size > 1
         assert nsrc >= 1
         self.size = int(size)  # No. of random variables
         self.nsrc = int(nsrc)  # No. of Sources
         self.constraints = constraints  # constraints
         
     def __repr__(self):
         if len(self.constraints) == 0:
             return "An empty network coding instance of size %s"%(self.size)
         else:
             return "A network coding instance of size %s with %s sources"%(
                     self.size, self.nsrc)
     def empty(self):
         self.constraints = [];
     def is_field_scalar_solvable(self,fieldsize):
         assert fieldsize.is_prime_power()
         return ratecertgen(self.constraints,[1]*self.size,self.nsrc,self.size,fieldsize)

     def is_achievable_rate_vector(self,rate_vector,fieldsize):
         assert fieldsize.is_prime_power()
         assert len(rate_vector)==self.size
         assert len([i for i in rate_vector if i==1])+len([i for i in rate_vector if i==0])==self.size
         return ratecertgen(self.constraints,rate_vector,self.nsrc,self.size,fieldsize)
        
         












    


