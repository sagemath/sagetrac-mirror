175,0
A,LatMod,1,GenusComp
S,IsPrimitiveLattice,"Check whether L is primitive at p, i.e., any Jordan decomposition consists of a 0-valuated component A and a 1-valuated component B (possibly empty), and rank(A) ge rank(B)",0,2,0,0,0,0,0,0,0,217,,0,0,LatMod,,36,-38,-38,-38,-38,-38
S,IsUnramified,check whether L is unimodular at p and Norm(L) equals 2*Scale(L) at p,0,2,0,0,0,0,0,0,0,217,,0,0,LatMod,,36,-38,-38,-38,-38,-38
S,IsUnramified,"check that IsUnramified(L, p) holds for every p over 2",0,1,0,0,0,0,0,0,0,LatMod,,36,-38,-38,-38,-38,-38
S,IsGoodBONG,"Returns true iff BONG is a good BONG at p, as defined by C. Beli",0,2,0,0,0,0,0,0,0,217,,0,0,82,,36,-38,-38,-38,-38,-38
S,IsMaximalNormSplitting,returns true iff the given list G of Gram matrices corresponds to a maximal norm splitting at p,0,2,0,0,0,0,0,0,0,217,,0,0,82,,36,-38,-38,-38,-38,-38
S,MaximalNormSplitting,"A maximal norm splitting of L at a dyadic prime p, as defined by C. Beli. Returns a Jordan decomposition into 1x1 and 2x2 components, and the corresponding list of basis vectors",0,2,0,0,0,0,0,0,0,217,,0,0,LatMod,,82,82,-38,-38,-38,-38
S,GoodBONG,"Return a good BONG of L at a dyadic prime p, as defined by C. Beli",0,2,0,0,0,0,0,0,0,217,,0,0,LatMod,,82,82,-38,-38,-38,-38
S,SpinorNorm,"The spinor norm of L at p, as calculated by C. Beli for dyadic p and by M. Kneser for non-dyadic p. Returns a subspace of LocalMultiplicativeGroupModSquares(p), and a boolean which is true iff the spinor norm is exactly the units",0,2,0,0,0,0,0,0,0,217,,0,0,LatMod,,159,175,36,-38,-38,-38
S,SpinorsNotExactlyTheUnits,those places of BaseRing(L) where the Spinor norm is not exactly the units,0,1,0,0,0,0,0,0,0,LatMod,,82,-38,-38,-38,-38,-38
S,MapIdeleIntoClassGroup,"map an idele into the class group associated to L. The parameter Idele must be a list of tuples <p, a_p>, where p is a prime of BaseRing(L), and a_p is an element of K^* (interpreted as an element of the completion K^*_p). The parameter AtInfinity can be a list of tuples <v, +1 or -1>, where v is an element of RealPlaces(NumberField(L)). All places, finite or infinite, which are unspecified are interpreted as 1",0,2,0,0,0,0,0,0,0,82,,0,0,LatMod,,108,-38,-38,-38,-38,-38
S,MapPrincipalIdeleIntoClassGroup,Map the principal idele defined by the element 'a' into the ray class group identified with the proper spinor genera in Genus(L),0,2,0,0,0,0,0,0,0,-45,,0,0,LatMod,,108,-38,-38,-38,-38,-38
S,PrepareClassGroup,internal use,0,1,0,0,1,0,0,0,0,LatMod,,-38,-38,-38,-38,-38,-38
S,IteratedNeighbours,The iterated neighbours of L at the prime p,0,2,0,0,0,0,0,0,0,217,,0,0,LatMod,,82,-38,-38,-38,-38,-38
S,SpinorGeneraInGenus,A sequence of lattices representing the spinor genera in the genus of L,0,1,0,0,0,0,0,0,0,LatMod,,82,-38,-38,-38,-38,-38
S,GenusRepresentatives,A list of representatives of the isometry classes in the genus of L,0,1,0,0,0,0,0,0,0,LatMod,,82,-38,-38,-38,-38,-38
S,IsOneClassGenus,"Returns (#GenusRepresentatives(L) eq 1), but is much faster than calling GenusRepresentatives",0,1,0,0,0,0,0,0,0,LatMod,,36,-38,-38,-38,-38,-38
