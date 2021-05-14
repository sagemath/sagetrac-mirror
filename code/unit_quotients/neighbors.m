//freeze;

/////////////////////////////////////////////////////////////////
//                                                             //
//      NEIGHBOR GRAPHS, SPINOR GENUS AND GENUS OF LATTICES    //
//              Based on Genus of a Lattice                    //
//              Bernd Souvignier, August 1997                  //  
//              Modified by David Kohel                        //
//              Last changed May 2000                          // 
//                                                             //
/////////////////////////////////////////////////////////////////

import "diagonal.m" : SpinorKernel;

forward SetDepth, Adjust, AdjoinNeighbor, ExtendNeighbors;
forward HasEvenNeighbor, Extend2Neighbors, BuildOddEven; 
forward IsNonsingularVector;

function MyLineOrbits(G)
  orb:=[ [g[2]]: g in OrbitsOfSpaces(G,1)];
  return orb;
end function;

/*******************************************************************************
				Basics
*******************************************************************************/

PREC := 12;
DUMP := true;
USE_AUT := true;

ZZ := Integers();
QQ := RationalField();

function MatrixR(X)
    if BaseRing(Parent(X)) cmpeq ZZ then
	return ChangeRing(Parent(X), QQ) ! X;
    else
	return X;
    end if;
end function;

/*******************************************************************************
			    AdjoinNeighbor
*******************************************************************************/

DTA := 0;

procedure AdjoinNeighbor(~TA, ~Lambda, L, dep)
    // Construct the neighbor of L with respect to v and 
    // append it to Lambda if it is not isometric to one 
    // of the lattices in Lambda

    N := LLL(L);
    if TA cmpne 0 then
	ts := ThetaSeries(L, PREC);
	if IsDefined(TA, ts) then
	    q := TA[ts];
	    if DUMP then vprintf Genus, 2: "\n***\n%o: %o\n", #q, ts; end if;

	    if USE_AUT then
		vprint Genus, 2: "Aut group:", #AutomorphismGroup(N);
		ao := #AutomorphismGroup(N);
	    end if;

vtime Genus, 2:
	    for i := 1 to #q do
		LL := q[i];

if USE_AUT and #AutomorphismGroup(LL) ne ao then
    vprintf Genus, 2: "Skip %o by aut\n", i, ao, #AutomorphismGroup(LL);
    continue;
end if;
//time "Centres:", #Centre(AutomorphismGroup(L)), #Centre(AutomorphismGroup(LL));
//time "Aut iso:", IsIsomorphic(AutomorphismGroup(L), AutomorphismGroup(LL));

vprintf Genus, 2: "Do IsIsometric on lattice %o/%o\n", i, #q;
vtime Genus, 2:
		if IsIsometric(N, LL: Depth := dep) then
		    vprint Genus, 2: "It IS Isometric";
		    return;
		end if;
	    end for;
	    TA[ts] := [];
	    N := CoordinateLattice(N);
	    Append(~q, N);
	    TA[ts] := q;
	else
	    if DUMP then vprintf Genus, 2: "0: %o\n", ts; end if;
	    N := CoordinateLattice(N);
	    TA[ts] := [N];
	end if;
	Append(~Lambda, N);
	vprintf Genus:
	    "New lattice number %o of minimum %o, tsl %o (memory usage %.1oM)\n", 
	    #Lambda, Minimum(N), #TA[ts], GetMemoryUsage()/10^6.0;
	return;
    end if;

vprintf Genus, 2: "Do isometry tests\n";
vtime Genus, 2:
    for i in [1..#Lambda] do
	if IsIsometric(Lambda[i], N : Depth := dep) then
	return;
    end if;
    end for;
    Append(~Lambda, CoordinateLattice(N));

    vprintf Genus: 
        "New lattice number %o of minimum %o\n", 
        #Lambda, Minimum(N); 
end procedure;

/*******************************************************************************
			    Genus funcs
*******************************************************************************/

function BinaryGenusRepresentatives(L)
   // Could do better by a factor of 2 by not computing full  
   // equivalence class group.
   c := Content(L);
   if c ne 1 then
      L := ScaledLattice(L,1/c);
   end if;
   if IsOdd(L) then
      L := ScaledLattice(L,2);
      c := c/2;
   end if;
   D := -Determinant(L); 
   Q := QuadraticForms(D);
   g := Q ! [ Norm(L.1) div 2, (L.1,L.2), Norm(L.2) div 2 ];
   G, h := ClassGroup(Q); 
   H := sub< G | [ 2*x : x in G ] >;
   S := { g * h(x) : x in H };
   return [ N : N in { LatticeWithGram( 
      c * MatrixAlgebra(RationalField(),2) !
         [ 2*f[1], Abs(f[2]), Abs(f[2]), 2*f[3] ]) : f in S } ];
end function;

function BinarySpinorRepresentatives(L)
   c := Content(L);
   if c ne 1 then
      L := ScaledLattice(L,1/c);
   end if;
   if IsOdd(L) then
      L := ScaledLattice(L,2);
      c := c/2;
   end if;
   D := -Determinant(L); 
   Q := QuadraticForms(D);
   g := Q ! [ Norm(L.1) div 2, (L.1,L.2), Norm(L.2) div 2 ];
   G, h := ClassGroup(Q); 
   H := sub< G | [ 4*x : x in G ] >;
   S := { g * h(x) : x in H };
   return [ N : N in { LatticeWithGram( 
      c * MatrixAlgebra(RationalField(),2) !
         [ 2*f[1], Abs(f[2]), Abs(f[2]), 2*f[3] ]) : f in S } ];
end function;

function CharacterValue(p,q)
   error if GCD(p,q) ne 1, "Error in character calculation.";
   if q eq -1 then
      return KroneckerSymbol(q,p);
   else 
      return KroneckerSymbol(p,q);
   end if;
end function;

function PruneCharacters(G,chars)
   // Remove the chars which are trivial.
   k := 1;
   while k le #chars do 
      I := [ 1 : i in [1..#chars] ];
      I[k] := -1;
      if I in G then
         Exclude(~chars,chars[k]);
         G := { Remove(g,k) : g in G };
      else 
         k +:= 1;
      end if;
   end while;
   return G, chars;      
end function;

function CokernelReduction(G,chars)
   repeat
      G1 := G join { [ a[i]*b[i] : i in [1..#chars] ] : a, b in G };
      G, chars := PruneCharacters(G1,chars);
   until G1 eq G;
   return G, chars;
end function;

function GenusGenerators(L, dep)
   // A small prime p for forming p-neighbor closure, and 
   // a sequence of primes (generally empty) generating the 
   // group of spinor operators, over p.

   D := Determinant(L);
   G, chars := SpinorKernel(L);
   if -1 in chars or 2 in chars then
      p := 3;
   elif D mod 2 eq 0 then
      p := 3;
   elif IsOdd(L) and not HasEvenNeighbor(L, dep) then // computes Aut
      p := 3;
   else 
      p := 2;
   end if;
   while D mod p eq 0 do
      p := NextPrime(p);
   end while;
   gens := [ p ];
   g := [ CharacterValue(p,q) : q in chars ];
   G, chars := CokernelReduction(G join { g },chars);
   while #chars gt 0 do
      p := NextPrime(p);
      while D mod p eq 0 do
         p := NextPrime(p);
      end while;
      g := [ CharacterValue(p,q) : q in chars ];
      if g notin G then
         Append(~gens,p);
         I := [ k : k in [1..#chars] | g[k] eq -1 ];
         if #I eq 1 then
            k := I[1];
            Remove(~chars,k);
            G := { Remove(g,k) : g in G };
         end if; 
         G, chars := CokernelReduction(G,chars);
      end if;
      while D mod p eq 0 do
         p := NextPrime(p);
      end while;
   end while;
   return gens;
end function;

function SpinorGenusPrime(L, dep)
   // A small prime p for forming p-neighbor closure, followed 
   // by 0 if in the spinor kernel and 1 if nontrivial spinor   
   // operator (in which case a bipartite partitioning algorithm 
   // must be used.

   D := Determinant(L);
   G, chars := SpinorKernel(L);
   if -1 in chars or 2 in chars then
      p := 3;
   elif D mod 2 eq 0 then
      p := 3;
   elif IsOdd(L) and not HasEvenNeighbor(L, dep) then // computes Aut
      p := 3;
   else 
      p := 2;
   end if;
   while D mod p eq 0 or 
      not [ Integers() | CharacterValue(p,q) : q in chars ] in G do
      p := NextPrime(p);
   end while;
   return p, 0;
   if [ Integers() | CharacterValue(p,q) : q in chars ] in G then
      return p, 0;
   end if;
   return p, 1;
end function;

intrinsic SpinorRepresentatives(L::Lat: Depth := -1, Bound := 2^32) -> SeqEnum
   {The genus of L as isometry class representatives for the
   p-neighbors of L.}

   // Compute the genus of L by exploring the neighboring graph 
   // The Depth parameter allows to specify the value of Depth 
   // for the automorphism group and isometry calculations.

   require IsExact(L) : "Argument 1 must be an exact lattice";

   L := CoordinateLattice(LLL(L));

   if Rank(L) eq 1 then
      return [ CoordinateLattice(LLL(L)) ];
   elif Rank(L) eq 2 then
      return BinarySpinorRepresentatives(L);
   end if;

   c := Content(L);
   if c ne 1 then
      L := ScaledLattice(L,1/c);
   end if;

   require Type(Depth) eq RngIntElt and Depth ge -1:
      "Parameter 'Depth' should be a non-negative integer.";
   dep := SetDepth(Rank(L),Depth);

   p, t := SpinorGenusPrime(L, dep);
   vprint Genus: "Using prime", p;
   require t eq 0 : "Choose a prime in spinor kernel.";

   if p^Rank(L) gt Bound then
      message := Sprintf("Requires computation of orbits of %o^%o points.\n"*
                         "Increase Bound parameter to proceed at your own risk.\n",
                         p, Rank(L));
      require false : message;
   elif p^Rank(L) gt (Bound div 2^8) then
      vprintf Genus : "Warning: may be slow and memory intensive\n" cat 
                      "(requires computation of orbits of %o^%o points)\n", 
                      p, Rank(L);
   end if;

TA := DTA;

   if IsOdd(L) and p eq 2 then
      return [ ScaledLattice(N,c) : 
         N in BuildOddEven(TA, [L],dep) ];
   end if;
   return [ ScaledLattice(N,c) : N in ExtendNeighbors(TA, [L], p, dep) ];
end intrinsic;

intrinsic MyGenusRepresentatives(L::Lat: Depth := -1, Bound := 2^50) -> SeqEnum
   {The genus of L as isometry class representatives for the
   p-neighbors of L.}

   // Compute the genus of L by exploring the neighboring graph 
   // The Depth parameter allows to specify the value of Depth 
   // for the automorphism group and isometry calculations.

   require IsExact(L) : "Argument 1 must be an exact lattice";

   L := CoordinateLattice(LLL(L));

   if Rank(L) eq 1 then
      return [ L ];
   elif Rank(L) eq 2 then
      return BinaryGenusRepresentatives(L);
   end if;

   c := Content(L);
   if c ne 1 then
      L := ScaledLattice(L,1/c);
   end if;

   require Type(Depth) eq RngIntElt and Depth ge -1:
      "Parameter 'Depth' should be a non-negative integer.";
   dep := SetDepth(Rank(L),Depth);

   gens := GenusGenerators(L, dep);
   p := gens[1];
   vprint Genus : "Using prime", p;
   vprint Genus : "Auxilliary generators", Exclude(gens,p);

   if p^Rank(L) gt Bound then
      message := Sprintf("Requires computation of orbits of %o^%o points.\n"*
                         "Increase Bound parameter to proceed at your own risk.\n",
                         p, Rank(L));
      require false : message;
   elif p^Rank(L) gt (Bound div 2^8) then
      vprintf Genus : "Warning: may be slow and memory intensive\n" cat 
                      "(requires computation of orbits of %o^%o points)\n", 
                      p, Rank(L);
   end if;

   TA := 0;
   TA := AssociativeArray();

   Lambda := [];
   AdjoinNeighbor(~TA, ~Lambda, L, dep);
   //Lambda := [ CoordinateLattice(LLL(L)) ];
   for i in [2..#gens] do
      n := #Lambda;
      for j in [1..n] do
         AdjoinNeighbor(~TA, ~Lambda, 
            CoordinateLattice( LLL(
               Neighbors(Lambda[j],gens[i])[1] ) ), dep );
      end for;
   end for;
   if IsOdd(L) and p eq 2 then
      // required to have an even neighbor
      return [ ScaledLattice(N,c) : N in BuildOddEven(TA, Lambda, dep) ];
   else
      return [ ScaledLattice(N,c) : N in ExtendNeighbors(TA, Lambda, p, dep) ];
   end if;
end intrinsic;

function ExtendNeighbors(TA, Lambda, p, dep)
   // Enumerate the genus of L by exploring the neighbor graph. 
   // Valid for all even L or odd p. 
   k := 1;
   while k le #Lambda do
      vprint Genus: "Candidate number", k, "of", #Lambda;
      L := Lambda[k];
      good := true;
      if (Determinant(L) mod p) eq 0 then
         good := false;
      end if; 
      G := ChangeRing( AutomorphismGroup(L : Depth := dep), GF(p));
      O := MyLineOrbits(G);
      vprint Genus: "Number of orbits:", #O;
      for o in O do
         v := Adjust(L!o[1].1,p);
         if not IsZero(v) then
            if good then
               AdjoinNeighbor(~TA, ~Lambda, Neighbor(L,v,p), dep);
            elif IsNonsingularVector(v,p) then
               AdjoinNeighbor(~TA, ~Lambda, Neighbor(L,v,p), dep);
            end if;
         end if;
      end for;
      k +:= 1;
   end while;
   return Lambda;
end function;

function Extend2Neighbors(Lambda, dep)
   // For odd L of odd determinant enumerate the genus of L.
   // If L has no even neighbor, the full neighboring graph 
   // of L is explored.
/*
   LO := CoordinateLattice(LLL(Lambda[1]));
   success, LE := HasEvenNeighbor(LO, dep); 
   if success then
      return BuildOddEven(TA, LO, LE, dep);
   end if;
   Lambda := [ LO ];
*/
   k := 1;
TA := 0;
   while k le #Lambda do
      vprint Genus: "Candidate number", k, "of", #Lambda;
      L := Lambda[k];
      G := ChangeRing( AutomorphismGroup(L : Depth := dep), GF(2));
      O := MyLineOrbits(G);
      vprint Genus: "Number of orbits:", #O;
      for o in O do
         v := L ! o[1].1;
         if Norm(v) mod 4 eq 0 and not IsZero(v) then
            // latter check to catch a bug in LineOrbit
            AdjoinNeighbor(~TA, ~Lambda, Neighbor(L,v,2), dep);
            B := [ b : b in Basis(L) | (v,b) mod 2 eq 1 ];
            if #B gt 0 then
               v +:= 2*B[1];
               AdjoinNeighbor(~TA, ~Lambda, Neighbor(L,v,2), dep);
            end if;
         end if;
      end for;
      k +:= 1;
   end while;
   return [ L : L in Lambda | IsOdd(L) ];
end function;

function HasEvenNeighbor(L, dep)
   G := ChangeRing( AutomorphismGroup(L : Depth := dep), GF(2));
   O := MyLineOrbits(G); // takes a long time
   for o in O do
      v := L ! o[1].1;
      if Norm(v) mod 4 eq 0 and not IsZero(v) then 
         // latter check to catch a bug in MyLineOrbits
         if not Norm(v) mod 8 eq 0 then
            B := [ b : b in Basis(L) | (v,b) mod 2 eq 1 ];
            if #B ne 0 then
               v +:= 2*B[1];
            end if;
         end if;
         if Norm(v) mod 8 eq 0 then
            N := LLL(Neighbor(L, v, 2));
            if IsEven(N) then
               return true, CoordinateLattice(N);
            end if;
         end if;   
      end if;
   end for;
   return false, L;
end function;

function BuildOddEven(TA, LOdd, dep)
   // For L not odd with odd determinant enumerate the genus of L
   // it is first checked whether L has an even neighbor LE.
   // In that case the genus of LE is computed and the genus of L
   // obtained from the edges of the neighboring graph of LE.
   // Otherwise, the function Odd2Genus is called which computes
   // the full neighboring graph of L.

   LEven := [ ];
   for L in LOdd do
      success, LE := HasEvenNeighbor(L, dep);
      error if not success, "No even neighbor in BuildOddEven.";
      Append(~LEven,LE);
   end for;
   k := 1;
TA := 0;
   while k le #LEven do
      vprint Genus: "Candidate number", k, "of", #LEven;
      L := LEven[k];
      G := ChangeRing( AutomorphismGroup(L : Depth := dep), GF(2));
      O := MyLineOrbits(G);
      vprint Genus: "Number of orbits:", #O ;
      for o in O do
         v := L ! o[1].1;
         if Norm(v) mod 4 eq 0 and not IsZero(v) then
            // latter check to catch a bug in MyLineOrbits
            if Norm(v) mod 8 ne 0 then
               v +:= 2*Rep({b : b in Basis(L) | (v,b) mod 2 eq 1});
            end if;
            AdjoinNeighbor(~TA, ~LEven, Neighbor(L,v,2), dep);
            v +:= 2*Rep({b : b in Basis(L) | (v,b) mod 2 eq 1});
            AdjoinNeighbor(~TA, ~LOdd, Neighbor(L,v,2), dep);
         end if;
      end for;
      k +:= 1;
   end while;
   return LOdd;
end function;

intrinsic NeighborClosure(L::Lat, p::RngIntElt : 
   Depth := -1, Bound := 2^32) -> SeqEnum
   {The closure of the lattice sequence under the p-neighbor relation.}

   // The Depth parameter allows to specify the value of Depth 
   // for the automorphism group and isometry calculations.

   require IsExact(L) : "Argument 1 must be an exact lattice";
   require IsPrime(p) : "Argument 2 is not a prime";
   p *:= Sign(p);

   L := CoordinateLattice(LLL(L));
   if Rank(L) eq 2 then
      return BinaryGenusRepresentatives(L,p);
   end if;

   c := Content(L);
   if c ne 1 then
      L := ScaledLattice(L,1/c);
   end if;

   require Type(Depth) eq RngIntElt and Depth ge -1:
      "Parameter 'Depth' should be a non-negative integer.";
   dep := SetDepth(Rank(L),Depth);

   if p^Rank(L) gt Bound then
      require false : 
         "Error: Requires computation of orbits of", p^Rank(L), "points\n"
         cat "Increase Bound parameter to proceed at your own risk.";
   elif p^Rank(L) gt (Bound div 2^8) then
      vprint Genus : "Warning: may be slow and memory intensive.";
   end if;

   if IsOdd(L) and p eq 2 then
      return [ ScaledLattice(N,c) : N in Extend2Neighbors([L], dep) ];
   end if;
   return [ ScaledLattice(N,c) : N in ExtendNeighbors([L], p, dep) ];
end intrinsic;

intrinsic NeighbourClosure(L::Lat,p::RngIntElt) -> Lat
   {The closure of L under the p-neighbour relation.}
   return NeighborClosure(L,p);
end intrinsic;

function BinaryNeighbors(L,p);
   if KroneckerSymbol(-Determinant(L),p) eq -1 then
      return [ Parent(L) | ];
   end if;
   c := Content(L);
   if c ne 1 then
      L := LatticeWithGram((1/c)*GramMatrix(L));
   end if;
   if IsOdd(L) then
      L := LatticeWithGram(2*GramMatrix(L));
      c := c/2;
   end if;
   D := -Determinant(L); 
   DK := FundamentalDiscriminant(D); 
   m := Isqrt(D div DK);
   error if GCD(p,m) ne 1, "Argument 2 must be prime to the conductor.";
   C := BinaryQuadraticForms(D);
   one := One(C);
   f := C ! [ Norm(L.1) div 2, (L.1,L.2), Norm(L.2) div 2 ];
   g := PrimeForm(C,p)^2;
   S := { C | f*g, f*g^-1 };
   return [ N : N in { LatticeWithGram( 
      c * MatrixAlgebra(Integers(),2) !
         [ 2*f[1], Abs(f[2]), Abs(f[2]), 2*f[3] ]) :
            f in S } ];
end function;

function TwoNeighbors(L, dep)
   Lambda := [ Parent(L) | ]; 
   L := CoordinateLattice(LLL(L)); 
   G := ChangeRing( AutomorphismGroup(L : Depth := dep), GF(2));
   O := MyLineOrbits(G);
   vprint Genus: "Number of orbits:", #O;
TA := 0;
   for o in O do
      v := L ! o[1].1;
      if Norm(v) mod 4 eq 0 and not IsZero(v) then
         // latter check to catch a bug in LineOrbit
         AdjoinNeighbor(~TA, ~Lambda, Neighbor(L,v,2), dep);
         B := [ b : b in Basis(L) | (v,b) mod 2 eq 1 ];
         if #B gt 0 then
            v +:= 2*B[1];
            AdjoinNeighbor(~TA, ~Lambda, Neighbor(L,v,2), dep);
         end if;
      end if;
   end for;
   return Lambda;
end function;

intrinsic Neighbours(L::Lat, p::RngIntElt : 
                    Depth := -1, Bound := 2^32) -> SeqEnum
   {The immediate p-neighbors of L.}

   require IsExact(L) : "Argument 1 must be an exact lattice";
   require IsPrime(p) : "Argument 2 is not a prime";
   p *:= Sign(p);

   if Rank(L) eq 2 then
      return BinaryNeighbors(L,p);
   end if;

   c := Content(L);
   if c ne 1 then
      L := LatticeWithGram((1/c)*GramMatrix(L));
   end if;
   require Type(Depth) eq RngIntElt and Depth ge -1:
      "Parameter 'Depth' should be a non-negative integer.";
   dep := SetDepth(Rank(L),Depth);

   if p^Rank(L) gt Bound then
      require false : 
         "Error: Requires computation of orbits of", p^Rank(L), "points\n"
         cat "Increase Bound parameter to proceed at your own risk.";
   elif p^Rank(L) gt (Bound div 2^8) then
      vprint Genus : "Warning: may be slow and memory intensive.";
   end if;

   if p eq 2 and IsOdd(L) then
      Lambda := TwoNeighbors(L, dep);
   else
      Lambda := [ Parent(L) | ]; 
      L := CoordinateLattice(LLL(L)); 
      G := ChangeRing( AutomorphismGroup(L : Depth := dep), GF(p));
      O := MyLineOrbits(G);
      good := true;
      if (Determinant(L) mod p) eq 0 then
         good := false;
      end if; 
      for o in O do
         v := Adjust(L!o[1].1,p);
         if not IsZero(v) then
            if good or IsNonsingularVector(v,p) then
               Append(~Lambda, Neighbor(L, v, p));
            end if;   
         end if;
      end for;
   end if;
   if c eq 1 then
      return Lambda;
   else
      return [ ScaledLattice(N,c) : N in Lambda ];
   end if;
end intrinsic;

function NeighborIndex(S, L, dep)
   for i in [1..#S] do
      if IsIsometric(S[i], L : Depth := dep) then
	 return i;       
      end if;
   end for;
   return 0;
end function;

intrinsic AdjacencyMatrix(G::SymGen,p::RngIntElt : 
   Depth := -1, Bound := 2^32) -> SeqEnum
   {The p-neighbor adjacency matrix with respect to the ordered 
   sequence of representatives of the genus.}

   S := Representatives(G);

   // Compute the genus of L by exploring the neighboring graph 
   // The Depth parameter allows to specify the value of Depth 
   // for the automorphism group and isometry calculations.

   if p eq 1 then 
      return One(MatrixAlgebra(Integers(),#S));
   end if;
   require IsPrime(p) : "Argument 2 must be prime";
   p *:= Sign(p);

   L := S[1]; 
   require Type(Depth) eq RngIntElt and Depth ge -1 :
      "Parameter 'Depth' should be a non-negative integer.";
   dep := SetDepth(Rank(L),Depth);

   if p^Rank(L) gt Bound then
      require false : 
         "Error: Requires computation of orbits of", p^Rank(L), "points\n"
         cat "Increase Bound parameter to proceed at your own risk.";
   elif p^Rank(L) gt (Bound div 2^8) then
      vprint Genus : "Warning: may be slow and memory intensive.";
   end if;

   A := Zero(MatrixAlgebra(Integers(),#S));
   if p eq 2 and IsOdd(L) then
      return A;
   end if;
   k := 1;
   while k le #S do
      L := S[k];
      vprint Genus: "Candidate number", k, "of", #S;
      G := ChangeRing( AutomorphismGroup(S[k] : Depth := dep), GF(p));
      O := MyLineOrbits(G);
      vprint Genus: "Number of orbits:", #O;
      good := true;
      if (Determinant(L) mod p) eq 0 then
         good := false;
      end if; 
      for o in O do
         v := Adjust(L!o[1].1,p);
         if not IsZero(v) then
            if good then
               N := LLL(Neighbor(L, v, p));
               j := NeighborIndex(S,N,dep);
               require j ne 0 : "Lattice sequence is not complete.";
               A[j,k] +:= #o;
            elif IsNonsingularVector(v,p) then  
               N := LLL(Neighbor(L, v, p));
               j := NeighborIndex(S,N,dep);
               require j ne 0 : "Lattice sequence is not complete.";
               A[j,k] +:= #o;
            end if;   
      end if;
      end for;
      k +:= 1;
   end while;
   return A;
end intrinsic;

function IsNonsingularVector(v,p)
   L := Parent(v);
   n := Rank(L);
   f := QuadraticForm(L);
   if p eq 2 and IsEven(L) then
      f := f div 2;
   end if;
   S := [ GF(p) | x : x in Eltseq(v) ];
   f := PolynomialRing(GF(p),n)!f;
   for i in [1..n] do
      if Evaluate(Derivative(f,i),S) ne 0 then
         return true;
      end if; 
   end for;
   return false;
end function;

intrinsic Neighbor(L::Lat, v::LatElt, p::RngIntElt) -> Lat
    {The p-neighbor of L with respect to v}

    require IsExact(L) : "Argument 1 must be an exact lattice";
    require IsIntegral(L) : "Argument 1 is not integral";
    require IsCompatible(L, Parent(v)) and v in L:
	"Argument 2 is not an element of argument 1";
    require IsPrime(p): "Argument 3 is not a prime";
    p *:= Sign(p);
    if IsDivisibleBy(Determinant(L), p) then
       require IsNonsingularVector(v,p) :
          "Argument 3 is singular point of quadratic form";
    end if;
    require IsDivisibleBy(ZZ!Norm(v), p^2):
	"Norm of argument 2 is not divisible by the square of argument 3";
    w := Coordinates(L!v);
    n := Degree(L);
    require not forall{i: i in [1 .. n] | IsDivisibleBy(w[i], p)}:
	"Argument 2 is in (argument 3) * (argument 1)";

    F := GramMatrix(L);
    m := Rank(L);
    n := Degree(L);

    R := GF(p);
    K := Kernel(MatrixRing(R, m)!F * RMatrixSpace(R, m, 1)!w);
    if Dimension(K) eq Degree(K) then
       return L;
    end if;
    u := Complement(Generic(K), K).1;
    C := MatrixRing(ZZ, m) !
	(
	    ChangeUniverse(Eltseq(BasisMatrix(K)), ZZ) cat
	    [p * ZZ!u[i]: i in [1 .. m]]
	);
    B := MatrixR(C) * MatrixR(BasisMatrix(L));
    B := VerticalJoin(B, 1/p * RSpace(QQ, n) ! Eltseq(v));
    return Lattice(B, MatrixR(InnerProductMatrix(L)));
end intrinsic;

intrinsic Neighbour(L::Lat, v::LatElt, p::RngIntElt) -> Lat
    {The p-neighbour of L with respect to v}
    require IsExact(L) : "Argument 1 must be an exact lattice";
    require IsIntegral(L) : "Argument 1 is not integral";
    require IsCompatible(L, Parent(v)) and v in L:
	"Argument 2 is not an element of argument 1";
    require IsPrime(p): "Argument 3 is not a prime";
    p *:= Sign(p);
    if IsDivisibleBy(Determinant(L), p) then
       require IsNonsingularVector(v,p) :
          "Argument 3 is singular point of quadratic form";
    end if;
    require IsDivisibleBy(ZZ!Norm(v), p^2):
	"Norm of argument 2 is not divisible by the square of argument 3";
    w := Coordinates(L!v);
    n := Degree(L);
    require not forall{i: i in [1 .. n] | IsDivisibleBy(w[i], p)}:
	"Argument 2 is in (argument 3) * (argument 1)";

    F := GramMatrix(L);
    m := Rank(L);
    n := Degree(L);

    R := GF(p);
    K := Kernel(MatrixRing(R, m)!F * RMatrixSpace(R, m, 1)!w);
    if Dimension(K) eq Degree(K) then
       return L;
    end if;
    u := Complement(Generic(K), K).1;
    C := MatrixRing(ZZ, m) !
	(
	    ChangeUniverse(Eltseq(BasisMatrix(K)), ZZ) cat
	    [p * ZZ!u[i]: i in [1 .. m]]
	);
    B := MatrixR(C) * MatrixR(BasisMatrix(L));
    B := VerticalJoin(B, 1/p * RSpace(QQ, n) ! Eltseq(v));
    return Lattice(B, MatrixR(InnerProductMatrix(L)));
end intrinsic;

function Adjust(v,p)
   // Adjust the vector v such that the neighbor with respect to v 
   // is integral (and even if L is even).
   // For p = 2 a vector of norm divisible by 4 is changed by a vector 
   // in 2*L such that the norm is divisible by 8.
   // For p > 2 a vector of norm divisible by p is changed by a vector 
   // in p*L such that the norm is divisible by p^2, if L is even the 
   // norm will be divisible by 2*p^2.
   // If v is not adjustable in this manner, then the zero vector 
   // is returned.

   n := Norm(v);
   if n mod p ne 0 then
      return 0*v;
   end if;
   if p eq 2 then
      if n mod 4 ne 0 then
         return 0*v;
      end if;
      if n mod 8 ne 0 then
         B := [ b : b in Basis(Parent(v)) | (v,b) mod 2 eq 1 ];
         if #B eq 0 then 
            return 0*v;
         else 
            v +:= 2 * B[1];
         end if;
      end if;
   else
      if n mod p^2 ne 0 then
         B := [ b : b in Basis(Parent(v)) | (v,b) mod p ne 0 ];
         if #B eq 0 then 
            return 0*v;
         else 
            v -:= (n*Modinv(2*(v,B[1]),p) mod p^2)*B[1]; 
         end if;
       end if;
   end if;
   return v;
end function;

function SetDepth(n,dep)
   if dep eq -1 then
      if n le 6 then
         dep := 0;
      elif n le 8 then
         dep := 1;
      elif n le 10 then
         dep := 2;
      elif n le 12 then
         dep := 3;
      else
         dep := 4;
      end if;
   elif dep gt 10 then
      dep := 10;
   end if;
   return dep;
end function;

