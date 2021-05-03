freeze;

////////////////////////////////////////////////////////////////////
//                                                                //
//         ORTHOGONAL DECOMPOSITION, P-ADIC DIAGONALIZATION       //
//          P-ADIC AUTOMORPHOUS NUMBERS, AND SPINOR KERNEL        //
//                        David Kohel                             //
//                  Created November 1999                         // 
//                  Last modified May 2000                        // 
//                                                                //
////////////////////////////////////////////////////////////////////

forward SpinorKernel, CokernelReduction, PruneCharacters;
forward BinarySpinorGenerators, UnitSquareClass, FieldSquareClass; 
forward pAdicDiagonals, FindRootVector; 
forward TwoAdicDiagonals, Find2RootVector, Jordan2Reduction;

intrinsic IsSpinorNorm(G::SymGen,r::RngIntElt) -> BoolElt
   {True if and only if r is coprime to 2 and the determinant
   and is the norm of an element of the spinor kernel.}
   xgens := SpinorCharacters(G);
   return &and[ Evaluate(x,r) eq 1 : x in xgens ];
end intrinsic; 

intrinsic SpinorGenerators(G::SymGen : Proper) -> SeqEnum
   {A sequence of primes which generate the spinor group.}
   if Rank(G) eq 2 then
      D:= -Determinant(G);
      if IsOdd( Representative(G) ) then D:= 4*D; end if;
      return BinarySpinorGenerators(D);
   end if;
   xgens := SpinorCharacters(G : Proper:= Proper); 
   X := Universe(xgens); 
   if #xgens eq 1 and xgens[1] eq X!1 then
      return [ Integers() | ];
   end if; 
   m := Modulus(Universe(xgens));
   d := Determinant(G);
   if m mod 2 eq 0 and d mod 2 eq 1 then d *:= 2; end if; 
   V := VectorSpace(GF(2),#xgens);
   U := sub< V | V!0 >;
   pgens := [ Integers() | ];   
   p := 2;
   while U ne V do
      if d mod p ne 0 then
         w := V![ (1-Evaluate(x,p)) div 2 : x in xgens ];
         if w notin U then
            U := sub< V | Basis(U) cat [ w ] >; 
            Append(~pgens,p);
         end if;
      end if;
      p := NextPrime(p);
   end while;
   return pgens;
end intrinsic

function BinarySpinorGenerators(D)
   // The goal is to compute primes p such that PrimeForm(Q,p)^2's
   // generate Cl(Q)^2/Cl(Q)^4, where Q is QuadraticForms(D).  We 
   // reduce to the 2-subgroup Cl2(Q) of Cl(Q), compute the order h2 
   // of the subgroup Cl2(Q)^2, and enumerate the set of elements S of
   // Cl2(Q)^4.  For each new generator found, we expand the set S.
   Q := QuadraticForms(D);
   // MW (Jan 2018): one can use Bosma-Stevenhagen here to get 2-part...
   A, f := ClassGroup(Q);
   r := #[ i : i in [1..Ngens(A)] | Order(A.i) mod 4 eq 0 ];
   gens := [ Integers() | ];
   if r eq 0 then return gens; end if;
   // Compute the order of Cl2(Q)^2.
   h2 := &*[ n div 2^Max(0,Valuation(n,2)-1) 
             where n is Order(A.i) : i in[1..Ngens(A)] ];
   // Compute the exponents for Cl2(Q)^4, and the abelian 
   // subgroup mapping to Cl2(Q)^4.
   e4 := [ n div 2^Max(0,Valuation(n,2)-2) 
           where n is Order(A.i) : i in[1..Ngens(A)] ];
   A4 := sub< A | [ e4[i]*A.i : i in [1..Ngens(A)] ] >;
   S := [ f(x) : x in A4 ];
   p := 2;
   c2 := #A div 2^(Valuation(#A,2)-1);
   while #gens lt r do   
      if KroneckerSymbol(D,p) eq 1 then
         f := PrimeForm(Q,p)^c2;
         if f notin S then
   	    Append(~gens,p); 
            if #gens lt r then
   	       S cat:= [ f*g : g in S ]; 
            end if;
         end if;
      end if;
      p := NextPrime(p);
   end while;
   return gens;
end function;

intrinsic SpinorCharacters(G::SymGen : Proper) -> SeqEnum
   {A sequence of Dirichlet characters dual to the spinor generators.}
   if Rank(G) eq 1 then return [ DirichletGroup(1,Integers()) | ]; end if;
   require Rank(G) gt 2 : "Argument can not have rank 2 -- \n" * 
      "spinor characters not given by Dirichlet characters";
   L := G`Representative;
   S, M := SpinorKernel(L : Proper:= Proper); 
   m := Abs(&*M);
   if 2 in M and -1 notin M then 
      M := [-1] cat M; 
      S := [ [1] cat s : s in S ];
   end if;
   if -1 in M then
      m *:= 4;
   end if; 
   X := DirichletGroup(m,Integers());
   V := VectorSpace(GF(2),Ngens(X));
   B := [ V | [ (1-s[i]) div 2 : i in [1..#s] ] : s in S ];
   W := Complement(V,sub< V | B >);
   chars := [ X | ];
   for e in Basis(W) do
      x := X!1;
      for i in [1..#M] do
         x *:= X.i^(Integers()!e[i]);
      end for;
      Append(~chars,x);
   end for;
   return chars;
end intrinsic;

function NormGeneratorAtp(L, p)
  n:= Rank(L);
  if p eq -1 then return L.1; end if;
  m:= Infinity();
  for i in [1..n] do
    val:= Valuation(Norm(L.i), p);
    if val lt m then v:= L.i; m:= val; end if;
  end for;
  for i in [1..n], j in [1..i-1] do
    val:= Valuation(Norm(L.i+L.j), p);
    if val lt m then v:= L.i+L.j; m:= val; end if;
  end for;
  return v;
end function;

// Compute a vector v in L such that (v,v) generates Norm(L_p)
// for all p in P.
function NormGenerator(L, P)
  n:= Rank(L);
  assert n gt 0;
  if #P eq 0 then return L.1; end if;
  N:= Norm(L);
  d:= 2 in P select 2 else P[1];
  v:= NormGeneratorAtp(L, d);
  nrm:= Norm(v);
  for p in P do
    if Valuation(nrm, p) ne Valuation(N, p) then
      w:= NormGeneratorAtp(L, p);
      // p is odd, so v+w or v-w does the trick.
      l:= Valuation(Norm(v+w), p) eq Valuation(N, p) select 1 else -1;
      l:= CRT([0, l], [d, p]);
      v:= v + l*w;
      nrm:= Norm(v);
    end if;
    d:= p*d;
  end for;
  assert forall{p: p in P | Valuation(nrm, p) eq Valuation(N, p)};
  return v;
end function;

function SpinorKernel(L : Proper)
   // {The image of the spinor kernel, followed by the sequence 
   // of primes p defining coordinatewise characters r |-> 
   // KroneckerSymbol(r,p) or if p = -1, then r |-> 
   // KroneckerSymbol(-1,r).}
   c := Content(L);
   if c ne 1 then
      L := ScaledLattice(L,1/c);
   end if;
   prms := [ p[1] : p in Factorization(2*Determinant(L)) ];
   chars := [ -1 ] cat prms;
   gens := { [ 1 : i in [1..#chars] ] };
   for p in prms do 
      A := AutomorphousClasses(L,p);
      // Remove "irrelevant" primes. 
      if p eq 2 and A eq {1,3,5,7} then
         k := Index(chars,-1);
         if k ne 0 then
            Remove(~chars,k);
            gens := { Remove(g,k) : g in gens };
         end if;
         k := Index(chars,2);
         if k ne 0 then
            Remove(~chars,k);
            gens := { Remove(g,k) : g in gens };
         end if;
      elif IsOdd(p) and #A eq 2 and 
         0 notin { a mod p : a in A } then
         k := Index(chars,p);
         if k ne 0 then
            Remove(~chars,k);
            gens := { Remove(g,k) : g in gens };
         end if;
      else
         // Remove "tractable" characters.
         if IsOdd(p) and #A eq 4 then
            k := Index(chars,p);
            if k ne 0 then
               Remove(~chars,k);
               gens := { Remove(g,k) : g in gens };
            end if;
         end if;
         for a in Exclude(A,1) do
            e := Valuation(a,p);
            u := a div p^e; 
            g := [ 1 : i in [1..#chars] ];
            for j in [1..#chars] do
               q := chars[j];
               if q eq -1 and p eq 2 then
                 g[j] := KroneckerSymbol(-1,a);
               elif q eq -1 then
                 g[j] := KroneckerSymbol(-1,p)^e;
               elif q ne p then
                  g[j] := KroneckerSymbol(p,q)^e;
               else
                  g[j] := KroneckerSymbol(u,q);
               end if;
            end for;
            I := [ k : k in [1..#chars] | g[k] eq -1 ];
            if #I eq 1 then
               k := I[1];
               Remove(~chars,k);
               gens := { Remove(g,k) : g in gens };
            elif #I ge 2 then
               Include(~gens,g);
            end if; 
         end for;
      end if; 
   end for;

   // Any reflection acts on the proper spinor genera.
   // The orbits are the full spinor genera.
   // We simply compute such a reflection v of norm coprime to 2*det(L).
   if IsEven(Rank(L)) and not Proper then
     v:= NormGenerator(L, prms);
     T:= Matrix([ x - 2*(x,v)/(v,v) * v : x in Basis(L) ]);
     LL:= Lattice( Matrix(Rationals(), T), Matrix(Rationals(), GramMatrix(L)) );
     r:= Index(L, L meet LL);
     assert GCD(r, 2*Determinant(L)) eq 1;
     Include(~gens, [ q eq -1 select KroneckerSymbol(-1, r) else KroneckerSymbol(r,q) : q in chars ]);
   end if;

   return CokernelReduction(gens,chars);
end function;

function CokernelReduction(G,chars)
   repeat
      G1 := G join { [ a[i]*b[i] : i in [1..#chars] ] : a, b in G };
      G, chars := PruneCharacters(G1,chars);
   until G1 eq G;
   return G, chars;
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

intrinsic AutomorphousClasses(G::SymGen,p::RngIntElt) -> SetEnum
   {A set of integer representatives of the p-adic automorphous 
   square classes for the genus G (see Conway & Sloane).}
   require IsPrime(p) : "Argument 2 is not a prime";
   return AutomorphousClasses(G`Representative,p);
end intrinsic;

intrinsic AutomorphousClasses(L::Lat,p::RngIntElt) -> SetEnum
   {A set of integer representatives of the p-adic automorphous 
   square classes for the lattice L (see Conway & Sloane).}
   // This follows Conway-Sloane p.390-392 (section 9.4).

   require IsExact(L) : "Argument 1 must be an exact lattice";
   require IsPrime(p) : "Argument 2 is not a prime";
   p *:= Sign(p);

   L := CoordinateLattice(LLL(L));
   L := ScaledLattice(L,1/Content(L));
   if IsOdd(p) then
      D := pAdicDiagonals(L,p);
      S := &join [ { FieldSquareClass(D[i]/D[j],p) 
         : j in [i+1..#D] } : i in [1..#D] ];
      if #{ Valuation(a,p) : a in D } lt Rank(L) then
         for u in [1..p] do
            if KroneckerSymbol(u,p) eq -1 then
               S join:= { 1, u };
               break u;
            end if;
         end for;  
      end if;
   else 
      if p eq 8 then p := 2; end if;
      D1, D2 := pAdicDiagonals(L,2);
      List1 := D1;
      V1 := [ Valuation(a,2) : a in List1 ];
      V2 := [ Valuation(GCD([s[i] : i in [1..3] ]),2)  : s in D2 ];
      List2 := [ u*2^(n+1) : n in V2, u in {1,3,5,7} ];
      TotList := List1 cat List2;
      // Ratio square classes of total list.
      S := { FieldSquareClass(a,2) : a in &cat[ [ TotList[j]/TotList[i] 
               : j in [i+1..#TotList] ] : i in [1..#TotList] ] };
      // Supplementary rules depend only on List1.
      R1 := { UnitSquareClass(a,2) : a in &cat[ [ List1[j]/List1[i] 
                  : j in [i+1..#List1] ] : i in [1..#List1] ] };
      U1 := [ Valuation(a,2) : a in R1 ]; 
      for i in [1..#V1-2] do
         if V1[i+2] le (V1[i] + 3) then
            S join:= { 1, 3, 5, 7 }; 
         end if;
      end for;    
      if 1 in R1 then
         S join:= { 2 }; 
      end if;
      if 5 in R1 then
         S join:= { 6 }; 
      end if;
      if &or[ a in R1 : a in [2,8,10,40] ] then
         S join:= { 3 }; 
      end if;
      if &or[ n in U1 : n in [0,2,4] ] then
         S join:= { 5 }; 
      end if;
      if &or[ a in R1 : a in [6,14,24,56] ] then
         S join:= { 7 }; 
      end if;
   end if;     
   return { FieldSquareClass(a*b,p) : a, b in S };
end intrinsic;

function UnitSquareClass(n,p)
   q := p^Valuation(n,p);
   u := Numerator(n/q)*Denominator(n/q);
   if p eq 2 then
     return (u mod 8)*q; 
   elif KroneckerSymbol(u,p) eq 1 then
     return q;
   end if;
   for a in [1..p] do
      if KroneckerSymbol(a,p) eq -1 then
         return a*q;
      end if;
   end for;   
end function;

function FieldSquareClass(n,p)
   e := Valuation(n,p); 
   n *:= 1/p^e;
   n := Numerator(n)*Denominator(n);
   if p eq 2 then
      return (n mod 8)*p^(e mod 2); 
   elif KroneckerSymbol(n,p) eq 1 then
      return p^(e mod 2);
   end if;
   for a in [1..p] do
      if KroneckerSymbol(a,p) eq -1 then
         return a*p^(e mod 2);
      end if;
   end for;   
end function;

intrinsic pAdicDiagonalization(L::Lat,p::RngIntElt) -> Lat
   {A lattice diagonalized lattice p-adically equivalent to L;
   if p is 2 then the diagonalized form may have Jordan blocks 
   of dimension 2. }
   require IsPrime(p) : "Argument 2 must be prime";
   if p eq 2 then
      D1, D2 := TwoAdicDiagonals(L);
      if #D1 eq 0 then
         return Jordan2Reduction( 
             DirectSum([ LatticeWithGram(
                MatrixAlgebra(Universe(f),2) 
                   ! [f[1],f[2],f[2],f[3]] ) : f in D2 ]) );  
      else 
         L1 := LatticeWithGram(
            DiagonalMatrix(MatrixAlgebra(Universe(D1),#D1),D1) );
      end if;
      if #D2 eq 0 then
         return Jordan2Reduction( L1 );
      else 
         L2 := DirectSum([ LatticeWithGram(
            MatrixAlgebra(Universe(f),2)
              ! [f[1],f[2], f[2],f[3]]) : f in D2 ]);  
      end if;
      return Jordan2Reduction( DirectSum(L1,L2) );
   end if;
   D := Sort( [ UnitSquareClass(a,p) : a in pAdicDiagonals(L,p) ] ); 
   for i in [2..#D] do
      if D[i-1] eq D[i] then
         q := p^Valuation(D[i],p); 
         D[i] := UnitSquareClass(D[i-1]*(D[i]/q),p);  
         D[i-1] := q;
      end if; 
   end for;
   return LatticeWithGram(
       DiagonalMatrix(MatrixAlgebra(Universe(D),#D),D) );
end intrinsic;

function pAdicDiagonals(L,p) 
   // Pre: L a lattice, p a prime.
   // Post: The diagonal entries of the diagonalization of L. 
   // For p = 2, the diagonals and quadratic forms <a,b,c>. 

   error if not IsPrime(p), "Argument 2 must be prime";
   if p eq 2 then
      return TwoAdicDiagonals(L);
   end if;
   if Rank(L) eq 1 then
      return [ Rationals() | Norm(L.1) ];
   end if;
   c := Content(L);
   L := LatticeWithGram((1/c)*GramMatrix(L));
   c := UnitSquareClass(c,p);
   v := L.1;
   n := (v,v);
   if (n mod p) ne 0 then
      B := [ n*L.i - (L.i,v)*v : i in [2..Rank(L)] ];
   else
      v := FindRootVector(L,p);
      G, f := quo< L | v >;
      n := (v,v);
      B := [ n*u - (u,v)*v : u in { g@@f : g in Generators(G) } ];
   end if;
   M := LatticeWithGram( MatrixAlgebra(RationalField(),Rank(L)-1) !
      [ (B[i],B[j]) : i, j in [1..Rank(L)-1] ] );
   return Sort([ RationalField() | UnitSquareClass(c*n,p) : 
            n in [ Norm(v) ] cat pAdicDiagonals(M,p) ]);
end function;

function TwoAdicDiagonals(L)
   K := BaseRing(L);
   if Rank(L) eq 1 then
      return [ K | UnitSquareClass((L.1,L.1),2) ], 
             [ Parent([ K | ]) | ]; 
   end if;
   c := Content(L);
   L := ScaledLattice(L,1/c);
   c := UnitSquareClass(c,2);
   q := 2^Valuation(c,2);   
   u := Find2RootVector(L);
   n := Norm(u);
   G, f := quo< L | u >;
   B := [ g@@f : g in Generators(G) ];
   if IsOdd(n) then
      B1 := [ ];
      B := [ n*v - (u,v)*u : v in B ];
      t := 1;
   else 
      m := n div 2;
      B1 := [ v : v in B | IsOdd((u,v)) ];
      if #B1 eq 0 then
         B := [ m*w - ((u,w) div 2)*u : w in B ];
         t := 1;
      else
         v := B1[1];
         B := [ m*(w - v) - ((u,w-v) div 2)*u : w in B1 | w ne v ] cat 
              [ m*w - ((u,w) div 2)*u : w in B | w notin B1 ]; 
         t := (u,v);
         d := t^2 - (u,u)*(v,v);
         B := [ d*w - (v,w)*((u,v)*u - (u,u)*v) : w in B ];
         t := 2;
      end if; 
   end if;
   if Rank(L) eq t then /* t equals 2 */ 
      D := -(((u,u)*(v,v)-(u,v)^2) mod 8); 
      if D eq -3 then
         return [ K | ], [ [ K | 2*q,q,2*q ] ]; 
      else
         return [ K | ], [ [ K | 2*q,q,4*q ] ]; 
      end if;
   end if;
   M := LatticeWithGram( c*MatrixAlgebra(K,Rank(L)-t) !
          [ (B[i],B[j]) : i, j in [1..Rank(L)-t] ] );
   D1, D2 := TwoAdicDiagonals(M);
   ChangeUniverse(~D1, K);
   if t eq 1 then
      Insert(~D1,1,UnitSquareClass(c*(u,u),2));
      return &cat[ Sort([ a : a in D1 | Valuation(a,2) eq i ]) : 
         i in {@ Valuation(a,2) : a in D1 @} ], D2;
   else /* t equals 2 */
      D := -(((u,u)*(v,v)-(u,v)^2) mod 8); 
      if D eq -3 then
         Insert(~D2,1,[ K | 2*q,q,2*q ]);
      else /* D eq -7 */
         Insert(~D2,1,[ K | 2*q,q,4*q ]);
      end if;
      D1 := &cat[ Sort([ a : a in D1 | Valuation(a,2) eq i ]) : 
         i in {@ Valuation(a,2) : a in D1 @} ];  
      return D1, D2; 
   end if;
end function;

function FindRootVector(L,p)
   for i in [1..Rank(L)] do
      if ((L.i,L.i) mod p) ne 0 then
         return L.i;
      end if;
      for j in [i+1..Rank(L)] do
         t := (L.i,L.j);
         if (t mod p) ne 0 then
            if ((2*t + (L.j,L.j)) mod p) ne 0 then 
               return L.i + L.j; 
            else 
               return L.i - L.j; 
            end if;
         end if;
      end for;
   end for;
end function;

function Find2RootVector(L)
   if IsOdd(L) then
      for i in [1..Rank(L)] do
         if (Norm(L.i) mod 2) eq 1 then
            return L.i;  
         end if;
      end for; 
   end if;
   for i in [1..Rank(L)] do
      if ((L.i,L.i) mod 4) ne 0 then
         return L.i;
      end if;
      for j in [i+1..Rank(L)] do
         if IsOdd((L.i,L.j)) then
            if (((L.j,L.j) + 2*(L.i,L.j)) mod 4) ne 0 then 
               return L.i + L.j; 
            else 
               return L.j; 
            end if;
         end if;
      end for;
   end for;
end function;

////////////////////////////////////////////////////////////////////
//                                                                //
//      Canonical Representatives of 2-Adic Jordan Forms          //
//                                                                //
////////////////////////////////////////////////////////////////////

forward 
   DiagonalizingReduction, NonDiag0Reduction, NonDiag1Reduction,
   Diag0Term2Reduction, Diag0Term4Reduction, 
   Diag1Term2Reduction, Diag1Term3Reduction, Diag1Term4Reduction,
   Diag2Term2Reduction, Diag2Term3Reduction;


Diags := [ [1], [3], [5], [7] ];
NSplt := [2,1,1,2];
Split := [2,1,1,4];

function Jordan2Reduction(L)
   // L is a lattice representing a 2-adic Jordan decomposition.    
   // D are integer sequences representing the Jordan blocks, 
   // Valid blocks are:
   //    [1], [3], [5], [7], [2,1,1,2], [2,1,1,4].
   // I are the scales, log base 2.
   M := GramMatrix(L);
   I := [ IntegerRing() | ]; 
   D := [ Parent(I) |  ];
   n := Degree(Parent(M));
   i := 1;
   while i le n do
      if i lt n and M[i,i+1] ne 0 then
         k := Valuation(M[i,i+1],2);
         q := 2^k;
         r := Max([0] cat [ j : j in [1..#I] | I[j] le k ]);
         Insert(~I,r+1,k);
         Insert(~D,r+1,[ IntegerRing() |
            M[i,i]/q, M[i,i+1]/q, M[i+1,i]/q, M[i+1,i+1]/q ] );
         i +:= 2; 
      else
         k := Valuation(M[i,i],2);
         q := 2^k;
         Append(~I,k);
         Append(~D,[ IntegerRing() | M[i,i]/q ]);
         i +:= 1; 
      end if;
   end while;
   D, I := DiagonalizingReduction(D,I);
   ScaleVals := {@ i : i in I @}; 
   repeat
      stable := true;
      for i0 in ScaleVals do
         I0 := [ j : j in [1..#I] | I[j] eq i0 ]; 
         I1 := [ j : j in [1..#I] | I[j] eq i0+1 ]; 
         D, changed := NonDiag0Reduction(D,I0);
         stable and:= not changed;
         D, changed := NonDiag1Reduction(D,I0,I1);
         stable and:= not changed;
      end for;
   until stable;
   repeat
      stable := true;
      for i0 in ScaleVals do
         I0 := [ j : j in [1..#I] | #D[j] eq 1 and I[j] eq i0 ]; 
         I1 := [ j : j in [1..#I] | #D[j] eq 1 and I[j] eq i0+1 ]; 
         I2 := [ j : j in [1..#I] | #D[j] eq 1 and I[j] eq i0+2 ]; 
         D, changed := Diag0Term2Reduction(D,I0);
         stable and:= not changed;
/*
         if changed then 
            print "Diag0Term2";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
         D, changed := Diag0Term4Reduction(D,I0);
         stable and:= not changed;
/*
         if changed then 
            print "Diag0Term4";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
         D, changed := Diag1Term2Reduction(D,I0,I1);
         stable and:= not changed;
/*
         if changed then 
            print "Diag1Term2";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
         D, changed := Diag1Term3Reduction(D,I0,I1);
         stable and:= not changed;
/*
         if changed then 
            print "Diag1Term3";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
         D, changed := Diag1Term4Reduction(D,I0,I1);
         stable and:= not changed;
/*
         if changed then 
            print "Diag1Term4";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
         D, changed := Diag2Term2Reduction(D,I0,I2);
         stable and:= not changed;
/*
         if changed then 
            print "Diag2Term2";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
         D, changed := Diag2Term3Reduction(D,I0,I1,I2);
         stable and:= not changed;
/*
         if changed then 
            print "Diag2Term3";
            print "D =", [ 2^I[i]*D[i][1] : i in [1..#D] ];
         end if;  
*/
      end for; 
   until stable;
   D := &cat[ 
         Sort([ D[j] : j in [1..#D] | I[j] eq i0 and #D[j] eq 1 ]) cat 
         Sort([ D[j] : j in [1..#D] | I[j] eq i0 and #D[j] eq 4 ]) 
      : i0 in ScaleVals ];
   // Since the zero lattice doesn't exist, must manually 
   // initialize a lattice.
   Mat1 := MatrixAlgebra(RationalField(),1);
   Mat2 := MatrixAlgebra(RationalField(),2);
   if #D[1] eq 1 then 
      L := LatticeWithGram(2^I[1]*Mat1!D[1]);
   else 
      L := LatticeWithGram(2^I[1]*Mat2!D[1]);
   end if;
   // Continue...
   for i in [2..#D] do
      if #D[i] eq 1 then 
         L := DirectSum(L, LatticeWithGram(2^I[i]*Mat1!D[i]));
      else 
         L := DirectSum(L, LatticeWithGram(2^I[i]*Mat2!D[i]));
      end if;
   end for;
   return L;
end function;

// BEGIN REDUCTION ALGORITHMS:

function DiagonalizingReduction(D,I)
   for i0 in I do
      I0 := [ j : j in [1..#I] | I[j] eq i0 ];
      E0 := [ D[i] : i in I0 ];
      for R in Diags do  
         while Split in E0 and R in E0 do
            k := I0[Index(E0,Split)];
            D[k] := R;
            Insert(~D,k,[ 7*R[1] mod 8]);
            Insert(~I,k,i0);
            I0 := [ j : j in [1..#I] | I[j] eq i0 ];
            E0 := [ D[j] : j in I0 ];
         end while;
         while NSplt in E0 and R in E0 do
            j := I0[Index(E0,R)];
            k := I0[Index(E0,NSplt)];
            D[j] := [ 3*R[1] mod 8];
            D[k] := [ 3*R[1] mod 8];
            Insert(~D,k,[ 3*R[1] mod 8]);
            Insert(~I,k,i0);
            I0 := [ j : j in [1..#I] | I[j] eq i0 ];
            E0 := [ D[j] : j in I0 ];
         end while;
      end for;
   end for;
   return D, I;
end function;

/*

// Diagonalizing relations (0-increment):
[
    u1*x7 - x3^2*x1,
    u1*x5 - x3*x1^2,
    u1*x3 - x5*x3^2,
    u1*x1 - x5*x3*x1,
    u0*x7 - x5*x1^2,
    u0*x5 - x7*x3^2,
    u0*x3 - x1^3,
    u0*x1 - x3^3,
];

*/

function NonDiag0Reduction(D,I0)
   Split := [2,1,1,4];
   NSplt := [2,1,1,2];
   IS := [ j : j in I0 | D[j] eq Split ];
   changed := false;
   while #IS ge 2 do
      j := IS[1]; Exclude(~IS,j); 
      k := IS[1]; Exclude(~IS,k); 
      D[j] := NSplt; 
      D[k] := NSplt;
      changed := true;
   end while;  
   return D, changed;
end function;

/*
// 0-increment non-diagonalizable relations:
[
    u1^2 - u0^2
];
*/

function Diag0Term2Reduction(D,I0)
   changed := false;
   if #I0 lt 2 then
      return D, changed;
   end if;
   for R in [ [5], [7] ] do
      J0 := [ j : j in I0 | D[j] eq R ];
      while #J0 ge 2 do
         j := J0[1]; Exclude(~J0,j); 
         k := J0[1]; Exclude(~J0,k); 
         D[j] := [(5*R[1]) mod 8];
         D[k] := [(5*R[1]) mod 8];
         changed := true;
      end while;
   end for;
   E0 := [ D[j] : j in I0 ];
   while [5] in E0 and [7] in E0 do
      j := I0[Index(E0,[5])]; 
      k := I0[Index(E0,[7])]; 
      D[j] := [1]; 
      D[k] := [3];
      changed := true;
      E0 := [ D[j] : j in I0 ];
   end while;
   while [1] in E0 and [7] in E0 do
      j := I0[Index(E0,[1])]; 
      k := I0[Index(E0,[7])]; 
      D[j] := [3]; 
      D[k] := [5];
      changed := true;
      E0 := [ D[j] : j in I0 ];
   end while;
   return D, changed;
end function;

/*
// 0-increment 2-term diagonalized relations:
[
    x7^2 - x3^2,
    x7*x5 - x3*x1,
    x7*x1 - x5*x3,
    x5^2 - x1^2,
];
*/

function Diag0Term4Reduction(D,I0)
   E0 := [ D[j] : j in I0 ];
   J3 := [ j : j in I0 | D[j] eq [3] ];
   changed := false;
   i7 := Index(E0,[7]);
   while #J3 ge 3 and i7 ne 0 do
      for i in [1..3] do  
         D[J3[1]] := [1];
         Remove(~J3,1); 
      end for;
      D[I0[i7]] := [5];
      changed := true;
      i7 := Index(E0,[7]);
   end while;  
   while #J3 ge 4 do
      for i in [1..4] do  
         D[J3[1]] := [1];
         Remove(~J3,1); 
      end for;
      changed := true;
   end while;
   return D, changed; 
end function;

/*
// 0-increment 4-term diagonalized relations:
[
    x7*x3^3 - x5*x1^3,
    x3^4 - x1^4
];
*/

function NonDiag1Reduction(D,I0,I1)
   Diags := [ [1], [3], [5], [7] ];
   NSplt := [2,1,1,2];
   Split := [2,1,1,4];
   E0 := [ D[j] : j in I0 ];
   E1 := [ D[j] : j in I1 ];
   changed := false;
   for R in Diags do
      while Split in E0 and R in E1 do
         j0 := Index(E0,Split);
         j1 := Index(E1,R);
         D[I0[j0]] := NSplt; 
         D[I1[j1]] := [5*R[1] mod 8];
         E0[j0] := NSplt; 
         E1[j1] := [5*R[1] mod 8];
         changed := true;
      end while;
      while R in E0 and Split in E1 do
         j0 := Index(E0,R); 
         j1 := Index(E1,Split); 
         D[I0[j0]] := [5*R[1] mod 8];
         D[I1[j1]] := NSplt;
         E0[j0] := [5*R[1] mod 8];
         E1[j1] := NSplt; 
         changed := true;
      end while;
   end for;
   return D, changed;
end function;

/*
// 1-increment non-diagonalizable relations:
[
    u1*y7 - u0*y3,
    u1*y5 - u0*y1,
    u1*y3 - u0*y7,
    u1*y1 - u0*y5,
    v1*x7 - v0*x3,
    v1*x5 - v0*x1,
    v1*x3 - v0*x7,
    v1*x1 - v0*x5
];
*/

function Diag1Term2Reduction(D,I0,I1)
   E0 := [ D[j] : j in I0 ];
   E1 := [ D[j] : j in I1 ];
   changed := false;
   Unreduced3 := [ [7,7], [7,3], [3,7], [3,3] ];    
   Unreduced7 := [ [7,5], [7,1], [5,7], [5,3] ]; 
   for P in Unreduced3 do
      u0, u1 := Explode(P);
      while [u0] in E0 and [u1] in E1 do
         i0 := Index(E0,[u0]); 
         i1 := Index(E1,[u1]); 
         D[I0[i0]] := [3*u0 mod 8];
         D[I1[i1]] := [3*u1 mod 8];
         E0[i0] := [3*u0 mod 8];
         E1[i1] := [3*u1 mod 8];
         changed := true;
      end while;
   end for;
   for P in Unreduced7 do
      u0, u1 := Explode(P);
      while [u0] in E0 and [u1] in E1 do
         i0 := Index(E0,[u0 mod 8]); 
         i1 := Index(E1,[u1 mod 8]); 
         D[I0[i0]] := [7*u0 mod 8];
         D[I1[i1]] := [7*u1 mod 8];
         E0[i0] := [7*u0 mod 8];
         E1[i1] := [7*u1 mod 8];
         changed := true;
      end while;
   end for;
   return D, changed;
end function;

/*
// 1-increment 2-term relations:
[
    x7*y7 - x5*y5,
    x7*y5 - x1*y3,
    x7*y3 - x5*y1,
    x7*y1 - x1*y7,
    x5*y7 - x3*y1,
    x5*y3 - x3*y5,
    x3*y7 - x1*y5,
    x3*y3 - x1*y1,
];
*/

function Diag1Term3Reduction(D,I0,I1)
   changed := false;
   if #I0 le 1 and #I1 le 1 then
      return D, changed;
   end if;
   E0 := [ D[j] : j in I0 ];
   E1 := [ D[j] : j in I1 ];
   Unreduced3 := [ [7,7], [7,3], [3,7], [3,3] ];    
   Unreduced7 := [ [7,5], [7,1], [5,7], [5,3] ]; 
   i5 := Index(E0,[5]);
   if i5 ne 0 then
      i1 := Index(E0,[1]);
      i3 := Index(E0,[3]);
      j1 := Index(E1,[1]); 
      j5 := Index(E1,[5]); 
      if i3 ne 0 and j5 ne 0 then
         D[I0[i5]] := [1]; 
         D[I0[i3]] := [1]; 
	 D[I1[j5]] := [3]; 
         changed := true;
      elif i3 ne 0 and j1 ne 0 then
         D[I0[i5]] := [1]; 
	 D[I0[i3]] := [1]; 
	 D[I1[j1]] := [7]; 
         changed := true;
      elif i1 ne 0 and j5 ne 0 then
         D[I0[i5]] := [3]; 
	 D[I0[i1]] := [3]; 
	 D[I1[j5]] := [1]; 
         changed := true;
      elif i1 ne 0 and j1 ne 0 then
         D[I0[i5]] := [3]; 
	 D[I0[i1]] := [3]; 
	 D[I1[j1]] := [5]; 
         changed := true;
      elif j5 ne 0 and j1 ne 0 then
         D[I0[i5]] := [1]; 
         D[I1[j1]] := [3]; 
	 D[I1[j5]] := [3]; 
         changed := true;
      elif j1 ne 0 then
         E1[j1] := [3]; 
         j2 := Index(E1,[1]);
         if j2 ne 0 then
            D[I0[i5]] := [1]; 
   	    D[I1[j1]] := [3];
   	    D[I1[j2]] := [7];
            changed := true;
         end if;
      end if;  
   end if;
   return D, changed;
end function;

/*
// 1-increment 3-term relations:
[
    x5*x3*y5 - x1^2*y3,
    x5*x3*y1 - x1^2*y7,
    x5*x1*y5 - x3^2*y1,
    x5*x1*y1 - x3^2*y5,
    x5*y5*y1 - x1*y3^2,
    x5*y1^2 - x1*y7*y3,
];
*/

function Diag1Term4Reduction(D,I0,I1)
   changed := false;
   if (#I0 le 2 and #I1 le 1) or (#I0 le 1 and #I1 le 2) then
      return D, changed;
   end if;
   E0 := [ D[j] : j in I0 ];
   E1 := [ D[j] : j in I1 ];
   I3 := [ j : j in [1..#E0] | E0[j] eq [3] ];
   if #I3 ge 3 then
      i1 := I3[1]; i2 := I3[2]; i3 := I3[3]; 
      j1 := Index(E1,[1]);  
      j5 := Index(E1,[5]);  
      if j5 ne 0 then
         D[I0[i1]] := [1];  
         D[I0[i2]] := [1];  
         D[I0[i3]] := [1];  
         D[I1[j5]] := [7];  
         changed := true;
      elif j1 ne 0 then
         D[I0[i1]] := [1];  
         D[I0[i2]] := [1];  
         D[I0[i3]] := [1];  
         D[I1[j1]] := [3];  
         changed := true;
      end if;
   elif #I3 ge 2 then
      i1 := I3[1]; i2 := I3[2];
      j1 := Index(E1,[1]);  
      j5 := Index(E1,[5]);  
      if j1 ne 0 and j5 ne 0 then
         D[I0[i1]] := [1];  
         D[I0[i2]] := [1];  
         D[I1[j1]] := [3];
         D[I1[j5]] := [7];
         changed := true;
      elif j1 ne 0 then
         E1[j1] := [3];
         j2 := Index(E1,[1]);         
         if j2 ne 0 then
            D[I0[i1]] := [1];
	    D[I0[i2]] := [1];
            D[I1[j1]] := [3];
            D[I1[j2]] := [3];
            changed := true;
         end if; 
      end if;  
   elif #I3 eq 1 then
      i1 := I3[1];
      J1 := [ j : j in [1..#I1] | E1[j] eq [1] ];
      if #J1 ge 3 then
         j1 := J1[1]; j2 := J1[2]; j3 := J1[3];
         D[I0[i1]] := [1];
         D[I1[j1]] := [3];
         D[I1[j2]] := [3];
         D[I1[j3]] := [3];
         changed := true;
      elif #J1 ge 2 then
         j1 := J1[1]; j2 := J1[2];
         j5 := Index(E1,[5]);  
         if j5 ne 0 then
            D[I0[i1]] := [1];
            D[I1[j1]] := [3];
            D[I1[j2]] := [3];
            D[I1[j5]] := [7];
            changed := true;
         end if;
      end if;  
   end if;
   return D, changed;
end function;

/*
// 1-increment 4-term relations:
[
    x3^3*y5 - x1^3*y7,
    x3^3*y1 - x1^3*y3,
    x3^2*y5*y1 - x1^2*y7*y3,
    x3^2*y1^2 - x1^2*y3^2,
    x3*y5*y1^2 - x1*y7*y3^2,
    x3*y1^3 - x1*y3^3,
];
*/

function Diag2Term2Reduction(D,I0,I2)
   changed := false;
   if #I0 eq 0 or #I2 eq 0 then
      return D, changed;
   end if;
   E0 := [ D[i] : i in I0 ];
   E2 := [ D[i] : i in I2 ];
   i5 := Index(E0,[5]);
   i7 := Index(E0,[7]);
   while i7 ne 0 or i5 ne 0 do
      if i7 ne 0 then
         E0[i7] := [3];
         D[I0[i7]] := [3]; 
         i7 := Index(E0,[7]);
      else 
         E0[i5] := [1];
         D[I0[i5]] := [1]; 
         i5 := Index(E0,[7]);
      end if;
      if [7] in E2 then
         D[I2[Index(E2,[7])]] := [3]; 
      elif [5] in E2 then
         D[I2[Index(E2,[5])]] := [1]; 
      elif [3] in E2 then
         D[I2[Index(E2,[3])]] := [7]; 
      elif [1] in E2 then
         D[I2[Index(E2,[1])]] := [5]; 
      end if; 
      E2 := [ D[i] : i in I2 ];
      changed := true; 
   end while;
   return D, changed;
end function;

/*
// 2-increment 2-term relations:
[
    x7*z7 - x3*z3,
    x7*z5 - x3*z1,
    x7*z3 - x3*z7,
    x7*z1 - x3*z5,
    x5*z7 - x1*z3,
    x5*z5 - x1*z1,
    x5*z3 - x1*z7,
    x5*z1 - x1*z5,
];
*/

function Diag2Term3Reduction(D,I0,I1,I2)
   changed := false;
   E0 := [ D[i] : i in I0 ];
   i3 := Index(E0,[3]);
   if i3 eq 0 then
      return D, changed;
   end if;
   E1 := [ D[i] : i in I1 ];
   E2 := [ D[i] : i in I2 ];
   j1 := Index(E1,[1]);   
   j5 := Index(E1,[5]);   
   k1 := Index(E2,[1]);   
   k5 := Index(E2,[5]);   
   while i3 ne 0 and j5 ne 0 and (k1 ne 0 or k5 ne 0) do
      E0[i3] := [1]; 
      E1[j5] := [3];
      D[I0[i3]] := [1];
      D[I1[j5]] := [3];
      if k5 ne 0 then
         E2[k5] := [1];
         D[I2[k5]] := [1];
         k5 := Index(E2,[5]);   
      else 
         E2[k1] := [5];
         D[I2[k1]] := [5];
         k1 := Index(E2,[1]);   
      end if;
      i3 := Index(E0,[3]);      
      j5 := Index(E1,[5]);   
      changed := true;
   end while;
   while i3 ne 0 and j1 ne 0 and #E2 ne 0 do // bugfix
      E0[i3] := [1]; 
      D[I0[i3]] := [1];
      if [7] in E2 then
         k7 := Index(E2,[7]); 
         E1[j1] := [5];
         E2[k7] := [1];
         D[I1[j1]] := [5];
         D[I2[k7]] := [1];
      elif [5] in E2 then
         k5 := Index(E2,[5]); 
         E1[j1] := [1];
         E2[k5] := [7];
         D[I1[j1]] := [1];
         D[I2[k5]] := [7];
      elif [3] in E2 then
         k3 := Index(E2,[3]); 
         E1[j1] := [5];
         E2[k3] := [5];
         D[I1[j1]] := [5];
         D[I2[k3]] := [5];
      elif [1] in E2 then
         k1 := Index(E2,[1]); 
         E1[j1] := [1];
         E2[k1] := [3];
         D[I1[j1]] := [1];
         D[I2[k1]] := [3];
      end if;
      i3 := Index(E0,[3]);      
      j1 := Index(E1,[1]);   
      changed := true;
   end while;
   return D, changed;
end function;

/*
// 2-increment 3-term relations:
[
    x3*y5*z5 - x1*y3*z1,
    x3*y5*z1 - x1*y3*z5,

    x3*y1*z7 - x1*y5*z1,
    x3*y1*z5 - x1*y1*z7,
    x3*y1*z3 - x1*y5*z5,
    x3*y1*z1 - x1*y1*z3
];
*/


