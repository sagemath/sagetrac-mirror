/********************************************
* Short pseudo bases for definite lattices. *
********************************************/

/* 
 The following functions implement the LLL-variant for definite lattices over number rings given in Section 4.3 of 
 A. Schiemann, Classification of Hermitian Forms with the Neighbour method, J. Symbolic Comp. (1998) 26, 487--508.
*/


/* Here is an example:
  > K:= QuadraticField(-3 );
  > L:= HermitianLattice( MatrixRing(K, 7) ! 1);
  > G:= GenusRepresentatives(L);
  > time X:= LLL(G[2] : k:= 4, q1:= 1, q:= 0.99);
  > GramMatrix(G[2], Matrix([ x[2] : x in X ] ));
*/

// Projection:= func< x, B, IP | IsEmpty(B) select x else x - &+[ IP(x, b) / IP(b,b) * b : b in B ] >;
Projection:= function(x, B, IP)
  if IsEmpty(B) then return x; end if;
  v:= Vector( [ IP(x,b) : b in B ] );
  G:= Matrix( [ [IP(b,c): c in B ]: b in B ] );
  l:= Solution(G, v);
  y:= x - &+[ l[i] * B[i] : i in [1..#B] ] ;
  assert forall{b: b in B | IP(y,b) eq 0};
  return y;
end function;

// Find a pseudo basis of Module(PB) having x as first basis vector.
//intrinsic MakePseudoBasis(PB :: [], x::.) -> .
//{}
function MakePseudoBasis(PB, x)
  n:= #PB; 
  K:= BaseRing(x);
  if ISA(Type(K), RngOrd) then 
    K:= FieldOfFractions(K);
  end if;
  x:= Vector(K, x);
  B:= Matrix(K, Matrix([ p[2] : p in PB]));
  ok, s:= IsConsistent(B, x); assert ok;
  d:= Depth(s);

  S:= Matrix( [ B[i]: i in [1..n] | i ne d ] cat [x] );
  T:= MatrixRing(K, n) ! 0;
  for i in [1..n-1] do
    if i lt d then T[i,i]:= 1; elif i ge d then T[i,i+1]:= 1; end if;
  end for;
  T[n]:= s;
  assert T*B eq S;

  H:= HermiteForm( PseudoMatrix( [ x[1]^-1 : x in PB ], Transpose(T) ));

  TT:= Transpose(Matrix(H));
  V:= Parent(PB[1,2]);
  X:= TT^-1 * S; // == TT^-1 * T * B
  C:= CoefficientIdeals(H);
  PA:= [ < C[i]^-1, V ! X[i] > : i in [n..1 by -1] ];
  a:= PA[1,2][d]/x[d] where d:= Depth(x); 
  if not IsOne(a) then PA[1,1] *:= a; PA[1,2] /:= a; end if;
  assert PA[1,2] eq x;
  assert Module(PB) eq Module(PA);
  return PA;
end function;
//end intrinsic;

procedure SRed(~PB, i, IP, aut, q1)
  for j in [i-1..1 by -1] do
    B:= AbsoluteBasis(PB[i,1]^-1 * PB[j,1]);
    C:= [ aut(b) : b in B ];
    yi:= Projection(PB[i,2], [ PB[l,2] : l in [1..j-1]], IP);
    yj:= Projection(PB[j,2], [ PB[l,2] : l in [1..j-1]], IP);
    p:= IP(yj, yi);
    q:= IP(yj, yj);

    // now solve a minimizing problem...
    G:= Matrix(Rationals(), #B, [ AbsoluteTrace(b*c*q) : b in B, c in C ]);
    assert IsSymmetric(G);
    v:= Vector(Rationals(), [ AbsoluteTrace(B[l]*p) : l in [1..#B]]);
    v:= Solution(G, -v);

    CP:= CartesianProduct( [ { Floor(x), Ceiling(x) } : x in Eltseq(v) ] );
    a:= 0;
    m:= q1 * AbsoluteTrace(IP(yi, yi));
    for c in CP do
      aa:= &+[ c[l] * B[l]: l in [1..#B] ];
      y:= yi+aa*yj;
      mm:= AbsoluteTrace(IP(y,y));
      //AbsoluteTrace(IP(yi, yi)), mm;
      if mm lt m then a:= aa; m:= mm; end if; 
    end for;
    if not IsZero(a) then
      PB[i,2] +:= a*PB[j,2];
    end if;
  end for;
end procedure;

procedure IRed(~PB, i, IP, aut, q)
  x:= Projection(PB[i,2], [ PB[j,2]: j in [1..i-1] ], IP);
  B:= AbsoluteBasis(PB[i,1]);
  C:= [ aut(b): b in B ];
  n:= #B;
  l:= IP(x,x);
  G:= Matrix(n, [ AbsoluteTrace(b * c * l) : b in B, c in C ] );
  L:= LatticeWithGram(G);
  m:= Minimum(L);
  if Minimum(PB[i,1]) ne 1 or q * AbsoluteTrace(l) gt m then
    a:= &+[ B[i] * sv[i]: i in [1..n] ] where sv:= ShortestVectors(L : Max:= 1)[1];
    PB[i,1] /:= a;
    PB[i,2] *:= a;
  end if;
end procedure;

intrinsic LLL(L::LatMod: q:= 0.75, q1:= 0.75, k:= 1) -> []
{A relatively short pseudo basis of L}
  require IsPositiveDefinite(L) : "The lattice must be positive definite";
  require 0 lt q  and q  lt 1: "The parameter q must be in the interval (0,1)";
  require 0 lt q1 and q1 le 1: "The parameter q1 must be in the interval (0,1]";
  requirerange k, 1, Rank(L);

  PB:= PseudoBasis(Module(L));
  // This requires to work with singular pseudo-matrices:
  //PB:= PseudoGenerators(Module(L));
  n:= #PB;
  aut:= Involution(L);
  G:= InnerProductMatrix(L);
  K:= BaseRing(G);
  IP:= func< x,y | ( Vector(K, x) * G, Vector(K, [ aut(e) : e in Eltseq(y) ]) ) >;
  wq:= Sqrt(q);

  for i in [1..n] do
    IRed(~PB, i, IP, aut, q);
  end for;

  i:= 1;
  while i lt n do
    Bi:= [ PB[j,2] : j in [1..i-1] ];
    g:= Min(n, i+k-1);
    yi:= Projection(PB[i,2], Bi, IP);
    normyi:= AbsoluteTrace( IP(yi,yi) );

    B:= [];
    C:= [];
    for l in [i..g] do
      yl:= Projection(PB[l,2], Bi, IP);
      B := B cat [ a * PB[l,2]: a in AbsoluteBasis(PB[l,1]) ];
      C := C cat [ a *      yl: a in AbsoluteBasis(PB[l,1]) ];
    end for;
    G:= Matrix(#C, [ AbsoluteTrace( IP(b,c) ) : b,c in C ]);
    LL:= LatticeWithGram(G: CheckPositive:= true);
    mm:= Minimum(LL);
    
    if q * normyi lt mm then
      SRed(~PB, i, IP, aut, q1);
      i +:= 1;
    else
      s:= ShortestVectors(LL : Max:= 1)[1];
      x:= &+[ B[i] * s[i] : i in [1..#B] ];

      C:= MakePseudoBasis( PB[i..g], x );
      for j in [i  ..g] do PB[j]:= C[j-i+1]; end for;
      for j in [i+1..g] do IRed(~PB, j, IP, aut, q ); end for;
      for j in [i  ..g] do SRed(~PB, j, IP, aut, q1); end for;
      i:= i eq 1 select 2 else Max(i-k+1, 1);
    end if;
  end while;
  SRed(~PB, n, IP, aut, q1);

/*
  for i in [1..n] do 
    ok, x:= IsPrincipal(PB[i,1]);
    if ok then
      PB[i]:= < ideal< BaseRing(L) | 1 >, x * PB[i,2] >;
    end if;
  end for;
*/

  assert Module(PB) eq Module(L);
  return [ < x[1], L ! x[2] > : x in PB ];
end intrinsic;
