declare type LatModHerm : LatMod;
//declare type LatModHermElt : LatModElt;
declare type LatModHerm[LatModElt];

declare attributes LatModHerm:
  Form,
  Module,
  Discriminant,
  Norm,
  Scale,
  Definite,
  Diagonal,
  GenusSymbols,
  ZForms,
  Auto;

//declare attributes LatModEltHerm: Parent, Element;

import "lat.m": SetupLatticeInSpace, GramMatrixOfBasis, GroupOrder, MySupport, FindQA, FixModule;


function MyWeakApproximation(P, S)
  O:= Order(P[1]);
  if not IsAbsoluteOrder(O) then
    A:= AbsoluteOrder(O);
    P:= [ A !! p: p in P ];
    x:= WeakApproximation(P, S);
    return NumberField(O) ! x;
  end if;
  return WeakApproximation(P, S);
end function;

get_disc:= func< G | Determinant(G) * (-1)^((m*(m-1)) div 2) where m:= Ncols(G) >;

// Let K be a number field with ring of integers R and let E be an ext. of K.
// Given a prime ideal p of R and an ideal a=x*Z_E for some fractional ideal 
// x of K, the function return the valuation v_p(x).
function MyValuation(a, p)
  if Type(p) eq RngInt then p:= Minimum(p); end if;
  if Type(p) eq RngIntElt then return Valuation(Minimum(a), p); end if;
  d:= Denominator(a);
  O:= Order(a);
  a:= O !! (a*d);
  R:= Order(p);
  return Valuation(a meet R, p) - Valuation(R ! d, p);
end function;

intrinsic IsLocalNorm(E::FldAlg, a::RngElt, p::RngOrdIdl) -> BoolElt
{Tests if a is a norm in the completion of p}
  require Degree(E) eq 2: "The first argument must be a field of degree 2";
  K:= BaseField(E);
  R:= Integers(K);
  if Order(p) cmpeq Integers(E) then
    if Type(R) eq RngInt then return IsLocalNorm(E, a, Minimum(p)); end if;
    p:= p meet R;
  else
    require Order(p) eq R:  "Incompatible arguments";
  end if;
  ok, a:= IsCoercible(K, a);
  require ok and IsPrime(p): "Incompatible arguments";
  x:= PrimitiveElement(E);
  return HilbertSymbol(a, K ! ((2*x - Trace(x))^2), p) eq 1;
end intrinsic;

intrinsic IsLocalNorm(E::FldAlg, a::RngElt, p::RngIntElt) -> BoolElt
{"}//"
  require Degree(E) eq 2 and Type(BaseField(E)) eq FldRat: "The first argument must be a field of degree 2 over the rationals";
  ok, a:= IsCoercible(Rationals(), a);
  require ok and IsPrime(p): "Incompatible arguments";
  x:= PrimitiveElement(E);
  return HilbertSymbol(a, Rationals() ! ((2*x - Trace(x))^2), p) eq 1;
end intrinsic;

function NormicDefect(E, a, p)
  R:= Integers(E);
  if IsZero(a) or IsLocalNorm(E, a, p) then return Infinity(); end if;
  return Valuation(a, p) + Valuation(Discriminant(R), p) - 1;
end function;

function PrimeAbove(L, p)
  if IsZero(p) or not IsPrime(p) then return false, "The ideal must be prime"; end if;
  R:= BaseRing(L);
  if Type(p) eq RngIntElt or Order(p) cmpne R then
    p:= Decomposition(R, p)[1,1];
  end if;
  if Order(p) cmpne R then return false, "Incompatible arguments"; end if;
  return true, p;
end function;

function PrimeBelow(L, p)
  R:= BaseRing(L);
  S:= BaseRing(R);
  if Type(p) eq PlcNumElt then
    inf, index:= IsInfinite(p);
    if inf then
      if NumberField(p) cmpeq NumberField(R) then
        if Type(index) eq RngIntElt then return true, Infinity(); end if;
	return true, InternalInfinitePlaceCreate(NumberField(S), #index eq 2 select index[1] else Remove(index, 1));
      elif NumberField(p) cmpeq NumberField(S) then
        return true, p;
      end if;
      return false, "Incompatible arguments";
    end if;
    p:= Ideal(p);
  elif Type(p) eq Infty then
    if Type(S) eq RngInt then return true, p; end if;
    return false, "Incompatible arguments";
  end if;
  if Order(p) cmpeq R then
     p:= Type(S) eq RngInt select Minimum(p) else p meet S;
  end if;
  if not IsZero(p) and IsPrime(p) and Order(p) cmpeq S then return true, p; end if;
  return false, "Incompatible arguments";
end function;


/* Constructors */

TestHerm:= func< F, a | F eq Transpose(Matrix(BaseRing(F), Ncols(F), [ a(x): x in Eltseq(F) ])) >;

intrinsic HermitianLattice(M::ModDed, F::Mtrx) -> LatModHerm
{The lattice consisting of the module M and the hermitian form F}
  R:= BaseRing(M);
  require Degree(R) eq 2 : "The base ring must be an extension of degree 2";
  require ISA(Type(R), RngOrd) : "The base ring must be an order in a number field";
  require Ncols(F) eq Degree(M) : "The arguments must have the same dimensions";
  FF:= FieldOfFractions(R);
  ok, F:= CanChangeRing(F, FF);
  require ok: "The first two arguments must have the same base ring";
  a:= Automorphisms(FF)[2];
  require TestHerm(F, a): "The form must be hermitian";

  L:= New(LatModHerm);
  L`Module:= FixModule(M);
  L`Form:= F;
//  L`aut:= a;
  return L;
end intrinsic;

intrinsic HermitianLattice(B::Mtrx, F::Mtrx) -> LatModHerm
{The hermitian standard lattice with basis B and form F}
  require Ncols(B) eq Ncols(F) : "The arguments must have the same dimensions";
  R:= BaseRing(F);
  T:= Type(R);
  if ISA(T, FldNum) then
    R:= FieldOfFractions(Integers(R));
    F:= ChangeRing(F, R);
  elif ISA(T, RngOrd) then
    R:= FieldOfFractions(R);
    F:= ChangeRing(F, R);
  else
    require ISA(T, RngOrd) or ISA(T, FldOrd): "Wrong base ring";
  end if;
  ok, B:= CanChangeRing(B, R);
  require ok : "Incompatible base rings";
  return HermitianLattice(Module(Rows(B)), F);
end intrinsic;

function IsHermitian(F)
  K:= BaseRing(F);
  if ISA(Type(K), RngOrd) then K:= NumberField(K); end if;
  if not ISA(Type(K), FldAlg) or Degree(K) ne 2 then
    return false, "The base field of the form must be a quadratic extension";
  end if;
  if not TestHerm(F, Automorphisms(K)[2]) then return false, "The matrix is not hermitian"; end if;
  return true, _;
end function;

intrinsic HermitianLattice(F::Mtrx) -> LatModHerm
{The hermitian standard lattice with form F}
  ok, str:= IsHermitian(F);
  require ok: str;
  M:= MatrixAlgebra(BaseRing(F), Ncols(F));
  return HermitianLattice(M ! 1, M ! F);
end intrinsic;

intrinsic HermitianLattice(L::Lat, h::Map) -> LatModHerm, Map
{Given a homomorphism h: E -> End(L), this return the resulting hermitian structure of L}

  E:= Domain(h);
  require ISA(Type(E), FldAlg) and Degree(E) eq 2 : "The second argument must be a map of a degree 2 extension over an absolute number field";
  K:= BaseField(E);
  require Type(K) eq FldRat or IsAbsoluteField(K) : "The second argument must be a map of a degree 2 extension over an absolute number field";

  n:= AbsoluteDegree(K);
  AS:= AmbientSpace(L);
  B:= ChangeUniverse(Basis(L), AS);
  A:= [ h(a): a in AbsoluteBasis(E) ];
  BE:= [];
  V:= sub< AS | >;
  i:= 1;
  while Dimension(V) ne Rank(L) do
    while B[i] in V do i +:= 1; end while;
    V:= sub< AS | V, [ B[i] * a : a in A ]>;
    BE:= BE cat [ B[i] * a : a in A ] ;
  end while;
  BE:= Matrix(Rationals(), Matrix(BE));
  S:= Matrix(Rationals(), Matrix(B)) * BE^-1;
  S:= Matrix( [ [ E | [K | p[1..n], p[n+1..2*n] ] : p in Partition(Eltseq(x), 2*n) ] : x in Rows(S) ] );
  F:= BE * Matrix(Rationals(), InnerProductMatrix(L)) * Transpose(BE) ;

  G:= [];
  T:= Matrix( Rationals(), 2*n, [ AbsoluteTrace(x*y) : x,y in AbsoluteBasis(E) ] );
  for i, j in [1..#B div (2*n) ] do
    ok, t:= IsConsistent( T, Vector(Rationals(), [ F[(i-1)*2*n+k, (j-1)*2*n+1] : k in [1..2*n]] ));
    assert ok;
    t:= Eltseq(t);
    Append(~G, E ! [ K | t[1..n], t[n+1..2*n] ] );
  end for;
  G:= Matrix(#B div (2*n), G);
  LL:= HermitianLattice(S, G);
  g:= map< L -> LL | x :-> LL ! Vector([ E | [K | p[1..n], p[n+1..2*n] ] : p in Partition(Eltseq(v), 2*n) ]) where v:= Solution(BE, AS ! x) >;

  return LL, g;
end intrinsic;

intrinsic HermitianLattice(L::Lat, x::AlgMatElt, w::AlgMatElt) -> LatModHerm, Map
{The Hermitian lattice over Q(x,w)/Q(x)}

  require x*w eq w*x: "The elements to not commute";
  K:= NumberField( MinimalPolynomial(x) );
  h:= hom< K -> Generic(Parent(x)) | x >;
  F:= [];
  for f in Factorization(ChangeRing(MinimalPolynomial(w), K)) do
    if Degree(f[1]) ne 2 or f[2] ne 1 then continue; end if;
    if w^2 + h(Coefficient(f[1], 1)) * w + h(Coefficient(f[1], 0)) eq 0 then Append(~F, f[1]); end if;
  end for;
  require #F eq 1: "Hermitian structure is not unique";
  E:= ext< K | F[1] >;
  B:= [ x^i: i in [0..Degree(K)-1] ];
  B:= B cat [ w*b: b in B ];
  h:= hom< E -> Generic(Parent(x)) | a :-> &+[ B[i] * e[i] : i in [1..#B] ] where e:= Flat(a) >;

  return HermitianLattice(L, h);
end intrinsic;

intrinsic HermitianLattice(L::Lat, w::Mtrx) -> LatModHerm, Map
{"} // "

  K:= NumberField( MinimalPolynomial(w) );
  return HermitianLattice( L, hom< K -> Generic(Parent(w)) | w > );
end intrinsic;

forward MyLocalBasis;
function IMI(L, p)
  R:= BaseRing(L);
  FF:= FieldOfFractions(R);
  D:= Decomposition(R, p);
  e:= Valuation(Discriminant(R), p);
  s:= e eq 0 select FF ! 1 else PrimitiveElement(D[1,1])^e;

  B:= MyLocalBasis(L, p : Type:= "Submodule");
  M:= Matrix(B);
  G:= M * L`Form * InternalConjugate(L, M);
  k, h:= ResidueClassField(D[1,1]);
  V:= KernelMatrix(Matrix(Ncols(G), [x @ h: x in Eltseq(s*G)]));
  if Nrows(V) eq 0 then return true, _; end if;

  PP:= ProjectiveLineProcess(k, Nrows(V));
  x:= Next(PP);
  while not IsZero(x) do
    ex:= [ e @@ h: e in Eltseq(x * V) ];
    v:= Vector(FF, ex);
//    valv:= Valuation( ((v*G) * Matrix(FF,1,ex))[1], D[1,1] );
    res:= L ! (v*M);
    valv:= Valuation(Norm(res), D[1,1]);
    assert valv ge 0;
    if valv ge 2 then
//      X:= [ Valuation( InnerProduct(res, g), D[1,1] ) : g in Generators(L) ];
//      assert Min(X) ge 1-e;
      return false, v*M;
    end if;
    x:= Next(PP);
  end while;
  return true, _;
end function;


/* Constructing maximal hermitian lattices */

// Check if L is max. integral at p. If not, return either:
// - a minimal integral overlattice at p (minimal flag set)
// - a maximal integral overlattice at p (minimal flag not set).
function MIL(L, p, Minimal)
  R:= BaseRing(L);
  // Already maximal?
  if Valuation(Norm(Volume(L)), p) eq 0 and not IsRamified(p, R) then return true, L; end if;

  B, G, S:= JordanDecomposition(L, p);
  D:= Decomposition(R, p);
  P:= D[1,1];
  IsMax:= true;

  if #D eq 2 then	// The split case is trivial.
    assert S[#S] ne 0;
    if Minimal then
      max:= 1;
      M:= P^-1 * B[#S,1];
    else
      max:= S[#S];
      M:= &+ [ P^-S[i] * b : b in Rows(B[i]), i in [1..#B] | S[i] ne 0 ];
    end if;
    M:= (ML + M) meet P^-max * ML where ML:= L`Module;
    return false, SetupLatticeInSpace(L, M);
  elif D[1,2] eq 1 then
    // The inert case, we first reduce to L squarefree.
    if S[#S] ge 2 then
      if Minimal then
        max:= 1;
        M:= P^-(S[#S] div 2) * B[#S, 1];
      else
        max:= S[#S];
        M:= &+ [ P^-(S[i] div 2) * b: b in Rows(B[i]), i in [1..#B] | S[i] ge 2 ];
      end if;
      M:= (ML + M) meet P^-max*ML where ML:= L`Module;
      L:= SetupLatticeInSpace(L, M);
      if Minimal then return false, L; end if;
      B, G, S:= JordanDecomposition(L, p);
      IsMax:= false;
    end if;
    // Now we seek for zeros of a*x^2+b*y^2;
    k,h:= ResidueClassField(D[1,1]);
    while &+[ S[i] * Nrows(B[i]) : i in [1..#B] ] gt 1 do
      ok:= exists(i){i: i in [1..#S] | S[i] eq 1};
      assert ok and Nrows(B[i]) ge 2;
      repeat
        r:= Random(k) @@ h;
      until Valuation(G[i,1,1] + G[i,2,2]*Norm(r), D[1,1]) ge 2;
      M:= P^-1 * (B[i,1] + r*B[i,2]);
      M:= (ML + M) meet P^-1*ML where ML:= L`Module;
      L:= SetupLatticeInSpace(L, M);
      if Minimal then return false, L; end if;
      IsMax:= false;
      B, G, S:= JordanDecomposition(L, p);
      assert S[1] ge 0;
    end while;
    return IsMax, L;
  else  // The ramified case. Again, we eliminate all blocks of scale >=2.
    if S[#S] ge 2 then
      if Minimal then
        max:= 1;
        M:= P^-(S[#S] div 2) * B[#S, 1];
      else
        max:= S[#S];
        M:= &+ [ P^-(S[i] div 2) * b : b in Rows(B[i]), i in [1..#B] | S[i] ge 2 ];
      end if;
      M:= (ML + M) meet P^-max*ML where ML:= L`Module;
      L:= SetupLatticeInSpace(L, M);
      if Minimal then return false, L; end if;
      B, G, S:= JordanDecomposition(L, p);
    end if;
    v:= Valuation(Volume(L), D[1,1]);
    ok, x:= IMI(L, p);
    while not ok do
      LL:= L;
      L:= SetupLatticeInSpace(L, Module(L) + P^-1 * x);
      v -:= 2;
      assert v eq Valuation(Volume(L), P);
      assert Valuation(Norm(L), P) ge 0;
      if Minimal then return false, L; end if;
      IsMax:= false;
      ok, x:= IMI(L, p);
    end while;
    // We know the val. of a max. lattice, so check
    // TODO:: Use valmax as stopping condition.
    assert IsEven(v); v div:= 2;
    m:= Rank(L);
    e:= Valuation(Discriminant(R), p);
    if IsOdd(m) then 
      valmax:= -(m-1)/2*e;
    else
      valmax:= -m/2 * e;
      disc:= get_disc(L`Form);
      if not IsLocalNorm(NumberField(R), disc, p) then
        valmax +:= 1;
      end if;
    end if;
//    if v ne valmax then "ERROR: Wrong volume at", p; end if; 
    assert v eq valmax;
    return IsMax, L;
  end if;
end function;

intrinsic IsMaximalIntegral(L::LatModHerm, p::RngOrdIdl) -> BoolElt, LatModHerm
{Checks whether L is p-maximal integral. If not, a minimal integral over-lattice at p is returned}
  ok, p:= PrimeBelow(L, p);
  require ok : p;
  require Valuation(Norm(Norm(L)), p) ge 0: "The lattice is not integral";
  return MIL(L, p, true);
end intrinsic;

intrinsic IsMaximalIntegral(L::LatModHerm, p::RngIntElt) -> BoolElt, LatModHerm
{"} //"
  ok, p:= PrimeBelow(L, p);
  require ok : p;
  require Valuation(Norm(Norm(L)), p) ge 0: "The lattice is not integral";
  return MIL(L, p, true);
end intrinsic;

intrinsic IsMaximalIntegral(L::LatModHerm) -> BoolElt, LatModHerm
{Checks whether L is maximal integral. If not, a minimal integral over-lattice is returned}
  require IsIntegral(Norm(L)) : "The lattice is not integral";
  S:= BaseRing(L);
  bad:= MySupport(Norm(Volume(L))) join MySupport(Discriminant(S));
  for p in bad do 
    ok, LL:= MIL(L, p, true);
    if not ok then return false, LL; end if;
  end for;
  return true, L;
end intrinsic;

intrinsic IsMaximal(L::LatModHerm, p::RngOrdIdl) -> BoolElt, LatModHerm
{Checks if L_p is Norm(L_p)-maximal}
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  require not IsZero(L): "The lattice must be non-zero";
  R:= BaseRing(BaseRing(L));
  v:= MyValuation(Norm(L), p);
  x:= (Type(p) eq RngIntElt select p else PrimitiveElement(p))^-v;
  ok, LL:= IsMaximalIntegral( Rescale(L, x), p);
  if ok then return true, L; end if;
  return false, SetupLatticeInSpace(L, Module(LL));
end intrinsic;

intrinsic IsMaximal(L::LatModHerm, p::RngIntElt) -> BoolElt, LatModHerm
{"} //"
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  return IsMaximal(L, BaseRing(L) * p);
end intrinsic;

intrinsic MaximalIntegralLattice(L::LatModHerm) -> LatModHerm
{A maximal integral lattice containing L}
  require IsIntegral(Norm(L)) : "The lattice is not integral";
  S:= BaseRing(L);
  bad:= MySupport(Norm(Volume(L))) join MySupport(Discriminant(S));
  for p in bad do
    ok, L:= MIL(L, p, false);
  end for;
  return L;
end intrinsic;

intrinsic MaximalIntegralLattice(L::LatModHerm, p::RngOrdIdl) -> LatMod
{A p-maximal integral lattice over L}
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  require MyValuation(Norm(L), p) ge 0 : "Lattices are not locally integral";
  ok, L:= MIL(L, p, false);
  return L;
end intrinsic;

intrinsic MaximalIntegralHermitianLattice(F::Mtrx) -> LatModHerm
{A maximal integral lattice with respect to the hermitian form F}
  ok, str:= IsHermitian(F);
  require ok: str;
  L:= HermitianLattice(F);
  FF:= Factorization(Scale(L));
  s:= 1;
  while #FF ne 0 do
    ok:= exists(i){i:i in [2..#FF] | Norm(FF[i,1]) eq nrm } where nrm:= Norm(FF[1,1]);
    if ok then
      assert FF[i,2] eq FF[1,2];
      s *:= FF[1,1]^FF[1,2];
      Remove(~FF, i);
    else
      s *:= FF[1,1]^-(FF[1,2] div 2);
    end if;
    Remove(~FF, 1);
  end while;
  if not IsOne(s) then L:= s*L; end if;
  return MaximalIntegralLattice(L);
end intrinsic;

intrinsic IsSquarefree(L::LatModHerm, p::RngIntElt) -> BoolElt
{Returns true iff L_p is squarefree}
  ok, p:= PrimeAbove(L, p);
  require ok: p;
  return IsSquarefree(L, p);
end intrinsic;

intrinsic IsSquarefree(L::LatModHerm, p::RngOrdIdl) -> BoolElt
{Returns true iff L_p is squarefree}
  ok, p:= PrimeAbove(L, p);
  require ok: p;
  if IsZero(L) then return true; end if;
  I:= IsSplit(p) select {0} else {0,1};
  return { g[2]: g in GenusSymbol(L, p) } subset I;
end intrinsic;


/* Basic functions */

function MyGramSchmidt(F, a)
  n:= Ncols(F);
  S:= MatrixRing(BaseRing(F), n) ! 1;
  ok, D:= IsDiagonal(F);
  if not ok then
    for i in [1..n] do
      if F[i,i] eq 0 then
        T:= Parent(S) ! 1;
        ok:= exists(j){j: j in [i+1..n] | F[j,j] ne 0 };
        if ok then
   	  T[i,i]:= 0;
  	  T[j,j]:= 0;
	  T[i,j]:= 1;
	  T[j,i]:= 1;
	else
          ok:= exists(j){j: j in [i+1..n] | F[i,j] ne 0 };
          error if not ok, "matrix is not of full rank";
	  T[i,j] := 1/(2*F[j,i]);
	end if;
        S:= T*S;
        F:= T * F * Transpose( Matrix(n, [ a(x): x in Eltseq(T) ]) );
      end if;
      T:= Parent(S) ! 1;
      for j in [i+1..n] do
        T[j,i]:= -F[j,i]/F[i,i];
       end for;
       F:= T * F * Transpose( Matrix(n, [ a(x): x in Eltseq(T) ]) );
       S:= T * S;
     end for;
     ok, D:= IsDiagonal(F);
     assert ok;
  end if;
  return D, S;
end function;

function IRE(L1, L2, p, Ambient)
  if Ambient then
    F1:= InnerProductMatrix(L1);
    F2:= InnerProductMatrix(L2);
  else
    F1:= GramMatrixOfBasis(L1);
    F2:= GramMatrixOfBasis(L2);
  end if;
  if F1 cmpeq F2 then 
    return true;
  elif Ncols(F1) ne Ncols(F2) then
    return false;
  elif ISA( Type(p), {RngIntElt, RngOrdIdl}) then
    return IsLocalNorm( NumberField(L1),  Determinant(F1) * Determinant(F2), p );
  elif Type(p) eq FldNumElt and IsComplex(p) then
    return true;
  end if;
  a:= Involution(L1);
  F:= FixedField(L1);
  I1:= #[ d: d in MyGramSchmidt(F1, a) | Evaluate(F ! d, p) lt 0 ];
  I2:= #[ d: d in MyGramSchmidt(F1, a) | Evaluate(F ! d, p) lt 0 ];
  return I1 eq I2; 
end function;

intrinsic IsRationallyEquivalent(L1::LatModHerm, L2::LatModHerm, p::RngOrdIdl : AmbientSpace:= false) -> BoolElt
{Tests if L1 and L2 are equivalent over the completion at p}
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible lattices";
  ok, p:= PrimeBelow(L1, p);
  require ok : p;
  return IRE(L1, L2, p, AmbientSpace);
end intrinsic;

intrinsic IsRationallyEquivalent(L1::LatModHerm, L2::LatModHerm, p::RngIntElt : AmbientSpace:= false) -> BoolElt
{"} //"
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible lattices";
  ok, p:= PrimeBelow(L1, p);
  require ok : p;
  return IRE(L1, L2, p, AmbientSpace);
end intrinsic;

intrinsic IsRationallyEquivalent(L1::LatModHerm, L2::LatModHerm, p::PlcNumElt : AmbientSpace:= false) -> BoolElt
{"} //"
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible lattices";
  ok, p:= PrimeBelow(L1, p);
  require ok : p;
  return IRE(L1, L2, p, AmbientSpace);
end intrinsic;

intrinsic IsRationallyEquivalent(L1::LatModHerm, L2::LatModHerm, p::Infty : AmbientSpace:= false) -> BoolElt
{"} //"
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible lattices";
  ok, p:= PrimeBelow(L1, p);
  require ok : p;
  return IRE(L1, L2, p, AmbientSpace);
end intrinsic;


intrinsic IsRationallyEquivalent(L1::LatModHerm, L2::LatModHerm : AmbientSpace:= false) -> BoolElt
{Tests if L1 and L2 are equivalent over their base field}
  if AmbientSpace then
    F1:= InnerProductMatrix(L1);
    F2:= InnerProductMatrix(L2);
  else
    F1:= GramMatrixOfBasis(L1);
    F2:= GramMatrixOfBasis(L2);
  end if;
  if F1 cmpeq F2 then return true; end if;
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible lattices";

  K:= FixedField(L1);
  E:= NumberField(BaseRing(L1));
  return Ncols(F1) eq Ncols(F2) and NormEquation(E, K ! (Determinant(F1) * Determinant(F2)) );
end intrinsic;

rat_gens:= function(S)
  if Type(Universe(S)) eq RngInt then
    return GCD(S);
  end if;
  d:= LCM({ Denominator(x) : x in S });
  return GCD({Integers() | x*d: x in S })/d;
end function;

ideal_gens:= function(I)
  O:= Order(I);
  if IsAbsoluteOrder(O) then return Basis(I); end if;
  F:= FieldOfFractions(O);
  G:= Generators(Module(I));
  G:= [ F | Eltseq(g) : g in G ];
  assert I eq ideal< O | G >;
  return G;
end function;

ideal_trace:= function(I)
  R:= Order(I);
//  Gens:= [ Trace(g) : g in ideal_gens(I) ];	// Klappt nicht!!!!
  Gens:= { NumberField(R) | Trace(g) : g in AbsoluteBasis(I) };
  return Type(Universe(Gens)) in {FldRat, RngInt} select rat_gens(Gens) else ideal< R | Gens >;
end function;

// Get a local unit at P such that T(u)=1.
// P is assumed to be even and inert.
/*
SpecialUnit:= function(P)
  R:= Order(P);
  E:= NumberField(R);
  assert Degree(E) eq 2;
  x:= PrimitiveElement(E);
  x:= x - Trace(x/2);
  v:= Valuation(x, P);
  if v ne 0 then
    x *:= PrimitiveElement(P)^-v;
  end if;
  e:= Valuation(R ! 2, P);
  k,h:= quo< R | P^e >;
  s:= ((x @ h)^-1) @@ h;
  u:= (1-s*x)/2;
  assert Valuation(u, P) eq 0;
  assert Trace(u) eq 1;
  return u;
end function;*/

function SpecialUnit(P)
  assert not IsRamified(P);
  R:= Order(P);
  E:= NumberField(R);
  assert Degree(E) eq 2;
  x:= PrimitiveElement(E);
  x -:= Trace(x)/2;
  K:= BaseField(E);
  a:= K ! (x^2);
  if Type(K) eq FldRat then
    p:= Minimum(P);
    pi:= p;
  else
    p:= P meet Integers(K);
    pi:= PrimitiveElement(p);
  end if;
  v:= Valuation(a, p);
  if v ne 0 then
    assert IsEven(v);
    a /:= pi^v;
    x /:= pi^(v div 2);
  end if;

  // x^2=a has quadratic defect (4).
  // Write a=t^2*(1+4d) where d is a unit.
  // Then u:=(1+x/t)/2 does the trick.
  k, h:= ResidueClassField(p);
  s:= SquareRoot(h(a));
  t:= s @@ h;
  a /:= t^2; x /:= t;
  w:= Valuation(a-1, p);
  e:= Valuation(Order(p) ! 2, p);
  while w lt 2*e do
    assert IsEven(w);
    s:= SquareRoot(h((a - 1) / pi^w));
    t:= 1+(s @@ h)*pi^(w div 2);
    a /:= t^2;
    x /:= t;
    w:= Valuation(a-1, p);
  end while;
  assert w eq 2*e;

  u:= (1+x)/2;
  assert Trace(u) eq 1;
  assert Valuation(u, P) eq 0;
  return u;
end function;

/*
intrinsic Scale(L::LatModHerm) -> RngOrdFracIdl
{The scale of the lattice L}
  if not assigned L`Scale then
    P:= PseudoBasis(Module(L));
    B:= Matrix( [ p[2] : p in P ] );
    F:= B * L`Form * InternalConjugate(L, B);
    a:= Involution(L);
    L`Scale:= &+{ F[i,j] * P[i,1] * a(P[j,1]) : i, j in [1..#P] | F[i,j] ne 0 };
    assert L`Scale eq a(L`Scale);
  end if;
  return L`Scale;
end intrinsic;
*/

intrinsic Norm(L::LatModHerm) -> RngOrdFracIdl
{The norm of the lattice L}
  if not assigned L`Norm then
    F, C:= GramMatrixOfBasis(L);
    a:= Involution(L);
    Gens:= { F[i,i] * C[i] * a(C[i]) : i in [1..#C] | F[i,i] ne 0 } 
      join { ideal_trace(C[i] * F[i,j] * a(C[j])) : i in [1..j-1], j in [1..#C] | F[i,j] ne 0 };
    L`Norm:= &+ Gens;
  //  L`Norm:= IsAbsoluteOrder(BaseRing(L)) select rat_gens(Gens) else &+Gens;
  end if;
  return L`Norm;
end intrinsic;

intrinsic NormGenerator(L::LatModHerm, p::RngOrdIdl) -> ModTupFldElt
{A generator of the norm of L_p}
  ok, p:= PrimeAbove(L, p);
  require ok: p;
  B:= LocalBasis(L, p);
  G:= GramMatrix(L, B);
  v:= Valuation(Norm(L), p);
  if exists(i){i: i in [1..#B] | Valuation( G[i,i], p) eq v} then return B[i]; end if;
  min:= < Infinity(), 0, 0 >;
  for i,j in [1..#B] do
    if i eq j then continue; end if;
    val:= Valuation(G[i,j], p);
    if val lt min[1] then min:= <val, i, j>; end if;
  end for;
  _, i, j:= Explode(min);
  I:= ideal< BaseRing(L) | G[i,j] >;
  ok:= exists(b){b: b in AbsoluteBasis(I) | Valuation(Trace(b*G[i,j]), p) eq v };
  assert ok;
  return b*B[i] + B[j];
end intrinsic;

intrinsic NormGenerator(L::LatModHerm, p::RngIntElt) -> ModTupFldElt
{"} //"
  ok, p:= PrimeAbove(L, p);
  require ok: p;
  return NormGenerator(L, p);
end intrinsic;

intrinsic IsModular(L::LatModHerm, p::RngOrdIdl) -> BoolElt, RngIntElt
{If L_p is p^v-modular, return true and v}
  ok, PP:= PrimeAbove(L, p);
  require ok : PP;
  v:= Valuation(Scale(L), PP);
//  G:= GenusSymbol(L, p);
  if v*Rank(L) eq Valuation(Volume(L), PP) then
//    assert #G eq 1 and G[1,2] eq v;
    return true, v;
  end if;
 // assert #G ne 1;
  return false, _;
end intrinsic;

intrinsic IsModular(L::LatModHerm, p::RngIntElt) -> BoolElt, RngIntElt
{"} //"
  ok, P:= PrimeAbove(L, p);
  require ok : P;
  return IsModular(L, P); 
end intrinsic;

intrinsic IsUnimodular(L::LatModHerm, p::RngIntElt) -> BoolElt
{Returns true iff L_p is unimodular}
  ok, P:= PrimeAbove(L, p);
  require ok : P;
  return IsUnimodular(L, P);
end intrinsic;

intrinsic BadPrimes(L::LatModHerm : Discriminant:= false) -> []
  {The set of places at with L is not unimodular}
  Bad:= MySupport(Norm(Scale(L))) join MySupport(Norm(Volume(L)));
  if Discriminant then 
    ok, Disc:= IsIntrinsic("Discriminant");	// Sigh
    assert ok;
    Bad join:= MySupport( Disc(BaseRing(L)) );
  end if;
  return Bad;
end intrinsic;


// A free submodule of Module(L), coinciding at p with L.
function MyLocalBasis(L, p: Type:= "")
  R:= BaseRing(L);
  if Order(p) cmpeq R then
    S:= IsSplit(p) select [ p, Involution(L)(p) ] else [p];
  else
    S:= [d[1]: d in Decomposition(R, p)];
  end if;
  B:= [ EmbeddingSpace(Module(L)) | ];
  for pp in PseudoBasis(Module(L)) do
    g:= Generators(pp[1]);
    if #g eq 1 then x:= g[1];
    elif Type eq "" then x:= MyWeakApproximation(S, [ Valuation(pp[1], P) : P in S ]);
    else
      Fact:= Factorization(pp[1]);
      Fact:= Type eq "Submodule" select [ f: f in Fact | f[1] in S or f[2] gt 0 ]
                                   else [ f: f in Fact | f[1] in S or f[2] lt 0 ];
      for p in S do
        if forall{ f: f in Fact | f[1] ne p } then Append(~Fact, <p, 0>); end if;
      end for;
      x:= MyWeakApproximation([ f[1] : f in Fact ], [ f[2] : f in Fact ]);
    end if;
    Append(~B, pp[2]*x);
  end for;

  return B;
end function;


function JordanSplitting(L, p)
  R:= BaseRing(L);
  aut:= Involution(L);
  even:= Minimum(p) eq 2;

  P:= MyLocalBasis(L, p);
  S:= Matrix(P);

  D:= Decomposition(R, p);
  split:= #D eq 2;
  ram:= D[1,2] eq 2;
  D:= [ d[1]: d in D ];
  n:= #P;

  P:= D[1];
  if split then
    pi:= MyWeakApproximation(D, [1,0]);
  elif ram then
    pi:= PrimitiveElement(P);
  else
    pi:= Type(p) eq RngIntElt select p else R ! PrimitiveElement(p);
    su:= even select SpecialUnit(P) else 1;
  end if;

  oldval:= Infinity();
  Blocks:= []; Exponents:= [];

  F:= L`Form;
  k:= 1;
  while k le n do
    G:= S * F * InternalConjugate(L, S);
    X:= [ Valuation(G[i,i], P): i in [k..n] ];

    m, ii:= Minimum( X ); ii +:= k-1;
    pair:= <ii, ii>;
    for i,j in [k..n] do
      tmp:= Valuation(G[i,j], P);
      if tmp lt m then m:= tmp; pair:= <i,j>; end if;
    end for;
    if m ne oldval then Append(~Blocks, k); oldval:= m; Append(~Exponents, m); end if;
    i,j:= Explode(pair);

    if (i ne j) and not (ram and (even or IsOdd(m))) then
      a:= G[i,j];
      if split then
        lambda:= Valuation( pi*a, P ) eq m select pi else aut(pi);
      elif ram then
        assert IsEven(m);
        lambda:= Norm(pi)^(m div 2) / a; 
      else
        lambda:= pi^m / a * su;
      end if;
      S[i] +:= aut(lambda) * S[j];
      G:= S * F * InternalConjugate(L, S);
      assert Valuation(G[i,i], P) eq m;
      j:= i;
    end if;

    if i ne j then
      assert i lt j;
      SwapRows(~S, i, k);
      SwapRows(~S, j, k+1);
      SF:= S*F;
      X1:= Eltseq(SF * InternalConjugate(L, S[k  ]));
      X2:= Eltseq(SF * InternalConjugate(L, S[k+1]));
      for l in [k+2..n] do
        den:= Norm(X2[k])-X1[k]*X2[k+1];
        S[l] -:= (X2[l]*X1[k+1]-X1[l]*X2[k+1])/den*S[k] + (X1[l]*X2[k] -X2[l]*X1[k])/den*S[k+1];
      end for;
      k+:= 2;
    else
      SwapRows(~S, i, k);
      X:= Eltseq(S * F * InternalConjugate(L, S[k]));
      for l in [k+1..n] do S[l] -:= X[l]/X[k] * S[k]; end for;
      k+:= 1;
    end if;
  end while;

  if not ram then
    G:= S * F * InternalConjugate(L, S);
    assert IsDiagonal(G);
  end if;

  Append(~Blocks, n+1);
  Matrices:= [* RowSubmatrix(S, Blocks[i], Blocks[i+1] - Blocks[i]) : i in [1..#Blocks-1] *];
  return Matrices, [* m*F*InternalConjugate(L, m): m in Matrices *], Exponents;
end function;

intrinsic JordanDecomposition(L::LatModHerm, p::RngOrdIdl) -> []
{A Jordan decomposition of the completion of L at p}
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  return JordanSplitting(L,p);
end intrinsic;

intrinsic JordanDecomposition(L::LatModHerm, p::RngIntElt) -> []
{"} //"
  R:= BaseRing(L);
  require Type(BaseRing(R)) eq RngInt and IsPrime(p): "Wrong arguments"; 
  return JordanSplitting(L, p);
end intrinsic;

function TraceIdeal(L, p)
  R:= Order(p);
  v:= Valuation(Different(R), p);
  V:= { Valuation(l, p) : l in L | l ne 0 };
  X:= [ 2*((l+v) div 2) : l in V ];
  return #X eq 0 select Infinity() else Min(X);
end function;

function GetNorm(G, P)
  n:= Ncols(G);
  T:= TraceIdeal([ G[i,j]: j in [i+1..n], i in [1..n] ], P);
  Diag:= Min([ Valuation(G[i,i], P) : i in [1..n] ]);
  return Type(T) eq Infty select Diag else Min( Diag, T );
end function;

function GS(L, p)
  if not assigned L`GenusSymbols then
    L`GenusSymbols:= AssociativeArray();
  end if;
  ok, sym:= IsDefined(L`GenusSymbols, p);
  if ok then return sym; end if; 

  B, G, S:= JordanSplitting(L, p);
  R:= BaseRing(L);
  E:= NumberField(R);
  K:= BaseField(E);
  if Minimum(p) ne 2 or not IsRamified(p, R) then
    sym:= [ < Nrows(B[i]), S[i], IsLocalNorm(E, K ! Determinant(G[i]), p) > : i in [1..#B] ];
  else
    P:= Decomposition(R, p)[1,1];
    pp:= Type(p) eq RngIntElt select p else PrimitiveElement(p);
    sym:= [];
    for i in [1..#B] do
      normal:= GetNorm(G[i], P) eq S[i];
      GG:= DiagonalJoin(< pp^Max(0, S[i]-S[j]) * G[j] : j in [1..#B] >);
      v:= GetNorm(GG, P);
      assert v eq Valuation(Norm(HermitianLattice(GG)), P);		// TODO::: Remove
      s:= < Nrows(B[i]), S[i], normal, v, K ! Determinant( DiagonalJoin(< G[j] : j in [1..i]>) ) >;
      Append(~sym, s);
    end for;
  end if;
  L`GenusSymbols[p]:= sym;
  return sym;
end function;

intrinsic GenusSymbol(L::LatModHerm, p::RngOrdIdl) -> []
{The genus symbol of L at p}
  ok, pp:= PrimeBelow(L, p);
  require ok : pp;
  return GS(L, pp);
end intrinsic;

intrinsic GenusSymbol(L::LatModHerm, p::RngIntElt) -> []
{"} //"
  ok, pp:= PrimeBelow(L, p);
  require ok : pp;
  return GS(L, pp);
end intrinsic;

function IsLocIso(L1, L2, p)
  R:= BaseRing(L1);
  E:= NumberField(R);
  S1:= GenusSymbol(L1, p);
  S2:= GenusSymbol(L2, p);
  if Minimum(p) ne 2 or not IsRamified(p, R) then 
    return S1 eq S2;
  end if;
  t:= #S1;
  // number of blocks:
  if t ne #S2 then return false; end if;
  // test ranks, scales, normal property and fundamental norm ideals
  if exists{<i, k>: i in [1..t], k in [1..4] | S1[i,k] ne S2[i,k] } then return false; end if;
  // Test if spaces are equivalent, i.e. Det matches
  if not IsLocalNorm(NumberField(R), S1[t,5] / S2[t,5], p) then return false; end if;
  for i in [1..t-1] do
    assert Valuation(S1[i,5], p) eq Valuation(S2[i,5], p);
    x:= (S1[i,5] / S2[i,5]);
    n:= 2*NormicDefect(E, x, p);
    if n lt (S1[i, 4] + S1[i+1, 4]) - 2*S1[i, 2] then return false; end if;
  end for;
  return true;
end function;

intrinsic IsLocallyIsometric(L1::LatModHerm, L2::LatModHerm, p::RngOrdIdl) -> BoolElt
{Returns true if and only if the p-adic completions of L1 and L2 are isometric}
  require BaseRing(L1) eq BaseRing(L2): "Incompatible lattices";
  ok, p:= PrimeBelow(L1, p);
  require ok : p;
  return IsLocIso(L1,L2,p);
end intrinsic;

intrinsic IsLocallyIsometric(L1::LatModHerm, L2::LatModHerm, p::RngIntElt) -> BoolElt
{"}//"
  require BaseRing(L1) eq BaseRing(L2): "Incompatible lattices";
  require IsPrime(p: Proof:= false) and Type(FixedField(L1)) eq FldRat : "Incompatible arguments";
  return IsLocIso(L1,L2,p);
end intrinsic;

intrinsic PartialDual(L::LatModHerm, p::RngOrdIdl) -> LatModHerm
{The dual lattice of L}
  require IsFull(L) : "The lattice must have full rank";
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  R:= BaseRing(L);
  D:= Dual(L);
  if IsRamified(p, BaseRing(L)) then 
    p:= Decomposition(R, p)[1,1]; 
    m:=  Valuation(Scale(L), p);
    M:= -Valuation(Scale(D), p);
  else
    m:=  MyValuation(Scale(L), p);
    M:= -MyValuation(Scale(D), p);
  end if;
  pp:= ideal< R | p >;
  X:= (D + pp^-m * L) meet pp^-M * L;
  d:= Volume(L);
  assert Volume(X) eq d * p^(-MyValuation(d, p)+MyValuation(Volume(D), p));
  return X;
end intrinsic;

intrinsic PartialDual(L::LatModHerm, p::RngIntElt) -> LatModHerm
{The dual lattice of L}
  require IsFull(L) : "The lattice must have full rank";
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  return PartialDual(L, BaseRing(L) * p);
end intrinsic;

// L=lattice, P=prime, C=Inv(P), B=local basis at P \cap C, G=Gem(B) mod C,
// split=(P ne C), x=global vector, h=epi mod C.
function Neighbour(L, B, xG, x, h, P, C, split)
  R:= Order(P);
  n:= #B;
  if C cmpeq 0 then 
    C:= split select Involution(L)(P) else P;
  end if;
  I:= [ i: i in [1..n] | xG[i] ne 0 ];
  assert #I ne 0;
  i:= I[1];
  a:= Involution(L);
  Gens:= [ < 1*R, B[j] > : j in {1..n} diff Set(I) ] cat [ < 1*R, B[j] - (xG[j]/xG[i]) @@ h @ a * B[i] > : j in I | j ne i ] cat [ < P, B[i] >, < C^-1, x > ];
  M:= Module(Gens) + (split select P*C else P)*L`Module;
  LL:= SetupLatticeInSpace(L, M);

  //test
  assert IsLocallyIsometric(L, LL, P);

  return LL;
end function;

function _Neighbours(L, P, Result, AutoOrbits, Max : CallBack:= false)
//  "Entering _Neighbours, #Result=", #Result, "Max=", Max;
  ok, scale:= IsModular(L, P);
  error if not ok, "non-modular lattice!";
  R:= BaseRing(L);
  F:= FieldOfFractions(R);
  a:= Involution(L);
  C:= a(P);
  B:= MyLocalBasis(L, P: Type:= "Submodule");
  T:= Matrix(B);
  k, h:= ResidueClassField(C);
  Form:= L`Form;
  special:= false;
  if scale ne 0 then
//    "rescaling lattice for neighbour method!";
    if IsRamified(P) then
      special:= IsOdd(scale);
      scale := (scale + 1) div 2;
    end if;
    Form:= PrimitiveElement(P meet BaseRing(R))^-scale * Form;
  end if;
  n:= Rank(L);
  W:= VectorSpace(k, n);

  if AutoOrbits cmpne false then
    if AutoOrbits cmpeq true then
      A:= AutomorphismGroup(L : NaturalAction:= false);
      BM:= BasisMatrix(L);
      S1:= Solution(BM, T);
      S2:= Solution(T, BM);
    else
      error if not IsFull(L), "Setting the automorphism group is not supported for lattices which are not full";
      A:= ChangeRing(AutoOrbits, F);
      S1:= T;
      S2:= T^-1;
    end if;
    A:= [Matrix(k, n, [ x @ h : x in Eltseq(S1 * A.i * S2) ]) : i in [1..Ngens(A)]];
    A:= [ a: a in A | not IsScalar(a) ];
    AutoOrbits:= #A ge 1;
  end if;
/*  if AutoOrbits cmpne false then
    A:= AutoOrbits cmpeq true select AutomorphismGroup(L) else ChangeRing(AutoOrbits, F);
    A:= [ Matrix(n, [ x @ h : x in Eltseq(T * A.i * TI) ]) : i in [1..Ngens(A)] ] where TI:= T^-1;
    // If A mod p is a scalar group, we should ignore it
    AutoOrbits:= exists{g: g in A | not IsScalar(g)};
  end if;*/
  if AutoOrbits then
    A:= sub< GL(n, k) | A >;
//    "computing orbits";
    LO:= IsPrimeField(k) select [ l[2].1: l in OrbitsOfSpaces(A, 1) ] else [ l[1].1: l in LineOrbits(A) ];
//    "done";
  else
    LO:= ProjectiveLineProcess(W);
  end if;

  //"Number of Orbits:", #LO;
  UseCallBack:= CallBack cmpne false;
  keep:= true; cont:= true;
  if P ne C then
    G:= Matrix(k, n, [ h(e) : e in Eltseq(T * Form * InternalConjugate(L, T)) ]);
    pi:= MyWeakApproximation([P, C], [1,0]);
    pih:= h(pi);
    for i in [1..#LO] do
      w:= AutoOrbits select W ! LO[i] else Next(LO);
      x:= &+[ B[i] * (w[i] @@ h): i in [1..n] | w[i] ne 0 ];
      LL:= Neighbour(L, B, pih * w * G, pi * x, h, P, C, true); 
      if UseCallBack then keep, cont:= CallBack(Result, LL); end if;
      if keep then Append(~Result, LL); end if;
      if not cont or #Result ge Max then break; end if;
    end for;
  elif special then
    pi:= PrimitiveElement(P);
    G:= Matrix(k, n, [ h(e) : e in Eltseq(pi * T * Form * InternalConjugate(L, T)) ]);
    for i in [1..#LO] do
      w:= AutoOrbits select W ! LO[i] else Next(LO);
      Gw:= G * Matrix(1, Eltseq(w));
      ok:= exists(d){d: d in [1..n] | Gw[d] ne 0 }; assert ok;
      x:= &+[ B[i] * (w[i] @@ h): i in [1..n] | w[i] ne 0 ];
      nrm:= (x*Form, v) where v:= Vector([ a(e): e in Eltseq(x) ]);
      b:= (-h(nrm) / (2*Gw[d,1])) @@ h;
      x:= x+b*pi*B[d];
      nrm:= (x*Form, v) where v:= Vector([ a(e): e in Eltseq(x) ]);
      LL:= Neighbour(L, B, w*G, x, h, P, P, false);
      if UseCallBack then keep, cont:= CallBack(Result, LL); end if;
      if keep then Append(~Result, LL); end if;
      if not cont or #Result ge Max then break; end if;
    end for;
  else
    G:= Matrix(k, n, [ h(e) : e in Eltseq(T * Form * InternalConjugate(L, T)) ]);
    ram:= IsRamified(P);
    if ram then
      pi:= PrimitiveElement(P);
      S:= [x @@ h : x in k];
    else
      p:= P meet BaseRing(R);
      pi:= PrimitiveElement(p);
      kp, hp:= ResidueClassField(p);
      alpha:= k.1 @@ h;
      // k/kp has basis [1,alpha] and T is the relative trace:
      T:= Matrix(kp, 1, [2, Trace(alpha) @ hp ]);
    end if;
    for i in [1..#LO] do
      w:= AutoOrbits select W ! LO[i] else Next(LO);
      x:= &+[ B[i] * (w[i] @@ h) : i in [1..n] | w[i] ne 0 ];
      nrm:= (x*Form, v) where v:= Vector([ a(e): e in Eltseq(x) ]);
      if Valuation(nrm, P) le 0 then continue; end if;
      wG:= w*G;
      ok:= exists(j){j: j in [1..n] | wG[j] ne 0};
      assert ok;
      NL:= [];
      if not ram then
        s, V:= Solution(T, Vector(kp, [hp(-nrm/pi)]));
        l:= a((1/wG[j]) @@ h);
        S:= [ l*( x[1] @@ h + (x[2] @@ h)*alpha ) where x:= s+v : v in V ];
      end if;
      for s in S do
        LL:= Neighbour(L, B, wG, x+pi*s*B[j], h, P, P, false);
        if UseCallBack then keep, cont:= CallBack(Result, LL); end if;
        if keep then Append(~Result, LL); end if;
        if not cont or #Result ge Max then break i; end if;
      end for;
    end for;
  end if;
  return Result;
end function;

intrinsic Neighbours(L::LatModHerm, P::RngOrdIdl: AutoOrbits:= true, Max:= Infinity()) -> []
{The immediate P-neighbours of L}
  require Order(P) eq BaseRing(L): "Arguments are incompatible";
  require IsPrime(P): "Second argument must be prime";
  require not IsRamified(P) or Minimum(P) ne 2: "Second argument cannot be a ramified prime over 2";
  require IsModular(L, P) : "The lattice must be locally modular";
  require Rank(L) ge 2: "The rank of the lattice must be at least 2";
  require IsIsotropic(L, P): "The lattice must be locally isotropic";

  if Max cmpeq 1 or (AutoOrbits cmpeq true and not IsDefinite(L)) then
    AutoOrbits:= false;
  end if;

  return _Neighbours(L, P, [], AutoOrbits, Max);
end intrinsic;

function StdCallBack(List, L)
  keep:= forall{LL: LL in List | not IsIsometric(LL, L) };
//  keep, #List;
  return keep, true;
end function;

intrinsic IteratedNeighbours(L::LatModHerm, P::RngOrdIdl: AutoOrbits, CallBack:= false, Max:= Infinity()) -> []
{The iterated P-neighbours of L}
  require Order(P) eq BaseRing(L): "Arguments are incompatible";
  require IsPrime(P): "Second argument must be prime";
  require not IsRamified(P) or Minimum(P) ne 2: "Second argument cannot be a ramified prime over 2";
  require IsModular(L, P) : "The lattice must be locally modular";
  require Rank(L) ge 2: "The rank of the lattice must be at least 2";
  require IsIsotropic(L, P): "The lattice must be locally isotropic";

  if AutoOrbits cmpeq true and not IsDefinite(L) then
    AutoOrbits:= false;
  end if;
  if CallBack cmpeq false and IsDefinite(L) then CallBack:= StdCallBack; end if;

  Result:= [ L ];
  i:= 1;
  while (#Result lt Max) and (i le #Result) do
    // _Neighbours and the callback only add new lattices if not isometric to known ones and stop if Max is reached.
    // So we have nothing to at all.
    Result:= _Neighbours(Result[i], P, Result, AutoOrbits, Max : CallBack:= CallBack);
/*    for X in N do
      if forall{Y: Y in Result | not IsIsometric(Y, X) } then 
        Append(~Result, X);
        if #Result ge Max then return Result; end if;
      end if;
    end for;*/
    i+:= 1;
  end while;

  return Result;
end intrinsic;

function SmallestNeighbourPrime(L)
  S:= BaseRing(L);
  R:= BaseRing(S);
  bad:= { p: p in BadPrimes(L) | not IsModular(L, p) } join { p : p in MySupport(Discriminant(S)) | Minimum(p) eq 2 };

  if not IsDefinite(L) then return 0, bad; end if;

  // get first candidate
  m:= Rank(L);
  p:= 1; P:= 1;
  repeat
    p:= NextPrime(p);
    P:= (Type(R) eq RngInt select {p} else {d[1]: d in Decomposition(R, p)}) diff bad;
    if m eq 2 then P:= { x: x in P | IsIsotropic(L, x) }; end if;
  until #P ne 0;
  P:= Decomposition(S, Rep(P))[1,1];

  // get smallest solution
  I:= IsEmpty(bad) select 1 else ideal< S | &*bad >;
  n:= AbsoluteNorm(P);
  if n gt 1000 then
    PP:= PrimesUpTo(1000, NumberField(S): coprime_to:= I);
    for P in PP do
      if IsIsotropic(L, P) then return P, bad; end if;
    end for;
  end if;
  PP:= PrimesUpTo(n, NumberField(S): coprime_to:= I);
  for P in PP do
    if IsIsotropic(L, P) then return P, bad; end if;
  end for;
  error "Impossible failure.";
end function;

function NeighboursWithPPower(L, P, e)
  C:= [];
  for i in [1..e] do
    if i eq 1 then
      L:= Neighbours(L, P: Max:= 1, AutoOrbits:= false)[1];
    else
      N:= Neighbours(L, P: Max:= 2, AutoOrbits:= false);
      L:= N[1] eq C[#C] select N[2] else N[1];
    end if;
    Append(~C, L);
  end for;
  return C;
end function;

function GenusGenerators(L)
  R:= BaseRing(L);
  D:= Different(R);
  P0, bad:= SmallestNeighbourPrime(L);
  a:= Involution(L);
  bad:= R !! (IsEmpty(bad) select 1 else &*bad);

  // First we get the ideals coming from the C/C0 quotient.
  A:= IsAbsoluteOrder(R) select R else AbsoluteOrder(R);
  C, h:= ClassGroup(A);
  RR:= BaseRing(R);
  C0:= MySupport(D) join { R !! h(g): g in Generators(C) } where C, h:= ClassGroup(RR);
  Q0,q0:= quo< C | { (A !! I) @@ h : I in C0 } >;
  q0:= q0^-1 * h;

  PP:= [];
  if IsEven(Rank(L)) then
    // Get the primes that matter:
    for f in Factorization(D) do
      P, e:= Explode(f);
      G:= GenusSymbol(L, P);
      if exists{g: g in G | IsOdd(g[1])} then
        continue;
      elif Minimum(P) ne 2 then
        if exists{g: g in G | IsEven(g[2])} then continue; end if;
      else
        if exists{g: g in G | IsOdd(e + g[2]) or g[2]+e ne g[4]} then continue; end if;
      end if;
      Append(~PP, P);
    end for;

    if not IsEmpty(PP) then
      U,f:= UnitGroup(A);
      UU,ff:= UnitGroup(RR);
      norm:= hom< U -> UU | [ Norm(f(U.i)) @@ ff : i in [1..Ngens(U)] ] >;
      V:= VectorSpace( GF(2), #PP );
      VD:= [ Valuation(D,P): P in PP ];
      W, w:= quo< V | [ V | [ Valuation(R ! f(u)-1, PP[i]) ge VD[i] select 0 else 1 : i in [1..#PP] ] : u in Generators(Kernel(norm)) ] >;
      if #W eq 1 then PP:= []; end if;
      // Now W is isomorphic to Epsilon(L)/R.
    end if;
  end if;

  // Search for a suitable generating set
  Gens:= [];
  if IsEmpty(PP) then
//    "case 1";
    S:= [];
    Q, q:= quo< Q0 | S >;
    Work:= IsDefinite(L) select [ P0 ] else [];
    p:= 2;
    while #Q ne 1 do
      while IsEmpty(Work) do
        p:= NextPrime(p);
	Work:= [ d[1]: d in Decomposition(R, p) | IsSplit(d[1]) and Valuation(bad, d[1]) eq 0 ];
      end while;
      P:= Work[1];
      Remove(~Work, 1);
      c:= P @@ q0;
      o:= Order(c @ q);
      if o ne 1 then 
        Append(~S, c);
        Q, q:= quo< Q0 | S >;
        Append(~Gens, <P, o>);
      end if;
    end while;
  else
//    "case 2";
    // A certain central product of C/C0 and W is isomorphic to the group we need.
    // So we first compute the corresponding cocycle.
    // It would be easier to use generic abelian groups, but ....

    cocycle:= [ [ (W ! 0)^^#Q0 ]^^#Q0 ];
    num:= NumberingMap(Q0);
    Ideals:= [ i @@ num @ q0 : i in [1..#Q0] ];
    assert Ideals[1] eq 1*R;
    for i in [1..#Q0], j in [1..i] do
      ij:= num((i @@ num) * (j @@ num));
      I:= Ideals[i]*Ideals[j] * Ideals[ij]^-1;
      J:= I * a(I)^-1;
      ok, x:= IsPrincipal(J);
      u:= (-((RR ! Norm(x)) @@ ff)) @@ norm @ f;
      x:= x*u;
      assert Norm(x) eq 1;
      y:= w(V ! [ Valuation(x-1, PP[i]) ge VD[i] select 0 else 1 : i in [1..#PP] ]);
      cocycle[i,j]:= y;
      cocycle[j,i]:= y;
    end for;

    S:= {@ < Q0 ! 0, W ! 0 > @};
    Work:= IsDefinite(L) select [ P0 ] else [];
    p:= 2;
    while #S ne #Q0*#W do
       while IsEmpty(Work) do
        p:= NextPrime(p);
        Work:= [ d[1]: d in Decomposition(R, p) | IsSplit(d[1]) and Valuation(bad, d[1]) eq 0 ];
      end while;
      P:= Work[1];
      Remove(~Work, 1);
      c:= P @@ q0;
      i:= c @ num;
      I:= P * Ideals[i]^-1;
      J:= I * a(I)^-1;
      ok, x:= IsPrincipal(J);
      u:= (-((RR ! Norm(x)) @@ ff)) @@ norm @ f;
      x:= x*u;
      assert Norm(x) eq 1;
      y:= V ! [ Valuation(x-1, PP[i]) ge VD[i] select 0 else 1 : i in [1..#PP] ];
      // TODO: Should be OK. P in PP means P is ramified and odd. Then the local det is \pi / \bar(\pi) = -1 and -1 is never in E_1.
      idx:= Index(PP, P);
      if idx ne 0 then y[idx]+:= 1; end if;
      elt:= < c, w(y) >;
      elt1:= elt;
      o:= 1;
      size:= #S;
      while elt1 notin S do
        j:= elt1[1] @ num;
	for l in [1..size] do
	  elt2:= S[l];
	  k:= elt2[1] @ num;
	  prod:= < elt1[1] * elt2[1], elt1[2] + elt2[2] + cocycle[j,k] >;
	  Include(~S, prod);
	end for;
        elt1:= < elt[1] * elt1[1], elt[2] + elt1[2] + cocycle[i,j] >;
	o +:= 1;
      end while;
      assert #S eq size * o;
      if o ne 1 then Append(~Gens, <P, o>); end if;
    end while;
  end if;

  if IsDefinite(L) then 
    return Gens, P0;
  else
    return Gens, _;
  end if;
end function;


intrinsic GenusRepresentatives(L::LatModHerm : Max:= Infinity(), AutoOrbits) -> []
{Returns a system of representatives of the isometry classes in the genus of L}
  require Rank(L) ge 2 : "The rank of the lattice must be at least 2";

  definite:= IsDefinite(L);
  Gens, P0:= GenusGenerators(L);
  a:= Involution(L);

  LL:= [ L ];
  for g in Gens do
    if definite and g[1] eq P0 then continue; end if;
    I:= g[1]^(g[2]-1);
    J:= a(I)^-1;
    N:= NeighboursWithPPower(L, g[1], g[2]-1);
    Inter:= [];
    for i in [2..#LL] do
      M:= Module(LL[i]);
      IM:= I*M; JM:=J*M;
      Inter:= Inter cat [ SetupLatticeInSpace(L, (IM + Module(x)) meet JM) : x in N ];
    end for;
    LL:= LL cat N cat Inter;
  end for;
  assert #LL eq &*[ Integers() | g[2]: g in Gens | not definite or g[1] ne P0 ];
  assert forall{ X: X in LL | IsSameGenus(X, L) };

  if definite then
    Result:= [];
    for L in LL do
      // Should never happen!
      assert forall{ X: X in Result | not IsIsometric(X, L) };
      Result := Result cat IteratedNeighbours(L, P0 : AutoOrbits:= AutoOrbits, Max:= Max);
      Max:= Max - #Result;
    end for;
    // assert forall{<i,j>: j in [1..i-1], i in [1..#Result] | not IsIsometric(Result[i], Result[j]) };
  else
    Result:= LL;
  end if;
  return Result;
end intrinsic;

intrinsic IsOneClassGenus(L::LatModHerm : AutoOrbits) -> BoolElt
{Returns (#GenusRepresentatives(L) eq 1), but is much faster than calling GenusRepresentatives}
  require Rank(L) ge 2 : "The rank of the lattice must be at least 2";

  Gens, P0:= GenusGenerators(L);
  if #Gens ne 0 then return false; end if;
  if not IsDefinite(L) then return true; end if;

  X:= _Neighbours(L, P0, [], AutoOrbits, 1 : CallBack:= func< Res, LL | not iso, iso where iso:= IsIsometric(L, LL) >);
  return #X eq 0;
end intrinsic;



/*   The mass of a lattice    */
gauss0:= func< m, q | &*[ Parent(q) | (1-q^i): i in [1..m] ] >;
gauss := func< m, k, q | Abs(gauss0(m,q) / gauss0(k,q) / gauss0(m-k, q)) >;

// local factor for ramified even primes.
function local_factor_even(L, p)
  S:= BaseRing(L);
  P:= Decomposition(S,p)[1,1];
  val:= Valuation(Scale(L), P) div 2;
  if val ne 0 then
    s:= Type(p) eq RngIntElt select p else PrimitiveElement(p);
    L:= Rescale(L, s^-val);
  end if;

  G:= GenusSymbol(L, p);
  if not ({g[2]: g in G} subset {0,1}) then return 0; end if;
  m0:= G[ 1,2] eq 0 select G[ 1,1] else 0;
  m1:= G[#G,2] eq 1 select G[#G,1] else 0;
  q:= Norm(p);
  m:= m0 + m1;
  assert IsEven(m1);
  k:= m div 2;
  k0:= m0 div 2;
  k1:= m1 div 2;

  if IsOdd(m) then return 1/2 * gauss(k, k1, q^2); end if;

  // The horror starts here.
  e:= Valuation(Discriminant(S), p);
  f2:= e div 2;
  d:= (-1)^((m-1)*k) * Determinant(L`Form);
  // b <=> L \isom M \perp H0^* \perp H1^* where M is unimodular.
  b:= m1 eq 0 or (m0 ne 0 and G[1,4] ne G[2,4]);
  lf:= 1/2 * gauss(k, k0, q^2);
  l:= (b select G[1,4] else G[#G,4]) div 2;

  if IsLocalNorm(NumberField(S), d, p) then     // the hyperbolic case
    h0:= m0 eq 0 or G[1,4] eq 2*f2;
    h1:= m1 eq 0 or G[#G,4] eq 2*((e+1) div 2);

    // L_p = H(0)^k0 \perp H(1)^k1
    if h0 and h1 then
      t:= IsEven(e) select k1 else k0;
      return (q^t+1) * lf;
    end if;

    if not b then
      if IsEven(e) then
        return q^(m*(f2-l)-k1) * (q^m1-1) * lf;
      else
        return q^(m*(f2-l)+k0) * (q^m1-1) * lf;
      end if;
    else
      if IsEven(e) then
        return q^(m*(f2-1-l)+k1) * (q^m0-1) * lf;
      else
        return q^(m*(f2-l)-k0)   * (q^m0-1) * lf;
      end if;
    end if;
  else                                          // non-hyperbolic case
    if IsOdd(e) then
      if b and l eq (e-1)/2 then
        return (q^k0-1) * lf;
      elif b then
        return q^(m*(f2-l)-k0) * (q^m0-1) * lf;
      else
        return q^(m*(f2-l)+k0) * (q^m1-1) * lf;
      end if;
    end if;

    if b then
      return q^(m*(f2-1-l)+k1) * (q^m0-1) * lf;
    elif l eq e/2 then
      return (q^k1-1) * lf;
    else
      return q^(m*(f2-l)-k1) * (q^m1-1) * lf;
    end if;
  end if;
  error "impossible failure";
end function;

// the mass of a maximal integral lattice in the space of L_p
function local_factor_max(L, p)
  S:= BaseRing(L);
  m:= Rank(L);
  if IsOdd(m) and IsRamified(p, S) then
      return 1/2;
  end if;
  G:= GramMatrixOfBasis(L);
  disc:= (-1)^(m div 2) * Determinant(G);
  if not IsLocalNorm( NumberField(S), disc, p) then
    q:= Norm(p);
    if IsRamified(p, S) then
      return (q^m-1)/(2*q+2);
    else
      return (q^m-(-1)^m)/(q+1);
    end if;
  end if;
  return 1;
end function;


// Compute the local factor by brute force.
function local_factor_generic(L, p)
  def, a:= IsDefinite(L);
  S:= BaseRing(L);
  P:= Decomposition(S,p)[1,1];
  val:= Valuation(Norm(L), P);
  if IsRamified(p, S) then val div:= 2; end if;
  if Type(p) eq RngIntElt then
    s:= p^-val;
  else
    s:= PrimitiveElement(p)^-val;
    if def then
      R:= BaseRing(S);
      s*:= CRT(p, [1..AbsoluteDegree(R)], R ! 1, Signature(NumberField(R) ! (s*a)) );
    end if;
  end if;
  L:= Rescale(L, s);
  Chain:= [L];
  ok, LL:= IsMaximalIntegral(L, p);
  while not ok do
    Append(~Chain, LL);
    ok, LL:= IsMaximalIntegral(LL, p);
  end while;
  f:= local_factor_max(L, p);
  for i in [1..#Chain-1] do
    M, E:= MaximalSublattices(Chain[i+1], P: AutoOrbits:= def);
    f *:= &+[ E[j]: j in [1..#M] | IsLocallyIsometric(Chain[i  ], M[j], p) ];
    M, E:= MinimalSuperlattices(Chain[i], P: AutoOrbits:= def);
    f /:= &+[ E[j]: j in [1..#M] | IsLocallyIsometric(Chain[i+1], M[j], p) ];
  end for;
  return f;
end function;

// General local factor driver
function local_factor(L, p)
  S:= BaseRing(L);
  q:= Norm(p);
  ram:= IsRamified(p, S);
  if ram and IsEven(q) then 
    if IsMaximal(L, p) then return local_factor_max(L, p); end if;
    lf:= local_factor_even(L, p);
    // TODO: Use Cho's recipe if p unramified over Z.
    if lf eq 0 then lf:= local_factor_generic(L, p); end if;
    return lf;
  end if;
  split:= not ram and IsSplit(p, S);
  _, G, s:= JordanDecomposition(L, p);
  if #s eq 1 and not ram then return 1; end if;

  m:= Rank(L);
  if ram then
    f := GroupOrder("Sp", 2*(m div 2), q);
  else
    f := GroupOrder(split select "GL" else "U", m, q);
  end if;
  d:= ram select 1 else 2;
  N:= 0;
  for i in [1..#s] do
    mi:= Ncols(G[i]);
    ri:= &+[ Integers() | Ncols(G[j]): j in [i+1..#s] ];
    if ram then
      N -:= (s[i] div 2) * mi^2;
      if IsOdd(s[i]) then
        N -:= (mi+1) * (mi div 2);
        f /:= GroupOrder("Sp", mi, q);
      else
        det:= FixedField(L) ! get_disc(G[i]);
        norm:= IsLocalNorm(NumberField(S), det, p);
        f /:= GroupOrder( norm select "O+" else "O-", mi, q);
      end if;
    else
      N -:= s[i] * mi^2;
      f /:= GroupOrder(split select "GL" else "U", mi, q);
    end if;
    N -:= d*s[i]*mi*ri;
    N +:= Integers() ! (d*s[i]*mi*m/2);               // coming from the volume!
  end for;

  // fix the difference coming from the discriminant
  if ram and IsEven(m) then
    N +:= m div 2;
  end if;
//    <q,N,1/f>;
  return q^N * f;
end function;

intrinsic LocalFactor(L::LatModHerm, p::RngOrdIdl) -> FldRatElt
{The local factor of L at the prime p in the Minkowsi-Siegel mass formula}
  ok, p:= PrimeBelow(L, p);
  require ok : p;
  return local_factor(L, p);
end intrinsic;

intrinsic LocalFactor(L::LatModHerm, p::RngIntElt) -> FldRatElt
{The local factor of L at the prime p in the Minkowsi-Siegel mass formula}
  ok, p:= PrimeBelow(L, p);
  require ok : p;
  return local_factor(L, p);
end intrinsic;

intrinsic Mass(L::LatModHerm : StandardMass:= 0) -> FldRatElt
{The mass of the definite lattice L}
  require IsDefinite(L): "The lattice must be definite";
  m:= Rank(L);
  if m eq 0 then return 1; end if;

  S:= BaseRing(L);
  E:= NumberField(S);
  K:= BaseField(E);
  n:= AbsoluteDegree(K);

  mass:= StandardMass ne 0 select StandardMass else
    2^(1-n*m) * Abs(&* [ Rationals() | DedekindZetaExact(K, 1-i) : i in [2..m by 2] ] * &* [ Rationals() | DedekindZetaExact(E, 1-i : Relative) : i in [1..m by 2] ]);

  // The primes that can possibly give a local factor different from 1.
  bad:= BadPrimes(L : Discriminant);
  for p in bad do
    mass *:= LocalFactor(L, p);
  end for;
  return mass;
end intrinsic;


/*  Isotropy tests  */

intrinsic IsIsotropic(L::LatModHerm, p::RngOrdIdl : AmbientSpace:= false) -> BoolElt
{Tests is L is isotropic at p}
  ok, p:= PrimeBelow(L, p);
  require ok : p;
  F:= AmbientSpace select InnerProductMatrix(L) else GramMatrixOfBasis(L);
  r:= Ncols(F);
  if r ge 3 then return true; end if;
  if r eq 0 then return false; end if;
  d:= Determinant(F);
  if r eq 1 then return d eq 0; end if;
  return IsLocalNorm(NumberField(BaseRing(L)), -d, p);
end intrinsic;

intrinsic IsIsotropic(L::LatModHerm, p::RngIntElt : AmbientSpace:= false) -> BoolElt
{"} //"
  ok, p:= PrimeBelow(L, p);
  require ok : p;
  F:= AmbientSpace select InnerProductMatrix(L) else GramMatrixOfBasis(L);
  r:= Ncols(F);
  if r ge 3 then return true; end if;
  if r eq 0 then return false; end if;
  d:= Determinant(F);
  if r eq 1 then return d eq 0; end if;
  return IsLocalNorm(NumberField(BaseRing(L)), -d, p);
end intrinsic;

/*
intrinsic IsIsotropic(L::LatModHerm) -> BoolElt
{Tests if L is isotropic}
  r:= Rank(L);
  if r ge 3 then return true; end if;
  if r eq 0 then return false; end if;
  M:= Matrix(Basis(L`Module));
  d:= Determinant(M * L`Form * InternalConjugate(L, M));
  if r eq 1 then return d eq 0; end if;
  return IsNorm(NumberField(BaseRing(L)), -d);
end intrinsic;
*/

intrinsic HermitianFormWithInvariants(E::FldAlg, dim::RngIntElt, P::Setq, N::[]) -> Mtrx
{The hermitian form of rank dim over E with determinant a non-norm only at the non-split primes in P and N[i] negative entries on the i-th real place}
  requirege dim, 1;
  require Degree(E) eq 2: "The base field must be a degree 2 extension";
  K:= BaseField(E);
  R:= Integers(K);
  P:= Setseq(Set(P));
  if not IsEmpty(P) then
    if Type(P[1]) eq RngInt then P:= [ Generator(x): x in P ]; end if;
    U:= Universe(P);
    if Type(K) eq FldRat then
      require Type(U) eq RngInt : "The third argument must be a sequence of primes";
    else
      require Type(U) eq PowIdl and Ring(U) cmpeq R : "Incompatible arguments";
    end if;
    require forall{p: p in P | IsPrime(p) } : "The third arguments must be a sequence of prime ideals";
  end if;

  require forall{n: n in N | n in {0..dim}}: "Number of negative entries is impossible";
  Inf:= [ i: i in RealPlaces(K) | #Decomposition(E, i) eq 1 ];
  require #N eq #Inf: "Wrong number of real places";

  // ignore split places
  S:= Integers(E);
  P:= { p: p in P | not IsSplit(p, S) };
  
  I:= [ Inf[i]: i in [1..#N] | IsOdd(N[i]) ];
  require IsEven(#I + #P): "The invariants do not satisfy the product formula";

  // now it boils down to the construction of a quaternion algebra with given ramification
  x:= 2*e-Trace(e) where e:= PrimitiveElement(E);
  b:= K ! (x^2);
  a:= FindQA(b, P, I);

  D:= [ E | ];
  for i in [1..dim-1] do
    if #Inf eq 0 then
      Append(~D, 1);
    elif Type(K) eq FldRat then
      Append(~D, N[1] ge i select -1 else 1);
    else
      Append(~D, RealWeakApproximation(Inf, [ N[j] ge i select -1 else +1 : j in [1..#Inf] ]) );
    end if;
  end for;
  Append(~D, a*&*D);
  D:= DiagonalMatrix(D);
  dim0, P0, N0:= HermitianFormInvariants(D);
  assert dim eq dim0 and P eq P0 and N eq N0;
  return D;
end intrinsic;

/*   Constructing spaces with given invariants   */
intrinsic HermitianFormInvariants(M::Mtrx[FldAlg]) -> RngIntElt, {}, []
{The invariants describing the hermitian form M}
  E:= BaseRing(M);
  require Degree(E) eq 2 : "The matrix is not hermitian";
  K:= BaseField(E);
  a:= Automorphisms(E)[2];
  require TestHerm(M, a) : "The marix is not hermitian";

  d:= K ! Determinant(M);
  P:= MySupport(d) join MySupport(Discriminant(Integers(E)));
  P:= { p: p in P | not IsLocalNorm(E, d, p) };

  D:= MyGramSchmidt(M, a);
  ChangeUniverse(~D, K);
  I:= [ #[ d: d in D  | Evaluate(d, i) lt 0 ] : i in RealPlaces(K) | #Decomposition(E, i) eq 1 ];
  return Ncols(M), P, I;
end intrinsic;


// P must be inert and odd
function NonSquareNorm(P)
  assert Minimum(P) ne 2 and IsInert(P);
  R:= Order(P);
  p:= P meet BaseRing(R);
  if Type(p) eq RngInt then p:= Generator(p); end if;
  k,h:= ResidueClassField(P);
  kp, hp:= ResidueClassField(p);
  repeat
    r:= Random(k);
    u:= r @@ h;
  until r ne 0 and not IsSquare( Norm(u) @ hp );
  return u;
end function;

intrinsic FindLattice(M::LatModHerm, L::LatModHerm, p::RngOrdIdl) -> LatModHerm
{Constructs a sublattice X of M such that X_p is isometrix to L_p and X_q = M_q for all other primes q}
  ok, p:= PrimeBelow(L, p);
  require ok: p;
  require BaseRing(L) cmpeq BaseRing(M): "Incompatible base rings";
  require IsRationallyEquivalent(M, L, p): "The ambient hermitian spaces are not locally equivalent";
  require IsMaximalIntegral(M, p): "The first lattice is not locally maximal";
  P:= Decomposition(BaseRing(L), p)[1,1];
  require IsIntegral(L, P): "The second lattice must be locally integral";

  // The split case is easy
  if IsSplit(P) then
    pi:= PrimitiveElement(P);
    BM,_,SM:= JordanSplitting(M, p);
    BL,_,SL:= JordanSplitting(L, p);
    SM:= &cat [ [SM[i]^^Nrows(BM[i])] : i in [1..#BM] ]; 
    SL:= &cat [ [SL[i]^^Nrows(BL[i])] : i in [1..#BL] ];
    BM:= &cat [ Rows(b): b in BM ];
    E:= [ SL[i]-SM[i] : i in [1..#BM] ];
    B:= [ BM[i] * pi^E[i] : i in [1..#BM] ];
    LL:= SetupLatticeInSpace( M, (Module(B) + P^Max(E)*MM) meet P^Min(E)*MM ) where MM:= Module(M);

  // The inert case is also easy
  elif IsInert(P) then
    GenL:= GenusSymbol(L, p);
    r0:= &+[ Integers() | g[1]: g in GenL | IsEven(g[2]) ];
    nsn:= Minimum(p) eq 2 select 0 else NonSquareNorm(P);
    B,G,S:= JordanSplitting(M, p);
    assert Set(S) subset {0,1};
    if S[1] eq 0 then
      BB:= Rows(B[1]);
      m:= (#BB-r0) div 2;
      k,h:= ResidueClassField(p);
      Y:= [ BB[i]: i in [2*m+1..#BB] ];
      for i in [1..m] do
        // transform <BB[2i-1], BB[2i]> into H(0). Then go from there.
        ok, s:= IsSquare((-G[1,2*i,2*i]/G[1,2*i-1,2*i-1]) @ h);
        if ok then
          Append(~Y, BB[2*i] + (s @@ h)*BB[2*i-1] );
        else
          ok, s:= IsSquare((-G[1,2*i,2*i]/G[1,2*i-1,2*i-1] * Norm(nsn)) @ h);
          assert ok;
          Append(~Y, nsn * BB[2*i] + (s @@ h) * BB[2*i-1]);
        end if;
      end for;
      if #B eq 2 then Y cat:= Rows(B[2]); end if;
      LL:= SetupLatticeInSpace(M, (Module(Y) + P*MM) meet MM) where MM:= Module(M);
    else	// only happens in rank 1.
      LL:= M;
    end if;
    Y:= &cat[ Rows(x): x in JordanSplitting(LL, p) ];

    // Now Y generates the Gerstein reduction of L_p locally.
    // So we simply rescale the generators in Y appropriately and assemble
    // the global solution. 
    pi:= Type(p) eq RngIntElt select p else PrimitiveElement(p);
    i:= 1; j:= r0+1;
    for g in GenL do
      s:=  pi^(g[2] div 2);
      if IsEven(g[2]) then 
        for k in [1..g[1]] do Y[i] *:= s; i +:= 1; end for; 
      else 
        for k in [1..g[1]] do Y[j] *:= s; j +:= 1; end for;
      end if;
    end for;
    max:= GenL[#GenL,2];
    LL:= SetupLatticeInSpace(M, (Module(Y) + P^max*MM) meet MM) where MM:= Module(M);

  // The odd ramified case
  elif Minimum(p) ne 2 then
    // C contains the genus symbols of all Gerstein reduced lattices above L_p.
    c:= GenusSymbol(L,p);
    C:= [ c ];
    while c[#c, 2] ge 2 do
      c0:= [ x: x in c | x[2] in {0,2} ];
      if #c0 eq 2 then 
        c0:= [ < c0[1,1] + c0[2,1], 0, not (c0[1,3] xor c0[2,3]) > ];
      elif #c0 eq 1 then
        c0[1,2]:= 0;
      end if;
      c1:= [ x: x in c | x[2] in {1,3} ];
      if #c1 eq 2 then
        c1:= [ < c1[1,1] + c1[2,1], 1, not (c1[1,3] xor c1[2,3]) > ];
      elif #c1 eq 1 then
        c1[1,2]:= 1;
      end if;
      c:= c0 cat c1 cat [ < x[1], x[2]-2, x[3] > : x in c | x[2] ge 4];
      Append(~C, c);
    end while;

    // Now we construct a square-free lattice with genus symbol c.
    B,G,S:= JordanDecomposition(M, p);
    assert S subset {-1,0};
    B0:= S[#S] eq  0 select Rows(B[#B]) else [];
    B1:= S[ 1] eq -1 select Rows(B[ 1]) else [];
    r0:= c[1,2] eq 0 select c[1,1] else 0;
    for i in [1..(r0-#B0) div 2] do
      Append(~B0, B1[2*i-1]);
    end for;
    if #B0 eq 0 then
      LL:= P*M;
    else
      LL:= SetupLatticeInSpace(M, (Module(B0) + P*MM) meet MM) where MM:= Module(M);
    end if;
    assert GenusSymbol(LL, p) eq c;

    K:= FixedField(M);
    k,h:= ResidueClassField(p);
    for j in [#C-1..1 by -1] do
      c:= C[j];
      // Rescaling does the trick
      if {x[2]: x in c} meet {0,1} eq {} then
        s:= (C[1,1,2] - Valuation(Scale(LL), P)) div 2;
        LL:= P^s*LL;
	break;
      end if;
      // We have to split L_0 and L_1 once more.
      // Splitting L_1 is trivial.
      B,G,S:= JordanDecomposition(LL, p);
      if exists(r){g[1]: g in c | g[2] eq 1} then
        i:= Index(S, 1); assert i ne 0;
        Y1:= Rows(B[i])[1..r];
      else
        Y1:= [];
      end if;
      // Now we split L_0.
      Y0:= [];
      if exists(r){g[1]: g in c | g[2] eq 0} then
        assert S[1] eq 0;
        B:= Rows(B[1]);
	n:= #B;
        G:= G[1];
	g:= c[1];
	NN:= [ i: i in [1..n] | not IsSquare(h(K!G[i,i])) ];
	if (#NN eq 0 and not g[3]) then
	  repeat
	    s:= Random(k) @@ h;
	  until not IsSquare(h(K ! G[n-1,n-1] + s^2 * G[n, n]));
          Y0:= B[1..r-1] cat [ B[n-1] + s*B[n] ];
	else
          SS:= [ i: i in [1..n] | i notin Bad ] where Bad:= Set(NN);
	  if g[3] then
	    Y0:= [];
	  else
	    Y0:= [ NN[1] ]; Remove(~NN, 1);
	  end if;
	  if IsOdd(r-#Y0) then
	    if #SS eq 0 then
	      repeat
	        s:= Random(k) @@ h;
                tmp:= h(K ! G[n-1,n-1] + s^2 * G[n, n]);
	      until tmp ne 0 and IsSquare(tmp);
              v:= B[n-1];
	      B[n-1] +:= s*B[n];
	      B[n] -:= s*G[n,n]/G[n-1,n-1] * v;
	      NN:= [i: i in NN | i lt n-1 ];
	      SS:= [n-1, n];
	    end if;
	    Append(~Y0, SS[1]); Remove(~SS, 1);
	  end if;
	  Y0 cat:= NN[1..2*(#NN div 2)] cat SS;
          Y0:= B[Y0[1..r]];
	end if;
      end if;
      LL:= SetupLatticeInSpace(M, (Module(Y0 cat Y1) + P*MM) meet MM) where MM:= Module(LL);
      assert GenusSymbol(LL, p) eq c;
    end for;

  // The even ramified case
  // The approach below is VERY lame.
  // What we should do is the following:
  // 1. Find an (suffiently good approximation of an) isometry between the ambient spaces of M_p and L_p.
  // 2. Let Y_p be the image of L_p under this map. If it is good enough, then Y_p \isom L_p.
  // 3. Contsruct a lattice X in the space of M such that X_p = Y_p and X_q = M_q for all q \ne p.
  else
    k,h:= ResidueClassField(P);
    FF:= FieldOfFractions(BaseRing(M));
    m:= Rank(M);
    auto:= IsDefinite(M) and (q^m-1)/(q-1) le 2^14 where q:= AbsoluteNorm(P);
    V:= VectorSpace(k, m);
    Chain:= [L];
    ok, LL:= IsMaximalIntegral(L, p);
    while not ok do
      Append(~Chain, LL);
      ok, LL:= IsMaximalIntegral(LL, p);
    end while;
    Remove(~Chain, #Chain);
    LL:= M;
    for X in Reverse(Chain) do
      if auto then
        MS:= MaximalSublattices(LL, P : Max:= 1, AutoOrbits, CallBack:= func< List, x | iso, not iso where iso:= IsLocallyIsometric(x, X, p)>);
        assert #MS eq 1;
        LL:= MS[1];
      else
        BM:= Matrix(LocalBasis(Module(LL), P : Type:= "Submodule"));
        pM:= P*Module(LL);
        repeat
          repeat v:= Random(V); until v ne 0;
          KM:= KernelMatrix(Matrix(1, Eltseq(v)));
          KM:= ChangeRing(ChangeRing(KM, h^-1), FF);
          LL:= SetupLatticeInSpace(M, Module(Rows(KM * BM)) + pM);
        until IsLocallyIsometric(X, LL, p);
      end if;
    end for;
  end if;

  assert IsLocallyIsometric(L, LL, p);
  return LL;
end intrinsic;

intrinsic FindLattice(M::LatModHerm, L::LatModHerm, p::RngIntElt) -> LatModHerm
{"} //"
  ok, P:= PrimeAbove(M, p);
  require ok: P;
  return FindLattice(M, L, P);
end intrinsic;
