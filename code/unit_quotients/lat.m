//freeze;

/***************************************
*                                      *
* Lattices over number fields          *
*                                      *
* 2012--2019 M. Kirschmer and D. Lorch *
*                                      *
***************************************/

declare type LatMod[LatModElt];

declare attributes LatMod:
  Form, 
  Module,
  Discriminant,
  Norm,
  Scale,
  Definite,
  Diagonal,
  GenusSymbols,
  ZForms,
  AutoIsomData,
  Auto;

declare attributes LatModElt: Parent, Element;

declare verbose LatMod, 1;

import "herm.m" : MyGramSchmidt;

AlwaysTryDecomposition:= true;

/* TODO:
 * Nonsense in misc
 * Tidy genus.m
 * WittInvariant for Mtrx and RngOrdIdl (why is this missing?).
 */

Conj:= func< T, a | Transpose(Matrix( Ncols(T), [ a(x): x in Eltseq(T) ] )) >;
get_disc:= func< F | (-1)^(((m-1)*m) div 2) * Determinant(F) where m:= Ncols(F) >;

// There seems to be the concept of integral modules over Dedekind rings.
// This breaks the code at various places (for example if one wants to extend a module by some vector).
// So we rewrite every module over its field of fractions ... it just never stops ...
function FixModule(M)
  P:= PseudoBasis(M);
  U:= Universe(P);
  // Is there a better way to detect this?
  if ISA(Type(BaseRing(U[2])), RngOrd) then
    V:= VectorSpace(FieldOfFractions(BaseRing(M)), Degree(M));
    UU:= CartesianProduct(U[1], V);
    ChangeUniverse(~P, UU);
    M:= Module(P);
  end if;
  return M;
end function;

function OrthDec(F)
  m:= Ncols(F[1]);
  S:= {1..m};
  I:= [];
  while #S ne 0 do
    II:= { Rep(S) };
    repeat 
      n:= #II;
      II join:= { i: i in S | exists{<j,f> : j in II, f in F | f[i,j] ne 0 } };
    until n eq #II;
    Append(~I, II); S diff:= II;
  end while;
  if #I eq 1 then
    return false, _;
  end if;
  vprint LatMod: "using orthogonal dec";
  B,G:= OrthogonalDecomposition(F);
  vprint LatMod: "Block sizes:", Sort([Nrows(b): b in B]);
  assert #G ge #I;
  for i in [1..#G] do 
    _, T:= LLLGram(G[i,1]: Delta:= 0.999, Eta:= 0.501);
    G[i]:= [ T*g*Transpose(T): g in G[i] ];
    B[i]:= T*B[i];
  end for;
  T:= RMatrixSpace(Integers(), 0, m) ! 0;
  n:= 1;
  A:= []; 
  M:= MatrixRing(Integers(), m);
  while #B ne 0 do
    idx:= Nrows(T)+1;
    c:= 1;
    AA:= AutomorphismGroup(G[1]);
    T:= VerticalJoin(T, B[1]);
    for i in [#G..2 by -1] do
      if Ncols(G[i,1]) ne Ncols(G[1,1]) then continue; end if;
      ok, TT:= IsIsometric(G[i], G[1]);
      if ok then T:= VerticalJoin(T, TT*B[i]); Remove(~G,i); Remove(~B, i); c+:= 1; end if;
    end for;
    Remove(~G, 1); Remove(~B, 1);
    n *:= (#AA)^c * Factorial(c); 
    W:= WreathProduct(AA, Sym(c));
    A cat:= [ InsertBlock(M!1, a, idx,idx) : a in Generators(W) ]; 
  end while;
  A:= [ TInv * a * T : a in A ] where TInv:= T^-1;
  assert forall{<f,a>: f in F, a in A | a*f*Transpose(a) eq f }; 
  A:= sub< GL(m, Integers()) | A >;
//  assert #A eq n;
  A`Order:= n;
  return true, A;
end function;

function MyAuto(F)
  ok, A:= OrthDec(F);
  if ok then return A; end if;
  return AutomorphismGroup(F);
end function;

// stuff that should be in Magma

intrinsic ZBasis(M::ModDed) -> []
{A Z-basis for the module M}
  P:= PseudoBasis(M);
  B:= [ b*p[2] : b in AbsoluteBasis(p[1]), p in P ];
  return B;
end intrinsic;

intrinsic LocalBasis(M::ModDed, p::RngOrdIdl : Type:= "") -> []
{A basis of a free module that coincides with M at the prime ideal p}
  require Order(p) cmpeq BaseRing(M) and IsPrime(p) : "Incompatible arguments";
  require Type in {"", "Submodule", "Supermodule"} : "Type must be \"Submodule\" or \"Supermodule\" when specified";
  if Type eq "" then
    pi:= PrimitiveElement(p);
  end if;
  B:= [ EmbeddingSpace(M) | ];
  for pp in PseudoBasis(M) do
    g:= Generators(pp[1]);
    if #g eq 1 then x:= g[1];
    elif Type eq "" then x:= pi^Valuation(pp[1], p);
    else
      Fact:= Factorization(pp[1]);
      Fact:= Type eq "Submodule" select [ f: f in Fact | f[1] eq p or f[2] gt 0 ]
                                   else [ f: f in Fact | f[1] eq p or f[2] lt 0 ];
      if forall{ f: f in Fact | f[1] ne p } then Append(~Fact, <p, 0>); end if;
      x:= WeakApproximation([ f[1] : f in Fact ], [ f[2] : f in Fact ]);
    end if;
    assert Valuation(x, p) eq Valuation(pp[1], p);
    Append(~B, pp[2]*x);
  end for;
  return B;
end intrinsic;


// Magmas FewGenerators for matrix group does exactly ... NOTHING!
// So try to be a bit more clever (code taken from bruce).
function order_mod(g, o, S)
  for d in Divisors(o) do
    if g^d in S then return d; end if;
  end for;
  error "Impossible failure";
end function;

function redgen(G)
  if Ngens(G) le 2 then return G; end if;
  gens := [G.i: i in [1 .. Ngens(G)]];
  ords := [Order(x): x in gens];
  S := sub<G | >;
  repeat
    _, p := Max(ords);
    S := sub<G | S, gens[p]>;
    ords := [order_mod(gens[i], ords[i], S): i in [1 .. #gens]];
    assert ords[p] eq 1;
  until #S eq #G;
//  "old", Ngens(G), "new", Ngens(S);
  return S;
end function;

function MySupport(I)
  T:= Type(I);
  if T eq FldRatElt then
    return { f[1]: f in FactorizationOfQuotient(I) };
  elif T eq RngIntElt then
    return { f[1]: f in Factorization(I) };
  elif ISA(T, FldAlgElt) then
    return Support( I * Integers(Parent(I)) );
  elif ISA(T, RngOrdElt) then
    return Support( I * Parent(I) );
  else
    return Support(I);
  end if;
end function;
// end stuff that should be in Magma

// The Element struff
intrinsic Print(v::LatModElt, Level::MonStgElt)
{internal}
  printf "%o", v`Element;
end intrinsic;

intrinsic Parent(v::LatModElt) -> LatMod
{internal}
  return v`Parent;
end intrinsic;

intrinsic Vector(v::LatModElt) -> ModDedElt
{internal}
  FF:= FieldOfFractions(BaseRing(v`Parent));
  return Vector(FF, Eltseq(v`Element));
end intrinsic;

intrinsic Element(v::LatModElt) -> ModDedElt
{internal}
  return v`Element;
end intrinsic;

intrinsic Eltseq(v::LatModElt) -> []
{internal}
  return Eltseq(v`Element);
end intrinsic;

function InitElement(L, x)
  v:= New(LatModElt);
  v`Parent:= L; v`Element:= L`Module ! x;
  return v;
end function;

intrinsic 'eq'(v::LatModElt, w::LatModElt) -> BoolElt
{internal}
  return v`Parent eq w`Parent and v`Element eq w`Element;
end intrinsic;

intrinsic Hash(x::LatModElt) -> RngIntElt
{internal}
  return Hash(x`Element);
end intrinsic;

intrinsic '+'(v::LatModElt, w::LatModElt) -> LatModElt
{internal}
  L:= v`Parent;
  require L eq w`Parent: "Incompatible elements";
  return InitElement(L, v`Element + w`Element);
end intrinsic;

intrinsic '-'(v::LatModElt, w::LatModElt) -> LatModElt
{internal}
  L:= v`Parent;
  require L eq w`Parent: "Incompatible elements";
  return InitElement(L, v`Element - w`Element);
end intrinsic;

function do_mult(v, a: inv:= false)
  x:= Vector(v);
  ok, a:= IsCoercible(BaseRing(x), a);
  if not ok then return ok, "Incompatible arguments"; end if;
  if inv then
    if IsZero(a) then return false, "Division by 0"; end if;
    a:= a^-1;
  end if;
  ok, x:= IsCoercible(v`Parent, x * a);
  if ok then return ok, x; end if;
  return false, "Element is not in the lattice";
end function;

intrinsic '*'(v::LatModElt, a::RngElt) -> LatModElt
{internal}
  ok, x:= do_mult(v, a);
  require ok : x;
  return x;
end intrinsic;

intrinsic '*'(a::RngElt, v::LatModElt) -> LatModElt
{internal}
  ok, x:= do_mult(v, a);
  require ok : x;
  return x;
end intrinsic;

intrinsic '/'(v::LatModElt, a::RngElt) -> LatModElt
{internal}
  ok, x:= do_mult(v, a: inv);
  require ok : x;
  return x;
end intrinsic;

intrinsic InnerProduct(v::LatModElt, w::LatModElt) -> RngElt
{The inner product (v,w)}
  L:= v`Parent;
  F:= InnerProductMatrix(L);
  o:= IsOrthogonal(L);
  require F cmpeq InnerProductMatrix(w`Parent) and o eq IsOrthogonal(w`Parent): "Incompatible arguments";
  y:= o select Vector(w) else (Vector( [a(e): e in Eltseq(w) ] ) where a:= Involution(w`Parent));
  return (Vector(v)*F, y);
end intrinsic;

intrinsic Norm(v::LatModElt) -> RngElt
{The norm of v}
  L:= v`Parent;
  x:= Vector(v);
  y:= IsOrthogonal(L) select x else (Vector( [ a(e) : e in Eltseq(x) ] ) where a:= Involution(L));
  return (x*InnerProductMatrix(L), y);
end intrinsic;

intrinsic IsZero(v::LatModElt) -> BoolElt
{Returns true if v is the zero vector}
  return IsZero(v`Element);
end intrinsic;

intrinsic Matrix(S::[LatModElt]) -> Mtrx
{The matrix of the sequence of lattice elements}
  if #S eq 0 then 
    L:= Universe(S);
    return Matrix(FieldOfFractions(BaseRing(L)), 0, Degree(L), []);
  end if;
  return Matrix([ Eltseq(s) : s in S ]);
end intrinsic;

// end of Element stuff


intrinsic Print(L::LatMod, Level::MonStgElt)
{internal}
  F:= InnerProductMatrix(L);
  if Level eq "Maximal" then
    printf "Lattice of dimension %o over %o\nGiven by %o\nand inner product matrix:\n%o", Rank(L), BaseRing(L), L`Module, F;
  elif Level eq "Magma" then
    form := Sprintf("%m", F);
    PB := PseudoBasis(L`Module);
    module := "Module([ ";
    for j in [1..#PB] do
      pb := PB[j];
      id := Sprintf("%m", pb[1]);
      vec := " Vector([ ";
      for k in [1..Ncols(pb[2])] do
        // the vector elements:
        vec cat:= Sprintf("%m%o", pb[2][k], k lt Ncols(pb[2]) select "," else "");
      end for;
      vec cat:= "])";
      module cat:= Sprintf("<%o , %o>", id, vec);
      module cat:= j lt #PB select "," else "";
    end for;
    module cat:=" ])";

    s := Sprintf("%o(%o, %o)", IsOrthogonal(L) select "LatticeModule" else "HermitianLattice", module, form);

    // find the text to replace:
    t:= Sprintf("%m", BaseRing(L));

    // StringReplace:
    repeat
      p := Position(s, t);
      if p eq 0 then break; end if;
      s := Substring(s, 1, p-1) cat "R" cat Substring(s, p+#t, #s - (p+#t) + 1);
    until false;

    printf "%o where R is MaximalOrder(%o)", s, t;
  else
    printf "Lattice of dimension %o over %o", Rank(L), BaseRing(L);
  end if;
end intrinsic;

intrinsic IsCoercible(L::LatMod, A::. ) -> BoolElt, .
{internal}
  T:= Type(A);
  if T eq LatModElt then	// A is lattice-module element
    if A`Parent eq L then return true, A; end if;
    A:= A`Element;		// A = 0 in some compatible ring
  elif ISA(T, RngElt) then
    R:= BaseRing(L);
    ok, x:= IsCoercible(R, A);
    if ok and IsZero(x) then A:= [ (R ! 0) ^^ Degree(L) ]; end if;
  end if;
  // Now A should be a sequence, vector or module element.
  ok, x:= IsCoercible(L`Module, A);
  if ok then
    return true, InitElement(L, x);
  end if;
  return false, "invalid coercion";
end intrinsic;

intrinsic 'in'(A::LatModElt, L::LatMod) -> BoolElt
{internal}
  ok := IsCoercible(L, A);
  return ok;
end intrinsic;

intrinsic 'in'(A::ModTupRngElt, L::LatMod) -> BoolElt
{internal}
  ok := IsCoercible(L, A);
  return ok;
end intrinsic;


/* The constructors */

intrinsic LatticeModule(M::ModDed, F::Mtrx) -> LatMod
{The lattice consisting of the module M and the bilinear form F}
  R:= BaseRing(M);
  require ISA(Type(R), RngOrd) : "The base ring must be an order in a number field";
  //require Ncols(F) eq Dimension(M) : "The arguments must have the same dimensions";
  require IsSymmetric(F): "The second argument must be symmetric";
  FF:= FieldOfFractions(R);
  ok, F:= CanChangeRing(F, FF);
  require ok: "The arguments must have the same base ring";
  L:= New(LatMod);
  L`Module:= FixModule(M);
  L`Form:= F;
  return L;
end intrinsic;

intrinsic LatticeModule(F::Mtrx) -> LatMod
{The free lattice with bilinear form F}
  require IsSymmetric(F): "The second argument must be symmetric";
  R:= BaseRing(F);
  T:= Type(R);
  if ISA(T, FldNum) then
    R:= Integers(R); 
  elif ISA(T, FldOrd) then
    R:= Order(R);
  else
    require ISA(T, RngOrd): "The form must be given over some number field";
  end if;
  return LatticeModule(Module(R, Ncols(F)), F);
end intrinsic;

intrinsic LatticeModule(B::Mtrx, F::Mtrx) -> LatMod
{The lattice generated by the rows of B with bilinear form F}
  require IsSymmetric(F): "The second argument must be symmetric";
  require Ncols(B) eq Ncols(F) : "The arguments must have the smae dimensions";
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
  return LatticeModule(Module(Rows(B)), F);
end intrinsic;

function SetupLatticeInSpace(L, M)
  LL:= New(Type(L));
  LL`Form:= L`Form;
  LL`Module:= FixModule(M);
  if assigned(L`Definite) and IsFull(L) and IsFull(LL) then LL`Definite:= L`Definite; end if;
  return LL;
end function;

/* end of constructors */

intrinsic GramMatrix(L::LatMod, S::Mtrx) -> Mtrx
{The Gram matrix of the sequence S (or the rows of S) in the ambient space of L}
  F:= InnerProductMatrix(L);
  ok, S:= CanChangeRing(S, BaseRing(F));
  require Ncols(S) eq Nrows(F) and ok : "Incompatible arguments";
  return S * F * InternalConjugate(L, S);
end intrinsic;

intrinsic GramMatrix(L::LatMod, S::[]) -> Mtrx
{"} //"
  F:= InnerProductMatrix(L);
  try
    S:= Matrix(S);
  catch e
    require false: "Incompatible arguments";
  end try;
  return GramMatrix(L, S);
end intrinsic;

intrinsic GramMatrix(L::LatMod) -> Mtrx, []
{The Gram matrix of a minimal set of generators of L}
  B:= Generators(L : Minimal);
  return GramMatrix(L, B), B;
end intrinsic;

//intrinsic GramMatrixOfBasis(L::LatMod) -> Mtrx, [], []
//{internal}
 function  GramMatrixOfBasis(L)
  P:= PseudoBasis(Module(L));
  U:= Universe(P);
  C:= [ U[1] | p[1]: p in P ];
  B:= [ U[2] | p[2]: p in P ] ;
  return GramMatrix( L, B ), C, B;
end function;
//end intrinsic;

function MyDiagonal(L, Ambient)
  if Ambient and not IsFull(L) then
    F:= InnerProductMatrix(L);
    return IsHermitian(L) select ChangeUniverse(MyGramSchmidt(F, Involution(L)), FixedField(L)) else Diagonal(OrthogonalizeGram(F));
  elif not assigned L`Diagonal then
    F:= IsFull(L) select InnerProductMatrix(L) else GramMatrixOfBasis(L);
    L`Diagonal:= IsHermitian(L) select ChangeUniverse(MyGramSchmidt(F, Involution(L)), FixedField(L)) else Diagonal(OrthogonalizeGram(F));
  end if;
  return L`Diagonal;
end function;

intrinsic IsOrthogonal(L::LatMod) -> BoolElt
{Returns true if the lattice L was created as an orthogonal lattice}
  return Type(L) eq LatMod;
end intrinsic;

intrinsic IsHermitian(L::LatMod) -> BoolElt
{Returns true if the lattice L was created as an hermitian lattice}
  return Type(L) eq LatModHerm;
end intrinsic;

intrinsic Involution(L::LatMod) -> Map
{The involution to the hermitian/ortohogonal space of L}
  return IsOrthogonal(L) select func< x | x > else Automorphisms(FieldOfFractions(BaseRing(L)))[2];
end intrinsic;

intrinsic InternalConjugate(L::LatMod, X::Mtrx) -> Mtrx
{internal}
  if IsHermitian(L) then 
    X:= Matrix(BaseRing(X), Ncols(X), [a(x): x in Eltseq(X) ]) where a:= Involution(L);
  end if;
  return Transpose(X);
end intrinsic;

intrinsic FixedField(L::LatMod) -> FldAlg
{The fixed field of the involution of L}
  return NumberField( IsOrthogonal(L) select BaseRing(L) else BaseRing(BaseRing(L)) );
end intrinsic;

intrinsic InnerProductMatrix(L::LatMod) -> Mtrx
{Returns the inner product matrix of the lattice L}
  return L`Form;
end intrinsic;

intrinsic Module(L::LatMod) -> ModDed
{Returns the module strucure of L}
  return L`Module;
end intrinsic;

intrinsic IsFree(L::LatMod) -> BoolElt
{Tests if the lattice L is a free module}
  return IsFree(L`Module);
end intrinsic;

intrinsic IsZero(L) -> BoolElt
{Tests if the lattice L is zero}
  return Rank(L) eq 0;
end intrinsic;

function DecideDefinite(L, AmbientSpace)
  K:= FixedField(L);
  R:= Integers(K);
  if not IsTotallyReal(K) or (IsHermitian(L) and not IsTotallyComplex(NumberField(L))) then return R ! 0; end if;
  D:= MyDiagonal(L, AmbientSpace);
  Signs:= [];
  for v in RealPlaces(K) do
    s:= { Sign( Evaluate(x, v) ) : x in D };
    if #s ne 1 then return 0; end if;
    Append(~Signs, Rep(s));
  end for;
  return #Set(Signs) eq 1 select R ! Signs[1] else CRT(1*R, [1..#Signs], R ! 1, Signs);
end function;

intrinsic IsDefinite(L::LatMod : AmbientSpace:= false) -> BoolElt, RngOrdElt
{Whether a number field lattice is definite. If so, an element
 such that the rescaled lattice is totally positive definite}

  if not AmbientSpace and assigned L`Definite then
    def:= L`Definite;
  else
    def:= DecideDefinite(L, AmbientSpace);
    if not AmbientSpace or IsFull(L) then L`Definite:= def; end if;
  end if;
  if def eq 0 then return false, _; end if;
  return true, def;
end intrinsic;

intrinsic IsPositiveDefinite(L::LatMod : AmbientSpace:= false) -> BoolElt
{Tests if the lattice L or its ambient space is (totally) positive definite}
  ok, a:= IsDefinite(L : AmbientSpace:= AmbientSpace);
  return ok and IsOne(a);
end intrinsic;

intrinsic IsTotallyPositiveDefinite(L::LatMod : AmbientSpace:= false) -> BoolElt
{"} // "
  ok, a:= IsDefinite(L : AmbientSpace:= AmbientSpace);
  return ok and IsOne(a);
end intrinsic;

intrinsic ChangeRing(L::LatMod, R::RngInt) -> Lat
{Return a lattice over the integers if L is over QNF()}
  require AbsoluteDegree(BaseRing(L)) eq 1 : "The module must be over the integers of QNF()";
  require IsOrthogonal(L) and IsPositiveDefinite(L) : "Over Z, only positive definite quadratic lattices are supported in Magma";
  B:= ChangeRing(Matrix( Generators(L: Minimal) ), Rationals());
  F:= ChangeRing(InnerProductMatrix(L), Rationals());
  return Lattice(B, F);
end intrinsic;

intrinsic CanChangeRing(L::LatMod, R::Rng) -> BoolElt, Lat
{Tests if the lattice L can be written over the ring R}
  if ISA(Type(R), RngOrd) then R:= FieldOfFractions(R); end if;
  try
    ok, FF:= CanChangeRing(InnerProductMatrix(L), R);
    if not ok then return false, _; end if;
    G:= [ Vector(R, Vector(x)): x in Generators(L) ];
    return true, LatticeModule(Matrix(G), FF);
  catch e
    return false, _;
  end try;
end intrinsic;

intrinsic ChangeRing(L::LatMod, R::Rng) -> Lat
{The lattice L over ring R}
  ok, L:= CanChangeRing(L, R);
  require ok: "The operation is not supported";
  return L;
end intrinsic;


intrinsic Discriminant(L::LatMod) -> RngOrdFracIdl
{The discriminant ideal of L}
  if not assigned L`Discriminant then
    G, C, B:= GramMatrixOfBasis(L);
    L`Discriminant:= Determinant(G) * I * Involution(L)(I) where I:= &* C;
  end if;
  return L`Discriminant;
end intrinsic;

intrinsic ZBasis(L::LatMod) -> []
{A Z-basis for the lattice L}
  return ChangeUniverse(ZBasis(L`Module), L);
end intrinsic;

intrinsic AbsoluteBasis(L::LatMod) -> []
{"} // "
  return ChangeUniverse(ZBasis(L`Module), L);
end intrinsic;

intrinsic Generators(L::LatMod : Minimal:= false) -> []
{A sequence of elements that generate L over its base ring}
  return ChangeUniverse(Generators(L`Module : Minimal:= Minimal), L);
end intrinsic;

intrinsic LocalBasis(L::LatMod, p::RngOrdIdl : Type:= "") -> []
{A basis of a free module that coincides with L at the prime ideal p}
  require Order(p) cmpeq BaseRing(L) : "Incompatible arguments";
  require Type in {"", "Submodule", "Supermodule"} : "Type must be \"Submodule\" or \"Supermodule\" when specified";
  return LocalBasis(Module(L), p: Type:= Type);
end intrinsic;

intrinsic BasisMatrix(L::LatMod) -> Mtrx
{Returns Matrix(Basis(Module(L)))}
  M:= Matrix(Basis(Module(L)));
  return ChangeRing(M, FieldOfFractions(BaseRing(L)));		// without this, the result sometimes lives over Z_K (shrug..)
end intrinsic;

intrinsic BaseRing(L::LatMod) -> RngOrd
{The base ring of the lattice L}
  return BaseRing(L`Module);
end intrinsic;

intrinsic NumberField(L::LatMod) -> FldNum
{The number field the lattice is over}
  return NumberField(BaseRing(L`Module));
end intrinsic;

intrinsic Rank(L::LatMod) -> RngIntElt
{The rank or dimension of the lattice L}
  return Dimension(L`Module);
end intrinsic;

intrinsic Dimension(L::LatMod) -> RngIntElt
{"} // "
  return Dimension(L`Module);
end intrinsic;

intrinsic Degree(L::LatMod) -> RngIntElt
{The degree of the lattice L}
  return Degree(L`Module);
end intrinsic;

intrinsic IsFull(L::LatMod) -> BoolElt
{Returns true if the lattice is of full rank}
  return Dimension(L) eq Degree(L);
end intrinsic;

intrinsic 'eq'(L1::LatMod, L2::LatMod) -> BoolElt
{internal}
  return InnerProductMatrix(L1) cmpeq InnerProductMatrix(L2) and L1`Module eq L2`Module;
end intrinsic;

intrinsic 'subset'(L1::LatMod, L2::LatMod) -> BoolElt
{internal}
  require InnerProductMatrix(L1) cmpeq InnerProductMatrix(L2) : "Incompatible arguments";
  return L1`Module subset L2`Module;
end intrinsic;

intrinsic Perp(L1::LatMod, L2::LatMod) -> LatMod
{The lattice \{x in L1 | (x,L2) = 0\}}
  require L2 subset L1 : "The second argument must be a sublattice of the first";

  // there must be a better way.
  FF:= FieldOfFractions(BaseRing(L1));
  Z:= ZBasis(L1);
  G:= Generators(L2);
  F:= Matrix(Rationals(), [ &cat [ Flat(InnerProduct(z, g)) : g in G ]  : z in Z ] );
  F:= Matrix(Integers(), Denominator(F) * F);
  K:= KernelMatrix(F);
  return SetupLatticeInSpace(L1, Module(Rows(Matrix(FF, K) * Matrix(Z))));
end intrinsic;

function GetGens(V, x)
  ok, y:= IsCoercible(V, x);
  if ok then return ok, [y]; end if;

  if ISA(Type(x), {LatMod, PMat}) then x:= Module(x); end if;
  if ISA(Type(x), ModDed) then 
    x:= Generators(x);
  elif ISA(Type(x), LatModElt) then
    x:= [ Vector(x) ];
  end if;

  if ISA(Type(x), {Setq, List, Tup}) then
    res:= [];
    for e in x do
      ok, y:= GetGens(V, e);
      if not ok then return ok, _; end if;
      res cat:= y;
    end for;
    return true, res;
  else
    return false, _;
  end if;
end function;

intrinsic SubConstructor(L::LatMod, RHS::Tup) -> LatMod
{The sublattice of L generated by RHS}
  ok, Gens:= GetGens(L, RHS);
  require ok: "Cannot coerce RHS into LHS";
  LL:= SetupLatticeInSpace(L, sub< Module(L) | [ x`Element : x in Gens ] >);
  return LL;
end intrinsic;

intrinsic ExtConstructor(L::LatMod, RHS::Tup) -> LatMod
{The lattice generated by L and RHS}
  M:= Module(L);
  ok, Gens:= GetGens(EmbeddingSpace(M), RHS);
  require ok: "Cannot coerce RHS into LHS";
  LL:= SetupLatticeInSpace(L, M + Module(Gens));
  return LL;
end intrinsic;

ExpandVector:= func< v | Vector(&cat [ Flat(x): x in Eltseq(v) ]) >;

intrinsic ZLattice(L::LatMod : TraceForm:= false) -> Lat, Map
{The lattice with a Z-basis of L}
  A:= AbsoluteBasis(BaseRing(L));
  B:= ZBasis(L);
  B:= Matrix(Rationals(), #B, Degree(L) * #A, [ ExpandVector(v) : v in B ] );
  if TraceForm then
    require IsPositiveDefinite(L) : "Magma only supports definite lattices over the integers";
    F:= InnerProductMatrix(L);
    AC:= [ I(x): x in A ] where I:= Involution(L);
    M:= Matrix( [ [ b*c: c in AC ] : b in A ] );
    F:= TensorProduct(F, Matrix(BaseRing(F), M)); 
    F:= Matrix(Ncols(F), [ AbsoluteTrace(x): x in Eltseq(F) ]);
    LL:= Lattice(B, F);
  else
    LL:= LatticeWithBasis(B);
  end if;
  h:= map< L -> LL | x :-> ExpandVector(x), y :-> L ! [ &+[ A[i] * p[i]: i in [1..#A] ] : p in  Partition(Eltseq(y), #A) ] >; 
  return LL, h;
end intrinsic;

intrinsic QuoConstructor(L::LatMod, RHS::Tup) -> GrpAb, Map
{The quotient of L by RHS}
  try
    LL:= SubConstructor(L, RHS);
  catch e;
    require false : "The right hand side does not yield a sublattice of the left hand side";
  end try;
  L1, h1:= ZLattice(L);
  L2:= ZLattice(LL);
  Q, h:= quo< L1 | L2 >;
  return Q, h1*h;
end intrinsic;

intrinsic '*'(a::RngElt, L::LatMod) -> LatMod
{internal}
  FF:= FieldOfFractions(BaseRing(L));
  ok, a:= IsCoercible(FF, a);
  require ok : "The arguments are not compatible";
  return SetupLatticeInSpace(L, a*Module(L));
end intrinsic;

intrinsic '*'(I::RngOrdFracIdl, L::LatMod) -> LatMod
{internal}
  require Order(I) cmpeq BaseRing(L) : "The arguments are not compatible";
  return SetupLatticeInSpace(L, I * Module(L));
end intrinsic;

intrinsic '*'(L::LatMod, I::RngOrdFracIdl) -> LatMod
{internal}
  require Order(I) cmpeq BaseRing(L) : "The arguments are not compatible";
  return SetupLatticeInSpace(L, I * Module(L));
end intrinsic;

intrinsic '*'(L::LatMod, a::Mtrx) -> LatMod
{internal}
  F:= FieldOfFractions(BaseRing(L));
  ok, a:= CanChangeRing(a, F);
  require ok and Nrows(a) eq Degree(L) : "The arguments are not compatible";
  P:= [ < p[1], p[2] * a> : p in PseudoBasis(Module(L)) ];
  return SetupLatticeInSpace(L, Module(P));
end intrinsic;


intrinsic '+'(L::LatMod, M::LatMod) -> LatMod
{internal}
  require InnerProductMatrix(L) cmpeq InnerProductMatrix(M) : "The arguments are not compatible";
  return SetupLatticeInSpace(L, Module(L) + Module(M));
end intrinsic;

intrinsic 'meet'(L::LatMod, M::LatMod) -> LatMod
{internal}
  require InnerProductMatrix(L) cmpeq InnerProductMatrix(M) : "The arguments are not compatible";
  return SetupLatticeInSpace(L, Module(L) meet Module(M));
end intrinsic;

intrinsic Rescale(L::LatMod, a::RngElt) -> LatMod
{Rescale the form of L with the scalar a}
  if IsOne(a) then return L; end if;
  ok, a:= IsCoercible(FixedField(L), a);
  require ok : "The arguments are not compatible";
  LL:= New(Type(L));
  LL`Form:= a*L`Form;
  LL`Module:= L`Module;
  // We could also copy norm, scale, volume etc.
  if assigned L`Auto then LL`Auto:= L`Auto; end if;
  if assigned L`Definite then
    if L`Definite eq 0 then
      def:= 0;
    else	// L definite, so FF totally real...
      a *:= L`Definite;
      def:= IsTotallyPositive(a) select 1 else a*Denominator(a);
    end if;
    LL`Definite:= Parent(L`Definite) ! def;
  end if;
  return LL;
end intrinsic;

intrinsic IsIntegral(L::LatMod) -> BoolElt
{Returns true if and only if L is integral (with respect to its bilinear or sesquilinear form)}
  return IsIntegral(Scale(L));
end intrinsic;

intrinsic IsIntegral(L::LatMod, p::RngOrdIdl) -> BoolElt
{"} //"
  require Order(p) cmpeq BaseRing(L) and IsPrime(p):
    "The second argument must be a prime ideal of the base ring of the lattice";
  return Valuation(Scale(L), p) ge 0;
end intrinsic;

intrinsic Dual(L::LatMod) -> LatMod
{The dual lattice of L}
  G, C, B:= GramMatrixOfBasis(L);
  if #B eq 0 then return L; end if;	// L==0
  GI:= G^-1;
//  BB:= [ &+[ GI[i,j] * B[j] : j in [1..#B] ] : i in [1..#B] ];
  BB:= GI * Matrix(B);
  a:= Involution(L);
  return SetupLatticeInSpace(L, Module( [ < a(C[i])^-1, BB[i] > : i in [1..#C] ] ));
end intrinsic;

// TODO:: Allow integral non-prime ideals?
intrinsic PartialDual(L::LatMod, p::RngOrdIdl) -> LatMod
{The dual lattice of L}
  require IsFull(L) : "The lattice must have full rank";
  require Order(p) cmpeq BaseRing(L) and IsPrime(p):
    "The second argument must be a prime ideal of the base ring of the lattice";
  D:= Dual(L);
  m:=  Valuation(Scale(L), p);
  M:= -Valuation(Scale(D), p);
  X:= (D + p^-m * L) meet p^-M * L;
  d:= Volume(L);
  assert Volume(X) eq d * p^(-Valuation(d, p)+Valuation(Volume(D), p));
  return X;
end intrinsic;

intrinsic Gerstein(L::LatMod, p::RngOrdIdl) -> LatMod
{"The Gerstein reduction of L_p"}
  require IsFull(L) : "The lattice must have full rank";
  require Order(p) cmpeq BaseRing(L) and IsPrime(p):
    "The second argument must be a prime ideal of the base ring of the lattice";
  if IsHermitian(L) and IsSplit(p) then
    return L + (p^-1*L meet Dual(L));
  end if;
  return L + (p^-1*L meet p*Dual(L));
end intrinsic;

intrinsic DirectSum(L1::LatMod, L2::LatMod) -> LatMod
{The direct (orthogonal) of the lattices L1 and L2}
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible base rings";
  require IsOrthogonal(L1) eq IsOrthogonal(L2) : "Incompatible lattices";
  LL:= New(Type(L1));
  LL`Module:= DirectSum(Module(L1), Module(L2));
  LL`Form:= DiagonalJoin(L1`Form, L2`Form);
  return LL;
end intrinsic;

declare attributes RngOrd: APE;
intrinsic AbsolutePrimitiveElement(R::RngOrd) -> RngOrdElt
{Returns an primitive element of R over the integers}
  if not assigned R`APE then
    R`APE:= R ! PrimitiveElement(AbsoluteOrder(R));
  end if;
  return R`APE;
end intrinsic;

function Forms(L : Generators := [])
  empty:= IsEmpty(Generators);
  if not assigned L`ZForms or empty then
    R:= BaseRing(L);
    if empty then Generators:= [1, AbsolutePrimitiveElement(R)]; end if;
    Form:= InnerProductMatrix(L);
    B:= Matrix(AbsoluteBasis(L`Module));
    B:= Matrix(BaseRing(Form), B);
    G:= B * Form * InternalConjugate(L, B);
    F:= [];
    d:= [];
    for x in Generators do
      X:= Matrix( Ncols(G), [ AbsoluteTrace(e) : e in Eltseq( x*G ) ] );
      X, den:= IntegralMatrix(X);
      c:= Content(X);
      Append(~F, X div c);
      Append(~d, c/den);
    end for;
    if not empty then return F, d, B; end if;
    G:= LLLGram(F[1]: Eta:= 0.501, Delta:= 0.999);
    quality:= Max(D)/Min(D) where D:= Diagonal(G);
    L`ZForms:= < F, d, B, Generators[2], quality >;
  end if;
  return Explode(L`ZForms);
end function;

intrinsic ElementsOfNorm(L::LatMod, a::RngElt) -> []
{Returns all elements of norm a in L up to sign}
  ok, a:= IsCoercible(FixedField(L), a);
  require ok : "The second argument must be in the field of fractions of the base ring";
  if IsZero(L) then return [ L | ]; end if;
  ok, s:= IsDefinite(L);
  require ok : "The lattice must be definite";

  R:= BaseRing(L);
  FF:= FieldOfFractions(R);
  if FF!a notin Norm(L) then return [ L | ];
  elif IsZero(a) then return [ L ! 0 ];
  elif not IsTotallyPositive(s*a) then return [ L | ]; end if;

  Z:= ZBasis(L); d:= #Z; n:= AbsoluteDegree(FF);
  B:= Matrix(FF, [ Eltseq(x) : x in Z ]);
  // We solve xGx^t = 1 where
  G:= a^-1 * GramMatrix(L, B);
  X:= LatticeWithGram( Matrix(Ncols(G), [ AbsoluteTrace(x) : x in Eltseq(G) ]) );

  // The easy case.
  if ideal< R | a > eq Norm(L) then
    S:= Minimum(X) eq n select Matrix(FF, ShortestVectorsMatrix(X)) else Matrix(FF, 0, 0, []);
  else
    S:= Matrix(FF, ShortVectorsMatrix(X, n, n));
    SG:= S * G;
    S:= RowSubmatrix(S, [ i: i in [1..Nrows(S)] | IsOne(InnerProduct(SG[i], S[i])) ]);
  end if;
  return Nrows(S) eq 0 select [ L | ] else ChangeUniverse(Rows( S * B ), L);
end intrinsic;


/* Automorphism group computation and isometry test via Nebe's lattice decomposition */

function IsDecomposableZ(L)
  if not assigned L`AutoIsomData then
    F:= Forms(L)[1];
    LL:= LatticeWithGram(F);

    n:= Ncols(F);
    B:= Matrix(Integers(), 0, n, []);
    Ranks:= [];
    SVCount:= [];
    Sat:= [* *];
    denom:= [];

    while Nrows(B) lt n do
      SV:= ShortestVectorsMatrix(LL); assert BaseRing(SV) eq Integers();
      Append(~SVCount, Nrows(SV));
      I:= BasisMatrix(Image(SV));
      Append(~Ranks, Nrows(I));
      B:= VerticalJoin(B, I);
      T:= Solution(Matrix(Rationals(), B), Matrix(Rationals(), Saturation(B)));
      d:= Denominator(T);
      Append(~denom, d);
      Append(~Sat, Matrix(Integers(), d*T));
      LL:= sub< LL | KernelMatrix(F * Transpose(B)) >;
    end while;
    assert Rank(B) eq n;
    if #SVCount eq 1 and Abs(Determinant(B)) eq 1 then 
      L`AutoIsomData:= false;
    else
      L`AutoIsomData:= <MatrixRing(Integers(), n) ! B, Ranks, SVCount, denom, Sat>;
    end if;
  end if;
  return L`AutoIsomData cmpne false;
end function;

function HasAutoData(L)
  if not assigned L`AutoIsomData or L`AutoIsomData cmpeq false then return false; end if;
  return #L`AutoIsomData eq 6; 
end function;

function GetAutoIsomData(L, Auto)
  assert assigned L`AutoIsomData and L`AutoIsomData cmpne false;
  if Auto and #L`AutoIsomData eq 5 then
    F:= Forms(L);
    B:= L`AutoIsomData[1];
    Ranks:= L`AutoIsomData[2];
    A:= [* *]; 
    n:= 1;
    for i in [1..#Ranks] do
      T:= RowSubmatrix(B, n, Ranks[i]);
      AA:= MyAuto( [ T * f * Transpose(T) : f in F] );
      Append(~A, < [ Matrix(g) : g in Generators(redgen(AA)) ], #AA > );	// Is the reduction of the number of generators worth the trouble?
      n +:= Ranks[i];
    end for;
    Append(~L`AutoIsomData, A);
  end if;
  return Explode(L`AutoIsomData);
end function;

procedure Swap(~X, ~Y)
  Z:= X; X:= Y; Y:= Z;
end procedure;

function AutoIsomDec(L1, L2)
  Iso:= L2 cmpne 0;	// Do we want an isometry or the automorphism group?
  assert Iso or not assigned L1`Auto;

  F1, d1, B1, a1:= Forms(L1);
  if Iso then
    F2, d2, B2, a2:= Forms(L2);
    error if a1 ne a2 , "Boo-boo: The primitive elements differ!";
    if d1 ne d2 then return true, false, _; end if;
  end if;
  dec:= IsDecomposableZ(L1);
  if Iso and (dec ne IsDecomposableZ(L2)) then return true, false, _; end if;
  if not dec then return false, _, _; end if;

  swap:= Iso and not HasAutoData(L1) and HasAutoData(L2);
  if swap then
    Swap(~L1, ~L2);
    Swap(~F1, ~F2);
    Swap(~B1, ~B2);
  end if;

  BB1, Ranks1, SV1, denom1, Sat1, AA:= GetAutoIsomData(L1, true );
  if Iso then
    BB2, Ranks2, SV2, denom2, Sat2  := GetAutoIsomData(L2, false);
    if Ranks1 ne Ranks2 or SV1 ne SV2 or denom1 ne denom2 then return true, false, _; end if;
  else
    Sat2:= Sat1;
    F2:= F1;
    BB2:= BB1;
  end if; 

  if Iso then
    // Isometries of small blocks:
    Isoms:= [**]; n:= 1;
    for i in [1..#Ranks1] do
      T1:= RowSubmatrix(BB1, n, Ranks1[i]);
      T2:= RowSubmatrix(BB2, n, Ranks1[i]);
      ok, T:= IsIsometric([ T1 * f * Transpose(T1) : f in F1], [ T2 * f * Transpose(T2) : f in F2]);
      if not ok then return true, false, _; end if;
      Append(~Isoms, T);
      n +:= Ranks1[i];
    end for;
  else
    Isoms:= [* MatrixRing(Integers(), r) ! 1: r in Ranks1 *];
  end if;

  for i in [1..#Isoms] do
    if i eq 1 then
      isomi:= Isoms[1];
      n:= Ranks1[1];
      G:= MatrixGroup< n, Integers() | AA[1,1] >;
      G`Order:= AA[1,2];
    else
      One1:= MatrixRing(Integers(), Ranks1[i]) ! 1;
      One2:= MatrixRing(Integers(), n) ! 1;
      Gens:= [DiagonalJoin(G.j, One1): j in [1..Ngens(G)]] cat
             [DiagonalJoin(One2, a) : a in AA[i,1] ];
      n +:= Ranks1[i];
      size:= #G * AA[i,2];
      G := MatrixGroup< n, Integers() | Gens >;
      G`Order:= size;
      isomi:= DiagonalJoin(isomi, Isoms[i]);
    end if;

    // Find stabilizer of the saturation
    if denom1[i] ne 1 then
      isommerk:= Matrix(Integers(), Sat2[i]*isomi);
      H:= ChangeBase(G,[RowSpace(Sat1[i])]);
      G:= BasicStabilizer(H,2);
      ok, trans:= IsInBasicOrbit(H,1,RowSpace(isommerk));
      if not ok then 
        assert Iso;
	//"sat failed in iter.", i;
	return true, false, _;
      end if;
      isomi:= isomi*trans;
    end if;

    // Find stabilizer of second form
    if i ne 1 then
      T1:= RowSubmatrix(BB1, n);
      T2:= RowSubmatrix(BB2, n);
      Gens:= [ G.i : i in [1..Ngens(G)]];
      GensT:= [ Transpose(g): g in Gens ];
      S:= sub< G | >;
      Orbit:= {@ T1 * F1[2] * Transpose(T1) @};
      Paths:= [ G ! 1 ];
      j:= 1;
      while j le #Orbit do
        for l in [1..#Gens] do
          f:= Gens[l] * Orbit[j] * GensT[l];
          idx:= Index(Orbit, f);
          if idx eq 0 then
            Include(~Orbit, f);
            Append(~Paths, Gens[l] * Paths[j]);
          else
            s:= Paths[idx]^-1 * Gens[l] * Paths[j];
            if s notin S then S:= sub< G | S, s >; end if;
          end if;
        end for;
        j +:= 1;
      end while;
      idx:= Index(Orbit, TT*F2[2]*Transpose(TT)) where TT:= Matrix(Integers(), isomi^-1 * T2);
      if idx eq 0 then 
        assert Iso;
      //  "2nd form failed";
        return true, false, _;
      end if;
      isomi:= isomi * Matrix(Paths[idx]);
      G:= redgen(S);
    end if;
  end for;

  FF:= FieldOfFractions(BaseRing(L1));

  // We can now set the automorphism group of L1:
  if not assigned L1`Auto then
    GG:= [ Inv * Matrix(G.i) * BB1 : i in [1..Ngens(G)] ] where Inv:= Matrix(Rationals(), BB1)^-1;
    assert forall{g: g in GG | [ Matrix(g) * f * Transpose(Matrix(g)) : f in F1 ] eq F1};

    BM:= BasisMatrix(L1);
    S1:= Solution(B1, BM);
    S2:= Solution(BM, B1);
    L1`Auto:= MatrixGroup< Rank(L1), FF | [ S1 * Matrix(FF, g) * S2 : g in GG ] >;
   // assert #L1`Auto eq #G;
    L1`Auto`Order:= #G;
    if not Iso then return true, _, _; end if;
  end if;

  T:= Matrix(Rationals(), BB2)^-1 * isomi * Matrix(Rationals(), BB1);
  assert [ T * f * Transpose(T) : f in F1 ] eq F2;

  S1:= Solution(B2, BasisMatrix(L2));
  S2:= Solution(BasisMatrix(L1), B1);
  S:= S1 * Matrix(FF, T) * S2;

  if swap then S:= S^-1; end if;  
  return true, true, S;
end function;

intrinsic OrthogonalDecomposition(L::LatMod) -> []
{Decomposes L into an orthogonal sum of indecomposables}
  ok, a:= IsDefinite(L);
  require ok : "The lattice must be definite";

  if not IsOne(a) then L:= Rescale(L, a); end if;
  F, _, B:= Forms(L : Generators:= AbsoluteBasis(BaseRing(L)));
  T:= OrthogonalDecomposition(F);

  if #T eq 1 then return [ L ]; end if;
  
  R:= BaseRing(B);
  return [ sub< L |  Rows(Matrix(R, t) * B) > : t in T ];
end intrinsic;

/* The functions to compute automorphism groups and isometries */

function SomePerp(L)
  if IsFull(L) then return Matrix( FieldOfFractions(BaseRing(L)), 0, Degree(L), []); end if;
  return KernelMatrix( InnerProductMatrix(L) * InternalConjugate(L, BasisMatrix(L)) );
end function;

intrinsic AutomorphismGroup(L::LatMod : NaturalAction:= true, Decomposition:= 0, Check:= false) -> GrpMat
{The automorphism group of the lattice L}
// Automorphism groups are internally stored with respect to BasisMatrix(L).
  if not assigned L`Auto then
    ok, s:= IsDefinite(L);
    require ok : "The lattice must be definite";
    LL:= IsOne(s) select L else Rescale(L, s);
    FF:= FieldOfFractions(BaseRing(L));

    if Decomposition cmpeq 0 then 
      if AlwaysTryDecomposition then
        Decomposition:= true; q:= -1;
      else
        _, _, _, _, q:= Forms(LL);
        Decomposition:= q ge 5;
      end if;
    end if;
    if Decomposition then
      vprint LatMod: "Decomposition attempted", q;
      vtime LatMod: Decomposition:= AutoIsomDec(LL, 0);
    end if;
    if Decomposition then
      vprint LatMod: "Decomposition succeeded";
      assert assigned LL`Auto;
      L`Auto:= LL`Auto;
// check
/*      AA:= LL`Auto;
      delete LL`Auto;
      A:= AutomorphismGroup(LL: Decomposition:= false, NaturalAction:= false);
      assert A eq AA;*/
    else
      vprint LatMod: "Will use no decomposition", q;
      F, _, B:= Forms(LL);
//      vprint LatMod: F;
//      vtime LatMod: A:= AutomorphismGroup(F);
      _, T:= LLLGram(F[1]: Delta:= 0.999, Eta:= 0.501);
      G:= [ T*f*Transpose(T) : f in F];
      vprint LatMod: G;
      vtime LatMod: A:= MyAuto(G);
      A:= A^(Generic(A)!T);
      assert forall{<a,f> : a in Generators(A), f in F | a * f * Transpose(a) eq f };

      BM:= BasisMatrix(L);
      S1:= Solution(B, BM);
      S2:= Solution(BM, B);
      L`Auto:= MatrixGroup< Rank(L), FF | [ S1 * Matrix(FF, A.i) * S2 : i in [1..Ngens(A)] ] >;
      if Check then assert #L`Auto eq #A; else L`Auto`Order:= #A; end if;
    end if;

    // the final check -- very expensive!
    if Check then 
      G, C, B:= GramMatrixOfBasis(L);
      M:= Module(L);
      for a in Generators(L`Auto) do 
        assert G eq Matrix(a) * G * InternalConjugate(L, a);
        assert M eq Module( [ < C[i], X[i] > : i in [1..#C] ] ) where X:= Matrix(a) * Matrix(B);
      end for;
    end if;
  end if;

  // Conjugate the group to the natural action, if desired.
  if NaturalAction then
    X:= BasisMatrix(L);
    FF:= BaseRing(X);
    n:= Nrows(X); m:= Ncols(X);
    if n eq m then
      A:= (L`Auto)^(GL(n, FF) ! X);
    else
      X:= VerticalJoin(X, SomePerp(L)); I:= X^-1;
      A:= sub< GL(m, BaseRing(X)) | [ I * DiagonalJoin(Matrix(L`Auto.i), IdentityMatrix(FF, m-n)) * X : i in [1..Ngens(L`Auto)] ] >;
    end if;
    if Check then
      assert #A eq #L`Auto;
      assert forall{ i: i in [1..Ngens(A)] | Matrix(A.i) * IP * InternalConjugate(L, A.i) eq IP } where IP:= InnerProductMatrix(L);
    end if;
    return A;
  end if;
  return L`Auto;
end intrinsic;

function ExtendBasis(M)
  m:= Ncols(M);
  if m eq Nrows(M) then return M; end if;
  V:= VectorSpace(BaseRing(M), m);
  D:= {1..m} diff { Depth(E[i]): i in [1..Nrows(M)] } where E:= EchelonForm(M);
  return VerticalJoin(M, Matrix( [ V.i: i in D ] ));
end function;

intrinsic IsIsometric(L1::LatMod, L2::LatMod : NaturalAction:= true, Check:= false, Decomposition:= 0) -> BoolElt, Mtrx
{Tests if the lattices L1 and L2 are isometric}
  require BaseRing(L1) cmpeq BaseRing(L2) : "The lattices must have identical base rings";
  require IsOrthogonal(L1) eq IsOrthogonal(L2) : "The lattices must both be orthogonal or both be hermitian";
  ok, s:= IsDefinite(L1);
  require ok : "Only implemented for definite lattices";
  if s ne 1 then
    L1:= Rescale(L1, s);
    L2:= Rescale(L2, s);
  end if;
  if not IsPositiveDefinite(L2) then return false, _; end if;

  if (assigned L1`Auto and assigned L2`Auto and #L1`Auto ne #L2`Auto) or Rank(L1) ne Rank(L2) or
       Scale(L1) ne Scale(L2) or Norm(L1) ne Norm(L2) or Discriminant(L1) ne Discriminant(L2) then 
    return false, _;
  end if;

  if Decomposition cmpeq 0 then 
    if AlwaysTryDecomposition then
      Decomposition:= true;
    else
      _, _, _, _, q1:= Forms(L1);
      _, _, _, _, q2:= Forms(L2);
      Decomposition:= Max(q1, q2) ge 5;
    end if;
    // "Trying decomposition? ", Decomposition;
  end if;

  if Decomposition then
    Decomposition, result, S:= AutoIsomDec(L1, L2);
    //if Decomposition and result then assert result eq IsIsometric(L1, L2 : Decomposition:= false); end if;
    if Decomposition and not result then return false, _; end if;
  end if;
  if not Decomposition then
    F1, d1, B1, a1:= Forms(L1);
    F2, d2, B2, a2:= Forms(L2);
    error if a1 ne a2 , "Boo-boo: The primitive elements differ!";
    if d1 ne d2 then return false, _; end if;

    G1, T1:= LLLGram(F1[1]: Delta:= 0.999, Eta:= 0.501);
    G2, T2:= LLLGram(F2[1]: Delta:= 0.999, Eta:= 0.501);
    H1:= T1*F1[2]*Transpose(T1);
    H2:= T2*F2[2]*Transpose(T2);

    if G1 eq G2 and H1 eq H2 then
      result:= true; T:= Parent(G1) ! 1;
    else
      result, T:= IsIsometric( [G1, H1], [G2, H2] );
      if not result then return false, _; end if;
    end if;
    T:= T2^-1 * T * T1;

    S1:= Solution(B2, BasisMatrix(L2));
    S2:= Solution(BasisMatrix(L1), B1);
    S:= S1 * Matrix(FieldOfFractions(BaseRing(L1)), T) * S2;
  end if;

  if Check then
    G1, C1, B1:= GramMatrixOfBasis(L1);
    G2, C2, B2:= GramMatrixOfBasis(L2);
    assert S * G1 * InternalConjugate(L1, S) eq G2;
    assert Module(L1) eq Module( [ < C2[i], X[i]> : i in [1..#C2] ] ) where X:= S * Matrix(B1);
  end if;

  // In general there is no chance that any x notin F*L1 is mapped to something useful.
  // So we simply extend BasisMatrix(Li) to bases of the ambient. 
  if NaturalAction then
    B1:= BasisMatrix(L1);
    B2:= BasisMatrix(L2);
    E1:= ExtendBasis(B1);
    E2:= ExtendBasis(B2);
    S:= E2^-1 * DiagonalJoin(S, IdentityMatrix(BaseRing(E1), Nrows(E1) - Nrows(B1))) * E1;

    if Check then
      G:= Generators(L2);
      I:= [ Vector(g) * S : g in G ];
      assert Module(L1) eq Module(I);
      assert GramMatrix(L2, G) eq GramMatrix(L1, I);
    end if;
  end if;

  return true, S;
end intrinsic;

intrinsic IsSimilar(L1::LatMod, L2::LatMod : Decomposition:= 0, NaturalAction:= true) -> BoolElt, Mtrx, FldAlgElt
{Tests if the lattices L1 and L2 are similar}
  require BaseRing(L1) cmpeq BaseRing(L2) : "The lattices must have identical base rings";
  m:= Rank(L1);
  if m ne Rank(L2) then return false, _, _; end if;
  orth:= IsOrthogonal(L1);
  require orth eq IsOrthogonal(L2) : "The lattices must both be orthogonal or both be hermitian";
  ok, s:= IsDefinite(L1);
  require ok : "Only implemented for definite lattices";
  if s ne 1 then
    L1:= Rescale(L1, s);
    L2:= Rescale(L2, s);
  end if;
  if not IsPositiveDefinite(L2) then return false, _, _; end if;

  I:= Scale(L1) / Scale(L2);
  if I ne Norm(L1) / Norm(L2) then return false, _, _; end if;
  if orth then
    ok, aa:= HasTotallyPositiveGenerator(I);
  else
    F:= FixedField(L1);
    if Type(F) eq FldRat then 
      aa:= [ Minimum(I) ];
      if I ne ideal< BaseRing(L1) | aa[1] > then return false, _, _; end if; 
      ok:= true;
    else
      J:= I;
      I:= I meet Integers(F);
      if J ne ideal< BaseRing(L1) | I > then return false, _, _; end if;
      ok, aa:= HasTotallyPositiveGenerator(I);
    end if;
  end if;
  if not ok then return false, _, _; end if;

  for a in aa do
    M:= Rescale(L2, a);
    assert Scale(L1) eq Scale(M) and Norm(L1) eq Norm(M);
    ok, x:= IsIsometric(L1, M : Decomposition:= Decomposition, NaturalAction:= NaturalAction);
    if ok then return true, x, a; end if;
  end for;
  return false, _, _;
end intrinsic;

intrinsic IsSimilarGenus(L1::LatMod, L2::LatMod) -> BoolElt
{Tests if the genera of L1 and L2 are similar}
  S:= BaseRing(L1);
  require S cmpeq BaseRing(L2) : "The lattices must have identical base rings";
  m:= Rank(L1);
  require m eq Rank(L2) : "The lattices must have the same rank";
  orth:= IsOrthogonal(L1);
  require orth eq IsOrthogonal(L2) : "The lattices must both be orthogonal or both be hermitian";

  K:= FixedField(L1);
  oo:= RealPlaces(K);
  I:= Scale(L1) / Scale(L2);
  if I ne Norm(L1) / Norm(L2) or I^m ne Volume(L1)/Volume(L2) then return false; end if;
  if orth then
    ok, x:= IsPrincipal(I);
    if not ok then return false; end if;
    D1:= MyDiagonal(L1, false);
    D2:= MyDiagonal(L2, false);
    I1:= [ #[ d: d in D1 | Evaluate(d, inf) lt 0 ] : inf in oo ];
    I2:= [ #[ d: d in D2 | Evaluate(d, inf) lt 0 ] : inf in oo ];
    Signs:= [];
    for i in [1..#I1] do 
      if I1[i] eq I2[i] then s:= 1;
      elif I1[i] eq m-I2[i] then s:= -1;
      else return false; 
      end if;
      Signs[i] := s * Sign(Evaluate(x, oo[i]));
    end for;
  else
    if Type(K) eq FldRat then
      x:= Minimum(I);
    else
      d:= Denominator(I);
      ok, x:= IsPrincipal((d*I) meet Integers(K));
      if not ok then return false; end if;
      x:= x/d;
    end if;
    if I ne ideal< S | x > then return false; end if;
    d1:= Determinant(L1`Form);
    d2:= Determinant(L2`Form);
    E:= NumberField(S);
    oo:= [ inf: inf in oo | #Decomposition(E, inf) eq 1 ];
    Signs:= [ Sign(Evaluate(tmp, inf)) : inf in oo ] where tmp:= K ! (d1*d2*x);
  end if;
  U:= UnitsWithSigns(K, oo, Signs);
  for u in U do
    M:= Rescale(L2, x*u);
    assert Scale(L1) eq Scale(M) and Norm(L1) eq Norm(M);
    if IsSameGenus(L1, M) then return true; end if;
  end for;
  return false;
end intrinsic;

// Norm with respect to x :-> (x,x)
intrinsic Norm(L::LatMod) -> RngOrdFracIdl
{The norm of the lattice L}
  assert IsOrthogonal(L);	// Hermitian forms have their own implementation
  if not assigned L`Norm then
    if IsZero(L) then return 0*BaseRing(L); end if;
    G, C:= GramMatrixOfBasis(L);
    L`Norm:= &+{ G[i,i] * C[i]^2 : i in [1..#C] } + 2*Scale(L);
  end if;
  return L`Norm;
end intrinsic;

intrinsic Scale(L::LatMod) -> RngOrdFracIdl
{The scale of the lattice L}
  if not assigned L`Scale then
    if IsZero(L) then return 0*BaseRing(L); end if;
    G, C:= GramMatrixOfBasis(L);
    a:= Involution(L);
    // sum over all i,j for hermitian lattices...
    Gens:= { G[i,j] * C[i] * a(C[j]) : i in [1..j], j in [1..#C] | G[i,j] ne 0 };
    if IsHermitian(L) then Gens join:= { a(g): g in Gens }; end if;
    L`Scale:= &+Gens;
  end if;
  return L`Scale;
end intrinsic;

intrinsic Volume(L::LatMod) -> RngOrdFracIdl
{The volume of the lattice L}
  return Discriminant(L);
end intrinsic;

intrinsic Determinant(L::LatMod) -> FldNumElt
{The determinant of the free lattice L}
  G:= GramMatrix(L);
  require Ncols(G) eq Rank(L) : "The lattice is not free";
  return Determinant(G);
end intrinsic;

intrinsic IsModular(L::LatMod) -> BoolElt, RngOrdFracIdl
{If L is modular, returns true and an ideal a such that a*Dual(L)=L}
  if IsZero(L) then return true, FieldOfFractions(BaseRing(L)) ! 1; end if;
  a:= Scale(L);
  if Volume(L) eq a^Rank(L) then
    return true, a;
  end if;
  return false, _;
end intrinsic;

intrinsic IsModular(L::LatMod, p::RngOrdIdl) -> BoolElt, RngIntElt
{If L_p is p^v-modular, return true and v}
  require Order(p) cmpeq BaseRing(L) and IsPrime(p): "The second argument must be a prime ideal of the base ring of the lattice";
  if IsZero(L) then return true, 0; end if;
  v:= Valuation(Scale(L), p);
  if v*Rank(L) eq Valuation(Volume(L), p) then
    return true, v;
  end if;
  return false, _;
end intrinsic;

intrinsic IsUnimodular(L::LatMod, p::RngOrdIdl) -> BoolElt
{Returns true iff L_p is unimodular}
  ok, v:= IsModular(L, p);
  return ok and v eq 0;
end intrinsic;

intrinsic IsUnimodular(L::LatMod) -> BoolElt
{Returns true iff L is unimodular}
  return forall{p : p in BadPrimes(L) | IsUnimodular(L, p)};
end intrinsic;

intrinsic IsSquarefree(L::LatMod, p::RngOrdIdl) -> BoolElt
{Returns true iff L_p is squarefree}
  require Order(p) cmpeq BaseRing(L) and IsPrime(p): "The second argument must be a prime ideal of the base ring of the lattice";
  if IsZero(L) then return true; end if;
  if Minimum(p) eq 2 then
    _, scales:= GenusSymbol(L,p);
    scales:= Set(scales);
  else
    G:= GenusSymbol(L,p);
    scales:= { g[2]: g in G };
  end if;
  return scales subset {0,1};
end intrinsic;

intrinsic IsSquarefree(L::LatMod) -> BoolElt
{Returns true iff L is squarefree}
  if IsZero(L) then return true; end if;
  for p in BadPrimes(L) do
    if not IsSquarefree(L, p) then return false; end if;
  end for;
  return true;
end intrinsic;

// Isometry testing over the field of fractions
Hasse:= func< D, p | &*[ Integers() | HilbertSymbol(K ! D[i], K ! (&*D[ [i+1..n] ]), p) : i in [1..n-1] ] where n:= #D where K:= NumberField(Order(p)) >;

intrinsic HasseInvariant(L::LatMod, p::RngOrdIdl : AmbientSpace:= false) -> RngIntElt
{The Hasse invariant of L at p}
  require Order(p) cmpeq BaseRing(L) and IsPrime(p):
    "The second argument must be a prime ideal of the base ring of the lattice";
  require IsOrthogonal(L): "The lattice must be quadratic";
  return Hasse(MyDiagonal(L, AmbientSpace), p);
end intrinsic;

intrinsic HasseInvariant(L::LatMod, p::PlcNumElt : AmbientSpace:= false) -> RngIntElt
{"} //"
  require IsOrthogonal(L): "The lattice must be quadratic";
  if IsFinite(p) then return HasseInvariant(L, Ideal(p)); end if;
  require NumberField(L) cmpeq NumberField(p): "Incompatible arguments";
  if IsComplex(p) then return 1; end if;
  n:= #[ d: d in MyDiagonal(L, AmbientSpace) | Evaluate(d, p) lt 0 ];
  return (n mod 4) ge 2 select -1 else 1;
end intrinsic;

intrinsic WittInvariant(L::LatMod, p::RngOrdIdl : AmbientSpace:= false) -> RngIntElt
{The Witt invariant of L at p}
  h:= HasseInvariant(L, p: AmbientSpace:= AmbientSpace);
  F:= AmbientSpace select InnerProductMatrix(L) else GramMatrixOfBasis(L);
  Det:= Determinant(F);
  K:= NumberField(L);
  c:= K ! case < Ncols(F) mod 8 | 3: -Det, 4: -Det, 5: -1, 6: -1, 7: Det, 0: Det, default : 1 >;
  return h * HilbertSymbol(K ! -1, c, p);
end intrinsic;

intrinsic WittInvariant(L::LatMod, p::PlcNumElt : AmbientSpace:= false) -> RngIntElt
{The Witt invariant of L at p}
  if IsFinite(p) then
    return WittInvariant(L, Ideal(p));
  elif IsComplex(p) then 
    return false;
  end if;
  h:= HasseInvariant(L, p : AmbientSpace:= AmbientSpace);
  F:= AmbientSpace select InnerProductMatrix(L) else GramMatrixOfBasis(L);
  Det:= Determinant(F);
  c:= NumberField(L) ! case < Ncols(F) mod 8 | 3: -Det, 4: -Det, 5: -1, 6: -1, 7: Det, 0: Det, default : 1 >;
  return Evaluate(c, p) lt 0 select -h else h;
end intrinsic;

intrinsic IsRationallyEquivalent(L1::LatMod, L2::LatMod, p::RngOrdIdl : AmbientSpace:= false) -> BoolElt
{Tests if L1 and L2 are equivalent over the completion at p}
  require IsOrthogonal(L1) and IsOrthogonal(L2) : "The lattices must both be hermitian or orthogonal.";
  // note: hermitian lattcies have their own implementation
  if AmbientSpace or (IsFull(L1) and IsFull(L2)) then
    F1:= InnerProductMatrix(L1);
    F2:= InnerProductMatrix(L2);
  else
    F1:= GramMatrixOfBasis(L1);
    F2:= GramMatrixOfBasis(L2);
  end if;
  // In most cases, we have full lattices in the same space anyway ...
  if F1 cmpeq F2 then return true; end if;

  R:= BaseRing(L1);
  require R cmpeq BaseRing(L2): "Incompatible arguments";
  require Order(p) cmpeq R and IsPrime(p):
    "The third argument must be a prime ideal of the base ring of the lattices";
  return Ncols(F1) eq Ncols(F2) and IsLocalSquare(Determinant(F1) * Determinant(F2), p) and 
         HasseInvariant(L1, p: AmbientSpace:= AmbientSpace) eq HasseInvariant(L2, p : AmbientSpace:= AmbientSpace);
end intrinsic;

intrinsic IsRationallyEquivalent(L1::LatMod, L2::LatMod, p::PlcNumElt : AmbientSpace:= false) -> BoolElt
{"} //"
  if IsFinite(p) then return IsRationallyEquivalent(L1, L2, Ideal(p): AmbientSpace:= AmbientSpace); end if;
  require BaseRing(L1) cmpeq BaseRing(L2) and NumberField(L1) cmpeq NumberField(p): "Incompatible arguments";
  D1:= MyDiagonal(L1, AmbientSpace);
  D2:= MyDiagonal(L2, AmbientSpace);
  return (#D1 eq #D2) and (IsComplex(p) or (#[ d: d in D1 | Evaluate(d, p) lt 0 ] eq #[ d: d in D2 | Evaluate(d, p) lt 0 ]));
end intrinsic;

intrinsic IsRationallyEquivalent(L1::LatMod, L2::LatMod : AmbientSpace:= false) -> BoolElt
{Tests if L1 and L2 are equivalent over their base field}
  require IsOrthogonal(L1) and IsOrthogonal(L2) : "The lattices must both be hermitian or orthogonal.";
  require BaseRing(L1) cmpeq BaseRing(L2): "Incompatible lattices";

  if AmbientSpace or (IsFull(L1) and IsFull(L2)) then
    F1:= InnerProductMatrix(L1);
    F2:= InnerProductMatrix(L2);
  else
    F1:= GramMatrixOfBasis(L1);
    F2:= GramMatrixOfBasis(L2);
  end if;
  if F1 cmpeq F2 then return true;
  elif Ncols(F1) ne Ncols(F2) then return false; end if;

  Det1, Finite1, I1:= QuadraticFormInvariants(F1);
  Det2, Finite2, I2:= QuadraticFormInvariants(F2);
  return I1 eq I2 and Finite1 eq Finite2 and IsSquare(Det1*Det2);
end intrinsic;


// The genera stuff

intrinsic BadPrimes(L::LatMod : Even:= false) -> {}
{The set of prime ideals at which L is not unimodular}
  Bad:= Support(Scale(L)) join Support(Volume(L));
  if Even then Bad join:= Support( 2*BaseRing(L) ); end if;
  return Bad;
end intrinsic;

intrinsic GenusSymbol(L::LatMod, p::RngOrdIdl : Uniformizer:= 0) -> []
{The genus symbol of L at p}

  // Use the cache and fix the odd-dimensional cases if needed
  if not assigned L`GenusSymbols then
    L`GenusSymbols:= AssociativeArray();
  end if;
  ok, sym:= IsDefined(L`GenusSymbols, p);
  if ok then
    if Minimum(p) eq 2 then
      return Explode(sym);
    elif Uniformizer cmpeq 0 then
      x, Uniformizer:= Explode(sym);
    else
      _, h:= ResidueClassField(p);
      x:= IsSquare( (Uniformizer / sym[2]) @ h ) select sym[1] else 
	[ < e[1], e[2], IsOdd(e[1]*e[2]) select -e[3] else e[3] > : e in sym[1] ];
    end if;
    return x, Uniformizer;
  end if;

  R:= BaseRing(L);
  require Order(p) cmpeq R and IsPrime(p) : "Incompatible arguments";
  J, G, E:= JordanDecomposition(L, p);
  F:= FieldOfFractions(R);
  if Uniformizer cmpeq 0 then
    Uniformizer:= F ! PrimitiveElement(p);
  else
    ok, Uniformizer := IsCoercible( F, Uniformizer );
    require ok and Valuation(Uniformizer, p) eq 1: "Wrong uniformizer";
  end if;
  if Minimum(p) ne 2 then
    _, h:= ResidueClassField(p);
    Gs:= [ h(&*Diagonal(G[i]) / Uniformizer^(E[i] * Nrows(J[i]))) :  i in [1..#J] ];
    assert 0 notin Gs;
    x:= [ < Nrows(J[i]), E[i], IsSquare(Gs[i]) select 1 else -1 >: i in [1..#J] ];
    L`GenusSymbols[p]:= <x, Uniformizer>;
    return x, Uniformizer;
  else
    t:= #G;
    sL := [Minimum({Valuation(g[i,j], p): i in [1..j], j in [1..Ncols(g)]}): g in G];
    assert sL eq E;
    e := RamificationIndex(p);

    aL:= []; uL:= []; wL:= [];
    for i in [1..t] do
      //TODO: Use Decomposition of L better, remove GG
      GG:= DiagonalJoin(< j lt i select Uniformizer^(2*(sL[i]-sL[j])) * G[j] else G[j] : j in [1..t] >);
      D:= Diagonal(GG);
      min, pos:= Minimum([ Valuation(d, p) : d in D ]);
      if e+sL[i] le min then
        Append(~aL, Uniformizer^(e+sL[i]));
      else
        Append(~aL, D[pos]);
      end if;
      Append(~uL, Valuation(aL[i], p));
      Append(~wL, Minimum({e+sL[i]} join {uL[i] + QuadraticDefect(d/aL[i], p): d in D}));
    end for;

    fL := [];
    for k in [1..t-1] do
      exp:= uL[k] + uL[k+1];
      Append(~fL, (IsOdd(exp) select exp else Minimum({QuadraticDefect(aL[k]*aL[k+1], p), uL[k]+wL[k+1], uL[k+1]+wL[k], e+(exp div 2) + sL[k]})) - 2*sL[k] );
    end for;
    L`GenusSymbols[p]:= < [Nrows(G[k]): k in [1..#G]], sL, wL, aL, fL, G >;
    return [Nrows(G[k]): k in [1..t]], sL, wL, aL, fL, G;
  end if;
end intrinsic;

my_little_helper := function(G1, G2, p)
  G1o := Diagonal(OrthogonalizeGram(G1));
  G2o := Diagonal(OrthogonalizeGram(G2));
  Append(~G1o, &*G1o * &*G2o); // G1o and G2o now have the same length
  return Hasse(G1o, p) eq Hasse(G2o, p);
end function;

intrinsic IsLocallyIsometric(L1::LatMod, L2::LatMod, p::RngOrdIdl: CheckIsometricSpaces) -> BoolElt
{Returns true if and only if the p-adic completions of L1 and L2 are isometric}
  R:= BaseRing(L1);
  require R cmpeq BaseRing(L2) : "Incompatible base rings";
  require R cmpeq Order(p) and IsPrime(p):
    "The third argument must be a prime ideal in the base ring of the lattices";

  d:= Dimension(L1);
  if d ne Dimension(L2) then return false;
  elif d eq 0 then return true; 
  end if;

  if Minimum(p) ne 2 then
    S, pi:= GenusSymbol(L1, p);
    return S eq GenusSymbol(L2, p : Uniformizer:= pi);
  end if;

  // Test the embedding spaces first
  if CheckIsometricSpaces and not IsRationallyEquivalent(L1, L2, p) then return false; end if;

  dimL1, sL1, wL1, aL1, fL1, G1 := GenusSymbol(L1, p);
  dimL2, sL2, wL2, aL2, fL2, G2 := GenusSymbol(L2, p);

  if (dimL1 ne dimL2) or (sL1 ne sL2) or (wL1 ne wL2) then return false; end if;
  uL1 := [Valuation(aL1[k], p): k in [1..#aL1]];
  uL2 := [Valuation(aL2[k], p): k in [1..#aL2]];
  if uL1 ne uL2 then return false; end if;

  bL := [aL1[k]/aL2[k]: k in [1..#aL1]];
  qL := [QuadraticDefect(bL[k], p): k in [1..#bL]];
 
  if &or[qL[k] lt wL1[k]-uL2[k]: k in [1..#qL]] then return false; end if;
 
  e := RamificationIndex(p);
  // G1, G2: Gram matrix
  d1 := [* DiagonalJoin(<G1[i]: i in [1..t]>): t in [1..#G1] *]; 
  d2 := [* DiagonalJoin(<G2[i]: i in [1..t]>): t in [1..#G2] *]; 

  for i in [1..#fL1] do   
    detquot := Determinant(d1[i])/Determinant(d2[i]);
    if Valuation(detquot, p) ne 0 then return false; end if;
    if QuadraticDefect(detquot, p) lt fL1[i] then return false; end if;
    if (fL1[i] gt 2*e + uL1[i+1] - wL2[i+1]) and 
      not my_little_helper(d1[i], DiagonalJoin(d2[i], MatrixRing(BaseRing(d2[i]), 1)![aL2[i+1]]), p) then return false;
    end if;
    if (fL1[i] gt 2*e + uL1[i  ] - wL2[i  ]) and
      not my_little_helper(d1[i], DiagonalJoin(d2[i], MatrixRing(BaseRing(d2[i]), 1)![aL2[i]  ]), p) then return false;
    end if;
  end for;

  return true;
end intrinsic;

intrinsic IsSameGenus(L1::LatMod, L2::LatMod) -> BoolElt
{Returns true if and only if the lattices L1 and L2 are in the same genus}
  require BaseRing(L1) cmpeq BaseRing(L2) : "Incompatible base rings";
  require IsHermitian(L1) eq IsHermitian(L2) : "Both lattices must be either hermitian or quadratic";
  if Rank(L1) ne Rank(L2) or Norm(L1) ne Norm(L2) or Discriminant(L1) ne Discriminant(L2) or Scale(L1) ne Scale(L2) or not IsRationallyEquivalent(L1, L2) then return false; end if;
  Bad:= IsHermitian(L1) select BadPrimes(L1) join { p: p in MySupport( Discriminant(BaseRing(L1)) ) | Minimum(p) eq 2 } 
                          else BadPrimes(L1 : Even);
  return forall{ p: p in Bad | IsLocallyIsometric(L1, L2, p) };
end intrinsic;

intrinsic JordanDecomposition(L::LatMod, p::RngOrdIdl) -> [], [], []
{A Jordan decomposition of the completion of L at p}
  require BaseRing(L) cmpeq Order(p): "The arguments are incompatible";
  require IsPrime(p) : "The second argument must be prime";

  F:= InnerProductMatrix(L);
  even:= Minimum(p) eq 2;
  if even then pi:= PrimitiveElement(p); end if;

  oldval:= Infinity();
  Blocks:= []; Exponents:= [];
  S:= Matrix( LocalBasis(L, p) );
  n:= Ncols(S);
  k:= 1;
  while k le n do
    G:= S*F*Transpose(S);

    // find an element <i, j> with minimal p-valuation.
    // Take it from the diagonal, if possible, and from the lowest-possible row number.
    X:= [ Valuation(G[i,i], p): i in [k..n] ];
    m, ii:= Minimum( X ); ii +:= k-1;
    pair:= <ii, ii>;

    for i in [k..n], j in [i+1..n] do
      tmp:= Valuation(G[i,j], p);
      if tmp lt m then m:= tmp; pair:= <i,j>; end if;
    end for;
    if m ne oldval then Append(~Blocks, k); oldval:= m; Append(~Exponents, m); end if;

    if even and pair[1] ne pair[2] then
      SwapRows(~S, pair[1],   k); // swap f_1 and e_i
      SwapRows(~S, pair[2], k+1); // swap f_2 and e_j

      T12:= (S[k] * F * Matrix(1, Eltseq(S[k+1])))[1];
      S[k] *:= pi^Valuation(T12, p)/T12;
      T := func<i,j|(S[i] * F * Matrix(1, Eltseq(S[j])))[1]>;
      T11 := T(k,k); T22 := T(k+1, k+1); T12:= T(k, k+1);

      d := T11*T22 - T12^2;
      for l in [k+2..n] do
        tl := T12*T(k+1,l)-T22*T(k  ,l); // t_k from step 4
        ul := T12*T(k  ,l)-T11*T(k+1,l); // u_k from step 4
        S[l] +:= (tl/d)*S[k] + (ul/d)*S[k+1];
      end for;
      k +:= 2;
    else
      if pair[1] eq pair[2] then
        SwapRows(~S, pair[1], k);
      else
        S[pair[1]] +:= S[pair[2]];
        SwapRows(~S, pair[1], k);
      end if;
      nrm:= (S[k] * F * Matrix(1, Eltseq(S[k])))[1];
      X:= S[k] * F * Transpose(S); // T(e_k, f_i), 1<=i<=n
      for l in [k+1..n] do S[l] -:= X[l]/nrm * S[k]; end for;
      k+:= 1;
    end if;
  end while;

  Append(~Blocks, n+1);
  Matrices:= [* RowSubmatrix(S, Blocks[i], Blocks[i+1] - Blocks[i]) : i in [1..#Blocks-1] *];
  return Matrices, [* m*F*Transpose(m): m in Matrices *], Exponents;
end intrinsic;

/*
function Reflection(x, F)
  V:= VectorSpace( BaseRing(F), Ncols(F) );
  norm:= (x*F, x);
  X:= [ b - 2*(b*F, x)/norm * x : b in Basis(V) ]; 
  return Matrix(X);
end function;

intrinsic ApproximateByIsometry(L1::LatMod, L2::LatMod, S::Setq : Precision:= 10) -> ModMatElt
{Returns a (proper) global isometry that maps L1_p to L2_p for all p in S}
  S:= Setseq(Set(S));
  require forall{p: p in S | IsLocallyIsometric(L1, L2, p)} : "The lattices must be locally isometric";
  F:= InnerProductMatrix(L2);
  if InnerProductMatrix(L1) ne F then
    ok, M:= AreEquivalentQuadraticForms( L1`Form, L2`Form : Isometry );
    require ok : "The lattices must be in isometric spaces for this to work";
    L:= SetupLatticeInSpace( L2, [ Vector(x) * M : x in Generators(L1) ] );
  else
    L:= L1;
    M:= Parent(F) ! 1;
  end if;

  p:= Rep(S);

  require Minimum(p) ne 2 : "Implement me";

  Precision:= Max(Precision, 1+2*Valuation(Order(p) ! 2, p));

  G, B:= JordanDecomposition(L, p);
  G:= DiagonalJoin(<g: g in G>);
  B:= VerticalJoin(<b: b in B>);
  D:= Diagonal(G);
  B2:= LocalBasis(L2, p);
  k, h:= ResidueClassField(p);
  Vectors:= [];
  for i in [1..#D] do
    g:= D[i];
    v:= Valuation(g, p);
    V:= VectorSpace(k, #B2);
    repeat
      repeat
        coords:= Random(V);
        x:= &+[ (coords[j] @@ h) * B2[j] : j in [1..#B2] ];
        norm:= (x*F, x);
      until Valuation(norm, p) eq v;
      ok, s:= IsSquare( (g/norm) @ h );
    until ok;
    x *:= (s @@ h);
    repeat
      norm:= (x*F, x);
      w:= Valuation(g/norm - 1, p);
      if w ge Precision then break; end if;
      x *:= 1 + (g/norm-1) / 2;
    until false;
    if B[i] notin {x, -x} then
      w:= (B[i]+x)/2;
      if (w*F, w) eq 0 then w:= (B[i]-x)/2; end if;
      Append(~Vectors, w);
    end if;
    B2:= [ B2[j] - (x*F, B2[j])/norm * x : j in [1..#B2] | j ne d] where d:= Depth(coords);
  end for;

  T:= M;
  for v in Vectors do
    T:= Reflection(v, F) * T;
  end for;
  LL:= SetupLatticeInSpace( L2, Module([ Vector(x) * T : x in Generators(L1) ]) );
  E:= ElementaryDivisors( Module(LL), Module(L2) );
  assert forall{ p: p in S | Valuation(E[1], p) eq 0 and Valuation(E[#E], p) eq 0};

  return T;
end intrinsic;
*/

// Magma's internal ProjectiveSpace point enumeration is way too slow.
// So we use our own.
function MyPSReps(k, n)
  assert n ge 0;
  V:= VectorSpace(k, n+1);
  v:= V.(n+1);
  L:= [ v ];
  S:= [ a*v : a in k ];
  size:= 1;
  for i in [1..n] do
    m:= #L;
    for j in [size..m] do
      w:= Rotate(L[j], -1);
      L cat:= [ w+s: s in S ];
    end for;
    size:= m+1;
  end for;
  return L;
end function;

intrinsic IsIsotropic(L::LatMod, p::RngOrdIdl : AmbientSpace:= false) -> BoolElt
{Tests if L is isotropic at p}
  require BaseRing(L) cmpeq Order(p) and IsPrime(p) : "The second argument must be a prime ideal of the base ring of the lattice";
 
  F:= AmbientSpace select InnerProductMatrix(L) else GramMatrixOfBasis(L);
  n:= Ncols(F);
  d:= Determinant(F);
  K:= NumberField(Order(p));
  if d eq 0 then return true;
  elif n le 1 then return false;
  elif n eq 2 then return IsLocalSquare(-d, p);
  elif n eq 3 then return HasseInvariant(F, p) eq HilbertSymbol(K ! -1, K! -1, p);
  elif n eq 4 then return not IsLocalSquare(d, p) or HasseInvariant(F, p) eq HilbertSymbol(K ! -1, K! -1, p);
  else return true;
  end if;
end intrinsic;

intrinsic IsIsotropic(L::LatMod, p::PlcNumElt : AmbientSpace:= false) -> BoolElt
{"} //"
  require FixedField(L) cmpeq NumberField(p) : "Incompatible arguments";
  if IsFinite(p) then return IsIsotropic(L, Ideal(p)); end if;

  F:= AmbientSpace select InnerProductMatrix(L) else GramMatrixOfBasis(L);
  n:= Ncols(F);
  d:= Determinant(F);
  K:= BaseRing(L);
  if d eq 0 then return true;
  elif n le 1 then return false;
  elif IsComplex(p) then return true;
  else
    D:= IsHermitian(L) select MyGramSchmidt(F, Involution(L)) else MyDiagonal(L, AmbientSpace);
    return #{ Sign( Evaluate(d, p) ) : d in D } eq 2;
  end if;
end intrinsic;

/*
intrinsic IsIsotropic(L::LatMod : AmbientSpace:= false) -> BoolElt
{Tests if L is isotropic}
  call IsIsotropic on a form ... 
end intrinsic;
*/

intrinsic Neighbours(L::LatMod, p::RngOrdIdl : AutoOrbits:= false, CallBack:= false, Max:= Infinity()) -> []
{The p-neighbours of L}
  R:= BaseRing(L);
  require R cmpeq Order(p): "The arguments are incompatible";
  require IsPrime(p): "The second argument must be a prime ideal";
  ok, rescale:= IsModular(L, p);
  require ok: "The lattice must be locally modular";
  require IsIsotropic(L, p) : "The lattice must be locally isotropic";
  assert IsOrthogonal(L); // hermitian lattices are different
  e:= Valuation(R ! 2, p);
  require e eq 0 or Valuation(Norm(L), p) ge e : "The lattice must be even";

  UseCallBack:= CallBack cmpne false;
  if Max eq 1 and not UseCallBack then AutoOrbits:= false; end if;

  B:= LocalBasis(L`Module, p : Type:= "Submodule");
  n:= #B;
  k, h:= ResidueClassField(p);
  pi:= PrimitiveElement(p);
  piinv:= WeakApproximation([p],[-1]);
  F:= FieldOfFractions(R);
  M:= GL(n, F) ! Matrix(B);
  Form:= GramMatrix(L, M);
  if rescale ne 0 then
    Form:= piinv^rescale * Form;
  end if;
  U:= VectorSpace(F, n, Form);
  pForm:= Matrix(k, n, [ x @ h : x in Eltseq(Form) ]);
  W:= VectorSpace(k, n, pForm);

  if AutoOrbits cmpne false then
    if AutoOrbits cmpeq true then
      A:= AutomorphismGroup(L : NaturalAction:= false);
      BM:= BasisMatrix(L);
      S1:= Solution(BM, M);
      S2:= Solution(M, BM);
    else
      require IsFull(L) : "Setting the automorphism group is not supported for lattices which are not full";
      A:= ChangeRing(AutoOrbits, F);
      S1:= M;
      S2:= M^-1;
    end if;
    A:= [Matrix(k, n, [ x @ h : x in Eltseq(S1 * A.i * S2) ]) : i in [1..Ngens(A)]];
    A:= [ a: a in A | not IsScalar(a) ];
    AutoOrbits:= #A ge 1;
  end if;
  if AutoOrbits then
    A:= sub< GL(n, k) | A >;
//    "start";
    LO:= [ l[2].1: l in MyOrbitsOfSpaces(A, 1) ];
//    "done";
  else
    LO:= ProjectiveLineProcess(W);
  end if;

  Result:= [];
  pM:= p*Module(L);
  keep:= true; cont:= true; found:= false;
//#LO, "orbits";
  for i in [1..#LO] do
//  if i mod 100000 eq 0 then i; end if;
    w:= AutoOrbits select W ! LO[i] else Next(LO);
    if (w,w) ne 0 then continue; end if;
    if IsZero(w*pForm) then continue; end if;	// for "bad primes"
    x:= U ! [e @@ h: e in Eltseq(w)];
    nrm:= (x,x);
    val:= Valuation(nrm, p);
    assert val gt 0;
    if val le e then continue; 	// can only happen if p even
    elif val eq e+1 then
      // make val > e+1
      ok:= exists(r){r : r in [1..n] | (w, W.r) ne 0}; assert ok;
      a:= R ! ((nrm / (2*pi * (x,U.r))) @ h @@ h);
      x := x - a * pi * U.r;
      assert Valuation( (x,x), p ) ge e+2;
    end if;
    found:= true;

    // normalize x
    ok:= exists(kk){ i : i in [1..n] | x[i] @ h ne 0 }; assert ok;
    x *:=  R ! (((x[kk] @ h)^-1) @@ h);
    if e ne 0 then
      x*:= 1 + ((((x[kk]-1)/pi) @ h) @@ h) * pi;
      assert Valuation(x[kk] - 1, p) ge 2;
    end if;

    xF:= x*Form;
    ok:= exists(mm){ r : r in [1..n] | r ne kk and Valuation(xF[r], p) eq 0 }; assert ok;

    VV:= Matrix([ x * piinv, 0 * pi*U.mm ] cat [ U.i - ((xF[i]/xF[mm]) @ h @@ h)*U.mm : i in [1..n] | i ne kk and i ne mm ]);
    V:= VV * M;
    LL:= SetupLatticeInSpace(L, pM + Module(Rows(V)));

    if UseCallBack then
      keep, cont:= CallBack(Result, LL);
    end if;
    if keep then Append(~Result, LL); end if;
    if not cont or #Result ge Max then break; end if;
  end for;
  if not found then "Warning: L_p/pL_p has no suitable isotropic vector!"; end if;

  // The lattices should all be in the same genus
  //assert forall{LL: LL in Result | IsLocallyIsometric(L, LL, p)};
  if (e eq 0 or Valuation(Norm(L), p) ge e) and assigned L`GenusSymbols then
    for i in [1..#Result] do Result[i]`GenusSymbols:= L`GenusSymbols; end for;
  end if;

  return Result;
end intrinsic;

function MaximalSubspaces(k, n)
  I:= MatrixRing(k, n) ! 1;
  L:= [];
  for i in [1..n] do
    X:= RemoveRow(I, i);
    V:= VectorSpace(k, i-1);
    for v in V do Append(~L, InsertBlock(X, Matrix(1, Eltseq(v)), 1, i)); end for;
  end for;
  return L;
end function;

intrinsic MaximalSublattices(L::LatMod, p::RngOrdIdl : AutoOrbits:= false, CallBack:= false, Max:= Infinity()) -> []
{The maximal sublattices of L that contain pL}
  require BaseRing(L) cmpeq Order(p): "The arguments are incompatible";
  require IsPrime(p): "The second argument must be a prime ideal";

  B:= LocalBasis(L`Module, p : Type:= "Submodule"); n:= #B;
  B:= Matrix(B); 
  F:= FieldOfFractions(BaseRing(L));
  k, h:= ResidueClassField(p); hinv:= h^-1;

  if AutoOrbits cmpne false then
    if AutoOrbits cmpeq true then
      A:= AutomorphismGroup(L : NaturalAction:= false);
      BM:= BasisMatrix(L);
      S1:= Solution(BM, B);
      S2:= Solution(B, BM);
    else
      require IsFull(L) : "Setting the automorphism group is not supported for lattices which are not full";
      A:= ChangeRing(AutoOrbits, F);
      S1:= B;
      S2:= B^-1;
    end if;
    AT:= [Transpose(Matrix(k, n, [ x @ h : x in Eltseq(S1 * A.i * S2) ])) : i in [1..Ngens(A)]];
    AT:= [ a: a in AT | not IsScalar(a) ];
    AutoOrbits:= #AT ge 1;
  end if;
  if AutoOrbits then
    AT:= sub< GL(n,k) | AT >;
    vprint LatMod: "calling orbits";
    vtime LatMod: Ls:= [ < o[2].1, o[1] > : o in MyOrbitsOfSpaces(AT, 1) ];
  else
    Ls:= MaximalSubspaces(k, n);
  end if;

  pML:= p*Module(L);
  Result:= []; keep:= true; cont:= true; E:= [];
  for i in [1..#Ls] do
    if AutoOrbits then
      m:= ChangeRing(KernelMatrix( Matrix(1, Eltseq(Ls[i,1])) ), hinv);
    else
      m:= ChangeRing(Ls[i], hinv);
    end if;
    m:= ChangeRing(m, F);
    LL:= SetupLatticeInSpace(L, Module(Rows(m * B)) + pML);
    if CallBack cmpne false then 
      keep, cont:= CallBack(Result, LL);
    end if;
    if keep then 
      Append(~Result, LL);
      Append(~E, AutoOrbits select Ls[i,2] else 1);
    end if;
    if not cont then break; end if;
    if #Result ge Max then break; end if;
  end for;

  return Result, E;
end intrinsic;

intrinsic Sublattices(L::LatMod, p::RngOrdIdl, i::RngIntElt : AutoOrbits:= false, CallBack:= false) -> [], []
{The sublattices of index p^i in L that contain pL}
  require BaseRing(L) cmpeq Order(p): "The arguments are incompatible";
  require IsPrime(p): "The second argument must be a prime ideal";

  n:= Rank(L);
  requirerange i, 0, n;
  case i:
    when   0: return [L], [1];
    when   n: return [p*L], [1];
    when   1: return MaximalSublattices(L, p : AutoOrbits:= AutoOrbits, CallBack:= CallBack);
    when n-1: return MinimalSuperlattices(p*L, p : AutoOrbits:= AutoOrbits, CallBack:= CallBack);
  end case;

  B:= LocalBasis(L`Module, p : Type:= "Submodule"); assert n eq #B;
  B:= Matrix(B);
  F:= FieldOfFractions(BaseRing(L));
  k, h:= ResidueClassField(p); hinv:= h^-1;

  if AutoOrbits cmpeq false then
    A:= sub< GL(n,k) | >;
  else
    if AutoOrbits cmpeq true then
      A:= AutomorphismGroup(L : NaturalAction:= false);
      BM:= BasisMatrix(L);
      S1:= Solution(BM, B);
      S2:= Solution(B, BM);
    else
      require IsFull(L) : "Setting the automorphism group is not supported for lattices which are not full";
      A:= ChangeRing(AutoOrbits, F);
      S1:= B;
      S2:= B^-1;
    end if;
    A:= sub< GL(n, k) | [ Matrix(k, n, [ x @ h : x in Eltseq(S1 * A.i * S2) ]) : i in [1..Ngens(A)]] >;
    A:= redgen(A);
  end if;
  Ls:= MyOrbitsOfSpaces(A, n-i);

  Gens:= Generators(p*L`Module);
  Result:= []; keep:= true; cont:= true;
  Lengths:= [];
  for s in Ls do
    m:= ChangeRing(BasisMatrix(s[2]), hinv);
    m:= ChangeRing(m, F);
    LL:= SetupLatticeInSpace(L, Module(Rows(m * B) cat Gens));
    if CallBack cmpne false then
      keep, cont:= CallBack(Result, LL);
    end if;
    if keep then Append(~Result, LL); Append(~Lengths, s[1]); end if;
    if not cont then break; end if;
  end for;

  return Result, Lengths;
end intrinsic;

intrinsic MinimalSuperlattices(L::LatMod, p::RngOrdIdl : AutoOrbits:= false, CallBack:= false, Max:= Infinity()) -> []
{The minimal superlattices of L contained in p^-1 * L}
  require BaseRing(L) cmpeq Order(p): "The arguments are incompatible";
  require IsPrime(p): "The second argument must be a prime ideal";

  B:= LocalBasis(L`Module, p : Type:= "Submodule"); n:= #B;
  B:= Matrix(B);
  F:= FieldOfFractions(BaseRing(L));
  k, h:= ResidueClassField(p); hinv:= h^-1;

  if AutoOrbits cmpne false then
    if AutoOrbits cmpeq true then
      A:= AutomorphismGroup(L : NaturalAction:= false);
      BM:= BasisMatrix(L);
      S1:= Solution(BM, B);
      S2:= Solution(B, BM);
    else
      require IsFull(L) : "Setting the automorphism group is not supported for lattices which are not full";
      A:= ChangeRing(AutoOrbits, F);
      S1:= B;
      S2:= B^-1;
    end if;
    AT:= [ Matrix(k, n, [ x @ h : x in Eltseq(S1 * A.i * S2) ]) : i in [1..Ngens(A)]];
    AT:= [ a: a in AT | not IsScalar(a) ];
    AutoOrbits:= #AT ge 1;
  end if;
  if AutoOrbits then
    AT:= sub< GL(n,k) | AT >;
    Ls:= MyOrbitsOfSpaces(AT, 1);
  else
    Ls:= ProjectiveLineProcess(k, n);
  end if;

  pinv:= p^-1;
  ML:= Module(L);
  Result:= []; keep:= true; cont:= true; E:= [];
  for i in [1..#Ls] do
    l:= AutoOrbits select Ls[i,2].1 else Next(Ls);
    m:= ChangeRing(ChangeRing(l, hinv), F);
    LL:= SetupLatticeInSpace(L, ML + Module([< pinv, m*B >]));
    if CallBack cmpne false then
      keep, cont:= CallBack(Result, LL);
    end if;
    if keep then
      Append(~Result, LL);
      Append(~E, AutoOrbits select Ls[i,1] else 1);
    end if;
    if not cont then break; end if;
    if #Result ge Max then break; end if;
  end for;

  return Result, E;
end intrinsic;


// Characteristic \ne 2
function IsIsotropicFinite(M)
  n:= Ncols(M);
  k:= BaseRing(M);
  assert IsField(k) and IsFinite(k) and Characteristic(k) ne 2;
  if n eq 0 then 
    ;
  elif n eq 1 then 
    if IsZero(M[1,1]) then return true, Vector(k, [1]); end if;
  else
    if n le 3 then
      G, T:= OrthogonalizeGram(M);
    else	// solution can be found within the first three vars.
      G, T:= OrthogonalizeGram(SubmatrixRange(M, 1,1,3,3));
      B:= RMatrixSpace(k, 3, n) ! 0;
      B[1,1]:= 1; B[2,2]:= 1; B[3,3]:= 1;
      T:= T * B;
    end if;
    if exists(i){i: i in [1..Ncols(G)] | G[i,i] eq 0} then return true, T[i];
    elif n eq 2 then
      ok, s := IsSquare( -G[1,1] / G[2,2] );
      if ok then return true, T[1] + s*T[2]; end if;
    else
      repeat
        x:= Random(k); y:= Random(k);
        ok, z:= IsSquare( (-x^2*G[1,1] - y^2*G[2,2])/G[3,3] );
       until ok and (x ne 0 or y ne 0);
       return true, x*T[1] + y*T[2] + z*T[3];
    end if;
  end if;
  
 return false, _;
end function;

function GuessMaxDet(L, p)
  m:= Rank(L); n:= m div 2;
  d:= Determinant(GramMatrixOfBasis(L));
  e:= 2*Valuation(Order(p) ! 2, p);
  if IsOdd(m) then
    v:= Valuation(d, p) mod 2;
    v:= WittInvariant(L, p) eq 1 select v-e*n else 2-v-e*n;
  else
    if IsOdd( (m*(m-1)) div 2 ) then d := -d; end if;
    qd:= QuadraticDefect(d, p);
    if Type(qd) eq Infty then
      v:= WittInvariant(L, p) eq 1 select -e*n else 2-e*n;
    else		// Wrong? Fix scaling \alpha
      vd:= Valuation(d, p);
      v:= vd - 2*(qd div 2) + e*(1-n);
/*K:= NumberField(BaseRing(L));
F:= ext< K | Polynomial([-d,0,1]) >;
ram:= IsRamified(p, Integers(F));
assert  (IsEven(vd) and qd eq vd+e) eq not ram;*/
      if IsEven(vd) and qd eq vd+e and WittInvariant(L,p) eq -1 then v:= -e*n+2; end if;
    end if;
  end if;
  return v;
end function;

intrinsic IsMaximalIntegral(L::LatMod, p::RngOrdIdl) -> BoolElt, LatMod
{Checks whether L is p-maximal integral. If not, a minimal integral over-lattice at p is returned}
  require Order(p) cmpeq BaseRing(L) and IsPrime(p): 
    "The second argument must be a prime ideal of the base ring of the lattice";
  if IsZero(L) then return true, L; end if;
  if Valuation(Norm(L), p) lt 0 then return false, _; end if;
  if GuessMaxDet(L, p) eq Valuation(Volume(L), p) then return true, L; end if;
 
  k, h:= ResidueClassField(p);
  BM:= Matrix(LocalBasis(Module(L), p: Type:= "Submodule"));
  FF:= BaseRing(BM);

  G:= 2*GramMatrix(L, BM);
  V:= KernelMatrix(Matrix(Nrows(BM), [x @ h: x in Eltseq(G)]));
  assert Nrows(V) ne 0;

  if Minimum(p) ne 2 then
    T:= ChangeRing(ChangeRing(V, h^-1), FF);
    H:= T * G * Transpose(T) / PrimitiveElement(p);
    ok, v:= IsIsotropicFinite(Matrix(Ncols(H), [x @ h : x in Eltseq(H)] ));
    assert ok;
    e:= Eltseq( v*V ) @@ h;
    v:= Vector(FF, e);
    valv:= Valuation( ((v*G) * Matrix(FF,1,e))[1], p );
    assert valv ge 2;
  else
    val2:= Valuation(BaseRing(L)!2, p);
    PP:= ProjectiveLineProcess(k, Nrows(V));
    repeat
      x:= Next(PP);
      assert not IsZero(x);
      e:= Eltseq(x * V) @@ h;
      v:= Vector(FF, e);
      valv:= Valuation( ((v*G) * Matrix(FF,1,e))[1], p );
      assert valv ge 1;
    until valv ge val2 + 2;
  end if;
  LL:= SetupLatticeInSpace(L, Module(L) + Module( [WeakApproximation([p], [-1]) *v*BM ] ) );
  assert Volume(L) eq Volume(LL) * p^2 and Valuation(Norm(LL), p) ge 0;
  return false, LL;
end intrinsic;

intrinsic IsMaximalIntegral(L::LatMod) -> BoolElt, LatMod
{Checks whether L is maximal integral. If not, a minimal integral over-lattice is returned}
  if IsZero(L) then return true, L; end if;
  if not IsIntegral(Norm(L)) then return false, _; end if;
  for p in BadPrimes(L : Even) do
    ok, LL:= IsMaximalIntegral(L, p);
    if not ok then return false, LL; end if;
  end for;
  return true, L;
end intrinsic;

intrinsic IsMaximal(L::LatMod, p::RngOrdIdl) -> BoolElt, LatMod
{Checks if L_p is Norm(L_p)-maximal}
  require Order(p) cmpeq BaseRing(L) and IsPrime(p):
    "The second argument must be a prime ideal of the base ring of the lattice";
  if IsZero(L) then return true, L; end if;
  v:= Valuation(Norm(L), p);
  x:= PrimitiveElement(p)^-v;
  ok, LL:= IsMaximalIntegral( Rescale(L, x), p);
  if ok then return true, L; end if;
  return false, Rescale(LL, 1/x);
end intrinsic;

intrinsic IsMaximal(L::LatMod) -> BoolElt, LatMod
{Checks if L is Norm(L)-maximal}
  if IsZero(L) then return true, L; end if;
  Bad:= IsHermitian(L) select BadPrimes(L : Discriminant) else BadPrimes(L : Even); 
  for p in Bad do
    ok, LL:= IsMaximal(L, p);
    if not ok then return false, LL; end if;
  end for;
  return true, L;
end intrinsic;

intrinsic MaximalIntegralLattice(L::LatMod, p::RngOrdIdl) -> LatMod
{A p-maximal integral lattice over L}
  require Order(p) eq BaseRing(L) and IsPrime(p): "The second argument must be a prime ideal of the base ring of L";
  require not IsZero(L) and Valuation(Norm(L), p) ge 0 : "The norm of the lattice must be locally integral";

  ok, LL:= IsMaximalIntegral(L, p);
  while not ok do 
    L:= LL;    
    ok, LL:= IsMaximalIntegral(L, p);
  end while;
  return L;
end intrinsic;

intrinsic MaximalIntegralLattice(L::LatMod) -> LatMod
{A maximal integral lattice containing L}
  require not IsZero(L) and IsIntegral(Norm(L)) : "The lattice must be integral and non-zero";

  for p in BadPrimes(L : Even) do
    L:= MaximalIntegralLattice(L, p);
  end for;
  return L;
end intrinsic;

intrinsic MaximalIntegralLattice(Q::Mtrx) -> LatMod
{A lattice which is maximal integral with respect to the quadratic form Q}
  require IsSymmetric(Q) : "The form must be symmetric";
  R:= BaseRing(Q); T:= Type(R);
  if ISA(T, RngOrd) then
    F:= FieldOfFractions(R);
    Q:= Matrix(F, Q);
  elif ISA(T, FldNum) then
    R:= Integers(R);
    F:= FieldOfFractions(R);
    Q:= Matrix(F, Q);
  else
    require ISA(T, FldOrd) : "The matrix must be over (the field of fractions of) an order";
    F:= R;
    R:= Integers(R);
  end if;
  n:= Nrows(Q);
  require Rank(Q) eq n : "The form must be non-degenerate";

  // We start with some integral lattice.
  L:= LatticeModule( Q );
  N:= Norm(L);
  if N ne 1*R then
    FN:= Factorization(N);
    d:= &*[ f[1]^(f[2] div 2) : f in Factorization(N) ];
    L:= d^-1*L;
    N:= Norm(L);
    assert IsIntegral(N);
  end if;

  return MaximalIntegralLattice(L);
end intrinsic;

// group actions
intrinsic '^'(L::LatMod, g::GrpMatElt) -> LatMod
{internal}
  require Degree(L) eq Degree(g) : "Incompatible arguments";
  F:= FieldOfFractions(BaseRing(L));
  ok, g:= CanChangeRing(g, F);
  require ok : "Incompatible arguments";

   P:= PseudoBasis(Module(L));
  return SetupLatticeInSpace(L, Module( [ < p[1], Vector(F, p[2]) * g > : p in P ] ) );
end intrinsic;


// The mass stuff.

// returns q^-dim(G) * #G(F_q)
function GroupOrder(G, m, q)
  if G[1] ne "O" then
    o:= &* case< G[1] |
          "G": [ 1-   q^-j : j in [1..m] ],             // GL_m
          "U": [ 1-(-q)^-j : j in [1..m] ],             //  U_m
      default: [ 1-   q^-j : j in [2..m by 2] ] >;      // Sp_m
  else // the orthogonal cases:
    if IsEven(m) then
      k:= m div 2;
      e:= G[2] eq "+" select 1 else -1;
      o:= 2 * (1-e * q^-k) * &* [ Rationals() | 1-q^-j: j in [2..(m-2) by 2] ];
    else
      o:= 2 * &* [ Rationals() | 1-q^-j: j in [2..m by 2] ];
    end if;
  end if;
  return o;
end function;

function WittToHasse(Dim, Det, Finite)
  K:= Parent(Det);
  c:= K ! case < Dim mod 8 | 3: -Det, 4: -Det, 5: -1, 6: -1, 7: Det, 0: Det, default : 1 >;
  return { x[1] : x in Finite | x[2] ne HilbertSymbol(K ! -1, c, x[1]) };
end function;

// Local factor of a maximal lattice in K*L following Shimura / Gan-Hanke-Yu.
function lf_max(L, p)
  m:= Rank(L); r:= m div 2;
  if m eq 1 then return 1; end if;
  w:= WittInvariant(L, p);
  d:= get_disc(L`Form);
  q:= Norm(p);
  if IsOdd(m) then
    if IsOdd(Valuation(d,p)) then return (q^r+w)/2;
    elif w eq -1 and IsEven(Valuation(d,p)) then return (q^(m-1)-1) / (q+1) / 2;
    end if;
  else
    if IsOdd(Valuation(d,p)) then return 1/2;           // ram.
    elif IsLocalSquare(d,p) then                        // split
      if w eq -1 then return (q^(r-1)-1)*(q^r-1)/(q+1)/2; end if;
    elif IsOdd(QuadraticDefect(d,p)) then return 1/2;   // ram.
    else                                                // inert
      if w eq -1 then return (q^(r-1)+1)*(q^r+1)/(q+1)/2; end if;
    end if;
  end if;
  return 1;
end function;

// Local factor of a unimodular lattice at an even prime p
GetEps:= func< c | (Type(c) eq Infty) select 1 else -1 >;
function lf_uni(L, p)
  d, s, w, a, _, G := GenusSymbol(L, p);
  assert s eq [0] and Minimum(p) eq 2;
  d:= d[1]; b:= w[1]; a:= Valuation(a[1], p);
  e:= RamificationIndex(p);
  q:= Norm(p);
  if IsEven(d) then
    assert b eq e or IsOdd(a+b);
    r:= d div 2;
    disc:= (-1)^r * Determinant(G[1]);
    c:= QuadraticDefect(disc, p);
    if d eq 2 then
      if a lt b and b eq e and c ge 2*e then
        lf:= q^( (e-a-1) div 2 ) * (q - GetEps(c));
      elif b eq e and a+e+1 le c and c lt 2*e then
        lf:= 2*q^((c-e-a) div 2);
      else
        lf:= 1;
      end if;
    elif IsEven(a+b) then
      if e eq a then
        lf:= 1;
      elif c ge 2*e then
        lf:= q^(Integers() ! ((e-a)*(r-1/2) - r)) * (q^r - GetEps(c));
      elif a+e+1 le c then
        lf:= 2*q^(Integers() ! ((c-e-a-1)*(r-1/2)));
      else
        lf:= q^((c-e-a-1)*(r-1));
      end if;
    else  // a+b odd
      if c eq a+b then
        lf:= 1;
      elif c ge 2*e then
        // We first compute the Hilbert symbol hs of (alpha, 1+gamma).
        hs:= (Type(c) eq Infty) or IsEven(a) select 1 else -1;
        // Compute c'
        if WittInvariant(L, p) eq hs then
          cc:= c;
        else
          cc:= Type(c) eq Infty select 2*e else Infinity();
        end if;

        if e eq b and Type(cc) eq Infty then
          ee:= GetEps(c);
          lf:= q^(Integers() ! ((e-a-1)*(r-1/2))) * (q^r-ee) * (q^(r-1)+ee) / 2;
        elif e eq b then
          lf:= q^(Integers() ! ((e-a-1)*(r-1/2))) * (q+1);
        elif Type(cc) eq Infty and c eq cc then
          ee:= GetEps(c);
          lf:= q^(Integers() ! ((2*e-b-a-3)*(r-1/2)+r)) * (q^r-ee) * (q^(2*r-2)-ee) / 2;
        elif Type(cc) eq Infty then
          lf:= q^(Integers() ! ((2*e-b-a-3)*(r-1/2)+r)) * (q+1);
        elif c eq cc then
          ee:= GetEps(c);
          lf:= q^(Integers() ! ((2*e-b-a-3)*(r-1/2)+r)) * (q^(r-1) + 1) * (q^r-ee) * (q^(r-1)+ee) / 2;
        else
          lf:= q^(Integers() ! ((2*e-b-a-3)*(r-1/2)+r)) * (q^(r-1) + 1) * (q+1);
        end if;
      elif b eq e then
        lf:= 2 * q^(Integers() ! ((c-e-a)*(r-1/2)));
      else
        lf:= (q^(2*r-2) - 1) * q^(Integers() ! ((c-a-b-2)*(r-1/2)+1));
      end if;
    end if;
  else  // Dim odd.
    assert a eq 0;
    r:= (d-3) div 2;
    Witt:= WittInvariant(L, p);

    if e eq b then
      lf:= 1;
    elif IsOdd(e) then
      lf:=(q^(r+1)-Witt) * q^((r+1)*(e-b-1));
    else
      lf:= (Witt eq 1 select (q^(2*r+2)-1)/2 else q+1) * q^((r+1)*(e-b-1));
    end if;
  end if;
  return lf_max(Rescale(L, 2), p) * lf;
end function;

// Local factor when p|2 is unramified following S. Cho.
function lf_Cho(L, p)
  assert Minimum(p) eq 2 and not IsRamified(p);
  m:= Rank(L);
  _, G, S:= JordanDecomposition(L, p);
  k,h:= ResidueClassField(p);
  V:= [**];

  for s in S do
    AG:= DiagonalJoin(< 2^(S[j] lt s select 2*(s- S[j]) else 0) * G[j] : j in [1..#G] >);
    B:= KernelMatrix( Matrix(k, 1, [ h(d/2^s) : d in Diagonal(AG) ]) );
    B:= Matrix(k, Ncols(B), [ Sqrt(x): x in Eltseq(B) ]);
    BK:= ChangeRing(ChangeRing(B, h^-1), BaseRing(AG));
    Q:= BK * AG * Transpose(BK) / 2^(s+1);
    for i,j in [1..Ncols(Q)] do
      if i gt j then Q[i,j]:= 0;
      elif i lt j then Q[i,j] *:= 2;
      end if;
    end for;
    Q:= Matrix(k, Ncols(Q), [ x @ h: x in Eltseq(Q) ]);

    VV:= VectorSpace(k, Nrows(B));
    BB:= Basis(VV);
    assert Matrix(k, #BB, [ (b*Q, c) + (c*Q, b): b, c in BB ]) eq Q + Transpose(Q);

    N:= KernelMatrix(Q+Transpose(Q));
    ok, D:= IsDiagonal(N * Q * Transpose(N));
    assert ok;
    Rad:= KernelMatrix(Matrix(k, 1, D)) * N;

    W, g:= quo< VV | RowSpace(Rad) >;
    if Dimension(W) eq 0 then
      Append(~V, < s, Matrix(k, 0, []) >);
    else
      BB:= Matrix([ b @@ g: b in Basis(W) ]);
      Append(~V, < s, MatrixRing(k, Nrows(BB)) ! (BB * Q * Transpose(BB)) > );
    end if;
  end for;
  M:= [ Ncols(g): g in G ];

  PT:= [ Valuation(Norm(LatticeModule(G[i])), p) eq S[i] : i in [1..#S] ];              // parity type I
  Alpha:= [];
  for i in [1..#G] do
    ok:= exists(v){v: v in V | v[1] eq S[i] }; assert ok;
    if Ncols(v[2]) ne 0 and (S[i]-1 notin S or not PT[i-1]) and (S[i]+1 notin S or not PT[i+1]) then Append(~Alpha, i); end if;
  end for;
  Beta:= [];
  for i in [1..#G] do
    if not PT[i] then continue; end if;
    idx:= Index(S, S[i]+2);
    if idx eq 0 or not PT[idx] then Append(~Beta, i); end if;
  end for;
  rk:= - &+[ (m*(m-1)) div 2 where m:= Ncols(Q[2]): Q in V ] + (m*(m-1)) div 2;
  assert rk ge 0;               // rk of max. red. quotient  

  q:= Norm(p);
  res:= 2^(#Alpha + #Beta);
  N:= Multiplicity(PT, true) - #[ i: i in [1..#S-1] | PT[i] and PT[i+1] and S[i+1] eq S[i]+1 ] + &+[ Integers() | M[i] : i in [1..#S] | not PT[i] ];
  for i in [1..#S] do
    N +:= S[i] * (M[i]*(M[i]+1)) div 2;
    N +:= S[i] * M[i] * &+[ Integers() | M[j]: j in [i+1..#S] ];
  end for;
  for v in V do
    mi:= Ncols(v[2]);
    if mi eq 0 then continue;
    elif IsOdd(mi) then
      QFT:= "O";
    elif mi gt 0 then
      QFT:= QuadraticFormType(v[2]) eq "orthogonalplus" select "O+" else "O-";
    end if;
//    v[1], QFT, mi;
    res *:= GroupOrder(QFT, mi, q)/2;
  end for;

  beta:= 1/2 * q^N * res;

  exp:= (m+1)/2 * (&+[ S[i] * M[i] : i in [1..#S] ]);           // from det.
  if IsOdd(m) then
    exp +:= (m+1)/2;
    H:= GroupOrder("O", m, q);
  else
    exp +:= m;
    d:= get_disc(L`Form);
    K:= NumberField(Order(p));
    if IsLocalSquare(d,p) then
      H:= GroupOrder("O+", m, q);
    else
      E:= ext< K | Polynomial([-d,0,1]) >;
      v:= Valuation(Discriminant(Integers(E)), p);
      if v eq 0 then
        H:= GroupOrder("O-", m , q);
      else
        H:= GroupOrder("O", m-1, q);
        exp +:= v * (1-m)/2;
      end if;
    end if;
  end if;
  assert IsIntegral(exp);
  exp:= Integers() ! exp;

  return q^exp * H/2 / beta;
end function;

// General local factor driver
intrinsic LocalFactor(L::LatMod, p::RngOrdIdl) -> FldRatElt
{The local factor of L at the prime p in the Minkowsi-Siegel mass formula}

  require Order(p) cmpeq BaseRing(L) and IsPrime(p):
    "The second argument must be a prime ideal of the base ring of the lattice";

  if Rank(L) le 1 then return 1; end if;

  if Minimum(p) eq 2 then
    if not IsRamified(p) then
      return lf_Cho(L,p);
    elif IsMaximal(L, p) then
      s:= PrimitiveElement(p)^-Valuation(Norm(L), p);
      return lf_max(Rescale(L, s), p);
    elif IsModular(L, p) then
      s:= PrimitiveElement(p)^-Valuation(Scale(L), p);
      return lf_uni(Rescale(L, s), p);
    else
      def, a:= IsDefinite(L);
      s:= PrimitiveElement(p)^-Valuation(Norm(L), p);
      if def then
        R:= BaseRing(L);
        s*:= CRT(p, [1..AbsoluteDegree(R)], R ! 1, Signature(NumberField(R) ! (s*a)) );
      end if;
      L:= Rescale(L, s);
      Chain:= [L];
      ok, LL:= IsMaximalIntegral(L, p);
      while not ok do
        Append(~Chain, LL);
        ok, LL:= IsMaximalIntegral(LL, p);
      end while;
      f:= lf_max(L, p);
      for i in [1..#Chain-1] do
        M, E:= MaximalSublattices(Chain[i+1], p: AutoOrbits:= def);
        f *:= &+[ E[j]: j in [1..#M] | IsLocallyIsometric(Chain[i  ], M[j], p) ];
        M, E:= MinimalSuperlattices(Chain[i], p: AutoOrbits:= def);
        f /:= &+[ E[j]: j in [1..#M] | IsLocallyIsometric(Chain[i+1], M[j], p) ];
      end for;
      return f;
    end if;
  end if;

  // The odd primes
   _, G, s:= JordanDecomposition(L, p);
  if #s eq 1 then return 1; end if;

  m:= Rank(L);
  q:= Norm(p);
  if IsOdd(m) then
    f := GroupOrder("O+", m, q);
  else
    d:= get_disc(L`Form);
    if IsOdd(Valuation(d, p)) then
      f := GroupOrder("O+", m-1, q);
    elif IsLocalSquare(d, p) then
      f := GroupOrder("O+", m, q);
    else
      f := GroupOrder("O-", m, q);
    end if;
  end if;
  N:= 0;

  for i in [1..#s] do
    mi:= Ncols(G[i]);
    ri:= &+[ Integers() | Ncols(G[j]): j in [i+1..#s] ];
    det:= get_disc(G[i]);
    sq:= IsLocalSquare(det, p);
    f /:= GroupOrder( sq select "O+" else "O-", mi, q);
    N -:= s[i] * mi*(mi+1)/2 + s[i]*mi*ri;
    N +:= s[i]*mi*(m+1)/2;              // coming from the volume!
  end for;

  // fix the wierd stuff
  if IsEven(m) and IsOdd(Valuation(d, p)) then
    N +:= (1-m)/2;
  end if;
  //    <q,N,1/f>;
  assert IsIntegral(N);
  N:= Integers() ! N;
  return q^N * f;
end intrinsic;

intrinsic Mass(L::LatMod : StandardMass:= 0) -> FldRatElt
{Returns the mass of L}
//  require IsDefinite(L) : "Only implemented for lattices in definite quadratic spaces";
  m:= Rank(L);
  if m eq 0 then return 1; end if;
  if m eq 1 then return 1/2; end if;

  R:= BaseRing(L);
  K:= NumberField(R);
  r:= m div 2;

  if StandardMass eq 0 then
    StandardMass:= 2^(-AbsoluteDegree(K) * r);
    if IsOdd(m) then
      StandardMass *:= &* [ DedekindZetaExact(K, -i) : i in [1..m-2 by 2] ];
    else
      StandardMass *:= &* [ Rationals() | DedekindZetaExact(K, -i) : i in [1..m-3 by 2] ];
      Form:= GramMatrixOfBasis(L);
      Disc:= get_disc(Form);
      if IsSquare(Disc) then
        StandardMass *:= DedekindZetaExact(K, 1-r);
      else
        E:= ext< K | Polynomial([-Disc,0,1]) >;
        StandardMass *:= DedekindZetaExact(E, 1-r : Relative);
      end if;
    end if;
  end if;

  mass:= Abs(StandardMass);
  for p in BadPrimes(L : Even) do
    mass *:= LocalFactor(L, p);
  end for;
  return mass;
end intrinsic;


// ------------------
// Markus' playground
// ------------------

intrinsic HasseInvariant(M::AlgMatElt, p::RngOrdIdl) -> RngIntElt
{The Hasse-Minkowski invariant of the quadratic form given by M at the place p}
  require IsPrime(p) : "The second argument must be prime";
  F:= FieldOfFractions(Order(p));
  ok, M:= CanChangeRing(M, F);
  require ok: "Incompatible base rings";
  require IsSymmetric(M) : "The form must be symmetric";
  D:= Diagonal(OrthogonalizeGram(M));
  return Hasse(D, p);
end intrinsic;

intrinsic HasseInvariant(M::AlgMatElt, p::RngIntElt) -> RngIntElt
{The Hasse-Minkowski invariant of the quadratic form given by M at the place p}
  require p eq -1 or IsPrime(p: Proof:= false) : "The second argument must be prime";
  ok, M:= CanChangeRing(M, Rationals());
  require ok: "Incompatible base rings";
  require IsSymmetric(M) : "The form must be symmetric";
  D:= Diagonal(OrthogonalizeGram(M));
  return Hasse(D, p);
end intrinsic;

intrinsic QuadraticFormInvariants(M::AlgMatElt[FldAlg]: Minimize:= true) -> FldAlgElt, SetEnum[RngOrdIdl], SeqEnum[RngIntElt]
{The invariants describing the quadratic form M}
  require IsSymmetric(M) and Rank(M) eq Ncols(M): "The form must be symmetric and regular";
  D:= Diagonal(OrthogonalizeGram(M));
  K:= BaseRing(M);
  R:= Integers(K);
  P:= Support(2*R);
  U:= Universe(P);
  for d in D do
    P join:= { f[1] : f in Factorization(d*R) | IsOdd(f[2]) };
  end for;
  F:= Minimize select {U | p : p in P | Hasse(D, p) eq -1 } else { <p, Hasse(D, p) > : p in P };
  I:= [ #[ d: d in D | Evaluate(d, f) le 0 ] : f in RealPlaces(K) ];
  return &* D, F, I;
end intrinsic;

intrinsic QuadraticFormInvariants(M::AlgMatElt[FldRat]: Minimize:= true) -> FldRatElt, SetEnum[RngIntElt], RngIntElt
{"} //"
  require IsSymmetric(M) and Rank(M) eq Ncols(M): "The form must be symmetric and regular";
  D:= Diagonal(OrthogonalizeGram(M));
  P:= { 2 };
  for d in D do
    P join:= { f[1] : f in FactorizationOfQuotient(d) | IsOdd(f[2]) };
  end for;
  det:= Rationals() ! Squarefree(Numerator(d) * Denominator(d)) where d:= &*D;
  F:= Minimize select {Integers() | p : p in P | Hasse(D, p) eq -1 } else { <p, Hasse(D, p) > : p in P };
  return det, F, #[ d: d in D | d le 0 ];
end intrinsic;

function MySquarefree(d)
  T:= Type(d);
  if T eq RngIntElt then 
    return Squarefree(d);
  elif T eq FldRatElt then
    return Rationals() ! Squarefree( Numerator(d) * Denominator(d) );
  end if;

  return d; // FIXME

  if ISA(T, FldElt) then
    _, x:= IsIntegral(d);
    d:= Integers(Parent(d)) ! (d*x^2);
  end if;
  R:= Parent(d); r:= R ! 1;
  if IsSquare(d) then return r; end if;
  F:= Reverse([ < f[1], f[2] div 2 >: f in Factorization(R*d) | f[2] ge 2 ]);
  for i in [1..#F] do
    for j in [1..F[i,2]] do
      ok, x:= IsPrincipal(F[i,1]);
      if ok then
	r *:= x^(F[i,2] div j);
        F[i,2] mod:= j;
	break;
      end if;
    end for;
  end for;

  F:= [f : f in F | f[2] ne 0 ];
  found := true;
  while #F ge 2 and found do
    found:= false;
    Prods:= [ 1*R ]; Vecs:= [ RSpace(Integers(), #F) ! 0 ];
    for i in [1..#F] do
      v:= Universe(Vecs).i;
      for j in [1..F[i,2]] do
	for k in [1..#Prods] do 
	  I:= Prods[k] * F[i,1];
          w:= Vecs[k] + v;
          ok, x:= IsPrincipal(I);
	  if ok then
	    found:= true;
	    for i in [1..#F] do F[i,2] -:= w[i]; end for;
	    r *:= x;
	    break i;
	  else
	    Append(~Prods, I);
	    Append(~Vecs, w);
	  end if;
	end for;
      end for;
    end for;
    F:= [f : f in F | f[2] ne 0 ];
  end while;
  return R ! (d / r^2);
end function;

// Returns all elements in R supported at the prime ideals in PP (up to squares).
function ElementsWithSupport(R, PP)
//  if ISA(Type(PP), SetEnum) then PP:= Setseq(PP); end if;
  if Type(R) eq RngInt then return [-1] cat PP, func<x|x>; end if;
  U, h:= UnitGroup(R);
  Result:= [ h(U.i): i in [1..Ngens(U)] ];
  Cl, h:= ClassGroup(R);
  if #Cl ne 1 then
    F:= FreeAbelianGroup(#PP);
    m:= hom< F -> Cl | [ p @@ h: p in PP ] >;
    K:= Kernel(m);
    for i in [1..Ngens(K)] do
      ok, x:= IsPrincipal(PowerProduct(PP, Eltseq(F ! K.i))); assert ok;
      Append(~Result, x);
    end for;
    _, hh:= quo< Cl | Image(m) >;
    f:= function(I)
      c:= (I @@ h);
      o:= Order(c @ hh);
      J:= PowerProduct(PP, Eltseq((-o * c) @@ m));
      ok, x:= IsPrincipal(I^o*J); assert ok;
      return x;
    end function;
  else
    for p in PP do
      ok, x := IsPrincipal(p); assert ok;
      Append(~Result, x);
    end for;
    f:= func< I | x where _, x:= IsPrincipal(I) >;
  end if;
  return Result, f;
end function;

// Given an element b in a number field K and sets of finite and infinite 
// places P and I of K, return an element a in K such that 
// { v: (a,b)_v = -1 } = P \cup I
// Note that the function coputes the unit and class groups of K!
function FindQA(b, P, I)
  assert IsEven(#I + #P);
  assert forall{ p: p in P | not IsLocalSquare(b, p) };
  assert forall{ i: i in I | Evaluate(b, i) lt 0 };
  R:= Parent(b);
  if ISA(Type(R), Fld) then R:= Integers(R); end if;
  K:= NumberField(R);
  n:= #P;
  m:= #I;

  P:= Setseq(P) cat Setseq((MySupport(R ! 2) join MySupport(b)) diff P); // the places at which we need +1
  I:= I cat [ i: i in RealPlaces(K) | i notin I and Evaluate(b, i) lt 0 ];

  V:= VectorSpace(GF(2), #P+#I);
  target:= V ! ([1^^n] cat [0^^(#P-n)] cat [1^^m] cat [0^^(#I-m)]);
  if IsZero(target) then
    return K ! 1;
  end if;

  found:= false;
  U, h:= UnitGroup(R);
  S:= sub< V | >; Base:= [ K | ]; M:= [];
  SignVector:= func< g | Vector(GF(2), [ (1-HilbertSymbol(K ! g, K ! b, p)) div 2 : p in P ] cat [ (1-Sign(Evaluate(g, i))) div 2 : i in I ]) >;

  // Get initial factor base
  L, f:= ElementsWithSupport(R, P);
  for l in L do
    x:= K ! l;
    v:= SignVector(x);
    if v in S then continue; end if;
    Append(~Base, x); Append(~M, v); S:= sub< V | S, v >;
    if target in S then
      found:= true; break;
    end if;
  end for;

  // Extend the factor base by one more prime ideal
  p:= 2;
  while not found do
    p:= NextPrime(p);
    for d in MySupport(R ! p) do
      if d in P then continue; end if;    // already there
      x:= K ! f(d);
      v:= SignVector(x);
      // target notin S thus we can only hope for target+v in S
      if (target+v) in S then
        Append(~M, v); Append(~Base, x);
        found:= true; break;
      end if;
    end for;
  end while;

  // solve for the exponents and assemble the solution.
  exp:= Solution(Matrix(M), target);
  a:= PowerProduct(Base, ChangeUniverse(Eltseq(exp), Integers()));
  return a;
//  return MySquarefree(a);
end function;

intrinsic MyQuaternionAlgebra(P::[RngOrdIdl], S::SeqEnum[PlcNumElt]) -> AlgQuat
{Returns a quaternion algebra only ramified at the places in P and S}
  ok1, U1 := HasUniverse(P);
  ok2, U2 := HasUniverse(S);
  require ok1 or ok2 : "At least one sequence must have a universe";

  if ok1 then
    K:= NumberField(Ring(U1));
    if ok2 then require K cmpeq NumberField(U2) : "Incompatible sequences"; end if;
  else
    K:= NumberField(U2);
  end if;
  if #P ne 0 then
    P:= Setseq(Set(P));
    require forall{ p: p in P | IsPrime(p) }: "The ideals must be prime";
    b:= WeakApproximation(P, [1^^#P]);
    D:= &* P;
  else
    b:= 1;
  end if;
  if #S ne 0 then
    S:= Setseq(Set(S));
    I:= [];
    for s in S do
      ok, i:= IsInfinite(s);
      require ok: "The second argument must be infinite places";
      if IsReal(s) then Append(~I, i); end if;
    end for;
    // try to keep everything small
    u:= RealWeakApproximation(S, [-Sign(Evaluate(b, s)) : s in S ] : UnitOnly);
    if u cmpne false then
      b:= b*u;
    else // nevermind
      I:= Sort(I);
      Inf:= InfinitePlaces(K);
      b:= CRT(D, I, Integers(K) ! 1, [-Sign(Evaluate(b, Inf[i])) : i in I ]) * b;
    end if;
  end if;
  a:= FindQA(b, Set(P), S);
  return QuaternionAlgebra< K | a, b >;
end intrinsic;



intrinsic QuadraticFormWithInvariants(Dim::RngIntElt, Det::RngIntElt, Finite::Setq[RngIntElt], Negative::RngIntElt) -> AlgMatElt
{Computes a quadratic form of dimension Dim and determinant Det that has Hasse invariants -1 at the primes in Finite.
 The number of negative entries of the real signature is given in Negative}
  requirege Dim, 1;
  require Det ne 0: "The determinant must be nonzero";
  requirerange Negative, 0, Dim;

  require Sign(Det) eq (-1)^(Negative mod 2) : "Real place information does not match the sign of the determinant";

  // Dimension 1
  if Dim eq 1 then
    require IsEmpty(Finite): "Impossible Hasse invariants";
    return Matrix(Rationals(), 1, [ Det ]);
  end if;
 
  Finite:= Set(Finite);
  require forall{p : p in Finite | IsPrime(p : Proof:= false)}:
    "The finite place invariants must be a set of primes";

  // Dimension 2 needs some more love
  if Dim eq 2 then
    ok:= forall(p){p: p in Finite | not IsLocalSquare(-Det, p)};
    require ok: Sprintf("A binary form with determinant %o must have Hasse invariant +1 at the prime:\n%o", Det, p);
  end if;

  // Product formula checking
  require IsEven((Negative mod 4 ge 2 select 1 else 0) + #Finite):
    "The number of places (finite or infinite) with Hasse invariant -1 must be even";

  // Reduce the number of bad primes!
  Det:= Squarefree(Det);

  // OK, a space with these invariants must exist.
  // For final testing, we store the invariants.
  Dim0:= Dim; Det0:= Det; Finite0:= Finite; Negative0:= Negative;

  // Pad with ones
  k:= Max(0, Dim - Max(3, Negative));
  D:= [ 1^^k ];
  Dim -:= k;

  // Pad with minus ones
  if Dim ge 4 then
    assert Dim eq Negative;
    k:= Dim - 3;
    d:= (-1)^k; f:= (k mod 4 ge 2) select {2} else {};
    PP:= Finite join { p[1] : p in Factorization(2*Det) };
    Finite:= { p: p in PP | HilbertSymbol(d/1, -Det/1, p) * (p in f select -1 else 1) * (p in Finite select -1 else 1) eq -1 };
    D cat:= [ (-1)^^k ];
    Det:= IsOdd(k) select -Det else Det;
    Dim := 3;
    Negative:= 3;
  end if;

  // The ternary case
  if Dim eq 3 then
    // The primes at which the form is anisotropic
    PP:= Finite join { f[1]: f in Factorization(2*Det) };
    PP:= [ p : p in PP | HilbertSymbol(-1/1, -Det/1, p) ne (p in Finite select -1 else 1) ];

    // Find some a such that for all p in PP: -a*Det is not a local square
    // TODO: Find some smaller a?! The approach below is very lame.
    a:= &*[ Integers() | Det mod p eq 0 select 1 else p : p in PP ];

    if Negative eq 3 then
      a*:= -1;
      Negative:= 2;
    end if;
    PP:= Finite join { p[1] : p in Factorization(2*Det*a) };
    Finite:= { p: p in PP | HilbertSymbol(a/1, -Det/1, p) * (p in Finite select -1 else 1) eq -1 };
    Det:= Squarefree(Det * a);
    Dim:= 2;
    Append(~D, a);
  end if;

  // The binary case
  a:= FindQA(-Det, Finite, Negative eq 2 select [Infinity()] else []);
  D cat:= [a, Squarefree(Integers() ! (Det*a))];

  ChangeUniverse(~D, Rationals());
  M:= DiagonalMatrix(D);

  // The final test...
  d, f, n:= QuadraticFormInvariants(M);
  assert Dim0 eq #D and d eq Det0 and f eq Finite0 and n eq Negative0;

  return M;
end intrinsic;

intrinsic QuadraticFormWithInvariants(Dim::RngIntElt, Det::FldRatElt, Finite::Setq[RngIntElt], Negative::RngIntElt) -> AlgMatElt
{"} //"
  return QuadraticFormWithInvariants(Dim, Numerator(Det)*Denominator(Det), Finite, Negative);
end intrinsic;

intrinsic QuadraticFormWithInvariants(Dim::RngIntElt, Det::FldAlgElt, Finite::Setq[RngOrdIdl], Negative::[RngIntElt]) -> AlgMatElt
{Computes a quadratic form of dimension Dim and determinant Det that has Hasse invariants -1 at the primes in Finite.
 The number of negative entries of the i-th real signature is given in Negative[i]}

  requirege Dim, 1;
  require Det ne 0: "The determinant must be nonzero";
  K:= Parent(Det);

  // Infinite places checking
  Inf:= RealPlaces(K);
  require #Negative eq #Inf: "The number of negative entries must be the number of real places";
  require forall(i){i : i in [1..#Inf] | Negative[i] in [0..Dim]}:
    Sprintf("Impossible negative entry at position %o", i);
  require forall(i){i : i in [1..#Inf] | Sign(Evaluate(Det, Inf[i])) eq (-1)^(Negative[i] mod 2)}:
    Sprintf("Information at the real place number %o does not match the sign of the determinant", i);

  // Dimension 1
  if Dim eq 1 then
    require IsEmpty(Finite): "Impossible Hasse invariants";
    return Matrix(1, [ Det ]);
  end if;

  // Finite places checking
  R:= Integers(K);
  PI:= PowerIdeal(R);
  require IsEmpty(Finite) or (Universe(Finite) cmpeq PI and forall{p : p in Finite | IsPrime(p)}):
    "The Invariants must be a set of prime ideals of the base field";
  Finite:= Set(Finite);

  // Dimension 2 needs some more love
  if Dim eq 2 then
    ok:= forall(p){p: p in Finite | not IsLocalSquare(-Det, p)};
    require ok: Sprintf("A binary form with determinant %o must have Hasse invariant +1 at the prime ideal:\n%o", Det, p);
  end if;

  // Product formula checking
  require IsEven(#[ n: n in Negative | n mod 4 ge 2 ] + #Finite):
    "The number of places (finite or infinite) with Hasse invariant -1 must be even";

  // OK, a space with these invariants must exist.
  // For final testing, we store the invariants.
  Dim0:= Dim; Det0:= Det; Finite0:= Finite; Negative0:= Negative;

  // Reduce Det
  Det:= K ! MySquarefree(Det);

  // Pad with ones
  k:= Max(0, Dim - Max(3, Max(Negative)));
  D:= [ (K ! 1)^^k ];
  Dim -:= k;

  if Dim ge 4 then
    // Pad with minus ones
    k:= Min(Dim-3, Min(Negative));
    D2 := [ (K ! -1)^^k ];
    Dim -:= k;
    Negative:= [ n-k: n in Negative ];

    // Pad with other entries
    while Dim ge 4 do
      V:= []; Signs:= [];
      for i in [1..#Negative] do
        if Negative[i] eq 0 then 
	  Append(~V, Inf[i]); Append(~Signs, +1);
	elif Negative[i] eq Dim then
	  Append(~V, Inf[i]); Append(~Signs, -1);
	end if;
      end for;
      x:= RealWeakApproximation(V, Signs : Al:= "Small");
      s:= RealSigns(x);
      k:= Min([Dim-3] cat [ s[i] eq 1 select (Dim - Negative[i]) else (Negative[i]) : i in [1..#Negative] ]);
      D2 cat:= [ (K ! x)^^k ];
      Dim -:= k;
      for i in [1..#Negative] do
        if s[i] eq -1 then Negative[i] -:= k; end if;
      end for;
    end while;

    // Adjust invariants: Dim and Negative are already ok.
    d, f:= QuadraticFormInvariants(DiagonalMatrix(D2));
    PP:= Finite join Support(2*R) join Support(Det*R);
    Finite:= { p: p in PP | HilbertSymbol(d, -Det, p) * (p in f select -1 else 1) * (p in Finite select -1 else 1) eq -1 };
    D cat:= D2;
    Det:= K ! MySquarefree(Det * d);
    delete d, f;
  end if;

  // The ternary case
  if Dim eq 3 then
//  "Dim 3", Det;
    // The primes at which the form is anisotropic
    PP:= Finite join Support(2*R) join Support(Det*R);
    PP:= [ PI | p : p in PP | HilbertSymbol(K ! -1, -Det, p) ne (p in Finite select -1 else 1) ];

    // Find some a such that for all p in PP: -a*Det is not a local square
    // TODO: Find some smaller a?! The approach below is very lame.
    // We simply make sure that a*Det has valuation 1 at each prime in PP....
    a:= K ! (#PP eq 0 select 1 else WeakApproximation(PP, [ (1+Valuation(Det, p)) mod 2 : p in PP ]));
    // Fix the signs of a if necessary.
    Idx:= [ i : i in [1..#Negative] | Negative[i] in {0,3} ];
    S:= [ Negative[i] eq 0 select s[i] else -s[i] : i in Idx ] where s:= RealSigns(a);
    a*:= K ! RealWeakApproximation(Inf[Idx], S: Al:= "Small", CoprimeTo:= &*PP);

    // Adjust invariants for the last time:
    s:= RealSigns(a);
    for i in [ i: i in [1..#Negative] | s[i] lt 0 ] do 
      Negative[i] -:= 1;
    end for;
    PP:= Finite join &join [ Support( R*x) : x in [K ! 2, Det, a] ];
    Finite:= { p: p in PP | HilbertSymbol(a, -Det, p) * (p in Finite select -1 else 1) eq -1 };
    Det:= K ! MySquarefree(Det * a);
    Append(~D, a);
  end if;

  // The binary case
  a:= FindQA(-Det, Finite, [ Inf[i]: i in [1..#Inf] | Negative[i] eq 2 ]);
  D cat:= [a, Det*a];

  M:= DiagonalMatrix(D);

  // The final test...
  d, f, n:= QuadraticFormInvariants(M);
  assert Dim0 eq #D and IsSquare(d*Det0) and f eq Finite0 and n eq Negative0;

  return M;
end intrinsic;

RealValuations:= func< V, x | [ GF(2) | Sign(Real(Evaluate(x, v))) eq -1 select 1 else 0 : v in V ] >;

// require 0 notin D
function LocalHasseInvariant(D, p)
  T:= Type(p);
  if (T eq Infty) or (p cmpeq -1) then return (-1)^((s*(s-1)) div 2) where s:= #[ d: d in D | d lt 0 ];
  elif T eq RngInt then p:= Generator(p);
  elif T eq PlcNumElt then
    if IsFinite(p) then p:= Ideal(p);
    elif IsReal(p) then return (-1)^((s*(s-1)) div 2) where s:= #[ d: d in D | Real(Evaluate(p, d)) lt 0 ];
    else return 1;
    end if;
  end if;
  return Hasse(D, p);
end function;

intrinsic IsLocallyIsotropic(D::., p::.) -> BoolElt
{}
  if exists{d: d in D | IsZero(d)} then 
    return true;
  elif #D eq 1 then
    return false;
  end if;

  T:= Type(p);
  if (T eq Infty) or (p cmpeq -1) then 
    return exists{i : i in [2..#D] | Sign(D[i]) ne s} where s:= Sign(D[1]); 
  elif T eq PlcNumElt then
    if IsFinite(p) then
      p:= Ideal(p);
    else
      return IsComplex(p) or exists{i : i in [2..#D] | Sign(Real(Evaluate(D[i], p))) ne s} where s:= Sign(Real(Evaluate(D[1], p)));
    end if;
  elif T eq RngInt then p:= Generator(p);
  end if;

  if #D eq 2 then
    return IsLocalSquare(-&*D, p);
  elif #D ge 5 then
    return true;
  end if;
  K:= Type(p) eq RngIntElt select Rationals() else NumberField(Order(p));
  if #D eq 3 then
    return HilbertSymbol(K!-1, -&*D, p) eq LocalHasseInvariant(D, p);
  else
    return HilbertSymbol(K!-1, K!-1, p) eq LocalHasseInvariant(D, p) or not IsLocalSquare(&*D, p);
  end if;
end intrinsic;

function MyRealWeakApproximation(i, Inf)
  K:= NumberField(i);
  a:= PrimitiveElement(K);

  repeat
    x:= BestApproximation( Evaluate(a, i), 4);
    a:= 2*(a - x);
  until forall{ j: j in Inf | i eq j or Abs(Evaluate(a, j)) ge 2 };
  assert Abs(Evaluate(a, i)) le 1/2;
  return a;
end function;

ToGF2:= func< a,b,p | HilbertSymbol(a,b,p) eq 1 select 0 else 1 >;
SignGF2:= func< x, p | Evaluate(x, p) lt 0 select 1 else 0 >;
MyFact:= func< R, d | Type(R) eq RngInt select FactorizationOfQuotient(Rationals() ! d) else Factorization(R*d) >;

intrinsic IsIsotropic(F::AlgMatElt : IsotropicVector:= false) -> BoolElt, ModTupFldElt
{Checks if the quadratic form F is isotropic. If IsotropicVector is set to true, an isotropic vector is also returned}
  if Ncols(F) eq 0 then return false, _; end if;
  require IsSymmetric(F) : "The form must be symmetric";
  K:= BaseRing(F); T:= Type(K);
  if T eq RngInt then
    K:= Rationals();
    F:= Matrix(K, F); 
  elif ISA(T, RngOrd) or ISA(T, FldOrd) then
    K:= NumberField(K);
    F:= Matrix(K, F); 
  end if;
  T:= Type(K); rat:= T eq FldRat;
  require rat or ISA(T, FldNum) : "The form must be over the rationals or a number field";

  D, T:= OrthogonalizeGram(F);
  D:= Diagonal(D);

  // The trivial cases
  if exists(i){i: i in [1..#D] | D[i] eq 0} then
    return true, T[i];
  elif #D le 1 then
    return false, _;
  elif exists(t){ <i,j> : i in [1..j-1], j in [2..#D] | D[i] eq -D[j] } then
    return true, T[t[1]] + T[t[2]];
  else
    ok, y:= IsSquare(-D[1]/D[2]);
    if ok then 
      return true, T[1] + y*T[2];
    elif #D eq 2 then
      return false, _;
    end if;
  end if;

  if #D eq 3 then
    // No need to reinvent the wheel.
    // Magma can solve norm equations and search for points on conics...
    ok, v:= HasRationalPoint(Conic(D));
    if ok then return ok, Vector(Eltseq(v)) * T; end if;
    return false, _;
  elif #D eq 4 then
    ok, v:= IsIsotropic(DiagonalMatrix(D[3..4]) : IsotropicVector);
    if ok then return true, v[1]*T[3] + v[2]*T[4]; end if;

    // Call the oracle and compute the bad primes.
    if exists{ v : v in RealPlaces(K) | not IsLocallyIsotropic(D, v) } then return false, _; end if;
    R:= Integers(K);
    P:= { };
    for d in Append(D, K!2) do
      for p in { p[1] : p in MyFact(R, d) | p[1] notin P } do
	if not IsLocallyIsotropic(D, p) then return false, p; end if;
        Include(~P, p);
      end for;
    end for;
    if not IsotropicVector then return true, _; end if;

    // Find x ne 0 such that <D[1],D[2]> and <-D[3],-D[4]> both represent x.
    P:= Setseq(P);
    target:= [ GF(2) | ToGF2(D[1], D[2], p) : p in P ] cat [ ToGF2(-D[3], -D[4], p) : p in P ];
    I:= [];
    for p in RealPlaces(K) do
      s:= Sign(Evaluate(D[1], p));
      if s eq Sign(Evaluate(D[2], p)) then 
        Append(~I, p);
	Append(~target, s eq -1 select 1 else 0 );
      else
        s:= Sign(Evaluate(-D[3], p));
        if s eq Sign(Evaluate(-D[4], p)) then 
          Append(~I, p);
  	  Append(~target, s eq -1 select 1 else 0 );
	end if;
      end if;
    end for;
    V:= VectorSpace(GF(2), #target);
    target:= V ! target;
    // We simply need to find some x such that target equals the vector
    // [ ToGF2(-D[1] * D[2], x, p) : p in P ] cat [ ToGF2(-D[3]*D[4], x, p) : p in P ] cat [ SignGF2(x, p): p in I ]

    if IsZero(target) then 
      x:= K ! 1;
    else
      found:= false;
      S:= sub< V | >; Basis:= []; Signs:= [];
      L, h:= ElementsWithSupport(R, P);
      for l in L do
        x:= K ! l;
        s:= V ! ([ ToGF2(-D[1] * D[2], x, p) : p in P ] cat [ ToGF2(-D[3]*D[4], x, p) : p in P ] cat [ SignGF2(x, p): p in I ]);
        if s notin S then 
          S:= sub< V | S, s >;
  	  Append(~Basis, x); Append(~Signs, s);
	  if target in S then found:= true; break; end if;
	end if;
      end for;
      p:= 2;
      while not found do
        p:= NextPrime(p);
        Dec:= rat select [ p ] else [ d[1] : d in Decomposition(R, p) ];
	for p in Dec do
	  if p in P then continue; end if;
          x:= K ! h(p);
          s:= V ! ([ ToGF2(-D[1] * D[2], x, p) : p in P ] cat [ ToGF2(-D[3]*D[4], x, p) : p in P ] cat [ SignGF2(x, p): p in I ]);
          if s+target in S then
	    Append(~Basis, x); Append(~Signs, s);
            found:= true; break;
	  end if;
	end for;
      end while;
      exp:= Solution(Matrix(Signs), target);
      x:= PowerProduct(Basis, ChangeUniverse(Eltseq(exp), Integers()));
    end if;

    ok, v:= IsIsotropic( DiagonalMatrix([ D[1], D[2], -x ]) : IsotropicVector); assert ok;
    ok, w:= IsIsotropic( DiagonalMatrix([ D[3], D[4],  x ]) : IsotropicVector); assert ok;
    v:= v/v[3];
    w:= w/w[3];
    v:= Vector([v[1], v[2], w[1], w[2]]) * T;
    assert InnerProduct(v*F, v) eq 0;
    return true, v;
  else 
    // Dim ge 5, here the real places are the only obstacles.
    ok:= forall{ v : v in RealPlaces(K) | IsLocallyIsotropic(D, v) };
    if not ok or not IsotropicVector then return ok, _; end if;

    // We need D[3..5] to yield both signs at every real place
    I:= RealPlaces(K);
    if exists(t){ <i,j> : j in [i+1..#D], i in [1..#D] | forall{ p: p in I | Sign(Evaluate(D[i], p)) ne Sign(Evaluate(D[j], p)) } } then
      S:= Sym(#D);
      g:= (( t[1] eq 3) select S!1 else S!(t[1],3)) * ((t[2] eq 4) select S!1 else S!(t[2],4));
      T:= PermutationMatrix(K, g) * T;
      ok, D:= IsDiagonal(T * F * Transpose(T) ); assert ok;
    else
      Fix:= [];
      Signs:= [];
      for i in [1..#I] do
        s:= Sign(Evaluate(D[3], I[i]));
        if s ne Sign(Evaluate(D[5], I[i])) then continue; end if;
        if s eq Sign(Evaluate(D[4], I[i])) then
          a:= 1/MyRealWeakApproximation(I[i], I[Fix]);
          ok:= exists(j){j: j in [1..#D] | Sign(Evaluate(D[j], I[i])) ne s}; assert ok;
          r:= 0;
          repeat
            r +:= 1;
            t := D[4] + a^(2*r)*D[j];
          until Sign(Evaluate(t, I[i])) ne s and forall{k: k in Fix | Sign(Evaluate(t, I[k])) eq Signs[k] };
          b:= -a^r * D[j] / D[4];
          v:= T[4];
          T[4] +:= a^r*T[j];
          T[j] +:= b*v;
        end if;
        Append(~Fix, i);
        Append(~Signs, -s);
        ok, D:= IsDiagonal(T * F * Transpose(T)); assert ok;
      end for;
    end if;

    ok, v:= IsIsotropic(DiagonalMatrix(D[3..5]) : IsotropicVector);
    if ok then return true, v[1]*T[3] + v[2]*T[4] + v[3]*T[5]; end if;

    R:= Integers(K);
    P:= [];
    X:= [];
    M:= [];
    for p in Setseq({ p[1] : p in MyFact(R, d), d in D cat [2] }) do
      if IsLocallyIsotropic(D[3..5], p) then continue; end if;

      if IsLocallyIsotropic([ D[3], D[4], D[5], D[1] ], p) then
        x:= 1; y:= 0;
      elif IsLocallyIsotropic([ D[3], D[4], D[5], D[2] ], p) then
        x:= 0; y:= 1;
      else
        // now D[1] and D[2] represent necessarily the same class
        // leaving this class is enough
        V1:= Valuation(D[1], p);
        V2:= Valuation(D[2], p);
        V:= Max(V1, V2);
        pi:= Type(p) eq RngIntElt select p else PrimitiveElement(p);
        k,h:= ResidueClassField(p);
        y:= pi^((V - V2) div 2);
        cnt:= 1;
        repeat
          cnt +:= 1;
          assert cnt le 1000;
          x:= (Random(k) @@ h) * pi^((V - V1) div 2);
        until IsLocallyIsotropic([ D[3], D[4], D[5], x^2*D[1] + y^2*D[2] ], p);
      end if;
      Append(~X, <R ! x, R ! y>);
      Append(~P, p);
      V:= Valuation( x^2*D[1] + y^2*D[2], p) + 1;
      if Minimum(p) eq 2 then
        V +:= Type(p) eq RngIntElt select 2 else 2*RamificationIndex(p);
      end if;
      Append(~M, p^V);
    end for;
    assert #P ne 0;

    x:= CRT( [ x[1]: x in X ], M ); 
    y:= CRT( [ x[2]: x in X ], M );
    t:= x^2*D[1] + y^2*D[2];
    ok, w:= IsIsotropic( DiagonalMatrix([ D[3], D[4], D[5], t ]) : IsotropicVector); assert ok;
    w:= w/w[4];
    v:= Vector([x, y, w[1], w[2], w[3]] cat [0^^(#D-5)]) * T;
    v:= Denominator(v) * v;
    assert InnerProduct(v*F, v) eq 0;
    return true, v;
  end if;
end intrinsic;

intrinsic QuadraticFormDecomposition(F::AlgMatElt) -> ModTupFld, ModTupFld, ModTupFld
{Decompose F into an anisotropic kernel, an hyperbolic space and its radical}
  require IsSymmetric(F) : "The form must be symmetric";
  R:= BaseRing(F);
  if ISA(Type(R), RngOrd) or Type(R) eq RngInt then
    F:= Matrix(NumberField(R), F);
    R:= BaseRing(F);
  end if;
  require Type(R) eq FldRat or ISA(Type(R), FldAlg): "Only implemented for forms over number fields";

  V:= VectorSpace(R, Ncols(F), F);
  R:= Radical(V);
  W:= R;

  H:= [ V | ];
  repeat
    C:= OrthogonalComplement(V, W);
    F:= GramMatrix(C);
    iso, v:= IsIsotropic(F : IsotropicVector);
    if iso then
      ok:= exists(b){ b: b in Basis(C) | InnerProduct(v,b) ne 0 };
      b:= b / InnerProduct(v,b);
      H:= H cat [ V | v, b ];
      W:= sub< V | W, v, b >;
    end if;
  until not iso;

  return C, sub< V | H >, R;
end intrinsic;

intrinsic MaximalIsotropicSubspace(F::AlgMatElt) -> ModTupFld
{Returns a maximal totally isotropic subspace of F}
  _, H, R:= QuadraticFormDecomposition(F);
  return sub< H | [ H.i : i in [1..Dimension(H) by 2] ] > + R;
end intrinsic;

function QFIsoField(F, G, Iso)
  if Ncols(F) ne Ncols(G) then return false, _; end if;
  d1, f1, i1:= QuadraticFormInvariants(F);
  d2, f2, i2:= QuadraticFormInvariants(G);
  if i1 ne i2 or f1 ne f2 or not IsSquare(d1*d2) then return false, _; end if;
  if not Iso then return true, _; end if;
  A1, H1, R1:= QuadraticFormDecomposition(F);
  A2, H2, R2:= QuadraticFormDecomposition(G);
  assert Dimension(H1) eq Dimension(H2) and Dimension(R1) eq Dimension(R2);
  V:= Generic(A1);
  W:= Generic(A2);
  X:= [V | ]; Y:= [W | ];
  while Dimension(A1) gt 0 do
    ok, v:= IsIsotropic( DiagonalJoin(GramMatrix(A1), -GramMatrix(A2)) : IsotropicVector);
    assert ok;
    e:= Eltseq(v);
    n:= #e div 2;
    x:= V ! (Vector(e[  1..  n]) * BasisMatrix(A1));
    y:= W ! (Vector(e[n+1..2*n]) * BasisMatrix(A2));
    Append(~X, x);
    Append(~Y, y);
    A1:= OrthogonalComplement(A1, sub< A1 | x >);
    A2:= OrthogonalComplement(A2, sub< A2 | y >);
//    A1:= OrthogonalComplement(V, sub< V | X >);
//    A2:= OrthogonalComplement(W, sub< W | Y >);
  end while;
  M:= Matrix( Y cat [H2.i: i in [1..Ngens(H2)] ] cat Basis(R2) )^-1 * Matrix( X cat [ H1.i: i in [1..Ngens(H1)] ] cat Basis(R1) );
  assert M * F * Transpose(M) eq G;
  return true, M;
end function;

intrinsic AreEquivalentQuadraticForms(F::AlgMatElt[FldRat], G::AlgMatElt[FldRat] : Isometry:= false) -> BoolElt, AlgMatElt
{"} //"
  return QFIsoField(F, G, Isometry);
end intrinsic;

intrinsic AreEquivalentQuadraticForms(F::AlgMatElt[FldAlg], G::AlgMatElt[FldAlg] : Isometry:= false) -> boolElt, AlgMatElt
{"} //"
  require BaseRing(F) eq BaseRing(G) : "The base rings are not the same";
  return QFIsoField(F, G, Isometry);
end intrinsic;
