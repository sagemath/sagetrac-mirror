//freeze;

/******************************************************************
*                                                                 *
* Spinor genera of lattices over number fields                    *
* BONGs of lattices over number fields                            *
* Genera of lattices over number fields                           *
*                                                                 *
* See D. Lorch, "Einklassige Geschlechter orthogonaler Gruppen",  * 
*                                                                 *
* David Lorch, 2015                                               *
*                                                                 *
******************************************************************/

import "lat.m": SetupLatticeInSpace;

// we fix some primes for the genus computations and add them as a global property of the field / lattice:

declare attributes LatMod:
  GenusComp;
  
GenusComputationFormat:=recformat<
  RCG:GrpAb, // a ray class group associated to some lattice L
  hom_RCG:Map, 

  MM:SeqEnum, // list of prime powers defining RCG
  M:RngOrdIdl, // M = &*MM
  RayPrimes, // like MM, but is a list of the prime ideals only (no exponents)
  Inf:SetIndx,                  // the infinite (real) places of the ray

  RCG_Subgroup:GrpAb, // a subgroup of RCG
  RCG_FactorGroup:GrpAb, // = RCG modulo RCG_Subgroup
  hom_RCG_FactorGroup:Map, 

  CriticalPrimes:SeqEnum, // a set of primes needed to list all the proper spinor genera in Genus(L)
  GroupElementsToPrimes: Assoc // maps the elements of RCG_Factorgroup to primes representing them
>; 

/*
intrinsic deletedata(L::LatMod)
  { deletes genus data stored in L }
  L`GenusComp := rec< GenusComputationFormat | GroupElementsToPrimes:=AssociativeArray() >;
  return;
end intrinsic;
  
// 0. Helper functions, and some basic things which might later be turned into intrinsics:

function GetExtra(L, Primes, Extra)
  // return more primes, the number of which is specified by Extra, which are coprime to Primes and BadPrimes(L).
  assert Extra ne 0;
  K:= NumberField(BaseRing(L));
  Coprime:= &*(Set(Primes) join BadPrimes(L: Even));

  Bnd:= 30;
  repeat
    Ps:= PrimesUpTo(Bnd, K : coprime_to:= 2*Coprime);
    Bnd +:= 10;
  until #Ps ge Extra;
  return Ps[1..Extra];
end function;
*/

function ScalesAndNorms(G, p, Uniformizer)
  // G is a splitting of a LatticeModule L into Gram matrices,
  // which can e.g. be calculated by JordanDecomposition(), but need not
  // be a Jordan decomposition.

  // We calculate the invariants of the entries of G, viewed as individual lattices.
  
  e := RamificationIndex(p);
  sL := [Minimum({Valuation(g[i,j], p): i in [1..j], j in [1..Ncols(g)]}): g in G];

  aL:= [FieldOfFractions(Order(p)) | ]; uL:= []; wL:= [];
  sL := [Minimum({Valuation(x, p): x in Eltseq(g)}): g in G];
  for i in [1..#G] do
    // similar, but not equal to, the code for obtaining a genus symbol
    // (for the genus symbol, we consider a scaling L^(s(L_i)) in O'Meara's notation)
    GG := G[i];

    D:= Diagonal(GG);
    if e+sL[i] le Minimum({ Valuation(d, p) : d in D }) then
      Append(~aL, Uniformizer^(e+sL[i]));
    else
      Append(~aL, D[b] where _, b is Minimum([Valuation(x,p): x in D]));
    end if;
    Append(~uL, Valuation(aL[i], p));
    assert uL[i] eq Valuation(Norm(LatticeModule(G[i])), p);
    Append(~wL, Minimum({e+sL[i]} join {uL[i] + QuadraticDefect(d/aL[i], p): d in D}));
  end for;

  return sL, aL, uL, wL; // scales, norm generators, norm valuations, weight valuations of a (Jordan) decomposition of L
end function;


function NormUpscaled(G, p)
  // calculate Norm(L^{scale(L_i)}), for all i.
  // for the definition of L^{a}, cf. § 82 I, p.237, in O'Meara
  
  // G is a splitting of a LatticeModule L into Gram matrices,
  // which can e.g. be calculated by JordanDecomposition(), but need not
  // be a Jordan decomposition.
    t:= #G;
    sL := [Minimum({Valuation(g[i,j], p): i in [1..j], j in [1..Ncols(g)]}): g in G];
    e := RamificationIndex(p);
    Uniformizer := PrimitiveElement(p);
    
    aL:= [FieldOfFractions(Order(p)) | ]; uL:= []; wL:= [];
    for i in [1..t] do
      // Die Norm ist 2*Scale + <die Ideale die von den Diagonalelementen erzeugt werden>, cf. § 94 O'Meara.
      GG:= DiagonalJoin(< j lt i select Uniformizer^(2*(sL[i]-sL[j])) * G[j] else G[j] : j in [1..t] >);
      D:= Diagonal(GG);
      // Is 2*s(L_i) eq n_i? (We always have scale \supseteq norm \supseteq 2*scale)
      vprintf Genus, 2: "val(2*scale) = %o, val(diagonal)=%o\n", e+sL[i], Minimum({ Valuation(d, p) : d in D }); 
      // um einen Normerzeuger zu finden, muss man einfach Erzeuger des größeren der Ideale 2*Scale bzw. <Diagonalelemente> nehmen:
      if e+sL[i] le Minimum({ Valuation(d, p) : d in D }) then
        // 20150505: note that "le", not "eq", is correct here (we are not comparing 2*s(L) with norm(L), but 2*s(L) with only the diagonal elements).
        // 2*scale ist das größere Ideal
        vprintf Genus, 2: "NU: case 1\n"; 
        Append(~aL, Uniformizer^(e+sL[i]));
      else
        // die Diagonalelemente bilden das größere Ideal:
        vprintf Genus, 2: "NU: case 2\n"; 
        Append(~aL, D[b] where _, b is Minimum([Valuation(x,p): x in D]));
      end if;
      Append(~uL, Valuation(aL[i], p));
      // Append(~wL, Minimum({e+sL[i]} join {uL[i] + QuadraticDefect(d/aL[i], p): d in D}));
    end for;

 // uL: the valuations of the upscaled norms
 // aL: the generators of the upscaled norms
 return uL, aL;
end function;


intrinsic IsPrimitiveLattice(L::LatMod, p::RngOrdIdl)->BoolElt
{ Check whether L is primitive at p, i.e., any Jordan decomposition consists of a 0-valuated component A and a 1-valuated component B (possibly empty), and rank(A) ge rank(B). }
  if Minimum(p) eq 2 then
    dims, scale := GenusSymbol(L, p);
  else
    sym := GenusSymbol(L, p);
    dims := [x[1]: x in sym];
    scale := [x[2]: x in sym];
  end if;
  
  if not (Set(scale) subset {0,1}) then return false; end if;
  
  b := exists(i){i: i in [1..#scale] | scale[i] eq 1};
  if b then dim1 := dims[i]; else dim1 := 0; end if;
  b := exists(i){i: i in [1..#scale] | scale[i] eq 0};
  if b then dim0 := dims[i]; else dim0 := 0; end if;
  
  return dim0 ge dim1;
end intrinsic;

intrinsic IsUnramified(L::LatMod, p::RngOrdIdl) -> BoolElt
{ check whether L is unimodular at p and Norm(L) equals 2*Scale(L) at p. }
  // see Pfeuffer, p.11
  return IsUnimodular(L, p) and (Valuation(Norm(L), p) eq Valuation(2*Scale(L), p));
end intrinsic;

intrinsic IsUnramified(L::LatMod) -> BoolElt
{ check that IsUnramified(L, p) holds for every p over 2 }
  if #[b: b in BadPrimes(L) | Minimum(b) ne 2] gt 0 then return false; end if;
  return &and[IsUnramified(L, p):  p in [x[1]: x in Factorization(2*BaseRing(L))]];
end intrinsic;

  
// 1. BONGs and good BONGs, as defined by C. Beli:

function IsInA(a, p) // needed for G_function
  e := RamificationIndex(p);
  // cf. Lemma 3.5 in Beli: A is the set of all a \in F^* for which either (-a\not\in (F^*)^2 and D(-a) \subseteq O_p) or (-a\in (F^*)^2 and ord(a)\geq -2e).
  def:= QuadraticDefect(-a, p);
  return ( (Type(def) ne Infty) and (def ge 0) ) or
         ( (Type(def) eq Infty) and (Valuation(a, p) ge -2*e) );
end function;

function HasPropertyA(L, p) // needed for SpinorNorm
  error if not Minimum(p) eq 2, "property A only applies to lattices over dyadic fields";
  rL, sL, wL, aL := GenusSymbol(L, p); 
  nL := [Valuation(aL[i], p): i in [1..#aL]];

  r:= Max(rL);
  if r gt 2 then
    vprintf Genus : "[at a prime over 2]: property A is violated because there is a %o-dimensional Jordan component\n", r;
    return false;
  end if;
  
  // GenusSymbol: rL, sL, wL, aL note that aL only contains valuations
  for i in [1..#sL] do for j in [i+1..#sL] do
    if not ((0 lt nL[j] - nL[i]) and (nL[j] - nL[i] lt 2*(sL[j] - sL[i]))) then
    vprintf Genus : "[at a prime over 2]: property A is violated at %o, %o (norm/scale valuations do not fit)\n", i, j;
    return false; 
    end if;
  end for; end for;
  return true;
end function;

function N_function(a, g, p)
  // g is the mapping obtained by LocalMultiplicativeGroupModSquares(p).
  // This is expensive to calculate, so we pass it as an argument.

  // project a into F^*/(F^*)^2:
  b := a@@g;
  V := Parent(b);  // Full F_2-vector space of some dimension over GF(2)

  // treat the squares separately:
  if IsZero(b) then return V; end if;

  // cf. paragraph before 1.2 in Beli 2003:
  // N(a) := N(F(\sqrt(a))/F), i.e. the values of the norm mapping of the quadratic extension
  // and N(a) is the orthogonal complement of (<a>(F^*)^2). 
  B := Basis(V);
  preimages := V!0;
  for i in [1..#B] do
    b := B[i];
    preim := b@g;
    // convert the Hasse symbol (as an element of the multiplicative group C_2)
    // to an element of the additive group F_2:
    // printf "Hilbert symbol with %o-th basis vector: %o\n", Index(B, b), HilbertSymbol(a, preim, p);
    preimages[i] := HilbertSymbol(a, preim, p) eq 1 select 0 else 1;
  end for;
  KernelGenerators := {B[i]: i in [1..#B] | preimages[i] eq 0};
  // (this fails if a is a square)
  j := Min({i: i in [1..#B] | preimages[i] eq 1});
  KernelGenerators join:= {B[i] + B[j]: i in [1..#B] | (preimages[i] eq 1) and (i ne j)};

  return sub<V | KernelGenerators>;
end function;


function OnePlusPowerOfP(k, V, g, p)
  // See Beli 2003, Def. 1. 
  // We expect V, g = LocalMultiplicativeGroupModSquares(p)
  
  S := {v@g: v in V};
  S := {s: s in S | RelativeQuadraticDefect(s, p) ge k};
  return sub<V | [s@@g: s in S]>;
end function;


function G_function(a, V, g, p)
  // use LocalMultiplicativeGroupModSquares to calculate in O*/(O^*)^2
  // (or F^*/(F^*)^2 -- the last component of the result is the p-valuation
  // mod 2, and  F^*/(F^*)^2 =  O*/(O^*)^2 \times C_2.)
  // Also we expect V, g = LocalMultiplicativeGroupModSquares(p)

  // cf. Beli 2003, Def. 4.
  
  e := RamificationIndex(p);
  R := Valuation(a, p);
  d := RelativeQuadraticDefect(-a, p);
 
  if not IsInA(a, p) then 
    vprintf Genus : "G_function, case F\n";
    return N_function(-a, g, p);
  elif Valuation(-4*a, p) eq 0 and IsLocalSquare(-4*a, p) then
    vprintf Genus : "G_function, case G\n";
    return sub<V|[V.i: i in [1..Dimension(V)-1]]>;
  elif (R gt 4*e) then 
    vprintf Genus : "G_function, case H\n";
    return sub<V| a@@g>;
  elif (2*e lt R) and (R le 4*e) then
    if (d le 2*e - R/2) then
      vprintf Genus : "G_function, case A\n";  
      return N_function(-a, g, p) meet (sub<V| (a@@g)> + OnePlusPowerOfP(R + d - 2*e, V, g, p));
    else
      vprintf Genus : "G_function, case B\n";
      assert R mod 2 eq 0;
      return (sub<V| (a@@g)>) + OnePlusPowerOfP(R div 2, V, g, p);
    end if;
  elif (R le 2*e) then 
    // printf "e = %o, R = %o, d = %o\n", e, R, d;
    if (d le e - R/2) then
      vprintf Genus : "G_function, case C\n";
      return N_function(-a, g, p);
    elif (e - R/2 lt d) and (d le 3*e/2 - R/4) then
      assert R mod 2 eq 0;
      vprintf Genus : "G_function, case D\n";
      return N_function(-a, g, p) meet OnePlusPowerOfP((R div 2) + d - e, V, g, p);
    else
      vprintf Genus : "G_function, case E\n";
      //assert R mod 4 eq 0 and e mod 2 eq 0;
      // attention! Use the Floor function wherever Beli writes stuff in square brackets. This is due to his citing Hsia's papers, which have this convention.
      return OnePlusPowerOfP(e - (Floor(e/2 - R/4)), V, g, p);
    end if;
  else error "this never happens.";
  end if;
end function;


intrinsic IsGoodBONG(BONG::[], p::RngOrdIdl) -> BoolElt
{Returns true iff BONG is a good BONG at p, as defined by C. Beli.}
  // Given: a BONG of a LatticeModule L at prime p.
  v := &and[Valuation(BONG[i], p) le Valuation(BONG[i+2], p): i in [1..#BONG-2]]; 
  return v;
end intrinsic;

procedure beli_correction(L, ~G, ~JJ, Steps, i, j, p);
  // this is helper procedure for GoodBONG().
  // Correct the Jordan components with number i and j of L from the Jordan decomposition given in G and JJ,
  // in such a way that the new component with number i has maximal norm. 
  assert #Steps[i] eq 2;

  /*
  // for debugging:
  NU := NormUpscaled(G, p);
  // assert orthogonality of the vectors in JJ[Steps[i]] and those in JJ[Steps[j]], i.e. those that make up G[i] and G[j]:
  if Nrows(G[j]) eq 2 then assert IsZero( C*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[i][1], Steps[i][2]]) where C is  RowSubmatrix(JJ, [Steps[j][2]])); end if;

  assert NU[i] lt Valuation(Norm(LatticeModule(G[i])), p); // otherwise, the norm condition is violated
  assert NU[j] le Valuation(Norm(LatticeModule(G[j])), p); // "<=" must always hold, otherwise something about the lattice is broken
  */

  // we will need the original Gram matrix for re-orthogonalization:
  GI := G[i]^(-1); // over the quotient field
  assert #Steps[j] in [1,2];
  
  
  // if G[j] is two-dimensional, make sure that a norm generator is on the [1,1] position:
  if (Ncols(G[j]) eq 2) and (Valuation(G[j][1,1], p) gt Valuation(G[j][2,2], p)) then
    temp := JJ[Steps[j][1]];
    JJ[Steps[j][1]] := JJ[Steps[j][2]];
    JJ[Steps[j][2]] := temp;
    G[j] := B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[j][1]..Steps[j][#Steps[j]]]);
  end if;

  JJ[Steps[i][1]] +:= JJ[Steps[j][1]];
  
  // update Gram matrix for component #i:
  G[i] := B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[i][1], Steps[i][2]]);
  x := JJ[Steps[i][1]];
  y := JJ[Steps[i][2]];
  // Ansatz: v = JJ[Steps[j][1]] + lambda*x + mu*y
  
  // re-orthogonalize the first basis vector of G[j]:
  v := Vector([-G[j][1,1], 0]) * GI; // G[j][1,1] is always non-zero (the lattice is definite)

  // assert integrality at p:
  assert &and[Valuation(v[k], p) ge 0: k in [1..Ncols(v)]];
  JJ[Steps[j][1]] +:= v[1]*JJ[Steps[i][1]] + v[2]*JJ[Steps[i][2]];
  
  // if applicable, re-orthogonalize the second basis vector of G[j]:
  if Ncols(G[j]) eq 2 then
    w := Vector([-G[j][1,2], 0]) * GI; // G[j][1,2] is always non-zero (or else the lattice would have been further decomposed into two 1x1-lattices here)
    assert &and[Valuation(v[k], p) ge 0: k in [1..Ncols(w)]];
    JJ[Steps[j][2]] +:= w[1]*JJ[Steps[i][1]] + w[2]*JJ[Steps[i][2]];
  end if;

  // update Gram matrix for component #j:
  G[j] := B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[j][1]..Steps[j][#Steps[j]]]);

  /*
  assert Valuation(Norm(LatticeModule(G[i])), p) eq NU[i];
  */
end procedure;

function IsMaximalNormSplittingInternal(GramMatrices, Scales, Norms, p)
  // Scales: list of valuations of scales
  // Norms: list of generators of norms
  // occurring in a Genus symbol of L at p, as calculated
  // by GenusSymbol(L, p). 

  Norms := [Valuation(Norms[i], p): i in [1..#Norms]];

  // test if each component is either unary or binary:
  fail := exists(i){i: i in [1..#GramMatrices] | Ncols(GramMatrices[i]) notin {1,2}};
  if fail then 
    vprintf Genus : "IsMaximalNormSplittingInternal failed: components are not all unary or binary\n";
    return false, -i, {}; 
  end if;

  // test if binary components are modular:
  fail := exists(i){i: i in [1..#GramMatrices] | Ncols(GramMatrices[i]) ne 1 and Valuation(Determinant(GramMatrices[i]), p) ne 2*Scales[i] };
  if fail then
    vprintf Genus : "IsMaximalNormSplittingInternal failed: at least one of the binary components is not modular\n";
    return false, -i, {}; 
  end if;
  

  // test if sL[1] \supseteq sL[2] \supseteq ... \supseteq sL[#sL]:
  for i in [1..#Scales-1] do if (Scales[i] gt Scales[i+1]) then 
    error Sprintf("scale condition at components %o/%o not met, perhaps something is wrong with your lattice?\nFor your reference, the scale valuations are %o", i, i+1, Scales);
    return false, 0, {}; end if; 
  end for;

  // test if nL[i] = n(L^{sL[i]}):
  NU, _ := NormUpscaled(GramMatrices, p);
  failset := {}; // for testing purposes, we count the components that do not have maximal norm
  for i in [1..#Scales] do
    assert NU[i] le Norms[i];
    if NU[i] lt Norms[i] then 
      vprintf Genus : "IsMaximalNormSplittingInternal failed: norm condition at component %o\n", i;
      Include(~failset, i);
    end if;
  end for;
  if #failset gt 0 then return false, SetToSequence(failset)[1], failset; end if;
  return true, 0, {};  
end function;

intrinsic IsMaximalNormSplitting(G::[], p::RngOrdIdl) -> BoolElt
{returns true iff the given list G of Gram matrices corresponds to a maximal norm splitting at p.}
  sL, aL, _, _ := ScalesAndNorms(G, p, PrimitiveElement(p));   
  return IsMaximalNormSplittingInternal(G, sL, aL, p);
end intrinsic;

intrinsic MaximalNormSplitting(L::LatMod, p::RngOrdIdl) -> [], []
{A maximal norm splitting of L at a dyadic prime p, as defined by C. Beli. 
 Returns a Jordan decomposition into 1x1 and 2x2 components, and the corresponding list of basis vectors.}
  // we follow Beli, "Representation of Integral Quadratic Forms over Dyadic Local Fields", section 7

  require BaseRing(L) cmpeq Order(p) : "Incompatible arguments";
  require Minimum(p) eq 2 and IsPrime(p): "p must be a dyadic prime.";

  e := RamificationIndex(p);
  Uniformizer := PrimitiveElement(p);

  J, G:= JordanDecomposition(L, p);

  // join J into one matrix of base vectors:
  JJ := VerticalJoin(<j: j in J>);
  // join the individual Gram matrices:
  A  := DiagonalJoin(<g: g in G>);

  // read the finer decomposition:
  G := [* *];
  Steps:= [];
  i:= 1; n:= Ncols(A);
  while i le n do
    size:= i ne n and A[i,i+1] ne 0 select 2 else 1;
    Append(~G, Submatrix(A, i, i, size, size));
    Append(~Steps, [i..i+size-1]);
    i +:= size;
  end while;

  // This handles the case that p is odd:
  if forall{s: s in Steps | #s eq 1} then
    return [[v]: v in Rows(JJ)], G;
  end if;

  // now turn this decomposition into unary/binary components into a maximal norm splitting:
  failset := {};
  while true do
    // Scale valuations sL, Norm generators aL, Norm valuations uL,
    // weight valuations wL of the decomposition:
    sL, aL, uL, wL := ScalesAndNorms(G, p, Uniformizer); 

    failset_old := failset;
    b, i, failset := IsMaximalNormSplittingInternal(G, sL, aL, p);
    assert (failset_old eq {}) or (#failset_old gt #failset);
    if b then break; end if; // maximal norm splitting reached!
    assert i gt 0; // if b is false and i=0, something is wrong with our lattice

    // here comes the interesting stuff:
    assert #Steps[i] eq 2; // unary components already have maximal norm.
  
    // the maximal norm splitting condition is violated at index i.
    find_j := [Valuation(aL[k], p): k in [i+1..#aL]] cat [2*(sL[i] - sL[k]) + Valuation(aL[k], p): k in [1..i-1]];
    assert #find_j gt 0;
    min, j := Minimum(find_j);

    if j le #aL - i then
      j := j + i;  // we want j to correspond to a Jordan component, not an index in find_j
      assert j gt i;

      // This is case 1 in Beli.
      beli_correction(L, ~G, ~JJ, Steps, i, j, p); 
    else 
      j := j - (#aL - i); // we want j to correspond to a Jordan component, not an index in find_j
      assert j lt i;
      
      // This is case 2 in Beli.
      // switch to a basis of the dual lattice L^#: 
      for k in [1..#Steps] do 
        for l in [1..#Steps[k]] do 
          JJ[Steps[k][l]] *:= Uniformizer^(-sL[k]); 
        end for; 
        assert Valuation(Scale(LatticeModule(B*L`Form*Transpose(B))), p) eq -sL[k] where B is RowSubmatrix(JJ, [Steps[k][1]..Steps[k][#Steps[k]]]);
      end for;

      // apply case 1 to the reversed orthogonal sum of the dual lattices:
      
      Steps := Reverse(Steps);
      // update Gram matrices for the reversed, dualized lattice:
      for k in [1..#G] do 
        G[k] :=  B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[k][1]..Steps[k][#Steps[k]]]);
        assert Nrows(G[k]) in {1,2};
      end for;
      
      assert #Steps[#aL - i + 1] eq 2; // component i is now at position #aL-i+1
      
      beli_correction(L, ~G, ~JJ, Steps, #aL - i + 1, #aL - j + 1, p);

      Steps := Reverse(Steps);
      G := Reverse(G); 

      // update norms/scales from the updated Gram matrices:
      sL, aL, uL, wL := ScalesAndNorms(G, p, Uniformizer); 
      // dualize again:
      for k in [1..#Steps] do for l in [1..#Steps[k]] do 
        JJ[Steps[k][l]] *:= Uniformizer^(-sL[k]); 
      end for; end for;

      // update Gram matrices:
      for k in [1..#G] do 
        G[k] :=  B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[k][1]..Steps[k][#Steps[k]]]);
        assert Nrows(G[k]) in {1,2};
      end for;
    end if;
  end while;

  assert &and[Nrows(G[k]) in {1,2}: k in [1..#G]];
  return [ JJ[x]: x in Steps ], G;
end intrinsic;

function MakeBONGdim2(A, p);
  // cf. Beli, Lemma 3.3
  // return a BONG of a 2-dimensional lattice represented by the Gram matrix A

  a:= Valuation(A[1,1], p);
  b:= Valuation(A[2,2], p);
  s:= Minimum( [ a, b, Valuation(A[1,2], p) ] );
  det:= Determinant(A);
  assert Valuation(det, p) eq 2*s;
  norm:= Minimum( [ a, b, RamificationIndex(p)+s ] );
  S := Parent(A) ! 1;

  // Make the (1,1) entry a norm generator
  if a ne norm then
    if b eq norm then
      SwapRows(~S, 1, 2);
    else
      S[1] +:= S[2];
    end if;
    A := S*A*Transpose(S);
    assert Valuation(A[1,1], p) eq norm;
  end if;
  S[2] +:= A[1,2]/A[1,1] * S[1];

  return S, [ A[1,1], det/A[1,1] ];
end function;

intrinsic GoodBONG(L::LatMod, p::RngOrdIdl) -> [], []
{Return a good BONG of L at a dyadic prime p, as defined by C. Beli.}
  // Any BONG obtained from a maximal norm splitting is good, cf. Corollary 4.2 in Beli.
  // If a maximal norm splitting is calculated, a good BONG can be read off by joining
  // together the BONGs of the 1- and 2-dimensional components.
  require Order(p) cmpeq BaseRing(L) : "Incompatible arguments";
  require Minimum(p) eq 2 : "p must be a dyadic prime.";

  JJ, G:= MaximalNormSplitting(L, p);
  BONG := []; Diag:= [];

  for i in [1..#G] do
    GG := G[i];
    assert Nrows(GG) le 2;
    if Nrows(GG) eq 2 then
      S, D:= MakeBONGdim2(GG, p);
      BONG cat:= Rows(S * Matrix(JJ[i]));
      Diag cat:= D;
    else
      Append(~BONG,  JJ[i,1]);
      Append(~Diag, GG[1,1]);
    end if;
  end for;
  return BONG, Diag;
end intrinsic;


// 2. Spinor norms of lattices over number fields:

intrinsic SpinorNorm(L::LatMod, p::RngOrdIdl) -> ModTupFld, Map, BoolElt
{The spinor norm of L at p, as calculated by C. Beli for dyadic p and by M. Kneser for non-dyadic p. 
 Returns a subspace of LocalMultiplicativeGroupModSquares(p), and a boolean which is true iff the spinor norm is exactly the units}
  require Order(p) cmpeq BaseRing(L) : "Incompatible arguments";
  V, g := LocalMultiplicativeGroupModSquares(p);  
  e:= Valuation(BaseRing(L)! 2, p);

  if Minimum(p) ne 2 then
    // cf. Satz 3 in Kneser, "Klassenzahlen indefiniter quadratischer Formen in drei oder mehr Veränderlichen":
    J, G, E := JordanDecomposition(L, p);
    if exists{g: g in G | Ncols(g) ge 2} then
      vprintf Genus : "This lattice has a 2-dimensional Jordan constituent, and p is odd. Spinor norm is either F^* or O_F^*(F^*)^2, i.e. we will find a vector space of dimension %o or %o.\n", Rank(V)-1, Rank(V);
      // well, which of the two is the case?
      if #{e mod 2: e in E} eq 1 then
        // the units only: 
        return sub<V | [V.i: i in [1..Dimension(V)-1]]>, g, true;
      else
        // the whole group:
        return V, g, false;
      end if;
    end if; 
    // this is the obscure case where p is odd and each
    // of its Jordan components has dimension 1. In this case, the units
    // need not be contained in the Spinor norm.
   
    // generators of the (principal) norm ideals of the Jordan
    // components: since p is odd, the norms (and the scales) are just the 
    // dimensions of the Jordan components
    NormGens := [g[1,1]: g in G];

    TwoNormGens := [];
    for i in [1..#NormGens] do for j in [i..#NormGens] do 
      Append(~TwoNormGens, NormGens[i]*NormGens[j]);
    end for; end for;
    TwoNormVectors := [x@@g: x in TwoNormGens];
    
    vprintf Genus : "odd p, norm generators of the %o Jordan components are: \n\n%o\n\n Products of two elements are \n\n%o\n\nTwoNormVectors:\n%o\n\n", #G, NormGens, TwoNormGens, TwoNormVectors; 
    
    // cf. Kneser 1956, Satz 3:
    SN := sub<V| TwoNormVectors>;
  else
    SN := sub<V| 0>;
    _, BONG := GoodBONG(L, p);
    vprintf Genus : "BONG = %o\n", BONG;
    // assert IsGoodBONG(BONG, p);
    if not HasPropertyA(L, p) then
      vprintf Genus : "This lattice does not have property A. Spinor norm is either F^* or O_F^*(F^*)^2, i.e. we will find a vector space of dimension %o or %o.\n", Rank(V)-1, Rank(V);
      // Using section 7, Thm. 3 in Beli 2002, one can decide which of the two cases applies. 
      // This is why we needed to compute a *good* BONG:
      for i in [1..#BONG-1] do
        BG := Basis(G_function(BONG[i+1]/BONG[i], V, g, p));
        for bg in BG do if bg[#Basis(V)] ne 0 then 
          // the whole group:
          return V, g, false;
        end if; end for;
      end for;
      for i in [1..#BONG-2] do
        if Valuation(BONG[i], p) eq Valuation(BONG[i+2], p) then 
          assert GF(2)!(Valuation(BONG[i+1], p) - Valuation(BONG[i], p)) eq 0;
          if GF(2)!((Valuation(BONG[i+1], p) - Valuation(BONG[i], p))/2) ne GF(2)!e then
            // the whole group:
            return V, g, false;
          end if;
        end if;
      end for;
       
      // if all checks have passed, the Spinor norm is exactly the units:
      return sub<V | [V.i: i in [1..Dimension(V)-1]]>, g, true;
    end if;

    // cf. Beli 2003, Thm. 1
    for i in [2..Rank(L)] do SN +:= G_function(BONG[i]/BONG[i-1], V, g, p); end for;
  
    // for why we can take the Minimum in what follows here, see the remark on p. 161 in Beli 2003:
    if exists{i: i in [1..Rank(L)-2]| GF(2)!(Valuation(BONG[i+2], p) - Valuation(BONG[i], p)) eq 0} 
    then 
      alpha := IntegerRing()!Minimum({(Valuation(BONG[i+2], p) - Valuation(BONG[i], p))/2 : i in [1..Rank(L)-2]|  GF(2)!(Valuation(BONG[i+2], p) - Valuation(BONG[i], p)) eq 0});
      // Uniformizer := PrimitiveElement(p);
      // SN +:= sub<V | (1+Uniformizer^alpha)@@g>;
      SN +:= OnePlusPowerOfP(alpha, V, g, p);
    end if;
  end if;

  return SN, g, SN eq sub<V | [V.i: i in [1..Dimension(V)-1]]>;
end intrinsic;

intrinsic SpinorsNotExactlyTheUnits(L::LatMod) -> SeqEnum
{those places of BaseRing(L) where the Spinor norm is not exactly the units.}
  assert Rank(L) ge 3;
  // Rank >= 2 necessary for SpinorNorm().
  // Rank >= 3 necessary for \Phi from O'Meara 102:7 to be surjective.
  PP := Factorization(Discriminant(L));
  S := []; // set of bad primes
  for i in [1..#PP] do
    P := PP[i][1];
    SN, _, b := SpinorNorm(L, P);
    /*
    // we double-check:
    ExactlyTheUnits := sub<V| [V.i: i in [1..Rank(V) - 1]]> where V is LocalMultiplicativeGroupModSquares(P);  

    assert b eq (SN eq ExactlyTheUnits);
    */
    if not b then Include(~S, [*P, SN*]); end if;
  end for;
  return S;
end intrinsic;



// 3. Enumerating genera of lattices over number fields:

intrinsic MapIdeleIntoClassGroup(L::LatMod, Idele::[]: AtInfinity:=false) -> GrpAbElt
{ map an idele into the class group associated to L. The parameter Idele must be a list of tuples <p, a_p>, where p is a prime of BaseRing(L), and a_p is an element of K^* (interpreted as an element of the completion K^*_p). The parameter AtInfinity can be a list of tuples <v, +1 or -1>, where v is an element of RealPlaces(NumberField(L)). All places, finite or infinite, which are unspecified are interpreted as 1.}

  PrepareClassGroup(L);

  // sanity checks: 
  F := NumberField(L);
  R := Integers(F);
  ok, Idele:= CanChangeUniverse(Idele, car< PowerIdeal(R), F >);
  require ok and forall{ i: i in Idele | IsPrime(i[1]) }: "Idele: must be a list of tuples <p, a_p> with p a prime ideal of BaseRing(L) and a_p an element of NumberField(L)";
  IP := L`GenusComp`Inf;
  the_idele_inf := [1: I in IP];
  if AtInfinity cmpne false then
    require Type(AtInfinity) eq SeqEnum: "AtInfinity: must be false, or a list";
    for i in AtInfinity do
      require (#i eq 2) and (i[2] in {-1,1}) and (Type(i[1]) eq PlcNumElt) and IsReal(i[1]): "AtInfinity: must be a list of tuples <v, +1/-1> with v in RealPlaces(NumberField(L))";
      idx:= Index(IP, i[1]);
      if idx ne 0 then the_idele_inf[idx]:= i[2]; end if;
    end for;
  end if;
  
  // The finite places we need to make congruent to 1 in order to be able to map into the class group:
  the_idele := [F!1: p in L`GenusComp`RayPrimes]; 
  for i in Idele do 
    j := Index(L`GenusComp`RayPrimes, i[1]);
    if j gt 0 then the_idele[j] := i[2]; end if;
  end for;
  
  // L`GenusComp`hom_RCG is the map from L`GenusComp`RCG into the divisors that Magma returns. 
  // So, only the ideals D that are coprime to the ray  lie in the codomain of the map and can be mapped 
  // into L`GenusComp`RCG via I@@L`GenusComp`hom_RCG.
  // The ideles we are considering here are considered to be representatives of classes of ideles (modulo (F^*_S)*J^L, where F^*_S is the set of principal ideles which are totally positive at the infinite places where the envelopping space V of L is anisotropic, and where J^L is the set of ideles which lie in the spinor norm of L at all finite places. 
  // So, we may modify the ideles in the following two ways without changing their class:
  // 1) multiply every component by an element of F^*
  // 2) multiply any individual component (at p) by an element of the spinor norms of L_p. In particular, any component is only relevant modulo squares.
  // Thus we can achieve that the idele in question is congruent to 1 modulo the ray M that defines L`GenusComp`RCG.
  // The idele is then interpreted as the ideal \prod_p{p^Valuation(idele_p, p)}, which can be mapped into the class group by Magma.

  // first, we need to be able to invert modulo the divisor L`GenusComp`M:
  if &or[Valuation(the_idele[i], L`GenusComp`RayPrimes[i]) ne 0: i in [1..#the_idele]] then
    // we need to correct with an element of F^* which is positive at the relevant infinite places.
    // Since InverseMod only works for finite places, we will first correct only the finite places ...
    s := WeakApproximation( IndexedSetToSequence(L`GenusComp`RayPrimes), [-Valuation(the_idele[i], L`GenusComp`RayPrimes[i]): i in [1..#L`GenusComp`RayPrimes]]); 
  else 
    s := F!1;
  end if;

  // Now s * Idele is a unit at the RayPrimes
  // Now rescale once more with some t in F^* such that t*s*Idele = 1 mod MM and has the correct signs.

  sIdele:= [ R ! quo< R | L`GenusComp`MM[i] > ! (s*the_idele[i]) : i in [1..#the_idele] ];
  x := (#L`GenusComp`RayPrimes eq 0) select R ! 1 else InverseMod(CRT(sIdele, L`GenusComp`MM), L`GenusComp`M);
  // Unfortunately, CRT wants SORTED integers ... 
  Inf:= [ idx where _, idx:= IsInfinite(v) : v in L`GenusComp`Inf ];
  Signs:= [ Sign(Evaluate(s,IP[j])) * the_idele_inf[j] : j in [1..#IP] ];
  ParallelSort(~Inf, ~Signs);
  t := CRT(L`GenusComp`M, Inf, x, Signs);

  // Now this element s should do the trick:
  s *:= t;

  // Check if everything is ok.
  assert &and[ IsOne(quo< R | L`GenusComp`MM[k] > ! (s * the_idele[k])): k in [1..#the_idele]];
  assert &and[ Evaluate(s * the_idele_inf[j], IP[j]) gt 0: j in [1..#IP]];

  // This idele can be sent into the ray class group L`GenusComp`RCG via the (inverse of the) homomorphism returned by Magma.
  // We first interpret it as the ideal which will actually have to be mapped:
  // i.e., we just collect the p-valuations at the noncritical places (p notin RayPrimes):
  temp := [ Idele[j][1]^Valuation(Idele[j][2], Idele[j][1]): j in [1..#Idele]];
  ideal:= #temp eq 0 select ideal< R | s > else s * &*temp;

  g:= ideal @@ L`GenusComp`hom_RCG;

  // The second return value is the element in the factor group of the ray class group actually corresponding to a proper spinor genus in Genus(L).
  // Note that this hom might not be set up yet, if this function gets called from PrepareClassGroup.
  if assigned L`GenusComp`hom_RCG_FactorGroup then
    return g, g@L`GenusComp`hom_RCG_FactorGroup;
  else
    return g, _ ;
  end if;
end intrinsic;


intrinsic MapPrincipalIdeleIntoClassGroup(L::LatMod, a::RngElt) -> GrpAbElt
{Map the principal idele defined by the element 'a' into the ray class group identified with the proper spinor genera in Genus(L)} 
  F := NumberField(L);
  ok, a:= IsCoercible(F, a);
  require ok : "Incompatible arguments";
  
  // Finite places:
  decomp := Factorization(a*Integers(F));
  the_idele := [ <x[1], a> : x in decomp ];
  // Real places:
  IP := RealPlaces(F);
  the_idele_infinite := [ <inf, Sign(Evaluate(a, inf))> : inf in IP ];
  
  return MapIdeleIntoClassGroup(L, the_idele: AtInfinity:=the_idele_infinite);
end intrinsic;

  
function GetCriticalPrimes(L: FullList:=false)
  // Set FullList:=true to calculate one prime for every element of L`GenusComp`RCG_FactorGroup,
  // instead of just one prime for every generator. This is needed for indefinite lattices.
  
  assert assigned L`GenusComp;
  
  CriticalPrimes := [];
  F := NumberField(BaseRing(L));
  MyRCG_FactorGroup := FullList select {Identity(L`GenusComp`RCG_FactorGroup) } else sub<L`GenusComp`RCG_FactorGroup | [] >; // set to the trivial subgroup and not to [], otherwise a trivial element might be picked as generator below
  
  maxnorm := 50; // this is really arbitrary -- we just grab some primes to start with and will calculate more later if needed.
  GoodPrimes := [p: p in PrimesUpTo(maxnorm, F: coprime_to:=&*BadPrimes(L : Even))];
  p_ind := 1;
  RCG_FactorGroupGens := []; 
  
  
  if FullList then vprintf Genus : "Looking for %o group elements.\n", #L`GenusComp`RCG; end if;

  while #MyRCG_FactorGroup lt #L`GenusComp`RCG_FactorGroup do
    while p_ind gt #GoodPrimes do
      // calculate more good primes!
      maxnorm := Floor(maxnorm * 1.2); 
      GoodPrimes := [p: p in PrimesUpTo(maxnorm, F: coprime_to:=&*BadPrimes(L : Even))];
    end while;
    p := GoodPrimes[p_ind]; 
    
    
    // map the (1,...,1,pi,1,...,1) idele into the class group, where pi is a uniformizing element of p:
    _, h := MapIdeleIntoClassGroup(L, [<p, F!PrimitiveElement(p)>]); // @L`GenusComp`hom_RCG_FactorGroup;
    
    if not(h in Keys(L`GenusComp`GroupElementsToPrimes)) then L`GenusComp`GroupElementsToPrimes[h] := p; end if;
    
    oldsize := #MyRCG_FactorGroup;
    Append(~RCG_FactorGroupGens, h);    
    if not FullList then 
      MyRCG_FactorGroup := sub<L`GenusComp`RCG_FactorGroup | RCG_FactorGroupGens>;
    else
      MyRCG_FactorGroup join:= {h};
    end if;
    if #MyRCG_FactorGroup gt oldsize then 
      vprintf Genus : " -- subgroup size increased from %o to %o.\n", oldsize, #MyRCG_FactorGroup;
      // this prime ideal is relevant for neighbour generation, store it:
      Append(~CriticalPrimes, p);
    end if;  
    p_ind +:= 1;
  end while;
  return CriticalPrimes;
end function;
  
  
  
intrinsic PrepareClassGroup(L::LatMod)
{internal use.}
  if assigned L`GenusComp then return; end if;
  F := NumberField(L);
  R := Integers(F);
  U := PowerIdeal(R);
  L`GenusComp := rec< GenusComputationFormat | M:= 1*R, MM:= [ U | ], RayPrimes:= {@ U | @}, GroupElementsToPrimes:=AssociativeArray() >;
  
  Gens := [];
  for p in BadPrimes(L : Even) do
    Spinors, g, ExactlyTheUnits := SpinorNorm(L, p);
    V:= Domain(g);
        
    // we only need to carry around those finite places where the Spinor norm is not exactly the units:
    if not ExactlyTheUnits then 
      vprintf Genus : "Found a prime over %o where the Spinor norm is not exactly the units of the order. dim(Spinors)=%o, dim(LocalMultGrpModSq)=%o\n", Minimum(p), Dimension(Spinors), Dimension(V);
      Include(~L`GenusComp`RayPrimes, p);
      // a basis of the Spinor norm of L at p, when viewed (modulo squares) as an F_2-vector space  (cf. LocalMultiplicativeGroupModSquares)
      b := Basis(Spinors); 
      Gens cat:= [<p, F!(a@g)>: a in b];
      assert &and[Valuation(gg, p) in {0,1}: gg in [a@g: a in b]];  // we rely on LocalMultiplicativeGroupModSquares mapping to 0- and 1-valued elements, and on the last base vector of V to represent a uniformizer of p

      I:= p^(1+2*Valuation(R!2, p));
      L`GenusComp`M *:= I;
      Append(~L`GenusComp`MM, I);
    end if;
  end for;
  // Now L`GenusComp`M contains a ray M and L`GenusComp`MM is the support of this ray.
  // We now compute the indefinite real places of L
  L`GenusComp`Inf:= {@ v: v in RealPlaces(F) | not IsIsotropic(L, v) @};
  // Now get the ray class group to M*I.
  I:= Sort([ i where _, i:= IsInfinite(v) : v in L`GenusComp`Inf ]);
  L`GenusComp`RCG, L`GenusComp`hom_RCG := RayClassGroup(L`GenusComp`M, I);
  
  // Step 1: map the generators into the class group to create the factor group.
  RCG_SubgroupGens := [ MapIdeleIntoClassGroup(L, [g]) : g in Gens ];
  L`GenusComp`RCG_Subgroup := sub<L`GenusComp`RCG | RCG_SubgroupGens, 2*L`GenusComp`RCG>;
  L`GenusComp`RCG_FactorGroup, L`GenusComp`hom_RCG_FactorGroup := quo<L`GenusComp`RCG | L`GenusComp`RCG_Subgroup >;

  vprintf Genus : "*** PrepareClassGroup: ray class group: size = %o\n", #L`GenusComp`RCG;
  vprintf Genus : "*** PrepareClassGroup: subgroup of ray class group: size = %o\n", #L`GenusComp`RCG_Subgroup;

  // Step 2: find good generators (that is: generators which are not dyadic, and not in BadPrimes(L) -- so that neighbour generation is always possible), 
  // of this factor group.
  // Only pick ideles of the form (1,...,1,p,1,...,1) so that:
  // a)  nothing has to be corrected before they can be mapped down into the factor group
  // b)  we can simply store the prime ideal and know that it corresponds to the (1,...,1,p,1,...,1) 
  // These prime ideals will be stored in L`GenusComp`CriticalPrimes.
  
  // the primes which generate the spinor genera in Genus(L):
  // if L is not definite, we need one prime for each factor group element.
  // if L is definite, we only need a generating set of primes for the factor group.
  
  L`GenusComp`CriticalPrimes := GetCriticalPrimes(L: FullList := true); // FullList := not IsDefinite(L));
  
  // vprintf Genus : "*** good primes over %o (together with the squares) generate the subgroup.\n", [Minimum(q): q in L`GenusComp`CriticalPrimes];
  // vprintf Genus : "*** PrepareClassGroup: factor group of ray class group: size = %o (this is the number of Spinor `+` genera in Genus(L))\n", #L`GenusComp`RCG_FactorGroup;
end intrinsic;

function SmallestNormGoodPrime(L)
  // return a smallest odd prime at which L is locally modular.
  // making neighbor generation possible (and hopefully fast!).
 
  K:= NumberField(L);
  Bad:= &*{ p : p in BadPrimes(L : Even) | Minimum(p) eq 2 or not IsModular(L, p) };
  limit:= 20;
  repeat
    PP:= PrimesUpTo(limit, K : coprime_to:= Bad);
    limit *:= 2;
  until #PP gt 0;
  
  return PP[1];
end function;

/*
// for use in for Genus enumeration:
intrinsic GetPrimes(L::LatMod) -> []
{A complete set of primes of BaseRing(L) needed to enumerate the genus of L.}
  require Rank(L) ge 3 : "The rank of L must be at least 3."; // otherwise we are unsure about the Spinor norms at unimodular primes.
  PrepareClassGroup(L);
  
  return L`GenusComp`CriticalPrimes;
end intrinsic;
*/

function SpinorGenerators(L : ModOut:= [])
  Bad:= { p: p in BadPrimes(L) | not IsModular(L, p) };
  S:= sub< L`GenusComp`RCG_FactorGroup | ModOut >;
  n:= #L`GenusComp`RCG_FactorGroup;
  Gens:= [];
  R:= BaseRing(L);
  p:= 2;		// omit even primes
  while #S ne n do
    p:= NextPrime(p);
    for P in [ d[1] : d in Decomposition(R, p) | d[1] notin Bad ] do
      _, h:= MapIdeleIntoClassGroup(L, [ < P, PrimitiveElement(P) > ]);
      if h notin S then
        S:= sub< L`GenusComp`RCG_FactorGroup | S, h >;
	Append(~Gens, P);
	if #S eq n then break; end if;
      end if;
    end for;
  end while;
  return Gens;
end function;

intrinsic IteratedNeighbours(L::LatMod, p::RngOrdIdl: UseAuto:= true, Max:= Infinity(), Mass:= 0) -> []
{The iterated neighbours of L at the prime p. }
  require IsDefinite(L) : "L must be definite"; // otherwise we cannot test for isometry

  Result := [ L ];
  if Mass gt 0 then
    Found:= 1/#AutomorphismGroup(L);
  end if;
  i:= 1;
  while (i le #Result) and (#Result lt Max) and (Mass le 0 or Found lt Mass) do
    // keep if not isometric, continue until the whole graph has been exhausted.
    call := func < Res, M | forall{ x: x in Res cat Result | not IsIsometric(M, x) }, true >;
    // The Max condition is checked here.
    // new CallBack function (February 2016).
    N:= Neighbours(Result[i], p: CallBack:=call, AutoOrbits:=UseAuto, Max:= Max - #Result);
    Result := Result cat N;
    if (Mass gt 0) and not IsEmpty(N) then
      Found +:= &+[ 1/#AutomorphismGroup(X) : X in N ];
      vprintf Genus : "Current batch: %o lattices, mass is %o out of %o (%o%%)\n",
        #Result, Found, Mass, RealField(Floor(Log(10,1+pc))+2)!pc where pc is 100*Found/Mass;
      assert Found le Mass;
    end if;
    i +:= 1;
  end while;
  return Result;
end intrinsic;

/*
function SpinorGeneraModuloSomething(L, g: CalculateCriticalPrimes:=true, GeneratorsOnly:=true);
  PrepareClassGroup(L);
  if ISA(Type(g), Setq) then
    error if #g ne 0 and Universe(g) cmpne L`GenusComp`RCG, "incompatible arguments, g is not an element of the class group";
  else
    error if not (Parent(g) eq L`GenusComp`RCG), "incompatible arguments, g is not an element of the class group";
  end if;
  
  // g is an element or a subset of L`GenusComp`RCG. Calculate the subgroup H generated by < L`GenusComp`RCG_Subgroup , g > .
  // Then the quotient Q := L`GenusComp`RCG / H is a quotient of L`GenusComp`RCG_FactorGroup. We return this quotient, and optionally return a set of critical primes for this quotient Q as well.

  NewRCG_FactorGroup, hom_NewRCG_FactorGroup := quo<L`GenusComp`RCG | L`GenusComp`RCG_Subgroup, g>;
  
  if not CalculateCriticalPrimes then return NewRCG_FactorGroup, hom_NewRCG_FactorGroup, _; end if;
  
  // TODO 20160111: oh god, redo this in a sane way ...
  // list of group elements which represent the factor group:
  to_find  := [(x@@hom_NewRCG_FactorGroup)@L`GenusComp`hom_RCG_FactorGroup : x in (GeneratorsOnly select Generators(NewRCG_FactorGroup) else NewRCG_FactorGroup)]; 
  Crit := [];

  // we use the mapping { Class group elements } -> { prime ideals } that we calculated in PrepareClassGroup:
  Crit := [L`GenusComp`GroupElementsToPrimes[g]: g in to_find];

  return NewRCG_FactorGroup, hom_NewRCG_FactorGroup, Crit;  
end function;
*/


intrinsic SpinorGeneraInGenus(L::LatMod : ModOut:= []) -> []
{A sequence of lattices representing the spinor genera in the genus of L}
  require Rank(L) ge 3 : "Only implemented for lattices of rank at least 3."; // otherwise the isomorphism to the class group fails, cf. §102 in O'Meara.

  // We have to find out whether ProperSpinorGenus == SpinorGenus.

  // 1) Find a nonzero element in the norm group of L.
  PrepareClassGroup(L); 
  G:= GramMatrix(L);
  if exists(d){d: d in Diagonal(G) | d ne 0} then
    spinornorm:= d;
  else
    ind:= Depth(G[1]);
    error if ind eq 0, "Lattice is degenerate!";
    // now Phi(v_1+v_ind, v_1+v_ind) = 2*Phi(v_1, v_ind)
    spinornorm := 2*G[1,ind];
  end if;

  // 2) at a place p where spinornorm does not generate norm(L_p)
  // at <p, spinornorm * normgenerator > to the Idele
  R:= BaseRing(L);
  F:= NumberField(L);
  norm:= Norm(L);
  Differ:= { p: p in Support(norm) join Support( ideal< R | spinornorm > ) | Valuation(norm, p) ne Valuation(spinornorm, p) };
  Idele:= [];
  for p in Differ join Set(L`GenusComp`RayPrimes) do
    if Minimum(p) eq 2 then
      _,_,_,a:= GenusSymbol(L, p);
      norm:= a[1];
    else
      G, pi:= GenusSymbol(L, p);
      norm:= pi^G[1,2];
      if G[1,1] eq 1 and G[1,3] eq -1 then
        k,h:= ResidueClassField(p);
        norm *:= Nonsquare(k) @@ h;
      end if;
    end if;
    Append(~Idele, <p, spinornorm * norm >);
  end for;
  
  // 3) Map the Idele into the IdeleClassGroup. 
  // Then h == 0 iff SpinorGenus == ProperSpinorGenus
  g, h := MapIdeleIntoClassGroup(L, Idele);

  // 4) Now there are 2^#Primes spinor genera, which can be reached by neighbours:
  Primes:= SpinorGenerators(L: ModOut:= [ h where _, h:= MapIdeleIntoClassGroup(L, [<p, PrimitiveElement(p)>]) : p in ModOut ] cat [h]);
  Result:= [ L ];
  for p in Primes do
    N:= Neighbours(L, p : Max:= 1, AutoOrbits:= false)[1];
    M:= Module(N);
    for i in [1..#Result] do
      Append(~Result, SetupLatticeInSpace(L, (p*MM + M) meet p^-1 * MM) where MM:= Module(Result[i]));
    end for;
  end for;
  return Result;
end intrinsic;


intrinsic GenusRepresentatives(L::LatMod : Max:= Infinity(), UseAuto, UseMass) -> []
{A list of representatives of the isometry classes in the genus of L}
  require Rank(L) ge 3 : "Only implemented for lattices of rank at least 3."; // otherwise the isomorphism to the class group fails, cf. §102 in O'Meara.
  require Max ge 1: "Max should be an integer or Infinity().";

  if not IsDefinite(L) then 
    return SpinorGeneraInGenus(L);
  end if;

  if UseMass then 
    UseAuto:= true;
    mass:= Mass(L);
  else
    mass:= -1;
  end if;

  Result:= [];
  p := SmallestNormGoodPrime(L);
  SpinorGenera:= SpinorGeneraInGenus(L: ModOut:= [p]);
  //if #SpinorGenera ne 1 then "here"; end if;

  for LL in SpinorGenera do
    assert forall{ X: X in Result | not IsIsometric(X, LL) };
    Result := Result cat IteratedNeighbours(LL, p: UseAuto:=UseAuto, Max:=Max - #Result, Mass:= mass / #SpinorGenera);
  end for;

  if Max gt #Result and UseMass then
    error if mass ne &+[1/#AutomorphismGroup(LL) : LL in Result ], "Something went wrong -- please report this example";
  end if;

  return Result; 
end intrinsic;

/*
intrinsic GenusRepresentatives(L::LatMod : Max:= Infinity(), UseAuto:= true, Extra:= 0, NoBreakOnError:=false, IndefiniteClassNumberOnly:=false) -> []
{A list of representatives of the isometry classes in the genus of L. Optional parameters: Max, the maximum number of isometry classes to calculate; Extra, the number of extra primes for safety}

  require Rank(L) ge 3 : "Only implemented for lattices of rank at least 3."; // otherwise the isomorphism to the class group fails, cf. §102 in O'Meara.
  requirege Extra, 0;
  require (Max ge 1): "Max should be an integer or Infinity().";
  
  
  // There are infinitely many primes generating our ray class group for L,
  // so there is no need to choose dyadic primes. Neighbour generation
  // for dyadic primes is a pain and should be avoided.
  
  F := NumberField(L);

  if IsDefinite(L) then
      // our favorite prime of small norm:
      good := SmallestNormGoodPrime(L);
      
      // Let's see if the favorite prime switches spinor genera too:
      g, h := MapIdeleIntoClassGroup(L, [<good, F!PrimitiveElement(good)>]);
      if not IsIdentity(h) then
        _,_,Primes := SpinorGeneraModuloSomething(L, g);
        // at this point, 'good' will not be an element of Primes.
        vprintf Genus : "We have %o critical primes (modulo our good prime).\n", #Primes;
      else
        Primes := {x: x in GetPrimes(L)};
        // here we are not sure whether 'good' is also an element of primes, and we can remove it to save some time:
        if good in Primes then
          /// error "why is this -- we thought h was the identity of our class group quotient?";
          // TODO 20160111: das liegt daran dass wir die CriticalPrimes mit FullList:=true berechnen und deswegen auch einen Repräsentanten für die 1 berechnen ._.
          Primes diff:= {good};
          vprintf Genus : "We have %o critical primes, after removing the good prime from this set.\n", #Primes; 
        else
          vprintf Genus : "We have %o critical primes.\n", #Primes;
        end if;
      end if;
      error if not &and[Minimum(q) gt 2: q in Primes], "PrepareClassGroup should not pick dyadic primes.";
      error if &or[p in BadPrimes(L): p in Primes], "PrepareClassGroup should not pick any primes where the lattice is not modular.";
    
      if Extra cmpne 0 then   
        if Type(Extra) eq RngIntElt then Extra:= GetExtra(L, Primes, Extra); end if;
        Mode:= [false, true];
      else
        Mode:= [false];
      end if;
      ExtraFail:= false;

      vprint Genus: "Calling iterated neighbours with prime over", Minimum(good), "of norm", Norm(good);
      Result:= IteratedNeighbours(L, good: UseAuto:=UseAuto, Max:=Max);
      vprint Genus: "Starting with", #Result, "lattices.";
      Indices:= [1];		// where new spinor genera start!
      for extra in Mode do
        i:= 1;
        while i le #Indices and #Result lt Max do
          for p in Primes do
            N:= Neighbours(Result[Indices[i]], p: Max:=1, AutoOrbits:= false)[1];
            vprint Genus: "Testing lattice", i, "with prime over", Minimum(p);
            if exists{L : L in Result | IsIsometric(L, N)} then continue; end if;
            if extra then
              ExtraFail:= true;
              error if not NoBreakOnError, "Oops, found new rep in the extra rounds! --- Please report";
            end if;
            vprint Genus: "New spinor genus found. Calling iterated neighbours again.";
            Append(~Indices, #Result + 1);
            Result cat:= IteratedNeighbours(L, good: UseAuto:=UseAuto, Max:=Max - #Result);
            if #Result ge Max then break; end if;
            vprint Genus: "New number of lattices:", #Result;
          end for;
          i +:= 1;
        end while;
      end for;
      if #Mode eq 2 then return Result, ExtraFail; end if;
      return Result;

  
  else
      // The indefinite case.
      // What makes this case easier: proper spinor genera and proper classes coincide (when rank(L) >= 3, cf. 104:5 in O'Meara). 
      // What makes this case harder: Aut(L) is infinite. We cannot test for isometry of global lattices.
      // How do we still achieve a system of representatives for the classes in Genus(L), and avoid hitting a class twice if it decomposes into two proper classes?: We explicitly calculate the group element corresponding to the improper automorphism in O(V) that switches two proper classes, and kick it out. The group element can be calculated by explicitly mapping the spinor norm of the reflection at any basis vector of L into the class group.
      
      G := GramMatrix(L);
      A, T := OrthogonalizeGram(G); // now T*G*Transpose(T) = A
      //M :=  DiagonalMatrix(BaseRing(T), [1: i in [1..Nrows(T)]]);
      // M[1,1] := -1; 

      ind := [i: i in [1..Nrows(G)] | not IsZero(G[i][i])];
      if #ind eq 0 then
        ind := [i: i in [2..Ncols(G)] | not IsZero(G[1][i])];
        error if #ind eq 0, "Lattice does not have full rank.";
        // now Phi(v_1+v_i, v_1+v_i) = 2*Phi(v_1, v_i), for any i in ind.
        spinornorm := 2*G[1][ind[1]];
      else
        spinornorm := G[ind[1],ind[1]];
      end if;
      // now spinornorm is a representative of the spinor norm of the reflection at x_i for some i, or at (x_1+x_i) for some i.
            
      
      h, g := MapPrincipalIdeleIntoClassGroup(L, F!spinornorm);
      // g := h@L`GenusComp`hom_RCG_FactorGroup;
      
      /// assert IsIdentity(g@L`GenusComp`hom_RCG_FactorGroup) eq (g in L`GenusComp`RCG_Subgroup);
      
     
      if not IsIdentity(g) then // equally: if not IsIdentity(g@L`GenusComp`hom_RCG_FactorGroup)
        // We factor out a group element that switches proper spinor genera (= proper classes in the indefinite case);
        // but only switches to a different proper spinor genus inside the same isometry class.
        //printf "Classes and proper classes do not coincide in Genus(L), i.e. there is no automorphism of L with determinant -1.\n";
        //printf "Adding an element to RCG_Subgroup ...";
        MyRCG_FactorGroup, Myhom_RCG_FactorGroup, Crit := SpinorGeneraModuloSomething(L, h: CalculateCriticalPrimes := not IndefiniteClassNumberOnly, GeneratorsOnly:=false);
        
        if IndefiniteClassNumberOnly then return #MyRCG_FactorGroup; end if;
      else
        //printf "Classes and proper classes coincide in Genus(L), i.e. there is an automorphism of L with determinant -1.\n";
        MySubgroup := L`GenusComp`RCG;
        MyRCG_FactorGroup := L`GenusComp`RCG_FactorGroup; Myhom_RCG_FactorGroup := L`GenusComp`hom_RCG_FactorGroup;
        if IndefiniteClassNumberOnly then return #MyRCG_FactorGroup; end if;
      end if;  

      // remove the representative of the trivial element:
      Crit := [p: p in L`GenusComp`CriticalPrimes | not IsIdentity(MapIdeleIntoClassGroup(L, [<p, F!PrimitiveElement(p)>])@L`GenusComp`hom_RCG_FactorGroup)];
      
      vprintf Genus : "Number of elements: RCG: %o, Subgroup: %o, RCG_FactorGroup: %o, MySubgroup: %o, MyRCG_FactorGroup: %o, Crit: %o\n", #L`GenusComp`RCG, #L`GenusComp`RCG_Subgroup, #L`GenusComp`RCG_FactorGroup, #MySubgroup, #MyRCG_FactorGroup, #Crit;
      
      //ClassGroupElements := [x: x in MyRCG_FactorGroup];
     
      Reps:= [];
      for p in Crit do
        // one neighbour each at all of the critical primes:
        // We now only need to calculate one neighbour for each critical
        // prime and have a complete set of representatives of the improper(!) isometry classes in Genus(L).
        /// call_one:= func< LL | true, false >;
        N:= Neighbours(L, p: Max:=1, AutoOrbits:= false);
        assert #N gt 0;
        Append(~Reps, N[1]);
      end for;
      return Reps;
  end if;
end intrinsic;
*/

intrinsic IsOneClassGenus(L::LatMod: UseAuto, UseMass) -> BoolElt
{Returns (#GenusRepresentatives(L) eq 1), but is much faster than calling GenusRepresentatives}
  if IsDefinite(L) then
    if UseMass then
      return Mass(L) eq 1/#AutomorphismGroup(L);
    end if;
    if #SpinorGeneraInGenus(L) ne 1 then return false; end if;
    p:= SmallestNormGoodPrime(L);
    N:= Neighbours(L, p : AutoOrbits:= UseAuto, CallBack:= func< Res, LL | not same, same where same:= IsIsometric(L, LL) >, Max:= 1);
    return IsEmpty(N);
  else
    return #SpinorGeneraInGenus(L) eq 1;
  end if;
end intrinsic;
