intrinsic get_quotient_split(E, P, p, i) -> Any
  {}
  OE := MaximalOrder(E);
  OK := Order(p);
  OEabs := AbsoluteOrder(OE);
  Eabs := NumberField(OEabs);

  RK, mRK := quo<OK | p^i>;
  RKu, mRKu := UnitGroup(RK);
  Pabs := OEabs!!P;
  REabs, mREabs := quo<OEabs | Pabs^i>;
  REabsu, mREabsu := UnitGroup(REabs);
  // We have to be a bit careful with denominators
  dlog := function(x)
    d := Denominator(x);
    xabs := d * Eabs!x;
    dabs := d;
    F := Factorization(Minimum(P) * OEabs);
    uni := [];
    for P in F do
      assert Valuation(Eabs!x, P[1]) ge 0;
      api := anti_uniformizer(P[1]);
      exp := Valuation(OEabs!d, P[1]);
      dabs := dabs * api^exp;
      xabs := xabs * api^exp;
    end for;

    xabs_image := ((OEabs!(xabs))@mREabs@@mREabsu);
    dabs_image := ((OEabs!(dabs))@mREabs@@mREabsu);

    ret := (xabs_image - dabs_image);

    if not (OE!d in P) then
      xd := d * x;
      xd_image := ((OEabs!(OE!xd))@mREabs@@mREabsu);
      d_image := ((OEabs!(OE!d))@mREabs@@mREabsu);
      ret2 :=  (xd_image - d_image);
      assert ret eq ret2;
    end if;

    return ret;
  end function;

  exp := function(k)
    x := E!(OE!(k@mREabsu@@mREabs));
    assert dlog(x) eq k;
    return x;
  end function;

  return REabsu, exp, dlog;
end intrinsic;

// The following functions give an abelian group K and two maps f : K -> E, g : E -> K, and K is isomorphic to E/E^i.
//
// Example:
//
// > Attach("quotients.m");
// > K<a> := NumberField(x^2 - 2); 
// > Qx<x> := PolynomialRing(Rationals());
// > K<a> := NumberField(x^2 - 2);
// > Kt<t> := PolynomialRing(K);
// > E<b> := NumberField(t^2 - a);
// > p := Decomposition(MaximalOrder(K), 2)[1][1];
// > P := Decomposition(MaximalOrder(E), p)[1][1];
// > KK, f, g := get_quotient_ramified(E, P, p, 5);
// > KK;    
// Abelian Group isomorphic to Z/2
// Defined on 1 generator
// Relations:
//     2*KK.1 = 0
// > f(KK.1);
// [[2, 1], [3, 3]]
// > g(f(KK.1));
// KK.1
// 

intrinsic get_quotient_inert(E, P, p, i) -> Any
  {}
  //p := Factorization(Norm(P))[1][1]; // Is there another way to do this?
  OE := Order(P);
  OK := Order(p);

  OEabs := AbsoluteOrder(OE);
  Eabs := NumberField(OEabs);

  RK, mRK := quo<OK | p^i>;
  RKu, mRKu := UnitGroup(RK);
  Pabs := OEabs!!P;
  REabs, mREabs := quo<OEabs | Pabs^i>;
  REabsu, mREabsu := UnitGroup(REabs);
  //f = hom<REabsu -> RKu | [ mRKu\(RK(OK(norm(EabstoE(elem_in_nf(mREabs\(mREabsu(REabsu[i])))))))) for i in 1:ngens(REabsu)])
  f := hom<REabsu -> RKu | [(Norm(OE!(REabsu.i@mREabsu@@mREabs)))@mRK@@mRKu : i in [1..Ngens(REabsu)]]>;
  K := Kernel(f);

  dlog := function(x)
    d := Denominator(x);
    xabs := d * Eabs!x;
    dabs := d;
    F := Factorization(Minimum(P) * OEabs);
    uni := [];
    for P in F do
      assert Valuation(Eabs!x, P[1]) ge 0;
      api := anti_uniformizer(P[1]);
      exp := Valuation(OEabs!d, P[1]);
      dabs := dabs * api^exp;
      xabs := xabs * api^exp;
    end for;

    xabs_image := ((OEabs!(xabs))@mREabs@@mREabsu);
    dabs_image := ((OEabs!(dabs))@mREabs@@mREabsu);

    ret := K!(xabs_image - dabs_image);

    if not (OE!d in P) then
      xd := d * x;
      xd_image := ((OEabs!(OE!xd))@mREabs@@mREabsu);
      d_image := ((OEabs!(OE!d))@mREabs@@mREabsu);
      ret2 :=  K!(xd_image - d_image);
      assert ret eq ret2;
    end if;

    return ret;
  end function;

  exp := function(k)
    x := E!(OE!((REabsu!k)@mREabsu@@mREabs));
    assert dlog(x) eq k;
    return x;
  end function;

  return K, exp, dlog;
end intrinsic;

intrinsic get_quotient_ramified(E, P, p, i) -> Any
  {}

  //p := Factorization(Norm(P))[1][1]; // Is there another way to do this?
  OE := Order(P);
  OK := Order(p);

  OEabs := AbsoluteOrder(OE);
  e := Valuation(Different(OE), P);
  t := e - 1;
  assert e ge 1;

  if i lt e then
    S := AbelianGroup([]);
    return S, func< k | One(E)>, func< x | Identity(S)>;
  end if;

  psi := func< x | t + 2 * (x - t)>;

  jj := (e - 1) + 1/2;
  while Ceiling(psi(jj)) ne i do
    jj := jj + 1/2;
  end while;
  j := Ceiling(jj);

  OEabs := AbsoluteOrder(OE);
  Eabs := NumberField(OEabs);

  RK, mRK := quo<OK | p^j>;
  RKu, mRKu := UnitGroup(RK);
  Pabs := OEabs!!P;
  REabs, mREabs := quo<OEabs | Pabs^i>;
  REabsu, mREabsu := UnitGroup(REabs);
  //f = hom<REabsu -> RKu | [ mRKu\(RK(OK(norm(EabstoE(elem_in_nf(mREabs\(mREabsu(REabsu[i])))))))) for i in 1:ngens(REabsu)])
  f := hom<REabsu -> RKu | [(Norm(OE!(REabsu.i@mREabsu@@mREabs)))@mRK@@mRKu : i in [1..Ngens(REabsu)]]>;
  K := Kernel(f);
  pi := UniformizingElement(P);
  dlog := function(x)
    d := Denominator(x);
    xabs := d * Eabs!x;
    dabs := d;
    F := Factorization(Minimum(P) * OEabs);
    uni := [];
    for P in F do
      assert Valuation(Eabs!x, P[1]) ge 0;
      api := anti_uniformizer(P[1]);
      exp := Valuation(OEabs!d, P[1]);
      dabs := dabs * api^exp;
      xabs := xabs * api^exp;
    end for;

    xabs_image := ((OEabs!(xabs))@mREabs@@mREabsu);
    dabs_image := ((OEabs!(dabs))@mREabs@@mREabsu);

    ret := K!(xabs_image - dabs_image);

    if not (OE!d in P) then
      xd := d * x;
      xd_image := ((OEabs!(OE!xd))@mREabs@@mREabsu);
      d_image := ((OEabs!(OE!d))@mREabs@@mREabsu);
      ret2 :=  K!(xd_image - d_image);
      assert ret eq ret2;
    end if;

    return ret;
  end function;

  if i eq e then
    assert #K eq 2;
  end if;

  exp := function(k)
    x := E!(OE!((REabsu!k)@mREabsu@@mREabs));
    assert dlog(x) eq k;
    return x;
  end function;

  for l in [1..10] do
    g := Random(K);
    assert dlog(exp(g)) eq g;
  end for;

  return K, exp, dlog;
end intrinsic;

intrinsic get_quotient(E, P, p, e) -> Any
  {}
  OE := MaximalOrder(E);
  l := Decomposition(OE, p);
  if #l eq 1 and l[1][2] eq 1 then // inert
    P := l[1][1];
    KK, f, g := get_quotient_inert(E, P, p, e);
  elif #l eq 1 and l[1][2] eq 2 then // ramified
    P := l[1][1];
    KK, f, g := get_quotient_ramified(E, P, p, e);
  else
    assert #l eq 2 and l[1][2] eq 1; // split
    P := l[1][1];
    KK, f, g := get_quotient_split(E, P, p, e);
  end if;
  return KK, f, g;
end intrinsic;

// Get \prod_p E_0,p/E_0^i,p given [<p, i>]

intrinsic get_big_quotient(E, Fac) -> Any
  {}
  OE := MaximalOrder(E);
  groups := [];
  exps := [];
  dlogs := [];
  Ps := []; 
  for i in [1..#Fac] do
    p := Fac[i][1];
    e := Fac[i][2];
    l := Decomposition(OE, p);
    if #l eq 1 and l[1][2] eq 1 then // inert
      P := l[1][1];
      KK, f, g := get_quotient_inert(E, P, p, e);
    elif #l eq 1 and l[1][2] eq 2 then // ramified
      P := l[1][1];
      KK, f, g := get_quotient_ramified(E, P, p, e);
    else
      assert #l eq 2 and l[1][2] eq 1; // split
      P := l[1][1];
      KK, f, g := get_quotient_split(E, P, p, e);
    end if;
    Append(~groups, KK);
    Append(~exps, f);
    Append(~dlogs, g);
    Append(~Ps, P);
  end for;

  if #groups eq 0 then
    A := AbelianGroup([]);
    dlog := func< x | Identity(A)>;
    exp := func< a | One(E)>;
    return A, dlog, exp;
  end if;

  if #groups eq 1 then
    return groups[1], dlogs[1], exps[1];
  end if;

  G, inj, proj := DirectSum(groups);

  dlog := function(x)
    return &+[x@(dlogs[i])@inj[i] : i in [1..#Fac]];
  end function;

  exp := function(k)
    v := [ OE!(k@(proj[i])@exps[i]) : i in [1..#Fac]];
    a := ChineseRemainderTheorem(v, [Ps[i]^(3*Fac[i][2]) : i in [1..#Fac]]);
    // There might be a problem with the precision here.
    assert dlog(a) eq k;
    return a;
  end function;

  for k in [1..10] do
    g := Random(G);
    assert dlog(exp(g)) eq g;
  end for;

  return G, dlog, exp;
end intrinsic;

intrinsic anti_uniformizer(P) -> Any
  {}
  iP := P^-1;
  O := Order(P);
  d := Denominator(iP);
  I := O!!(d * iP);
  while true do
    x := Random(iP, 10);
    if not IsIntegral(x) then
      return x;
    end if;
  end while;
end intrinsic;

intrinsic get_the_group(L) -> Any
  {}
  // First get the ideal [D^-1 L^# / L ]
  LD := Dual(L);
  OE := BaseRing(L);
  E := NumberField(OE);
  K := BaseField(E);
  OK := MaximalOrder(K);
  if IsAbsoluteField(E) then
    D := Different(MaximalOrder(E));
  else
    D := Different(MaximalOrder(E)) * Different(MaximalOrder(K));
  end if;
  Drel := Different(MaximalOrder(E));
  I := &*ElementaryDivisors(Module(D^-1*LD), Module(L));
  Fac := Factorization(I);
  for i in [1..#Fac] do
    if Fac[i][2] lt 0 then
      error("L <= A^-1 L^# violated");
    end if;
  end for;
  S := Set([Factorization(Norm(Fac[i][1]))[1][1] : i in [1..#Fac]]);
  NL := Norm(L);
  // Compute the invariants to determine F^#(L)
  Fsharpdata := [];
  for p in S do
    lp := Decomposition(OE, p);
    P := lp[1][1];
    n := Valuation(NL, P);
    e := Valuation(D, P);
    // The det group is E_0^(v_P(n(L)) + e)
    Append(~Fsharpdata, <p, n + e>);
    //assert Valuation(OE!!(p^n), P) eq Valuation(NL, P);
  end for;

  RmodFsharp, Fsharplog, Fsharpexp := get_big_quotient(E, Fsharpdata);

  Fdata := [];
  for p in S do
    lp := Decomposition(OE, p);
    if not is_special(L, p) then
      continue;
    end if;
    e := Valuation(Drel, P);
    // The det group E_0^(2*n + e)
    Append(~Fdata, <p, e>);
  end for;

  RmodF, Flog, Fexp := get_big_quotient(E, Fdata);

  // Compute f : R/F^#(L) -> R/F(L)

  A := [ Flog(Fsharpexp(RmodFsharp.i)) : i in [1..Ngens(RmodFsharp)]];
  f := hom<RmodFsharp -> RmodF | A>;
  K := Kernel(f);
  FmodFsharp := K;
  
  // "R/F^#(L) of order", #RmodFsharp;
  // "R/F(L) of order", #RmodF;
  // "F(L)/F^#(L) of order", #K;

  // Now Compute F(O_E):

  OEabs := AbsoluteOrder(OE);
  OEabsU, mOEabsU := UnitGroup(OEabs);
  OKU, mOKU := UnitGroup(OK);

  f := hom<OEabsU -> OKU | [ (Norm(OE!((OKU.i@mOKU))))@@mOKU : i in [1..Ngens(OKU)]]>;
  K := Kernel(f);

  gens_norm_one := [ OE!(K.i@mOEabsU) : i in [1..Ngens(K)]];

  FOEmodFsharpL := sub<RmodFsharp | [Fsharplog(x) : x in gens_norm_one]>;

  I := FOEmodFsharpL meet FmodFsharp;
  Q, mQ := quo< FmodFsharp | I>;
  S := InvariantRepresentation(Q);
  assert IsDiagonal(RelationMatrix(S));
  IF := InvariantFactors(S);
  
  dlog := function(x)
    return Eltseq(S!((FmodFsharp!Fsharplog(x))@mQ));
  end function;

  return Q, IF, dlog;
end intrinsic;

intrinsic is_special(L, p) -> Any
  {}
  OE := BaseRing(L);
  E := NumberField(OE);
  lp := Decomposition(OE, p);
  if lp[1][2] ne 2 or IsOdd(Rank(L)) then
    // unramified or odd rank
    return false;
  end if;
  _,_R,S := JordanDecomposition(L, p);
  R := [Nrows(x) : x in _R]; 
  for r in R do
    if r ne 2 then
      return false;
    end if;
  end for;
  P := lp[1][1];
  e := Valuation(Different(OE), P);
  for i in [1..#S] do
    if (e - S[i]) mod 2 ne 0 then
      return false;
    end if;
  end for;
  // Now L has splitting type [2,...,2], [s1,...,sr] and si = e mod 2
  pi := UniformizingElement(P);
  A := Automorphisms(E);
  assert A[2](E.1) ne E.1;
  pico := A[2](pi);
  H := DiagonalJoin([Matrix(E, 2, 2, [0, pi^S[i], pico^S[i], 0]) : i in [1..#S]]);
  return IsLocallyIsometric(L, HermitianLattice(H), p);
end intrinsic;
