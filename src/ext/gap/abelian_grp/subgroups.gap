
OnSubgroups:=function(x,g)
# indeed this is a right action.
return Image(g,x);
end;


ElemAbelSubgrpOrbsStabs:=function(G, aut, gens_aut, gens_act)
  # Return Orbits and Stabilizers of the subgroups of G modulo aut.
  # INPUT:
  # G -- an elementary abelian group
  # aut -- a group
  # gens_aut -- a list of generators
  # gens_act -- a list of elements of Aut(G) defining the action
  #             of aut on G.
  # gens_aut and gens_act must have the same length and define a group homomorphism
  # no input checks
  # OUTPUT:
  # A list consisting of records with entries
  # repr  -- a subgroup of G
  # stab  -- a subgroup of aut

  local subs,i,k, n, p, gens_mats, MatrixRepresentation, V, orbits, Gmat,
  subgrp, stab, subgroups, pcgs, orb;

  pcgs:= Pcgs(G); # a basis
  p:= Order(pcgs[1]); # a prime number
  # we compute with respect to a given pcgs
  MatrixRepresentation:= function(f, pcgs)
    return List(pcgs, i -> ExponentsOfPcElement(pcgs, Image(f,i)))*One(GF(p));
  end;

  gens_mats:= List(gens_act, f->MatrixRepresentation(f, pcgs));

  Gmat:= Group(gens_mats);
  n:= Size(pcgs);
  V:= GF(p)^n;
  subgroups:= [rec(repr:=Subgroup(G,[]), stab:=aut)];
  for k in [1..n] do
    subs:= List(Subspaces(V, k), i->Basis(i)); # is this the fastest way?
    orbits:= ExternalOrbitsStabilisers(aut, subs, gens_aut, gens_mats, OnSubspacesByCanonicalBasis);
    # transform orbit reps back to subgroups of G
    for orb in orbits do
      subgrp:= Subgroup(G, List(Representative(orb),i -> PcElementByExponents(pcgs,i)));
      stab:= StabilizerOfExternalSet(orb);
      Add(subgroups, rec(repr:=subgrp, stab:=stab) );
   od;
  od;
  return subgroups;
end;


SubgroupReps1:=function(epi, aut, gens_aut, gens_act)
  # INPUT:
  # epi: G0 --> G1 an epimorphism of abelian groups
  # aut < Aut(G0) a subgroup
  # the action is defined by the homomorphism
  # gens_aut --> gens_act where gens_act lies in in Aut(G1)
  #
  # OUTPUT:
  # A list consisting of records with entries
  # repr  -- a subgroup of G
  # stab  -- a subgroup of aut
  # We follow an algorithm of Hulpke and (mostly) his notation in Lemma 3.1
  local N, G, GmodN, subgrps_rep, complement_reps, aut_on_GmodN,
      Alist, Blist, U, A, B, i, invs, min, UmodBList,epi_mod,epi_new, G0, G1,epi_G0_B,
      stabA, StabB, AmodB, NmodB, epiB, C_AB, Blist_A , act,gens_C_AB,act_C_AB;

  # gens_aut:= GeneratorsOfGroup(aut);
  G0:=Source(epi);
  G1:=Image(epi);

  N:=InvariantElementaryAbelianSeries(G1, gens_act, G1, true);  # true means the series is fine
  N:=N[Size(N)-1];   # an invariant elementary abelian group
  # is it minimal? Not necessarily.
  # remedy this here
  invs:=InvariantSubgroupsElementaryAbelianGroup(N, gens_act);  # brute force
  invs:=Filtered(invs, IsNonTrivial);
  min:=Minimum(List(invs, i->Size(i)));
  for N in invs do
    if Size(N)=min then
      break;
    fi;
  od;
  # now N is minimal
  Blist:=ElemAbelSubgrpOrbsStabs(N, aut, gens_aut, gens_act);
  subgrps_rep:= Blist;

  if Size(N) = Size(G1) then
    return subgrps_rep;
  fi;

  epi_mod:=NaturalHomomorphismByNormalSubgroupNC(G1, N);
  # LockNaturalHomomorphismsPool(G,epi_mod); #what does this do?
  GmodN:=Image(epi_mod);
  aut_on_GmodN:= List(gens_act, i->InducedAutomorphism(epi_mod, i));
  epi_new:=CompositionMapping(epi_mod,epi);
  Alist:= SubgroupReps1(epi_new, aut, gens_aut, aut_on_GmodN);
  for A in Alist do
    A.repr:= PreImage(epi_mod, A.repr);
  od;
  # remove N
  Assert(0, Size(Alist[1].repr)=Size(N));
  Remove(Alist, 1);
  subgrps_rep:= Concatenation(subgrps_rep, Alist);

  for A in Alist do
    stabA:= A.stab;
    A:= A.repr;
    # could this be handled quicker by Blist and fusion? How?
    act:= List(GeneratorsOfGroup(stabA), i->InducedAutomorphism(epi,i));
    Blist_A:=ElemAbelSubgrpOrbsStabs(N, stabA, GeneratorsOfGroup(stabA), act );
    # we want only proper subgroups -> remove N
    Assert(0, Size(Blist_A[Size(Blist_A)].repr)=Size(N)); # asser the last entry is N
    Remove(Blist_A, Size(Blist_A));
    for B in Blist_A do
      C_AB:= B.stab;
      B:= B.repr;
      epiB:= NaturalHomomorphismByNormalSubgroupNC(A, B);
      AmodB:= Image(epiB);
      NmodB:= Image(epiB, N);
      # TODO we only want classes of complements up to the action
      # can we work directly with H^1 first and
      # take orbits there?
      UmodBList:= ComplementClassesRepresentatives(AmodB, NmodB);
      gens_C_AB:=GeneratorsOfGroup(C_AB);
      act_C_AB:= List(gens_C_AB, i->InducedAutomorphism(epi, i));
      act_C_AB:= List(act_C_AB, i->InducedAutomorphism(epiB, i));
      complement_reps:=ExternalOrbitsStabilisers(C_AB,UmodBList, gens_C_AB, act_C_AB, OnSubgroups);
      complement_reps:=List(complement_reps, i->rec(repr:=PreImage(epiB, Representative(i)), stab:=StabilizerOfExternalSet(i)));
      subgrps_rep:= Concatenation(subgrps_rep, complement_reps);
    od;
  od;
  # sort subgroups by their order
  Sort(subgrps_rep, function(a,b) return Size(a.repr) < Size(b.repr); end);
  return subgrps_rep;
end;

SubgroupReps:=function(G, aut)
  # Compute representatives and stabilizers of all subgroups
  # {B < G} / aut
  # up to the action of aut
  # sorted by their order
  #
  # INPUT:
  # - G -- an abelian group
  # - aut -- a subgroup of AutomorphismGroup(G)
  #
  # OUTPUT:
  # A list consisting of records with entries
  # repr  -- a subgroup of G
  # stab  -- a subgroup of aut
  #
  local gens_aut;
  gens_aut:=GeneratorsOfGroup(aut);
  return SubgroupReps1(gens_aut[1]^0, aut, gens_aut, gens_aut);
end;
