using Hecke

# Not clear how we represent E_p
function get_quotient_split(E::Hecke.NfRel, p, i::Int)
  Eabs, EabstoE = absolute_field(E, cached = false)
  pabs = (EabstoE\p)
  OEabs = maximal_order(Ebas)
  R, mR = quo(Eabs, pabs^i)
  # not done yet
end

function get_quotient_unramified(E::Hecke.NfRel, P, i::Int)
  K = base_field(E)
  OE = order(P)
  Eabs, EabstoE = absolute_field(E, cached = false)
  p = minimum(P)
  OK = order(p)
  RK, mRK = quo(OK, p^i)
  RKu, mRKu = unit_group(RK)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  REabs, mREabs = quo(OEabs, Pabs^i)
  REabsu, mREabsu = unit_group(REabs)
  f = hom(REabsu, RKu, [ mRKu\(RK(OK(norm(EabstoE(elem_in_nf(mREabs\(mREabsu(REabsu[i])))))))) for i in 1:ngens(REabsu)])

  K, mK = kernel(f)

  S, mS = snf(K)

  # exp
  exp = function(s)
    @assert parent(s) == S
    return EabstoE(elem_in_nf(mREabs\(mREabsu(mK(mS(s))))))
  end

  log = function(x)
    @assert parent(x) == E
    @assert isintegral(x)
    return mS\(mK\(mREabsu\(mREabs(OEabs(EabstoE\x)))))
  end


  return S, exp, log
end

function _hilbert_90(a)
  @assert isone(norm(a))
  E = parent(a)
  A = automorphisms(E)
  s = A[1](gen(E)) == gen(E) ? A[2] : A[1]
  B = 10
  theta = rand(E, -B:B)
  while true
    b = theta + a * s(theta) 
    if !iszero(b)
      break
    end
  end
  @assert a == b//s(b)
  return b
end

function get_quotient_ramified(E::Hecke.NfRel, P, i::Int)
  e = valuation(different(maximal_order(E)), P)
  p = minimum(P)
  if i < e
    S = abelian_group(Int[])
    return S, x -> one(E), x -> id(S)
  end

  t = e - 1

  psi(x) = x <= t ? x : t + 2(x - t)

  pi = uniformizer(P)

  #if i == e
  #  S = abelian_group([2])

  #  exp = function(s)
  #    error("Not implemented yet")
  #  end

  #  log = function(x)
  #    b = _hilbert_90(a)
  #    v = valuation(x, P)
  #    return isodd(v) ? S[1] : id(S)
  #  end

  #  return S, exp, log
  #end

  # We find jj with ceil(psi(jj)) = i
  # then j = ceil(jj)
  jj = (e - 1) + 1//2
  while ceil(psi(jj)) != i
    jj += 1//2
  end
  j = Int(ceil(jj))


  Pi = P^i

  OE = order(P)
  Eabs, EabstoE = absolute_field(E, cached = false)
  p = minimum(P)
  OK = order(p)
  RK, mRK = quo(OK, p^j)
  RKu, mRKu = unit_group(RK)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  REabs, mREabs = quo(OEabs, Pabs^i)
  REabsu, mREabsu = unit_group(REabs)

  f = hom(REabsu, RKu, [ mRKu\(RK(OK(norm(EabstoE(elem_in_nf(mREabs\(mREabsu(REabsu[i])))))))) for i in 1:ngens(REabsu)])

  K, mK = kernel(f)

  S, mS = snf(K)

  # exp
  exp = function(s)
    @assert parent(s) == S
    return EabstoE(elem_in_nf(mREabs\(mREabsu(mK(mS(s))))))
  end

  log = function(x)
    @assert parent(x) == E
    @assert isintegral(x)
    return mS\(mK\(mREabsu\(mREabs(OEabs(EabstoE\x)))))
  end


  return S, exp, log
end
