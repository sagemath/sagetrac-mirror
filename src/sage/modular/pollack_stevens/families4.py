class WeightKAction_fam(WeightKAction_vector):
    """
    """
    def __init__(self, Dk, character, tuplegen, on_left):
        #Only difference is that it adds a cache for automorphy factors.
        self._autfactors = {}
        WeightKAction_vector.__init___(self, Dk, character, tuplegen, on_left)
    
    def clear_cache(self):
        #Only difference is that it clears the cache for automorphy factors.
        self._actmat = {}
        self._maxprecs = {}
        self._autfactors = {}
    
    def _compute_aut_factor_matrix(self, g, M):
        #compute the power series
        D = self.underlying_set()
        p_prec, var_prec = D.precision_cap()
        return automorphy_factor_matrix(D.prime(), g[0, 0], g[1,0], self._k, self._character, M, var_prec, D.base_ring())
    
    def get_action_matrices(self, g, M):
        g = M2Z(g)
        g.set_immutable()
        if not self._maxprecs.has_key(g):
            AF = self._compute_aut_factor_matrix(g, M)
            self._autfactors[g] = {M : AF}
        else:
            auts = self._autfactors[g]
            if auts.has_key(M):
                AF = auts[M]
            else:
                maxprec = self._maxprecs[g]
                if M < maxprec:
                    AF = auts[maxprec][:M,:M]
                    auts[M] = AF
                else:
                    if M < 2 * maxprec:
                        maxprec = 2 * maxprec
                    else:
                        maxprec = M
                    auts[maxprec] = self._compute_aut_factor_matrix(g, maxprec)
                    if M == maxprec:
                        AF = auts[maxprec]
                    else:
                        AF = auts[maxprec][:M,:M]
                        auts[M] = AF
        
        A = self.acting_matrix(g, M)
        return [AF, A]
    
    def _call_(self, v, g):
        AF, A = get_action_matrices(self, g, len(v.moments))
        return v.parent()((v.moments * AF) * A)
