from __future__ import absolute_import




def isom_fqf_brute_force(A, B):
    r"""
    """
    assert A.invariants() == B.invariants()
    na = len(A.smith_form_gens())
    nb = len(B.smith_form_gens())

    b_cand = [[b for b in B if b.q()==a.q() and b.order() == a.order()] for a in A.smith_form_gens()]

    res = []
    waiting = [[]]
    while len(waiting) > 0:
        f = waiting.pop()
        i = len(f)
        if i == na:
            res.append(f)
            continue
        a = A.smith_form_gens()[i]
        for b in b_cand[i]:
            if all(b.b(f[k])==a.b(A.smith_form_gens()[k]) for k in range(i)):
                waiting.append(f + [b])
    return res










