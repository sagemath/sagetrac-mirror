# the last one may be too hard
L = [i for i in range(1,43) if fixed[i].rank()>=4]
L.reverse()
rw = 'w'
############### primes #############

for p in [2,3,5,7,11]:
    for e<=4:
        if euler_phi(p^e) >= 9:
            continue
        file_name = 'order%s.txt'%p^e
        classify_ord_pe(L, p, e, file_name,rw):
