import enumerate_short_vectors_boost_python

def short_vectors(lattice, lower_bound, upper_bound, up_to_sign = False):
    return enumerate_short_vectors_boost_python \
        .enumerate_short_vectors_boost_python([[int(e) for e in r] for r in lattice],
                                              int(lower_bound), int(upper_bound),
                                              bool(up_to_sign))
