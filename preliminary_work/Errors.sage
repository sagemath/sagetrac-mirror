def random_error_pos(n, noErrs):
    errPos = []
    i = 0
    while i<n and noErrs>0:
        if random() < noErrs/(n-i):
            errPos.append(i)
            noErrs -= 1
        i += 1
    return errPos

def random_error_vect(n, F, errPos):
    vect=[F.zero()]*n
    for i in errPos:
        while vect[i].is_zero():
            vect[i] = F.random_element()
    return vector(vect)

