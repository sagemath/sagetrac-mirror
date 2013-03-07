




def _jacobi_theta(self, Vec, Max=10):
    """
    Compute the theta series as a power series in the variable given
    in var_str (which defaults to '`q`'), up to the specified precision
    `O(q^max)`.

    """
    from sage.modular.jacobi.jacobi_form import Jacobi_Form
    JF = Jacobi_Form(self,Vec,Max)
    return JF.Fourier_expansion()


def __jacobi_theta(self,Vec, Max):
    return _jacobi_theta(self,Vec,Max)
