def _jacobi_theta(self, vec, max=10):
    """
    Compute the theta series as a power series in the variable given
    in ``var_str`` (which defaults to `q`), up to the specified precision
    `O(q^max)`.

    EXAMPLES:
    """
    from sage.modular.jacobi.jacobi_form import Jacobi_Form
    return Jacobi_Form(self, vec, max).Fourier_expansion()


__jacobi_theta = _jacobi_theta
