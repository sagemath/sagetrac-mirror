#include <fmpz_poly.h>

/* Based on code by Sebastian Pancratz from the FLINT repository.
   TODO: Compare performance with methods based on root isolation (VCA).
    */

#define SWAP(f, g)  do { t = f; f = g; g = t; } while (0)

/*
    Assumes that:
        - {poly, n} is a normalized vector with n >= 2
        - {w, 3 * n + 8} is scratch space.

    Return values:
        1: poly has all real roots
        j <= 0: poly does not have all real roots, and this remains
                true for any choice of coefficients in degrees less than -j.
*/
int _fmpz_poly_all_real_roots(fmpz *poly, slong n, fmpz *w)
{
    fmpz *f0     = w + 0 * n;
    fmpz *f1     = w + 1 * n;
    fmpz *f2     = w + 2 * n;
    fmpz *c      = w + 3 * n;
    fmpz *d      = w + 3 * n + 1;

    fmpz *l0;
    fmpz *l1;
    fmpz *t;

    int sgn0_l;
    int sgn1_l;
    int i;
    int j;
    slong n0 = n-1;

    if (n == 1)
        return 1;

    _fmpz_vec_set(f0, poly, n);
    _fmpz_poly_derivative(f1, f0, n);
    n--;
    sgn0_l = fmpz_sgn(f0+n);
    
    for ( ; ; )
      {
        /* Invariant:  n = len(f0) - 1, len(f1) <= n */

        l0 = f0 + n;
        l1 = f1 + n-1;
	sgn1_l = fmpz_sgn(l1);
        /* If we miss any one sign change, we cannot have enough */
	if (sgn1_l == 0) return(0);
	if (sgn1_l != sgn0_l) {
	  j = 2*n - n0+1; 
	  if (j>0) return(-j); /* Independent of terms of degree <j */
	  return 0;
	}

        /* 
            Explicitly compute the pseudoremainder f2 of f0 modulo f1
            f2 := l0 * x * f1 - l1 * f0 
            f2 := l1 * f2 - f2[n-1] * f1
         */
        fmpz_zero(f2);
        _fmpz_vec_scalar_mul_fmpz(f2 + 1, f1, n-1, l0);
        _fmpz_vec_scalar_submul_fmpz(f2, f0, n, l1);
	fmpz_set(c, f2+n-1);
        _fmpz_vec_scalar_mul_fmpz(f2, f2, n-1, l1);
        _fmpz_vec_scalar_submul_fmpz(f2, f1, n-1, c); 

        if (_fmpz_vec_is_zero(f2, n - 1))
            return 1;

        n--; // len(f2) = n

        /* Extract content from f2; in practice, this seems to do better than
	 an explicit subresultant computation. */
        _fmpz_vec_content(d, f2, n);
        _fmpz_vec_scalar_divexact_fmpz(f0, f2, n, d);

        /* Rotate the polynomials */
        SWAP(f0, f1);
      }

    return 1;
}
