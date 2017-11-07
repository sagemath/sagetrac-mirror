/*
  Low-level code to exhaust over trees of Weil polynomials

*/

#include <flint.h>
#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>
#include <arith.h>

/*
    Use a subresultant (Sturm-Habicht) sequence to test whether a given
    polynomial has all real roots. Note that this test has an early abort
    mechanism: having all real roots means that the sign sequence has
    the maximal number of sign changes, so the test can be aborted as soon
    as a sign change is missed. Moreover, if this miss is early enough, it
    does not even depend on the lowest-order coefficients of the polynomial.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2
        - {w, 3 * n + 8} is scratch space.

    Return values:
        1: poly has all real roots
        j <= 0: poly does not have all real roots, and this remains
                true for any choice of coefficients in degrees less than -j.

   Based on code by Sebastian Pancratz from the FLINT repository.
   TODO: Compare performance with methods based on root isolation (VCA).
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
	  return(0);
	}

        /* 
            Explicitly compute the pseudoremainder f2 of f0 modulo f1:

            f2 := l0*x*f1 - l1*f0 
            f2 := l1*f2 - f2[n-1]*f1
         */
        fmpz_zero(f2);
        _fmpz_vec_scalar_mul_fmpz(f2 + 1, f1, n-1, l0);
        _fmpz_vec_scalar_submul_fmpz(f2, f0, n, l1);
	fmpz_set(c, f2+n-1);
        _fmpz_vec_scalar_mul_fmpz(f2, f2, n-1, l1);
        _fmpz_vec_scalar_submul_fmpz(f2, f1, n-1, c); 

        if (_fmpz_vec_is_zero(f2, n - 1)) /* We made it! */
	  return(1);

        n--; // len(f2) = n

        /* Extract content from f2; in practice, this seems to do better than
	 an explicit subresultant computation. */
        _fmpz_vec_content(d, f2, n);
        _fmpz_vec_scalar_divexact_fmpz(f0, f2, n, d);

        /* Rotate the polynomials */
	t = f0; f0 = f1; f1 = t;
      }
}

/* Primary data structures for tree exhaustion.
 */

typedef struct ps_static_data {
  int d, sign;
  long node_limit;
  fmpz_t a, b, lead, q;
  fmpz_mat_t binom_mat;
  fmpz *cofactor;
  fmpz *modlist;
  fmpq_mat_t *sum_mats;
  fmpq_t *f;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, ascend, flag;
  long node_count;
  fmpq_mat_t sum_col, sum_prod;
  fmpz *pol, *sympol, *upper;

  /* Scratch space */
  fmpz *w;
  int wlen; /* = 4*d+12 */
  fmpq *w2;
  int w2len; /* = 6 */
} ps_dynamic_data_t;

/* Set res to floor(a). */
void fmpq_floor(fmpz_t res, const fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Set res to ceil(a). */
void fmpq_ceil(fmpz_t res, const fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpz_sqrt_f(fmpz_t res, const fmpz_t a) {
  fmpz_sqrt(res, a);
}

void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  int s = fmpz_is_square(a);
  fmpz_sqrt(res, a);
  if (!s) fmpz_add_ui(res, res, 1);
}

/* Set res to floor(a + b sqrt(q)). 
   For efficiency, we do not assume a and b are canonical;
   we must thus be careful about signs. */
void fmpq_floor_quad(fmpz_t res, const fmpq_t a,
		     const fmpq_t b, const fmpz_t q) {
  if (b==NULL) fmpq_floor(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    int anum_s = fmpz_sgn(anum);
    fmpz *aden = fmpq_denref(a);
    int aden_s = fmpz_sgn(aden);
    fmpz *bnum = fmpq_numref(b);
    int bnum_s = fmpz_sgn(bnum);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);
    
    fmpz_mul(res, aden, bnum);
    fmpz_mul(res, res, res);
    fmpz_mul(res, res, q);
    if (bnum_s*bden_s >= 0) fmpz_sqrt_f(res, res);
    else {
      fmpz_sqrt_c(res, res);
      fmpz_neg(res, res);
    }
    fmpz_mul_si(res, res, aden_s*bden_s);
    fmpz_addmul(res, anum, bden);
    if (bden_s > 0) fmpz_fdiv_q(res, res, aden);
    else fmpz_cdiv_q(res, res, aden);
    fmpz_fdiv_q(res, res, bden);
  }
}

/* Set res to ceil(a + b sqrt(q)). */
void fmpq_ceil_quad(fmpz_t res, const fmpq_t a,
		     const fmpq_t b, const fmpz_t q) {
  if (b==NULL) fmpq_ceil(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    int anum_s = fmpz_sgn(anum);
    fmpz *aden = fmpq_denref(a);
    int aden_s = fmpz_sgn(aden);
    fmpz *bnum = fmpq_numref(b);
    int bnum_s = fmpz_sgn(bnum);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);

    fmpz_mul(res, aden, bnum);
    fmpz_mul(res, res, res);
    fmpz_mul(res, res, q);
    if (bnum_s*bden_s >= 0) fmpz_sqrt_c(res, res);
    else {
      fmpz_sqrt_f(res, res);
      fmpz_neg(res, res);
    }
    fmpz_mul_si(res, res, aden_s*bden_s);
    fmpz_addmul(res, anum, bden);
    if (bden_s > 0) fmpz_cdiv_q(res, res, aden);
    else fmpz_fdiv_q(res, res, aden);
    fmpz_cdiv_q(res, res, bden);
  }
}

/* Memory allocation and initialization. */
ps_static_data_t *ps_static_init(int d, fmpz_t q, int coeffsign, fmpz_t lead,
				 int cofactor, 
				 fmpz *modlist,
				 long node_limit) {
  int i, j;
  ps_static_data_t *st_data;
  fmpz_poly_t pol;
  fmpz_t m, const1;
  fmpq *k1;

  fmpz_poly_init(pol);
  fmpz_init(m);
  fmpz_init_set_ui(const1, 1);

  st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  st_data->sign = coeffsign;
  fmpz_init(st_data->q);
  fmpz_set(st_data->q, q);
  st_data->node_limit = node_limit;

  fmpz_init(st_data->a);
  fmpz_init(st_data->b);
  if (fmpz_is_one(st_data->q)) {
    fmpz_set_si(st_data->a, -2);
    fmpz_set_si(st_data->b, 2);
  } else {
    fmpz_set_si(st_data->a, 0);
    fmpz_mul_si(st_data->b, st_data->q, 4);    
  }
  fmpz_init(st_data->lead);
  fmpz_set(st_data->lead, lead);

  st_data->cofactor = _fmpz_vec_init(3);
  switch (cofactor) {
  case 0: /* Cofactor 1 */
    fmpz_set_si(st_data->cofactor, 1);
    fmpz_set_si(st_data->cofactor+1, 0);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 1: /* Cofactor x+sqrt(q) */
    fmpz_set(st_data->cofactor, st_data->q);
    fmpz_sqrt(st_data->cofactor, st_data->cofactor);
    fmpz_set_si(st_data->cofactor+1, 1);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 2:  /* Cofactor x-sqrt(q) */
    fmpz_set(st_data->cofactor, st_data->q);
    fmpz_sqrt(st_data->cofactor, st_data->cofactor);
    fmpz_neg(st_data->cofactor, st_data->cofactor);
    fmpz_set_si(st_data->cofactor+1, 1);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 3: /* Cofactor x^2-q */
    fmpz_neg(st_data->cofactor, st_data->q);
    fmpz_set_si(st_data->cofactor+1, 0);
    fmpz_set_si(st_data->cofactor+2, 1);
    break;
  }

  st_data->modlist = _fmpz_vec_init(d+1);
  st_data->f = _fmpq_vec_init(d+1);
  for (i=0; i<=d; i++) {
    fmpz_set(st_data->modlist+i, modlist+d-i);
    fmpq_set_si(st_data->f+i, d-i, 1);
    fmpq_div_fmpz(st_data->f+i, st_data->f+i, st_data->lead);
    /* In order to apply power sums and Descartes' rule of signs
       when the modulus is 0, we must pretend that the modulus is 1. */
    if (!fmpz_is_zero(st_data->modlist+i))
      fmpq_mul_fmpz(st_data->f+i, st_data->f+i, st_data->modlist+i);
  }

  fmpz_mat_init(st_data->binom_mat, d+1, d+1);
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(fmpz_mat_entry(st_data->binom_mat, i, j), i, j);
  
  st_data->sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  for (i=0; i<=d; i++) {

    fmpq_mat_init(st_data->sum_mats[i], 9, d+1);
    fmpq_mat_zero(st_data->sum_mats[i]);

    arith_chebyshev_t_polynomial(pol, i);
    for (j=0; j<=d; j++) {
      
      /* Row 0: coeffs of 2*(i-th Chebyshev polynomial)(x/2). 
         If q != 1, the coeff of x^j is multiplied by q^{floor(i-j)/2}. */
      if (j <= i) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 0, j);
	fmpq_set_fmpz_frac(k1, fmpz_poly_get_coeff_ptr(pol, j), const1);
	fmpz_mul_2exp(m, const1, j);
	fmpq_div_fmpz(k1, k1, m);
	fmpz_set_ui(m, 2);
	fmpq_mul_fmpz(k1, k1, m);
	if (!fmpz_is_one(st_data->q) && i%2==j%2) {
	  fmpz_set(m, st_data->q);
	  fmpz_pow_ui(m, m, (i-j)/2); 
	  fmpq_mul_fmpz(k1, k1, m);
	}
      }

      /* The other rows are currently used only when q==1. */
      
      /* Row 1: coeffs of row 0 from matrix i-2, multiplied by -2. */
      if (i >= 2) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 1, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 0, j));
	fmpz_set_si(m, -2);
	fmpq_mul_fmpz(k1, k1, m);
      }

      /* Row 2: coeffs of row 0 from matrix i-2, shifted by 2. */
      if (i>= 2 && j >= 2) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 2, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 0, j-2));
      }

      /* Row 3: coeffs of (2+x)^i. */
      if (j<= i) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 3, j);
	fmpq_set_fmpz_frac(k1, fmpz_mat_entry(st_data->binom_mat, i, j), const1);
	fmpq_mul_2exp(k1, k1, i-j);
      }
      
      /* Row 4: coeffs of (2+x)^(i-1). */
      if (i >= 1) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 4, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-1], 3, j));	
      }

      /* Row 5: coeffs of (2+x)^(i-2). */
      if (i>=2)	{
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 5, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 3, j));	
      }

      /* Row 6: coeffs of (-2+x)^i. */
      k1 = fmpq_mat_entry(st_data->sum_mats[i], 6, j);
      fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i], 3, j));
      if ((i-j)%2==1) fmpq_neg(k1, k1);

      /* Row 7: coeffs of (-2+x)^(i-1). */
      if (i >= 1) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 7, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-1], 6, j));	
      }

      /* Row 8: coeffs of (-2+x)^(i-2). */
      if (i >= 2) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 8, j);
	fmpq_set(k1, fmpq_mat_entry(st_data->sum_mats[i-2], 6, j));
      }

    }
  }
  
  fmpz_poly_clear(pol);
  fmpz_clear(m);
  fmpz_clear(const1);

  return(st_data);
}

ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data;
  int i;

  dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));
  dy_data->d = d;

  /* Initialize mutable quantities */
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+3);
  if (coefflist != NULL)
    for (i=0; i<=d; i++) 
      fmpz_set(dy_data->pol+i, coefflist+i);
  
  fmpq_mat_init(dy_data->sum_col, d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->sum_col, 0, 0), d, 1);

  dy_data->upper = _fmpz_vec_init(d+1);

  /* Allocate scratch space */
  fmpq_mat_init(dy_data->sum_prod, 9, 1);
  dy_data->wlen = 4*d+12;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);
  dy_data->w2len = 6;
  dy_data->w2 = _fmpq_vec_init(dy_data->w2len);
  return(dy_data);
}

ps_dynamic_data_t *ps_dynamic_clone(ps_dynamic_data_t *dy_data) {
  ps_dynamic_data_t *dy_data2;
  int i, d = dy_data->d;

  dy_data2 = ps_dynamic_init(d, NULL);
  dy_data2->n = dy_data->n;
  dy_data2->node_count = dy_data->node_count;
  dy_data2->ascend = dy_data->ascend;
  _fmpz_vec_set(dy_data2->pol, dy_data->pol, d+1);
  _fmpz_vec_set(dy_data2->upper, dy_data->upper, d+1);
  fmpq_mat_set(dy_data2->sum_col, dy_data->sum_col);
  return(dy_data2);
}

ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data) {
  if (dy_data==NULL) return(NULL);

  ps_dynamic_data_t *dy_data2;
  int i, d = dy_data->d, n = dy_data->n;

  for (i=d; i>n+1; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) <0) {
      dy_data2 = ps_dynamic_clone(dy_data);
      fmpz_set(dy_data->upper+i, dy_data->pol+i);
      dy_data2->n = i-1;
      dy_data2->ascend = 1;
      dy_data2->node_count = 0;
      return(dy_data2);
  }
  return(NULL);
}

/* Memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return(NULL);
  int i, d = st_data->d;
  fmpz_clear(st_data->a);
  fmpz_clear(st_data->b);
  fmpz_clear(st_data->lead);
  fmpz_clear(st_data->q);
  _fmpz_vec_clear(st_data->cofactor, 3);
  fmpz_mat_clear(st_data->binom_mat);
  _fmpq_vec_clear(st_data->f, d+1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  for (i=0; i<=d; i++) 
    fmpq_mat_clear(st_data->sum_mats[i]);
  free(st_data->sum_mats);
  free(st_data);
}

void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  if (dy_data == NULL) return(NULL);
  int d = dy_data->d;
  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->sympol, 2*d+3);
  _fmpz_vec_clear(dy_data->upper, d+1);
  fmpq_mat_clear(dy_data->sum_col);
  fmpq_mat_clear(dy_data->sum_prod);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  _fmpq_vec_clear(dy_data->w2, dy_data->w2len);
  free(dy_data);
}

/* Return values: 
   -r, r<0: if the n-th truncated polynomial does not have roots in the
       interval, and likewise for all choices of the bottom r-1 coefficients
   1: if lower <= upper
   0: otherwise. 

   The case n=0 is allowed. In this case, we return 1 if the polynomial is
   admissible and 0 otherwise.
*/
int set_range_from_power_sums(ps_static_data_t *st_data,
			  ps_dynamic_data_t *dy_data) {
  int i, j, r, r1, r2, s;
  int d = st_data->d;
  int n = dy_data->n;
  int k = d+1-n;
  fmpz *modulus = st_data->modlist + n-1;
  fmpz *pol = dy_data->pol;
  fmpz *q = st_data->q;
  fmpq *f;
    
  /* Allocate temporary variables from persistent scratch space. */
  fmpz *tpol = dy_data->w;
  fmpz *tpol2 = dy_data->w + d+1;
  fmpz *tpol3 = dy_data->w + 2*d+2;

  fmpz *t0z = dy_data->w + 3*d+3;
  fmpz *t1z = dy_data->w + 3*d + 4;
  fmpz *t2z = dy_data->w + 3*d + 5;
  fmpz *lower = dy_data->w + 3*d + 6;
  fmpz *upper = dy_data->w + 3*d + 7;
  
  fmpq *t0q = dy_data->w2;
  fmpq *t1q = dy_data->w2 + 1;
  fmpq *t2q = dy_data->w2 + 2;
  fmpq *t3q = dy_data->w2 + 3;
  fmpq *t4q = dy_data->w2 + 4;
  fmpq *t5q = dy_data->w2 + 5;

  /* Subroutines to adjust lower and upper bounds. 
   These use t0q and t4q as persistent scratch space. */

  void set_lower(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_ceil(lower, t0q);
  }
  
  void set_upper(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_floor(upper, t0q);
  }

  void set_lower_quad(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_ceil(lower, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_ceil_quad(lower, t0q, t4q, q);
    }
  }
  
  void set_upper_quad(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_floor(upper, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_floor_quad(upper, t0q, t4q, q);
    }
  }

  void change_lower(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_ceil(t0z, t0q);
    if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
  }
  
  void change_upper(const fmpq_t val) {
    fmpq_div(t0q, val, f);
    fmpq_floor(t0z, t0q);
    if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
  }

  /* Here, we allow lower and upper bounds in the quadratic field
     Q(sqrt(q)). */
  void change_lower_quad(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_ceil(t0z, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_ceil_quad(t0z, t0q, t4q, q);
    }
    if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
  }
  
  void change_upper_quad(const fmpq_t val1, const fmpq_t val2) {
    fmpq_div(t0q, val1, f);
    if (val2==NULL) fmpq_floor(t0z, t0q);
    else {
      fmpq_div(t4q, val2, f);
      fmpq_floor_quad(t0z, t0q, t4q, q);
    }
    if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
  }
    
  /* Compute the divided n-th derivative of pol. */
  for (i=0; i<=k-1; i++)
    fmpz_mul(tpol+i, fmpz_mat_entry(st_data->binom_mat, n+i, n), pol+n+i);

  /* TODO: try using real root isolation instead of Sturm sequences. */
  r = _fmpz_poly_all_real_roots(tpol, k, dy_data->w+d+1);
  if (r<=0) return(r-1);
    
  /* If r=1 and k>d, no further coefficients to bound. */
  if (k>d) return(1);

  /* Compute the k-th power sum. */
  f = fmpq_mat_entry(dy_data->sum_col, k, 0);
  fmpq_set_si(f, -k, 1);
  fmpq_mul_fmpz(f, f, pol+d-k);
  fmpq_div_fmpz(f, f, pol+d);
  for (i=1; i<k; i++) {
    fmpq_set_fmpz_frac(t0q, pol+d-i, pol+d);
    fmpq_neg(t0q, t0q);
    fmpq_addmul(f, t0q, fmpq_mat_entry(dy_data->sum_col, k-i, 0));
  }
  
  /* First, set bounds using symmetrized power sums. */
  f = st_data->f + n-1;
  fmpq_mat_mul(dy_data->sum_prod, st_data->sum_mats[k], dy_data->sum_col);
  
  if (q == 1) {
    fmpq_set_si(t1q, 2*d, 1);
    fmpq_sub(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
    set_lower(t0q);
    fmpq_add(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
    set_upper(t0q);
  }
  else if (k%2==0) {
    fmpq_set_si(t1q, 2*d, 1);
    fmpz_set(t0z, q);
    fmpz_pow_ui(t0z, t0z, k/2);
    fmpq_mul_fmpz(t1q, t1q, t0z);
    fmpq_sub(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
    set_lower(t0q);
    fmpq_add(t0q, fmpq_mat_entry(dy_data->sum_prod, 0, 0), t1q);
    set_upper(t0q);
  } else {
    fmpq_zero(t1q); 
    fmpq_set_si(t2q, 2*d, 1);
    fmpz_set(t0z, q);
    fmpz_pow_ui(t0z, t0z, k/2);
    fmpq_mul_fmpz(t2q, t2q, t0z);
    set_upper_quad(fmpq_mat_entry(dy_data->sum_prod, 0, 0), t2q);
    fmpq_neg(t2q, t2q);
    set_lower_quad(fmpq_mat_entry(dy_data->sum_prod, 0, 0), t2q);
  }

  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Second, apply Descartes' rule of signs at -2*sqrt(q), +2*sqrt(q);
     this enforces the roots being in the correct interval (if real). */

  fmpq_set_si(t3q, -k, 1);
  fmpq_div_fmpz(t3q, t3q, pol+d);

  /* Currently tpol is the divided n-th derivative of pol.
     Undo one derivative, then evaluate at the endpoints. */
  for (i=k; i>=1; i--) {
    fmpz_mul_si(tpol+i, tpol+i-1, n);
    fmpz_divexact_si(tpol+i, tpol+i, i);
  }
  fmpz_set(tpol, pol+d-k);
  
  if (q == 1) {
    _fmpz_poly_evaluate_fmpz(t0z, tpol, k+1, st_data->a);    
    fmpq_mul_fmpz(t1q, t3q, t0z);
    if (k%2==1) change_upper(t1q);
    else change_lower(t1q);
    
    _fmpz_poly_evaluate_fmpz(t0z, tpol, k+1, st_data->b);
    fmpq_mul_fmpz(t1q, t3q, t0z);
    change_lower(t1q);
  } else {
    for (i=0; 2*i <= k; i++)
      fmpz_set(tpol2+i, tpol+2*i);
    for (i=0; 2*i+1 <= k; i++)
      fmpz_set(tpol3+i, tpol+2*i+1);
    fmpz_mul_si(t2z, q, 4);
    _fmpz_poly_evaluate_fmpz(t0z, tpol2, (k+2) / 2, t2z);
    _fmpz_poly_evaluate_fmpz(t1z, tpol3, (k+1) / 2, t2z);
    fmpz_mul_si(t1z, t1z, 2);
    fmpq_mul_fmpz(t1q, t3q, t0z);
    fmpq_mul_fmpz(t2q, t3q, t1z);

    change_lower_quad(t1q, t2q);

    fmpq_neg(t2q, t2q);
    if (k%2==1) change_upper_quad(t1q, t2q);
    else change_lower_quad(t1q, t2q);
  }

  /* If modulus=0, then return 1 if [lower, upper] contains 0
     and 0 otherwise. After this, we may assume modulus>0.
   */
  if (fmpz_is_zero(modulus)) {
    if ((fmpz_sgn(lower) > 0) || (fmpz_sgn(upper) < 0)) return(0);
      fmpz_zero(lower);
      fmpz_zero(upper);
      return(1);
  }

  /* Third, if q=1, compute additional bounds using power sums. 
  */
  if (q==1 && (fmpz_cmp(lower, upper) <= 0) && k >= 2) {
    /* The k=2 case requires separate attention; this corrects a bug
       in the 2008 implementation.
    */
    fmpq_add(t1q, fmpq_mat_entry(dy_data->sum_prod, 1, 0),
	     fmpq_mat_entry(dy_data->sum_prod, 2, 0));
    fmpq_set_si(t2q, 4*d, 1);
    fmpq_sub(t0q, t1q, t2q);
    // fmpq_sub_si(t3q, t1q, 4*d);
    if (k==2) fmpq_div_fmpz(t0q, t0q, st_data->b); // b=2
    change_lower(t0q);
    fmpq_add(t0q, t1q, t2q);
    // fmpq_add_si(t3q, t1q, 4*d);
    if (k==2) fmpq_div_fmpz(t0q, t0q, st_data->b); // b=2
    change_upper(t0q);
    
    /* t1q, t2q, t3q are no longer needed, so can be reassigned. */
    t1q = fmpq_mat_entry(dy_data->sum_prod, 3, 0);
    t2q = fmpq_mat_entry(dy_data->sum_prod, 4, 0);
    t3q = fmpq_mat_entry(dy_data->sum_prod, 5, 0);
    if (fmpq_sgn(t3q) > 0) { // t0q <- t1q - t2q^2/t3q
      fmpq_mul(t0q, t2q, t2q);
      fmpq_div(t0q, t0q, t3q);
      fmpq_sub(t0q, t1q, t0q);
      change_upper(t0q);
    }
    fmpq_set_si(t3q, -4, 1);
    fmpq_mul(t0q, t3q, t2q);
    // fmpq_mul_si(t0q, t2q, -4);
    fmpq_add(t0q, t0q, t1q);
    change_lower(t0q);
    
    t1q = fmpq_mat_entry(dy_data->sum_prod, 6, 0);
    t2q = fmpq_mat_entry(dy_data->sum_prod, 7, 0);
    t3q = fmpq_mat_entry(dy_data->sum_prod, 8, 0);
    if ((k%2 == 0) && (fmpq_sgn(t3q) > 0)) {
      fmpq_mul(t0q, t2q, t2q);
      fmpq_div(t0q, t0q, t3q);
      fmpq_sub(t0q, t1q, t0q);
      change_upper(t0q);
    } else if ((k%2 == 1) && (fmpq_sgn(t3q) < 0)) {
      fmpq_mul(t0q, t2q, t2q);
      fmpq_div(t0q, t0q, t3q);
      fmpq_sub(t0q, t1q, t0q);
      change_lower(t0q);
    }
    fmpq_set_si(t0q, 4, 1);
    fmpq_mul(t0q, t0q, t2q);
    // fmpq_mul_si(t0q, t2q, 4);
    fmpq_add(t0q, t0q, t1q);
    if (k%2 == 0) change_lower(t0q);
    else change_upper(t0q);
    
    if (k%2 == 0) {
	fmpq_set_si(t0q, -4, 1);
	fmpq_mul(t0q, t0q, fmpq_mat_entry(dy_data->sum_col, k-2, 0));
      // fmpq_mul_si(t0q, fmpq_mat_entry(dy_data->sum_col, k-2, 0), -4);
	fmpq_add(t0q, t0q, fmpq_mat_entry(dy_data->sum_col, k, 0));
	change_lower(t0q);	
      }
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);
    
  /* Set the new upper bound. Note that modulus>0 at this point. */
  fmpz_mul(upper, upper, modulus);
  fmpz_add(dy_data->upper+n-1, pol+n-1, upper);
  /* Correct the k-th power sum. */
  t1q = fmpq_mat_entry(dy_data->sum_col, k, 0);
  fmpq_mul_fmpz(t0q, f, lower);
  fmpq_sub(t1q, t1q, t0q);
  /* Set the new polynomial value. */
  fmpz_mul(lower, lower, modulus);
  fmpz_add(pol+n-1, pol+n-1, lower);

  return(1);

}

/* Return value sent back in dy_data->flag:
    1: if a solution has been found
    0: if the tree has been exhausted
   -1: if the maximum number of nodes has been reached
*/

void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data) {
  if (dy_data==NULL) return(0);

  int d = st_data->d;
  int node_limit = st_data->node_limit;
  fmpz *modlist = st_data->modlist;

  int ascend = dy_data->ascend;
  int n = dy_data->n;
  int count = dy_data->node_count;
  fmpz *upper = dy_data->upper;
  fmpz *pol = dy_data->pol;
  fmpz *sympol = dy_data->sympol;

  int i, j, t, r;
  fmpq *tq;

  if (n>d) return(0);
  while (1) {
    if (ascend > 0) {
      n += 1;
      if (n>d) {
	/* We have exhausted the entire tree. */
	t=0;
	break;
      }
    } else {
      i = dy_data->n;
      dy_data->n = n;
      r = set_range_from_power_sums(st_data, dy_data);
      if (r > 0) {
	n -= 1;
	if (n<0) { 
	  t=1;
	  /* We have found a solution! Convert it back into reciprocal form for output. */
	  _fmpz_vec_zero(sympol, 2*d+3);
	  fmpz *temp = dy_data->w;
	  for (i=0; i<=d; i++) {
	    fmpz_one(temp);
	    for (j=0; j<=i; j++) {
	      fmpz_addmul(sympol+2*d-(d-i+2*j), pol+i, temp);
	      if (j<i) {
		fmpz_mul(temp, temp, st_data->q);
		fmpz_mul_si(temp, temp, i-j);
		fmpz_divexact_si(temp, temp, j+1);
	      }
	    }
	  }
	  _fmpz_vec_scalar_mul_si(sympol, sympol, 2*d+1, st_data->sign);
	  _fmpz_poly_mul_KS(sympol,sympol, 2*d+1, st_data->cofactor, 3);
	  break; 
	}
	continue;
      } else {
	count += 1;
	if (node_limit != -1 && count >= node_limit) { t= -1; break; }
	if (r<-1) {
	  /* Early abort: Sturm test failed on a coefficient determined at 
	     a previous level. */
	  ascend = -r-1;
	  continue;
	} else if (r==-1 && i<n) { 
	/* Early abort: given the previous coefficient, the set of values for
	   a given coefficient giving the right position of real roots for
	   the corresponding derivative is always an interval. */
	ascend = 1;
	continue;
	}
      }
    }
    if (ascend>1) ascend -= 1;
    else if (fmpq_is_zero(modlist+n)) ascend = 1;
    else {
      fmpz_add(pol+n, pol+n, modlist+n);
      if (fmpz_cmp(pol+n, upper+n) > 0) ascend = 1;
      else {
	ascend = 0;
	/* Update the (d-n)-th power sum. */
	tq = fmpq_mat_entry(dy_data->sum_col, d-n, 0);
	fmpq_sub(tq, tq, st_data->f+n);
      }
    }
  }
  dy_data->ascend = (n<0);
  dy_data->n = n;
  dy_data->node_count = count;
  dy_data->flag = t;
}
