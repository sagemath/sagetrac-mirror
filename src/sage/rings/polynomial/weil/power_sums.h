#include <fmpz_poly.h>
#include <fmpq.h>
#include <fmpq_mat.h>

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
  slong wlen; /* = 4*d+12 */
  fmpq *w2;
  slong w2len; /* = 5 */
} ps_dynamic_data_t;

ps_static_data_t *ps_static_init(int d, fmpz_t q, int coeffsign, fmpz_t lead,
				 int cofactor, fmpz *modlist, long node_limit);
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist);
void ps_static_clear(ps_static_data_t *st_data);
void ps_dynamic_clear(ps_dynamic_data_t *dy_data);
void extract_pol(int *Q, ps_dynamic_data_t *dy_data);
// void extract_symmetrized_pol(int *Q, ps_dynamic_data_t *dy_data);
// long extract_count(ps_dynamic_data_t *dy_data);
ps_dynamic_data_t *ps_dynamic_clone(ps_dynamic_data_t *dy_data);
ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data);
void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data);

