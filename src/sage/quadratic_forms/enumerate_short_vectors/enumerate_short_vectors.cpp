/**
 *
 * Copyright (C) 2013 Martin Raum
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
 * Implement Fincke-Pohst algorithm with multi precision floating points.
 */

#include "enumerate_short_vectors_internal.hpp"

#include "enumerate_short_vectors.hpp"

using namespace std;


void
enumerate_short_vectors
(
 const vector<vector<int>> &qfmatrix,
 unsigned int lower_bound,
 unsigned int upper_bound,
 map<unsigned int, vector<vector<int>>> &result
 )
{
  // check dimensions of the quadratic form matrix
  size_t m = qfmatrix.size();
  for ( auto &row : qfmatrix )
    if ( row.size() != m )
      return;

  if (lower_bound > upper_bound )
    return;

  // standard precision
  mp_prec_t precision = 53;

  for ( size_t ind = lower_bound; ind <= upper_bound; ind += 2 )
    result[ind] = vector<vector<int>>();


  // create varibles
  mpfr_ptr mpfr_tmp = new __mpfr_struct;
  mpfi_ptr mpfi_tmp = new __mpfi_struct;
  mpfi_ptr mpfi_tmp2 = new __mpfi_struct;

  vector<vector<mpfi_ptr>> rmatrix;
  vector<mpfi_ptr> rdiag_sqrt;

  vector<mpfi_ptr> vec_Ti;
  vector<mpfi_ptr> vec_Ui;
  vector<vector<mpfi_ptr>> vec_Uij;

  mpfi_ptr C = new __mpfi_struct;
  mpfi_ptr Z = new __mpfi_struct;

  auto vec_UB = vector<int>(m, 0);
  auto vec_LB = vector<int>(m, 0);
  int IB_lower, IB_upper;

  vector<int> vec_x = vector<int>(m, 0);
  int int_tmp;
  bool x_is_zero;

  // we use zero based indices, so we start with i = m - 1
  size_t i{ m - 1 };

  init( m, vec_Ti, vec_Ui, vec_Uij, C, Z, rmatrix, rdiag_sqrt, mpfr_tmp, mpfi_tmp, mpfi_tmp2 );

  recompute( i, m, vec_x, lower_bound, upper_bound,
	     vec_LB, vec_UB, IB_lower, IB_upper,
	     vec_Ti, vec_Ui, vec_Uij,
	     C, Z,
	     qfmatrix, rmatrix, rdiag_sqrt,
	     mpfi_tmp, mpfi_tmp2, mpfr_tmp,
	     precision );
  vec_x[i] = vec_LB[i] - 1;
  for ( size_t j = 0; j < i; ++j )
    mpfi_mul_si( vec_Uij[j][i], rmatrix[j][i], vec_x[i] );

  while ( true )
    {
      // step 3
      ++vec_x[i];

      // in order to implement the lower bound, we have this extra condition
      if ( i == 0 && vec_x[0] == IB_lower )
	{
	  vec_x[0] = IB_upper;
	  mpfi_mul_si( vec_Uij[0][0], rmatrix[0][0], IB_upper);
	}
      else
	for ( size_t j = 0; j <= i; ++j )
	  mpfi_add( vec_Uij[j][i], vec_Uij[j][i], rmatrix[j][i] );

      if ( vec_x[i] > vec_UB[i] ) // goto step 5
        {
          ++i;
          continue; // goto step 3
        }
      else
        {
          if ( i == 0 ) // goto step 6
            {
              x_is_zero = true;
	      for ( auto e : vec_x )
		if ( e != 0 )
		  {
		    x_is_zero = false;
		    break;
		  }
              if ( x_is_zero )
                break;

              // Q(x) = C - T_1 + q_{1, 1} * (x_1 + U_1)^2
              mpfi_add_si( mpfi_tmp, vec_Ui[0], vec_x[0] );
              mpfi_sqr( mpfi_tmp, mpfi_tmp );
              mpfi_mul( mpfi_tmp, rmatrix[0][0], mpfi_tmp );
              mpfi_sub( mpfi_tmp, vec_Ti[0], mpfi_tmp );
	      if ( !mpfi_get_unique_si( int_tmp, mpfi_tmp, mpfr_tmp ) )
		{
		  bool IB_jump = ( vec_x[0] == IB_upper );
		  recompute( i, m, vec_x, lower_bound, upper_bound,
			     vec_LB, vec_UB, IB_lower, IB_upper,
			     vec_Ti, vec_Ui, vec_Uij,
			     C, Z,
			     qfmatrix, rmatrix, rdiag_sqrt,
			     mpfi_tmp, mpfi_tmp2, mpfr_tmp,
			     precision );

		  if ( IB_jump )
		    vec_x[0] = IB_lower - 1;
		  else
		    --vec_x[0];

		  continue;
		}

	      // These conditions can only occur, if 
	      // any of the adaption steps for sqrt's, div's, floor's
	      // or ceil's was not correct. In this case, increasing
	      // the precision would not necessarily lead to
	      // correct bounds.
	      if ( int_tmp < 0 || lower_bound > ( int_tmp = upper_bound - int_tmp ) )
		continue;

	      result[int_tmp].push_back( vec_x );
            }
          else // step 5
            {
              --i;

              // U_i = \sum_{j = i + 1}^m q_ij x_j
              mpfi_set_si( vec_Ui[i], 0 );
              for ( size_t j = i + 1; j < m; ++j )
		mpfi_add( vec_Ui[i], vec_Ui[i], vec_Uij[i][j] );

              // T_i = T_{i + 1} - q_{i + 1, i + 1} * (x_{i + 1} + U_{i + 1})^2
              mpfi_add_si( mpfi_tmp, vec_Ui[i + 1], vec_x[i + 1] );
              mpfi_sqr( mpfi_tmp, mpfi_tmp );
              mpfi_mul( mpfi_tmp, rmatrix[i + 1][i + 1], mpfi_tmp );
              mpfi_sub( vec_Ti[i], vec_Ti[i + 1], mpfi_tmp );

	      // step 2
	      if ( !step_2( i, vec_x, lower_bound, upper_bound,
			    vec_LB, vec_UB, IB_lower, IB_upper, Z,
			    vec_Ti, vec_Ui, vec_Uij, rmatrix, rdiag_sqrt,
			    mpfi_tmp, mpfi_tmp2, mpfr_tmp ) )
		{
		  ++i;
		  --vec_x[i];
		  recompute( i, m, vec_x, lower_bound, upper_bound,
			     vec_LB, vec_UB, IB_lower, IB_upper,
			     vec_Ti, vec_Ui, vec_Uij,
			     C, Z,
			     qfmatrix, rmatrix, rdiag_sqrt,
			     mpfi_tmp, mpfi_tmp2, mpfr_tmp,
			     precision );
		  continue;
		}
	      else
		// We don't have to recompute vec_Uij, because
		// it was set to the right value in step_2(...),
		// which predicts that vec_x[i] is later set to
		// vec_LB[i] - 1.
		vec_x[i] = vec_LB[i] - 1;
            }
        }
    }

  // clear variables
  clear( vec_Ti, vec_Ui, vec_Uij, C, Z, rmatrix, rdiag_sqrt, mpfr_tmp, mpfi_tmp, mpfi_tmp2 );
}
