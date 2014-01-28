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

#ifndef __ENUMERATE_SHORT_VECTORS_INTERNAL_H
#define __ENUMERATE_SHORT_VECTORS_INTERNAL_H

#include <vector>

#include "mpfr.h"
#include "mpfi.h"


static const size_t precision_increment = 3;


inline
bool
cholesky_decomposition( const std::vector<std::vector<int>>&,
			std::vector<std::vector<mpfi_ptr>>, std::vector<mpfi_ptr>,
			mpfi_ptr );

inline
bool
step_2( size_t, std::vector<int>&, int, int,
	std::vector<int>&, std::vector<int>&, int&, int&, mpfi_ptr&,
	std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
	std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
	mpfi_ptr&, mpfi_ptr&, mpfr_ptr& );

inline
void
recompute( size_t&, size_t, std::vector<int>&, unsigned int, unsigned int,
	   std::vector<int>&, std::vector<int>&, int&, int&,
	   std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
	   mpfi_ptr&, mpfi_ptr&,
	   const std::vector<std::vector<int>>&, std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
	   mpfi_ptr&, mpfi_ptr&, mpfr_ptr&,
	   mp_prec_t );

inline
void
init( size_t, std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
      mpfi_ptr&, mpfi_ptr&, std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
      mpfr_ptr&, mpfi_ptr&, mpfi_ptr& );

inline
void
clear( std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
       mpfi_ptr&, mpfi_ptr&, std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
       mpfr_ptr&, mpfi_ptr&, mpfi_ptr& );

inline bool mpfi_get_unique_si( int&, mpfi_ptr, mpfr_ptr );
inline bool mpfi_get_unique_floor_si( int&, mpfi_ptr, mpfr_ptr, bool );
inline bool mpfi_get_unique_ceil_si( int&, mpfi_ptr, mpfr_ptr, bool );

inline void mpfi_init_vector( std::vector<mpfi_ptr>&, size_t );
inline void mpfi_init_matrix( std::vector<std::vector<mpfi_ptr>>&, size_t );
inline void mpfi_set_prec_vector( std::vector<mpfi_ptr>&, mp_prec_t );
inline void mpfi_set_prec_matrix( std::vector<std::vector<mpfi_ptr>>&, mp_prec_t );
inline void mpfi_clear_vector( std::vector<mpfi_ptr>& );
inline void mpfi_clear_matrix( std::vector<std::vector<mpfi_ptr>>& );


/*
 * Implementations
 */

using namespace std;


bool
cholesky_decomposition
(
 const vector<vector<int>>& qfmatrix,
 vector<vector<mpfi_ptr>> rmatrix,
 vector<mpfi_ptr> rdiag_sqrt,
 mpfi_ptr mpfi_tmp
 )
{
  size_t m = qfmatrix.size();

  for ( size_t i = 0; i < m; ++i )
    {
      for ( size_t j = 0; j < m; ++j )
	mpfi_set_si( rmatrix[i][j], qfmatrix[i][j] );

      // If the input is too large, then the intervals associated
      // to the diagonal entries could contain zero.
      if ( mpfr_cmp_si( &rmatrix[i][i]->left, 0 ) <= 0 )
	return false;
    }
    
  for ( size_t i = 0; i < m; ++i )
    {
      // step 1 in Fincke-Pohst
      // q_ji <- q_ij, q_ij <- q_ij / q_ii

      for ( size_t j = i + 1; j < m; ++j )
	{
	  mpfi_set( rmatrix[j][i], rmatrix[i][j] );
	  mpfi_div( rmatrix[i][j], rmatrix[i][j], rmatrix[i][i] );
	}

      // step 2 in Fincke-Pohst
      // q_kl <- q_kl - q_ki q_il

      for ( size_t k = i + 1; k < m; ++k )
	for ( size_t l = k; l < m; ++l )
	  {
	    mpfi_mul( mpfi_tmp, rmatrix[k][i], rmatrix[i][l] );
	    mpfi_sub( rmatrix[k][l], rmatrix[k][l], mpfi_tmp );
	  }
    }

  // we later need \sqrt{q_ii}, which we precompute here.
  for ( size_t i = 0; i < m; ++i )
    {
      if ( mpfr_cmp_si( &rmatrix[i][i]->left, 0 ) < 0 )
	return false;

      mpfi_sqrt( rdiag_sqrt[i], rmatrix[i][i] );
    }

  return true;
}

inline
bool
step_2
(
 size_t i,
 vector<int> &vec_x,
 int lower_bound,
 int upper_bound,
 vector<int> &vec_LB,
 vector<int> &vec_UB,
 int &IB_lower,
 int &IB_upper,
 mpfi_ptr &Z,
 vector<mpfi_ptr> &vec_Ti,
 vector<mpfi_ptr> &vec_Ui,
 vector<vector<mpfi_ptr>> &vec_Uij,
 vector<vector<mpfi_ptr>> &rmatrix,
 vector<mpfi_ptr> &rdiag_sqrt,
 mpfi_ptr &mpfi_tmp,
 mpfi_ptr &mpfi_tmp2,
 mpfr_ptr &mpfr_tmp
 )
{
  // Z = (T_i / q_ii)^(1/2)

  // By the following we potentially increase the number of vectors
  // that are computed.
  if ( mpfr_cmp_si( &vec_Ti[i]->left, 0 ) < 0 )
    {
      if ( mpfr_cmp_si( &vec_Ti[i]->right, 0 ) < 0 )
	return false;

      mpfr_set_si( &vec_Ti[i]->left, 0, MPFR_RNDZ );
    }
  mpfi_sqrt( mpfi_tmp, vec_Ti[i] );
  mpfi_div( Z, mpfi_tmp, rdiag_sqrt[i] );

  // UB_i = floor( Z - U_i )
  mpfi_sub( mpfi_tmp, Z, vec_Ui[i] );
  if ( !mpfi_get_unique_floor_si( vec_UB[i], mpfi_tmp, mpfr_tmp, false ) )
    return false;

  // set_xi is only false if constants are recomuted, which happens rarely
  // this justifies to compute LB in all cases so initialization
  // can later access LB

  // x_i = ceil( - Z - U_i ) - 1
  mpfi_add( mpfi_tmp, Z, vec_Ui[i] );
  mpfi_neg( mpfi_tmp, mpfi_tmp );
  if ( !mpfi_get_unique_ceil_si( vec_LB[i], mpfi_tmp, mpfr_tmp, false ) )
    return false;

  // U_ij = q_ij x_j
  // We compute the new values as if vec_x was set to vec_LB[i] - 1
  // In case we step_2(...) is called during a recomputation (which happens
  // rarely), then vec_Uij need to be computed a second time.
  for ( size_t j = 0; j < i; ++j )
    mpfi_mul_si( vec_Uij[j][i], rmatrix[j][i], vec_LB[i] - 1);

  // compute intermediate bounds corresponding to lower_bound
  if ( i == 0 )
    {
      mpfi_sub_ui( mpfi_tmp, vec_Ti[0], upper_bound - lower_bound );

      if ( mpfr_cmp_si( &mpfi_tmp->left, 0 ) < 0 )
	{
	  if ( mpfr_cmp_si( &mpfi_tmp->right, 0 ) < 0 )
	    {
	      // independently of vec_x[0], the lower bound will never be attained.
	      IB_lower = vec_LB[i] - 1;
	      return true;
	    }

	  // adjust the bound, so that we can compute the square root
	  // by this we only enlarge the number of vector enumerated
	  mpfr_set_si( &mpfi_tmp->left, 0, MPFR_RNDZ );
	}

      mpfi_sqrt( mpfi_tmp, mpfi_tmp );
      mpfi_div( mpfi_tmp, mpfi_tmp, rdiag_sqrt[0] );
      mpfi_set( mpfi_tmp2, mpfi_tmp );
		  
      mpfi_add( mpfi_tmp, mpfi_tmp, vec_Ui[0] );
      mpfi_neg( mpfi_tmp, mpfi_tmp );
      mpfi_add_si( mpfi_tmp, mpfi_tmp, 1 );

      if ( !mpfi_get_unique_floor_si( IB_lower, mpfi_tmp, mpfr_tmp, false ) )
	return false;

      mpfi_sub( mpfi_tmp, mpfi_tmp2, vec_Ui[0] );
      if ( !mpfi_get_unique_ceil_si( IB_upper, mpfi_tmp, mpfr_tmp, false ) )
	return false;
      
      // This can happen because of intervals representing integers that
      // give wrong answers when being rounded.
      if ( IB_upper < IB_lower )
	IB_upper = IB_lower;

      // We must not prevent the algorithm from terminating.  If
      // all but the first entry vanish, we therefore set
      // IB_upper to 0.
      bool x_is_zero = true;
      for ( auto it = vec_x.begin() + 1, it_end = vec_x.end();
	    it != it_end; ++it )
	if ( *it != 0 )
	  {
	    x_is_zero = false;
	    break;
	  }
      if ( x_is_zero )
	IB_upper = 0;
    }

  return true;
}

/*
 * After changing the precision, we recompute all quantities. 
 */

inline void
recompute
(
 size_t &i_current,
 size_t m,
 vector<int> &vec_x,
 unsigned int lower_bound,
 unsigned int upper_bound,
 vector<int> &vec_LB,
 vector<int> &vec_UB,
 int &IB_lower,
 int &IB_upper,
 vector<mpfi_ptr> &vec_Ti,
 vector<mpfi_ptr> &vec_Ui,
 vector<vector<mpfi_ptr>> &vec_Uij,
 mpfi_ptr &C,
 mpfi_ptr &Z,
 const vector<vector<int>> &qfmatrix,
 vector<vector<mpfi_ptr>> &rmatrix,
 vector<mpfi_ptr> &rdiag_sqrt,
 mpfi_ptr &mpfi_tmp,
 mpfi_ptr &mpfi_tmp2,
 mpfr_ptr &mpfr_tmp,
 mp_prec_t precision
 )
{
  while ( true )
    {
      precision += precision_increment;

      // set precisions
      mpfr_set_prec( mpfr_tmp, precision );
      mpfi_set_prec( mpfi_tmp, precision );

      mpfi_set_prec_matrix( rmatrix, precision );
      mpfi_set_prec_vector( rdiag_sqrt, precision );

      mpfi_set_prec_vector( vec_Ti, precision );
      mpfi_set_prec_vector( vec_Ui, precision );
      mpfi_set_prec_matrix( vec_Uij, precision );

      mpfi_set_prec( C, precision );
      mpfi_set_prec( Z, precision );


      if ( !cholesky_decomposition( qfmatrix, rmatrix, rdiag_sqrt, mpfi_tmp ) )
	continue;

      // step 1
      mpfi_set_ui( vec_Ti[m - 1], upper_bound );
      mpfi_set_ui( vec_Ui[m - 1], 0 );

      // step 2
      if( !step_2( m - 1, vec_x, lower_bound, upper_bound,
		   vec_LB, vec_UB, IB_lower, IB_upper, Z,
		   vec_Ti, vec_Ui, vec_Uij, rmatrix, rdiag_sqrt,
		   mpfi_tmp, mpfi_tmp2, mpfr_tmp ) )
	continue;
      else if ( vec_x[m - 1] < vec_LB[m - 1] )
	{
	  vec_x[m - 1] = vec_LB[m - 1] - 1;
	  i_current = m - 1;

	  // recompute vec_Uij
	  for ( size_t j = 0; j < m - 1; ++j )
	    mpfi_mul_si( vec_Uij[j][m - 1], rmatrix[j][m - 1], vec_x[m - 1] );

	  return;
	}
      // recompute vec_Uij
      for ( size_t j = 0; j < m - 1; ++j )
	mpfi_mul_si( vec_Uij[j][m - 1], rmatrix[j][m - 1], vec_x[m - 1] );

      // we make use of overflow in the for loop
      size_t i;
      for ( i = m - 2; i >= i_current && i < m - 1; --i )
	{
	  // step 5
      
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
	  if( !step_2( i, vec_x, lower_bound, upper_bound,
		       vec_LB, vec_UB, IB_lower, IB_upper, Z,
		       vec_Ti, vec_Ui, vec_Uij, rmatrix, rdiag_sqrt,
		       mpfi_tmp, mpfi_tmp2, mpfr_tmp ) )
	    {
	      i = m - 1;
	      break;
	    }
	  // This can happen after correcting wrong sqrt, div, floor, or ceil
	  // operations
	  else if ( vec_x[i] < vec_LB[i] )
	    {
	      vec_x[i] = vec_LB[i] - 1;
	      i_current = i;

	      // recompute vec_Uij
	      for ( size_t j = 0; j < i; ++j )
		mpfi_mul_si( vec_Uij[j][i], rmatrix[j][i], vec_x[i] );

	      return;
	    }
	  // recompute vec_Uij, which was set wrongly in step_2(...)
	  for ( size_t j = 0; j < i; ++j )
	    mpfi_mul_si( vec_Uij[j][i], rmatrix[j][i], vec_x[i] );
	}
      if ( i == m - 1 )
	continue;

      break;
    }
}

/*************************************************************
 * initialization and clearing for the algorithm
 *************************************************************/

inline void
init
(
 size_t m,
 vector<mpfi_ptr> &vec_Ti,
 vector<mpfi_ptr> &vec_Ui,
 vector<vector<mpfi_ptr>> &vec_Uij,
 mpfi_ptr &C,
 mpfi_ptr &Z,
 vector<vector<mpfi_ptr>> &rmatrix,
 vector<mpfi_ptr> &rdiag_sqrt,
 mpfr_ptr &mpfr_tmp,
 mpfi_ptr &mpfi_tmp,
 mpfi_ptr &mpfi_tmp2
 )
{
  mpfr_init( mpfr_tmp );
  mpfi_init( mpfi_tmp );
  mpfi_init( mpfi_tmp2 );

  mpfi_init_matrix( rmatrix, m );
  mpfi_init_vector( rdiag_sqrt, m );

  mpfi_init_vector( vec_Ti, m );
  mpfi_init_vector( vec_Ui, m );
  mpfi_init_matrix( vec_Uij, m );

  mpfi_init( C );
  mpfi_init( Z );
}

inline void
clear
(
 vector<mpfi_ptr> &vec_Ti,
 vector<mpfi_ptr> &vec_Ui,
 vector<vector<mpfi_ptr>> &vec_Uij,
 mpfi_ptr &C,
 mpfi_ptr &Z,
 vector<vector<mpfi_ptr>> &rmatrix,
 vector<mpfi_ptr> &rdiag_sqrt,
 mpfr_ptr &mpfr_tmp,
 mpfi_ptr &mpfi_tmp,
 mpfi_ptr &mpfi_tmp2
 )
{
  mpfr_clear( mpfr_tmp );
  mpfi_clear( mpfi_tmp );
  mpfi_clear( mpfi_tmp2 );

  mpfi_clear_matrix( rmatrix );
  mpfi_clear_vector( rdiag_sqrt );

  mpfi_clear_vector( vec_Ti );
  mpfi_clear_vector( vec_Ui );
  mpfi_clear_matrix( vec_Uij );

  mpfi_clear( C );
  mpfi_clear( Z );
}

/*************************************************************
 * arithmetic helper functions
 *************************************************************/

inline
bool
mpfi_get_unique_si
(
 int &si,
 mpfi_ptr srcptr,
 mpfr_ptr mpfr_tmp 
 )
{
  mpfi_get_left( mpfr_tmp, srcptr );
  mpfr_ceil( mpfr_tmp, mpfr_tmp );
  si = mpfr_get_si( mpfr_tmp, MPFR_RNDD );

  mpfi_get_right( mpfr_tmp, srcptr );
  mpfr_floor( mpfr_tmp, mpfr_tmp );
  if ( si != mpfr_get_si( mpfr_tmp, MPFR_RNDU ) )
    return false;

  return true;
}

/*
 * This is almost the floor of the lower and the upper bound.
 * However consider the following case: The exact result is an
 * integer. The numerical calculations, no matter how exact, will
 * produce a small interval around the exact value, so that the
 * floor of the left and right boundary are never the same.
 */

inline
bool
mpfi_get_unique_floor_si
(
 int &si,
 mpfi_ptr srcptr,
 mpfr_ptr mpfr_tmp,
 bool round_down
 )
{
  mpfi_get_left( mpfr_tmp, srcptr );
  mpfr_floor( mpfr_tmp, mpfr_tmp );
  int si_lower = mpfr_get_si( mpfr_tmp, MPFR_RNDD );

  mpfi_get_right( mpfr_tmp, srcptr );
  mpfr_floor( mpfr_tmp, mpfr_tmp );
  int si_upper = mpfr_get_si( mpfr_tmp, MPFR_RNDU );
  if ( si_lower == si_upper )
    {
      si = si_lower;
      return true;
    }
  else if ( si_lower + 1 == si_upper )
    {
      mpfi_diam_abs( mpfr_tmp, srcptr );
      if ( mpfr_get_exp( mpfr_tmp ) < -10 ) // standard precision is 53
	{
	  si = round_down ? si_lower : si_upper;
	  return true;
	}
    }
  return false;
}

inline
bool
mpfi_get_unique_ceil_si
(
 int &si,
 mpfi_ptr srcptr,
 mpfr_ptr mpfr_tmp,
 bool round_up
 )
{
  mpfi_get_left( mpfr_tmp, srcptr );
  mpfr_ceil( mpfr_tmp, mpfr_tmp );
  int si_lower = mpfr_get_si( mpfr_tmp, MPFR_RNDD );

  mpfi_get_right( mpfr_tmp, srcptr );
  mpfr_ceil( mpfr_tmp, mpfr_tmp );
  int si_upper = mpfr_get_si( mpfr_tmp, MPFR_RNDU );
  if ( si_lower == si_upper )
    {
      si = si_lower;
      return true;
    }
  else if ( si_lower + 1 == si_upper )
    {
      mpfi_diam_abs( mpfr_tmp, srcptr );
      if ( mpfr_get_exp( mpfr_tmp ) < -10 ) // standard precision is 53
	{
	  si = round_up ? si_upper : si_lower;
	  return true;
	}
    }
  return false;
}

/*************************************************************
 * initialization and clearing of mpfi objects
 *************************************************************/

inline
void
mpfi_init_vector
(
 vector<mpfi_ptr> &vec,
 size_t length )
{
  for ( size_t i = 0; i < length; ++i )
    {
      mpfi_ptr tmp = new __mpfi_struct;
      mpfi_init( tmp );
      vec.push_back( tmp );
    }
}

inline
void
mpfi_init_matrix
(
 vector<vector<mpfi_ptr>> &mat,
 size_t size )
{
  for ( size_t i = 0; i < size; ++i )
    {
      auto row = vector<mpfi_ptr>();
      mpfi_init_vector( row, size );
      mat.push_back( row );
    }
}

inline
void
mpfi_clear_vector
(
 vector<mpfi_ptr>& vec
)
{
  for ( auto &ptr : vec )
    {
      mpfi_clear( ptr );
      ptr = nullptr;
    }
}

inline
void
mpfi_clear_matrix
(
 vector<vector<mpfi_ptr>>& mat
)
{
  for ( auto &row : mat )
    mpfi_clear_vector( row );
}

inline
void
mpfi_set_prec_vector
(
 vector<mpfi_ptr> &vec,
 mp_prec_t prec )
{
  for ( auto &ptr : vec )
    mpfi_set_prec( ptr, prec );
}

inline
void
mpfi_set_prec_matrix
(
 vector<vector<mpfi_ptr>>& mat,
 mp_prec_t prec
)
{
  for ( auto &row : mat )
    mpfi_set_prec_vector( row, prec );
}

#endif
