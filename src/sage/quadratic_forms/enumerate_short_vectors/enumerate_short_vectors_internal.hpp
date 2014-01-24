#ifndef __ENUMERATE_SHORT_VECTORS_INTERNAL_H
#define __ENUMERATE_SHORT_VECTORS_INTERNAL_H

#include "mpfr.h"
#include "mpfi.h"

bool
cholesky_decomposition( const std::vector<std::vector<int>>&,
			std::vector<std::vector<mpfi_ptr>>, std::vector<mpfi_ptr>,
			mpfi_ptr );

bool
step_2( size_t, std::vector<int>&, int, int,
	std::vector<int>&, std::vector<int>&, int&, int&, mpfi_ptr&,
	std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
	std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
	mpfi_ptr&, mpfi_ptr&, mpfr_ptr& );

void
recompute( size_t&, size_t, std::vector<int>&, unsigned int, unsigned int,
	   std::vector<int>&, std::vector<int>&, int&, int&,
	   std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
	   mpfi_ptr&, mpfi_ptr&,
	   const std::vector<std::vector<int>>&, std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
	   mpfi_ptr&, mpfi_ptr&, mpfr_ptr&,
	   mp_prec_t );

void
init( size_t, std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
      mpfi_ptr&, mpfi_ptr&, std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
      mpfr_ptr&, mpfi_ptr&, mpfi_ptr& );

void
clear( std::vector<mpfi_ptr>&, std::vector<mpfi_ptr>&, std::vector<std::vector<mpfi_ptr>>&,
       mpfi_ptr&, mpfi_ptr&, std::vector<std::vector<mpfi_ptr>>&, std::vector<mpfi_ptr>&,
       mpfr_ptr&, mpfi_ptr&, mpfi_ptr& );

bool mpfi_get_unique_si( int&, mpfi_ptr, mpfr_ptr );
bool mpfi_get_unique_floor_si( int&, mpfi_ptr, mpfr_ptr, bool );
bool mpfi_get_unique_ceil_si( int&, mpfi_ptr, mpfr_ptr, bool );

void mpfi_init_vector( std::vector<mpfi_ptr>&, size_t );
void mpfi_init_matrix( std::vector<std::vector<mpfi_ptr>>&, size_t );
void mpfi_set_prec_vector( std::vector<mpfi_ptr>&, mp_prec_t );
void mpfi_set_prec_matrix( std::vector<std::vector<mpfi_ptr>>&, mp_prec_t );
void mpfi_clear_vector( std::vector<mpfi_ptr>& );
void mpfi_clear_matrix( std::vector<std::vector<mpfi_ptr>>& );

#endif
