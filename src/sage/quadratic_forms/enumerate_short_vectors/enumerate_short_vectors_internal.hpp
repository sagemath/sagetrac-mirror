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
