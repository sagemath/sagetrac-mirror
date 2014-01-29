/*
 * This file is part of decoding
 * Copyright 2012  Guillaume Quintin, INRIA
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <decoding/decoding.h>

/**
The following functions allow one to compute the parameters to
give the Sudan and Guruswami-Sudan algorithms. They apply to
all Reed-Solomon codes over any kind of ring. That's the reason
they take as arguments only the length and the dimension of
the Reed-Solomon code and, for some, the number of errors to
be corrected.
*/

/* Sudan stuff here */

/**
> `int sudan_max_errors(int n,int k);`
  Return the maximum number of errors that can be corrected by the
  Sudan algorithm for a Reed-Solomon code of length `n` and
  dimension `k`. (Lemma 8, page 6 of Decoding Reed-Solomon codes
  beyond the Error-Correction Diameter, Madhu Sudan.)
*/

int sudan_max_errors(int n,int k) {
  double nn,kk,rkn,tau,num,denom;
  nn = (double)n;
  kk = (double)k;
  rkn = floor(sqrt(((2.0 * (nn + 1.0)) / (kk - 1.0)) + .25) - .5);
  num = 2.0 * nn - (kk - 1) * rkn * (rkn + 1);
  denom = 2.0 * (rkn + 1.0);
  tau = nn - (kk - 1.0) * rkn - ceil(num / denom) - 2.0;
  return (int)tau;
}

/**
> `int sudan_list_size(int n,int k,int tau);`
  Return the list-size for the Sudan algorithm given the length
  `n` and the dimension `k` of the Reed-Solomon code, and the
  number of errors `tau` that one wants to correct.
*/

int sudan_list_size(int n,int k,int tau) {
  double nn,kk,tt;
  nn = (double)n;
  kk = (double)k;
  tt = (double)tau;
  return (int)floor((nn - tt - 1.0) / (kk - 1.0));
}

/* Guruswami-Sudan here */

/**
> `int guruswami_sudan_max_errors(int n,int k);`
  Return the maximum number of errors that can be corrected by the
  Guruswami-Sudan algorithm for a Reed-Solomon code of length `n` and
  dimension `k`. (Theorem 8, page 8 of Improved Decoding of
  Reed-Solomon and Algebraic-Geometry Codes, Venkatesan Guruswami
  and Madhu Sudan.)
*/

int guruswami_sudan_max_errors(int n,int k) {
  double nn,kk;
  nn = (double)n;
  kk = (double)k;
  return (int)(ceil(nn - sqrt(nn * (kk - 1.0))) - 1.0);
}

/**
> `int guru_sudan_mult(int n,int k,int tau);`
  Return the multitplicity of the interpolating curve
  of the Guruswami-Sudan algorithm for a Reed-Solomon code
  of length `n` and dimensions `k` where `tau` errors have
  to be corrected. (Lemma 7, page 8 of Improved Decoding of
  Reed-Solomon and Algebraic-Geometry Codes, Venkatesan Guruswami
  and Madhu Sudan.)
*/

int guruswami_sudan_mult(int n,int k,int tau) {
  double nn,kk,kknn,tt,rr,num,denom;
  nn = (double)n;
  kk = (double)k;
  tt = nn - (double)tau;
  kknn = (kk - 1.0) * nn;
  num = kknn + sqrt((kknn * kknn + 4.0 * (tt * tt - kknn)));
  denom = 2.0 * (tt * tt - kknn);
  rr = 1.0 + floor(num / denom);
  return (int)rr;
}

/**
> `int guruswami_sudan_list_size(int n,int k,int tau,int m);`
  Return the list-size for the Guruswami-Sudan algorithm given the
  length `n` and the dimension `k` of the Reed-Solomon code, and the
  number of errors `tau` that one wants to correct with multiplicity
  `m`.
*/

int guruswami_sudan_list_size(int n,int k,int tau,int m) {
  double nn,kk,tt,mm;
  nn = (double)n;
  kk = (double)k;
  tt = (double)tau;
  mm = (double)m;
  return (int)(floor((mm * (nn - tt) - 1.0) / (kk - 1)));
}

/**
> `int guruswami_sudan_errors(int n,int k,int m);`
  Return the maximum number of errors that can be corrected
  by the Guruswami-Sudan algorithm given a Reed-Solomon code
  of length `n`, dimension `k` and the multiplicity `m`.
*/

int guruswami_sudan_errors(int n,int k,int m) {
  double nn,kk,mm,num;
  nn = (double)n;
  kk = (double)k;
  mm = (double)m;
  num = (kk - 1.0) * nn * (mm + 1.0) + 1.0;
  return (int)(ceil((nn - sqrt(num / mm)) - 1.0));
}

