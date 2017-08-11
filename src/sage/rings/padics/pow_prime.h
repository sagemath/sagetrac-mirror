#ifndef _POW_PRIME_H
#define _POW_PRIME_H

#include <stdio.h>
#include "gmp.h"


class PowPrime {
    mpz_t *powers;
    mpz_t current_power;
    mpz_t tmp;

    int cache_loglength;
    long cache_length;

    void update_current_power(const long exponent);

public:
    PowPrime();
    ~PowPrime();

    void init(mpz_t prime, const int cache_loglength);

    mpz_t* power(const long);

    // Use if you need
    mpz_t* temp1;
    mpz_t* temp2;

};


#endif
