#include <stdio.h>
#include <stdlib.h>
#include "pow_prime.h"
#include "gmp.h"


// Constructor

PowPrime::PowPrime() {
    puts("constructor powprime");
    mpz_init(this->current_power);
    mpz_init(this->tmp);
    this->temp1 = (mpz_t*) (malloc(2*sizeof(mpz_t)));
    this->temp2 = this->temp1 + 1;
    this->cache_length = -1;
}


// Initialization

void PowPrime::init(mpz_t prime, const int cache_loglength) {
    puts("initialization powprime");
    if (this->cache_length >= 0) {
        abort(); // already initialized
    }

    this->cache_loglength = cache_loglength;
    this->cache_length = 1 << cache_loglength;

    this->powers = (mpz_t*) (malloc(((this->cache_length) + 1) * sizeof(mpz_t)));
    mpz_init_set_ui(powers[0], 1);
    for (long i = 1; i <= this->cache_length; i++) {
        mpz_init(this->powers[i]);
        mpz_mul(this->powers[i], this->powers[i-1], prime);
    }
    printf("cache length: %d\n", (int)(this->cache_length));
}


// Destructor

PowPrime::~PowPrime() {
    mpz_clear(this->current_power);
    mpz_clear(tmp);
    mpz_clear(*(this->temp1));
    mpz_clear(*(this->temp2));
    for (long i = 0; i <= this->cache_length; i++) {
        mpz_clear(this->powers[i]);
    }
    free(this->temp1);
    if (this->cache_length >= 0) {
        free(this->powers);
    }
}


// Compute powers of p

void PowPrime::update_current_power(const long exponent) {
    long init_exponent = exponent & (this->cache_length - 1);
    long exp = exponent >> this->cache_loglength;

    mpz_set(this->current_power, this->powers[init_exponent]);
    mpz_set(this->tmp, this->powers[this->cache_length]);

    while (1) {
        if ((exp & 1) == 1) {
            mpz_mul(this->current_power, this->current_power, this->tmp);
        }
        exp >>= 1;
        if (exp == 0) break;
        mpz_mul(this->tmp, this->tmp, this->tmp);
    }
}

mpz_t* PowPrime::power(const long exponent) {
    if (exponent <= this->cache_length) {
        return this->powers + exponent;
    } else {
        this->update_current_power(exponent);
        return &(this->current_power);
    }
}
