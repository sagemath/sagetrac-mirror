#ifndef _ZP_APPROXIMATED_ELEMENT_H
#define _ZP_APPROXIMATED_ELEMENT_H


#include "padic_approximated_element.h"
#include "pow_prime.h"
#include "gmp.h"


class ZpApproximatedElement : public pAdicRingApproximatedElement<ZpApproximatedElement> {
    mpz_t value;
    mpz_t cached_unit;
    long cached_valuation;
    bool is_valuation_cached;
    bool is_unit_cached;
    PowPrime *pow;

    void val_unit();

public:

    ZpApproximatedElement();
    ZpApproximatedElement(const ZpApproximatedElement&);
    void operator=(const ZpApproximatedElement&);
    ZpApproximatedElement(PowPrime*);
    ZpApproximatedElement(const mpz_t v, PowPrime*);
    ZpApproximatedElement(const mpz_t v, const long, bool, PowPrime*);
    inline ZpApproximatedElement cast() const;

    ~ZpApproximatedElement();

    char* repr();

    inline void set_value(const mpz_t);
    inline void set_zero();
    inline void set_one();
    inline void set_valuation(const long);

    inline bool is_zero();
    bool is_zero_at_prec(const long) const;
    inline bool operator==(const ZpApproximatedElement&) const;
    bool is_equal_to_at_prec(const ZpApproximatedElement&, const long) const;

    long valuation();

    void ireduce(const long, bool);
    void operator+=(ZpApproximatedElement&);
    inline void ineg();
    void operator-=(ZpApproximatedElement&);
    void operator*=(ZpApproximatedElement&);
    void iinvert(const long, bool);
    void idiv(ZpApproximatedElement&, const long, bool);
    inline void operator <<=(const long n);
    inline void operator >>=(const long n);
    inline void iunit();

};


#endif
