#ifndef _QP_APPROXIMATED_ELEMENT_H
#define _QP_APPROXIMATED_ELEMENT_H


#include "padic_approximated_element.h"
#include "pow_prime.h"
#include "gmp.h"


class QpApproximatedElement : public pAdicFieldApproximatedElement<QpApproximatedElement> {
    mpz_t value;
    long exponent;
    bool is_normalized;

    PowPrime *pow;

public:

    QpApproximatedElement();
    QpApproximatedElement(const QpApproximatedElement&);
    QpApproximatedElement(PowPrime*);
    QpApproximatedElement(const mpz_t v, PowPrime*);
    QpApproximatedElement(const mpz_t u, const long exp, PowPrime*);
    // QpApproximatedElement(const mpq_t v, const long prec, PowPrime&);
    inline QpApproximatedElement cast() const;

    ~QpApproximatedElement();

    char* repr();

    void normalize();

    void operator=(const QpApproximatedElement &);
    inline void set_value(const mpz_t v);
    inline void set_value(const mpz_t u, const long exp);
    inline void set_zero();
    inline void set_one();

    inline bool is_zero();
    bool is_zero_at_prec(const long) const;
    inline bool operator==(const QpApproximatedElement&) const;
    bool is_equal_to_at_prec(const QpApproximatedElement&, const long) const;

    long valuation();

    void ireduce(const long, bool);
    void operator+=(QpApproximatedElement&);
    void iadd(QpApproximatedElement&, const long, bool);
    inline void ineg();
    void operator-=(QpApproximatedElement&);
    void isub(QpApproximatedElement&, const long, bool);
    void operator*=(QpApproximatedElement&);
    void iinvert(const long, bool);
    inline void operator <<=(const long n);
    inline void operator >>=(const long n);
    inline void iunit();

};


#endif
