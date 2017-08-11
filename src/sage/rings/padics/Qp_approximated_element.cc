#include <stdio.h>
#include <string.h>
#include "Qp_approximated_element.h"
#include "gmp.h"

#define INFTY 10000

// Constructors

QpApproximatedElement::QpApproximatedElement() {
    this->pow = 0;
    mpz_init(this->value);
    this->is_normalized = 0;
}

QpApproximatedElement::QpApproximatedElement(const QpApproximatedElement &other) {
    this->pow = other.pow;
    mpz_init_set(this->value, other.value);
    this->exponent = other.exponent;
    this->is_normalized = other.is_normalized;
}

QpApproximatedElement::QpApproximatedElement(PowPrime *p) {
    this->pow = p;
    mpz_init(this->value);
    this->is_normalized = 0;
}

QpApproximatedElement::QpApproximatedElement(const mpz_t v, PowPrime *p) {
    this->pow = p;
    mpz_init(this->value);
    this->set_value(v);
}

QpApproximatedElement::QpApproximatedElement(const mpz_t u, const long exp, PowPrime *p) {
    this->pow = p;
    mpz_init(this->value);
    this->set_value(u, exp);
}

inline QpApproximatedElement QpApproximatedElement::cast() const {
    return *this;
}


// Destructor

QpApproximatedElement::~QpApproximatedElement() {
    mpz_clear(this->value);
}


// Printing functions

char* QpApproximatedElement::repr() {
    char *svalue;
    char *sprime;
    char *ans;

    if (mpz_cmp_ui(this->value, 0) == 0) {
        return "0";
    }

    mpz_t *p = this->pow->power(1);
    if (mpz_cmp_ui(*p, 10) == -1) {
        svalue = mpz_get_str(NULL, mpz_get_ui(*p), this->value);
    } else {
        svalue = mpz_get_str(NULL, 10, this->value);
    }

    if (this->exponent == 0) {
        return svalue;
    }

    sprime = mpz_get_str(NULL, 10, *p);

    ans = (char*) malloc((strlen(sprime) + strlen(svalue) + 30) * sizeof(char));
    sprintf(ans, "%s^%d * %s", sprime, (int)this->exponent, svalue);
    return ans;
}


// Normalization

void QpApproximatedElement::normalize() {
    if (this->is_normalized || mpz_cmp_ui(this->value, 0) == 0) {
        return;
    }
    this->exponent += mpz_remove(this->value, this->value, *(this->pow->power(1)));
    this->is_normalized = 1;
}


// Set functions

void QpApproximatedElement::operator=(const QpApproximatedElement &other) {
    this->pow = other.pow;
    mpz_set(this->value, other.value);
    this->exponent = other.exponent;
    this->is_normalized = other.is_normalized;
    this->pow = other.pow;
}

inline void QpApproximatedElement::set_value(const mpz_t v) {
    mpz_set(this->value, v);
    this->exponent = 0;
    this->is_normalized = 0;
}

inline void QpApproximatedElement::set_value(const mpz_t u, const long exp) {
    mpz_set(this->value, u);
    this->exponent = exp;
    this->is_normalized = 0;
}

inline void QpApproximatedElement::set_zero() {
    mpz_set_ui(this->value, 0);
    this->exponent = INFTY;
    this->is_normalized = 1;
}

inline void QpApproximatedElement::set_one() {
    mpz_set_ui(this->value, 1);
    this->exponent = 0;
    this->is_normalized = 1;
}


// Comparison functions

inline bool QpApproximatedElement::is_zero() {
    return mpz_cmp_ui(this->value, 0) == 0;
}


// Valuation

long QpApproximatedElement::valuation() {
    this->normalize();
    return this->exponent;
}

void QpApproximatedElement::ireduce(const long prec, bool absolute) {
    mpz_t *modulo;
    this->normalize();
    if (absolute) {
        modulo = this->pow->power(prec - this->exponent);
    } else {
        modulo = this->pow->power(prec);
    }
    mpz_fdiv_r(this->value, this->value, *modulo);
}

void QpApproximatedElement::operator+=(QpApproximatedElement &other) {
    mpz_t* tmp;
    long diffval = this->exponent - other.exponent;
    if (diffval == 0) {
        mpz_add(this->value, this->value, other.value);
        this->is_normalized = 0;
    } else if (diffval > 0) {
        mpz_mul(this->value, this->value, *(this->pow->power(diffval)));
        mpz_add(this->value, this->value, other.value);
        this->exponent = other.exponent;
        this->is_normalized = other.is_normalized;
    } else {
        tmp = this->pow->temp1;
        mpz_mul(*tmp, other.value, *(this->pow->power(-diffval)));
        mpz_add(this->value, this->value, *tmp);
    }
}

void QpApproximatedElement::iadd(QpApproximatedElement &other, const long prec, bool absolute) {
    mpz_t* tmp;
    long diffval = this->exponent - other.exponent;
    long powred;

    if (diffval == 0) {
        mpz_add(this->value, this->value, other.value);
        this->is_normalized = 0;
    } else if (diffval > 0) {
        if (absolute) {
            powred = prec - this->exponent;
        } else {
            powred = prec - this->exponent + other.exponent;
        }
        if (powred > 0) {
            mpz_fdiv_r(this->value, this->value, *(this->pow->power(powred)));
            mpz_mul(this->value, this->value, *(this->pow->power(diffval)));
            mpz_add(this->value, this->value, other.value);
        } else {
            mpz_set(this->value, other.value);
        }
        this->exponent = other.exponent;
        this->is_normalized = other.is_normalized;
    } else {
        if (absolute) {
            powred = prec - other.exponent;
        } else {
            powred = prec - other.exponent + this->exponent;
        }
        if (powred > 0) {
            tmp = this->pow->temp1;
            mpz_fdiv_r(*tmp, other.value, *(this->pow->power(powred)));
            mpz_mul(*tmp, *tmp, *(this->pow->power(-diffval)));
            mpz_add(this->value, this->value, *tmp);
        }
    }
    this->ireduce(prec, absolute);
}

void QpApproximatedElement::operator-=(QpApproximatedElement &other) {
    mpz_t* tmp;
    long diffval = this->exponent - other.exponent;
    if (diffval == 0) {
        mpz_sub(this->value, this->value, other.value);
        this->is_normalized = 0;
    } else if (diffval > 0) {
        mpz_mul(this->value, this->value, *(this->pow->power(diffval)));
        mpz_sub(this->value, this->value, other.value);
        this->exponent = other.exponent;
        this->is_normalized = other.is_normalized;
    } else {
        tmp = this->pow->temp1;
        mpz_mul(*tmp, other.value, *(this->pow->power(-diffval)));
        mpz_sub(this->value, this->value, *tmp);
    }
}

inline void QpApproximatedElement::ineg() {
    mpz_neg(this->value, this->value);
}

void QpApproximatedElement::isub(QpApproximatedElement &other, const long prec, bool absolute) {
    mpz_t* tmp;
    long diffval = this->exponent - other.exponent;
    long powred;

    if (diffval == 0) {
        mpz_sub(this->value, this->value, other.value);
        this->is_normalized = 0;
    } else if (diffval > 0) {
        if (absolute) {
            powred = prec - this->exponent;
        } else {
            powred = prec - this->exponent + other.exponent;
        }
        if (powred > 0) {
            mpz_fdiv_r(this->value, this->value, *(this->pow->power(powred)));
            mpz_mul(this->value, this->value, *(this->pow->power(diffval)));
            mpz_sub(this->value, this->value, other.value);
        } else {
            mpz_set(this->value, other.value);
        }
        this->exponent = other.exponent;
        this->is_normalized = other.is_normalized;
    } else {
        if (absolute) {
            powred = prec - other.exponent;
        } else {
            powred = prec - other.exponent + this->exponent;
        }
        if (powred > 0) {
            tmp = this->pow->temp1;
            mpz_fdiv_r(*tmp, other.value, *(this->pow->power(powred)));
            mpz_mul(*tmp, *tmp, *(this->pow->power(-diffval)));
            mpz_sub(this->value, this->value, *tmp);
        }
    }
    this->ireduce(prec, absolute);
}

void QpApproximatedElement::operator*=(QpApproximatedElement &other) {
    mpz_mul(this->value, this->value, other.value);
    this->exponent += other.exponent;
    this->is_normalized &= other.is_normalized;
}

void QpApproximatedElement::iinvert(const long prec, bool absolute) {
    mpz_t gcd, *modulo;
    this->normalize();
    if (absolute) {
        modulo = this->pow->power(prec - 2*this->exponent);
    } else {
        modulo = this->pow->power(prec);
    }
    mpz_init(gcd);
    mpz_gcdext(gcd, this->value, NULL, this->value, *modulo);
    mpz_clear(gcd);
    this->exponent = -this->exponent;
}

inline void QpApproximatedElement::operator <<=(const long n) {
    this->exponent += n;
}

inline void QpApproximatedElement::operator >>=(const long n) {
    this->exponent -= n;
}

inline void QpApproximatedElement::iunit() {
    this->normalize();
    this->exponent = 0;
}
