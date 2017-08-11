#include "Zp_approximated_element.h"
#include "pow_prime.h"
#include "gmp.h"

#define INFTY 10000



// Constructors

ZpApproximatedElement::ZpApproximatedElement() {
    mpz_init(this->value);
    this->pow = 0;
    this->is_valuation_cached = 0;
    this->is_unit_cached = 0;
}

ZpApproximatedElement::ZpApproximatedElement(const ZpApproximatedElement &other) {
    mpz_init_set(this->value, other.value);
    this->cached_valuation = other.cached_valuation;
    this->is_valuation_cached = other.is_valuation_cached;
    this->is_unit_cached = 0;  // Observe that we do not copy the unit
    this->pow = other.pow;
}

void ZpApproximatedElement::operator=(const ZpApproximatedElement &other) {
    mpz_init_set(this->value, other.value);
    this->cached_valuation = other.cached_valuation;
    this->is_valuation_cached = other.is_valuation_cached;
    this->is_unit_cached = 0;  // Observe that we do not copy the unit
    this->pow = other.pow;
}

ZpApproximatedElement::ZpApproximatedElement(PowPrime *p) {
    this->pow = p;
    mpz_init(this->value);
    this->is_valuation_cached = 0;
    this->is_unit_cached = 0;
}

ZpApproximatedElement::ZpApproximatedElement(const mpz_t v, PowPrime *p) {
    this->pow = p;
    mpz_init(this->value);
    this->set_value(v);
    this->is_valuation_cached = 0;
    this->is_unit_cached = 0;
}

ZpApproximatedElement::ZpApproximatedElement(const mpz_t v, const long prec, bool absolute, PowPrime *p) {
    this->pow = p;
    mpz_init(this->value);
    this->set_value(v);
    this->ireduce(prec, absolute);
}

inline ZpApproximatedElement ZpApproximatedElement::cast() const {
    return *this;
}


// Destructor

ZpApproximatedElement::~ZpApproximatedElement() {
    mpz_clear(this->value);
    if (this->is_unit_cached) {
        mpz_clear(this->cached_unit);
    }
}


// Printing functions

char* ZpApproximatedElement::repr() {
    mpz_t *p = this->pow->power(1);
    if (mpz_cmp_ui(*p, 10) == -1) {
        return mpz_get_str(NULL, mpz_get_ui(*p), this->value);
    } else {
        return mpz_get_str(NULL, 10, this->value);
    }
}


// Set functions

inline void ZpApproximatedElement::set_value(const mpz_t v) {
    mpz_set(this->value, v);
    this->is_valuation_cached = 0;
    this->is_unit_cached = 0;
}

inline void ZpApproximatedElement::set_zero() {
    mpz_set_ui(this->value, 0);
    this->cached_valuation = INFTY;
    this->is_valuation_cached = 1;
    this->is_unit_cached = 0;
}

inline void ZpApproximatedElement::set_one() {
    mpz_set_ui(this->value, 1);
    this->cached_valuation = 0;
    this->is_valuation_cached = 1;
    this->is_unit_cached = 0;
}

inline void ZpApproximatedElement::set_valuation(const long val) {  // We trust the caller
    this->cached_valuation = val;
    this->is_valuation_cached = 1;
}


// Comparison functions

inline bool ZpApproximatedElement::is_zero() {
    return mpz_cmp_ui(this->value, 0);
}

bool ZpApproximatedElement::is_zero_at_prec(const long prec) const {
    if (this->is_valuation_cached) {
        return (this->cached_valuation >= prec);
    }
    mpz_t *tmp = this->pow->temp1;
    mpz_mod(*tmp, this->value, *(this->pow->power(prec)));
    bool ans = mpz_cmp_ui(*tmp, 0);
    return ans;
}

inline bool ZpApproximatedElement::operator==(const ZpApproximatedElement &other) const {
    return mpz_cmp(this->value, other.value);
}

bool ZpApproximatedElement::is_equal_to_at_prec(const ZpApproximatedElement &other, const long prec) const {
    mpz_t *tmp = this->pow->temp1;
    mpz_sub(*tmp, this->value, other.value);
    mpz_mod(*tmp, *tmp, *(this->pow->power(prec)));
    bool ans = mpz_cmp_ui(*tmp, 0);
    return ans;
}


// Valuation

void ZpApproximatedElement::val_unit() {
    // WARNING: do not call this function if cached_unit is already initialized!
    mpz_init(this->cached_unit);
    if (mpz_cmp_ui(this->value, 0) == 0) {
        this->cached_valuation = INFTY;
        mpz_set_ui(this->cached_unit, 0);
    } else {
        // Should we rewrite this function in order to take advantage of pow?
        this->cached_valuation = (long)mpz_remove(this->cached_unit, this->value, *(pow->power(1)));
    }
    this->is_valuation_cached = 1;
    this->is_unit_cached = 1;
}

long ZpApproximatedElement::valuation() {
    if (! this->is_valuation_cached) {
        this->val_unit();
    }
    return this->cached_valuation;
}


// Inplace operations

void ZpApproximatedElement::ireduce(const long prec, bool absolute) {
    if (mpz_cmp_ui(this->value, 0) == 0) {
        return;
    }
    if (absolute) {
        mpz_mod(this->value, this->value, *(pow->power(prec)));
        if (this->cached_valuation > prec) {
            this->cached_valuation = INFTY;
        }
    } else {
        mpz_mod(this->value, this->value, *(pow->power(prec + this->valuation())));
    }
}

void ZpApproximatedElement::operator+=(ZpApproximatedElement &other) {
    mpz_add(this->value, this->value, other.value);
    this->is_valuation_cached = 0;
    this->is_unit_cached = 0;
}

inline void ZpApproximatedElement::ineg() {
    mpz_neg(this->value, this->value);
}

void ZpApproximatedElement::operator-=(ZpApproximatedElement &other) {
    mpz_sub(this->value, this->value, other.value);
    this->is_valuation_cached = 0;
    this->is_unit_cached = 0;
}

void ZpApproximatedElement::operator*=(ZpApproximatedElement &other) {
    mpz_mul(this->value, this->value, other.value);
    this->cached_valuation += other.cached_valuation;
    this->is_valuation_cached &= other.is_valuation_cached;
    this->is_unit_cached = 0;
}

void ZpApproximatedElement::iinvert(const long prec, bool absolute) {
    // TODO: implement Newton iteration
    mpz_t *gcd = this->pow->temp1;
    mpz_t *modulo;

    modulo = pow->power(prec);
    mpz_gcdext(*gcd, this->value, 0, this->value, *modulo);
    if (mpz_cmp_ui(*gcd,1) != 0) {
        abort(); // inverse of non-unit
    }
    this->cached_valuation = 0;
    this->is_valuation_cached = 1;
    this->is_unit_cached = 0;
}

void ZpApproximatedElement::idiv(ZpApproximatedElement &other, const long prec, bool absolute) {
    long valother = other.valuation();
    mpz_t *tmp = this->pow->temp1;
    mpz_t *inverse = this->pow->temp2;
    mpz_t *modulo;

    if (valother > 0) {
        mpz_fdiv_qr(this->value, *tmp, this->value, *(pow->power(valother)));
        if (mpz_cmp_ui(*tmp,0) != 0) {
            abort(); // quotient is not in Zp
        }
    }
    if (absolute) {
        modulo = pow->power(prec+valother);
    } else {
        modulo = pow->power(prec);
    }
    mpz_gcdext(*tmp, *inverse, 0, other.cached_unit, *modulo);
    /*
    if (mpz_cmp_ui(inverse, 0) == -1) {
        mpz_add(inverse, inverse, *modulo);
    }
    */
    mpz_mul(this->value, this->value, *inverse);
    this->ireduce(prec, absolute);
    this->cached_valuation -= valother;
    this->is_unit_cached = 0;
}

inline void ZpApproximatedElement::operator <<=(const long n) {
    mpz_mul(this->value, this->value, *(pow->power(n)));
    this->cached_valuation += n;
}

inline void ZpApproximatedElement::operator >>=(const long n) {
    mpz_fdiv_q(this->value, this->value, *(pow->power(n)));
    if (this->is_valuation_cached) {
        if (this->cached_valuation < n) {
            this->cached_valuation = 0;
        } else {
            this->cached_valuation -= n;
        }
    }
}

inline void ZpApproximatedElement::iunit() {
    if (mpz_cmp_ui(this->value, 0) != 0) {
        if (! this->is_unit_cached) {
            this->val_unit();
        }
        mpz_set(this->value, this->cached_unit);
        this->cached_valuation = 0;
    }
}
