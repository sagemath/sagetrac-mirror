#ifndef _PADIC_APPROXIMATED_ELEMENT_H
#define _PADIC_APPROXIMATED_ELEMENT_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

template<class DerivedClass>
class pAdicApproximatedElement {

public:
    // Constructors
    // to be implemented in derived classes

    pAdicApproximatedElement() { };
    virtual DerivedClass cast() const = 0;

    // Set functions

    virtual void set_zero() = 0;
    virtual void set_one() = 0;


    // Print functions

    virtual char* repr() = 0;


    // Comparison functions

    virtual bool is_zero() = 0;      // We are not sure that this will be const in all subclasses

    bool is_zero_at_prec(const long prec) const {
        return this->reduce(prec).is_zero();
    }

    bool operator==(const DerivedClass &other) const {
        return (*this-other).is_zero();
    }

    bool is_equal_to_at_prec(const DerivedClass &other, const long prec) const {
        pAdicApproximatedElement<DerivedClass> tmp = *this - other;
        tmp.ireduce(tmp);
        return tmp.is_zero();
    }


    // Valuation

    virtual long valuation() = 0;    // We are not sure that this will be const in all subclasses


    // Inplace operations

    virtual void ireduce(const long, bool) = 0;
    virtual void operator+=(DerivedClass&) = 0;
    virtual void ineg() = 0;
    virtual void operator-=(DerivedClass&) = 0;
    virtual void operator*=(DerivedClass&) = 0;
    virtual void iinvert(const long, bool) = 0;
    virtual void idiv(DerivedClass&, const long, bool) = 0;
    virtual void operator<<=(const long) = 0;
    virtual void operator>>=(const long) = 0;

    inline void iadd(DerivedClass &other, const long prec, bool absolute) {
        *this += other;
        this->ireduce(prec, absolute);
    }

    inline void isub(DerivedClass &other, const long prec, bool absolute) {
        *this -= other;
        this->ireduce(prec, absolute);
    }

    inline void imul(DerivedClass &other, const long prec, bool absolute) {
        *this *= other;
        this->ireduce(prec, absolute);
    }

    inline void iunit() {
        *this >>= this->valuation();
    }


    // Non inplace operations

    DerivedClass reduce(const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.ireduce(prec, absolute);
        return ans;
    }

    DerivedClass operator+(DerivedClass &other) {
        DerivedClass ans = this->cast();
        ans += other;
        return ans;
    }

    DerivedClass add(DerivedClass &other, const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.iadd(other, prec, absolute);
        return ans;
    }

    DerivedClass neg() {
        DerivedClass ans = this->cast();
        ans.ineg();
        return ans;
    }

    DerivedClass neg(const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.ineg();
        ans.ireduce(prec, absolute);
        return ans;
    }

    DerivedClass operator-(DerivedClass &other) {
        DerivedClass ans = this->cast();
        ans -= other;
        return ans;
    }

    DerivedClass sub(DerivedClass &other, const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.isub(other, prec, absolute);
        return ans;
    }

    DerivedClass operator*(DerivedClass &other) {
        DerivedClass ans = this->cast();
        ans *= other;
        return ans;
    }

    DerivedClass mul(DerivedClass &other, const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.imul(other, prec, absolute);
        return ans;
    }

    DerivedClass invert(const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.iinvert(prec, absolute);
        return ans;
    }

    DerivedClass div(DerivedClass &other, const long prec, bool absolute) {
        DerivedClass ans = this->cast();
        ans.idiv(other, prec, absolute);
        return ans;
    }

    DerivedClass operator<<(const long n) {
        DerivedClass ans = this->cast();
        ans <<= n;
        return ans;
    }

    DerivedClass operator>>(const long n) {
        DerivedClass ans = this->cast();
        ans >>= n;
        return ans;
    }

    DerivedClass unit() {
        DerivedClass ans = this->cast();
        ans.iunit();
        return ans;
    }

};


template<class DerivedClass>
class pAdicRingApproximatedElement : public pAdicApproximatedElement<DerivedClass> {

public:

    void idiv(DerivedClass &other, const long prec, bool absolute) {
        long valthis = this->valuation();
        long valother = other.valuation();
        if (valother > valthis) {
            abort();   // quotient is not integral
        }
        DerivedClass unit = other.unit();
        if (absolute) {
            unit.iinvert(prec-valthis+valother, 0);
        } else {
            unit.iinvert(prec, 0);
        }
        *this *= unit;
        *this >>= valother;
        this->reduce(prec, absolute);
    }

};


template<class DerivedClass>
class pAdicFieldApproximatedElement : public pAdicApproximatedElement<DerivedClass> {

public:

    void idiv(DerivedClass &other, const long prec, bool absolute) {
        DerivedClass inverse = other.cast();
        if (absolute) {
            inverse.iinvert(prec - this->valuation() + other.valuation(), 1);
        } else {
            inverse.iinvert(prec, 0);
        }
        this->imul(inverse, prec, absolute);
    }

};





#endif
