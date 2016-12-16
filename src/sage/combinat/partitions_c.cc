/*
 *      Author:     Jonathan Bober
 *      Version:    .6
 *
 *      This program computes p(n), the number of integer partitions of n, using Rademacher's
 *      formula. (See Hans Rademacher, On the Partition Function p(n),
 *      Proceedings of the London Mathematical Society 1938 s2-43(4):241-254; doi:10.1112/plms/s2-43.4.241,
 *      currently at
 *
 *      http://plms.oxfordjournals.org/cgi/content/citation/s2-43/4/241
 *
 *      if you have access.)
 *
 *      We use the following notation:
 *
 *      p(n) = lim_{n --> oo} t(n,N)
 *
 *      where
 *
 *      t(n,N) = sum_{k=1}^N a(n,k) f_n(k),
 *
 *      where
 *
 *      a(n,k) = sum_{h=1, (h,k) = 1}^k exp(\pi i s(h,k) - 2 \pi i h  n / k)
 *
 *      and
 *
 *      f_n(k) = \pi sqrt{2} cosh(A_n/(sqrt{3}*k))/(B_n*k) - sinh(C_n/k)/D_n;
 *
 *      where
 *
 *      s(h,k) = \sum_{j=1}^{k-1}(j/k)((hj/k))
 *
 *      A_n = sqrt{2} \pi * sqrt{n - 1/24}
 *      B_n = 2 * sqrt{3} * (n - 1/24)
 *      C_n = sqrt{2} * \pi * sqrt{n - 1.0/24.0} / sqrt{3}
 *      D_n = 2 * (n - 1/24) * sqrt{n - 1.0/24.0}
 *
 *      and, finally, where ((x)) is the sawtooth function ((x)) = {x} - 1/2 if x is not an integer, 0 otherwise.
 *
 *      Some clever tricks are used in the computation of s(h,k), and perhaps at least
 *      some of these are due to Apostol. (I don't know a reference for this.)
 *
 *      TODO:
 *
 *          -- Search source code for other TODO comments.
 *
 *      OTHER CREDITS:
 *
 *      I looked source code written by Ralf Stephan, currently available at
 *
 *              http://www.ark.in-berlin.de/part.c
 *
 *      while writing this code, but didn't really make use of it, except for the
 *      function cospi(), currently included in the source (slightly modified), but not
 *      used.
 *
 *      More useful were notes currently available at
 *
 *              http://www.ark.in-berlin.de/part.pdf
 *
 *      and others at
 *
 *              http://www.math.uwaterloo.ca/~dmjackso/CO630/ptnform1.pdf
 *
 *      Also, it is worth noting that the code for GCD(), while trivial,
 *      was directly copied from the NTL source code.
 *
 *      Also, Bill Hart made some comments about ways to speed up this computation on the SAGE
 *      mailing list.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#if defined(__sun) || defined(__CYGWIN__)
extern "C" long double fabsl (long double);
extern "C" long double sinl (long double);
extern "C" long double cosl (long double);
extern "C" long double sqrtl (long double);
extern "C" long double coshl (long double);
extern "C" long double sinhl (long double);
#endif

#include <stdio.h>
#include <cfloat>
#include <time.h>

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <iomanip>
#include <limits>

#include <mpfr.h>
#include <gmp.h>

using namespace std;

using std::cout;
using std::endl;

/*****************************************************************************
 *
 *      We begin be declaring all the constant global variables.
 *
 ****************************************************************************/

// First, some variables that can be set to have the program
// output some information that might be useful for debugging.

const bool debug = false;                                       // If true, output some random stuff

const bool debug_precision = false;                             // If true, output information that might be useful for
                                                                // debugging the precision setting code.

const bool debugs = false;                                      // If true, output information that might be useful for
                                                                // debugging the various s() functions.
const bool debugf = false;                                      // Same for the f() functions.
const bool debuga = false;                                      // Same for the a() functions.
const bool debugt = false;                                      // Same for the t() functions.

#if FLT_RADIX != 2
#error I don't know what to do when the float radix is not 2.
#endif
// Note - it might be unreasonable to not support a float radix other than
// 2, but apparently gcc doesn't either. See http://gcc.gnu.org/ml/fortran/2006-12/msg00032.html.


const unsigned int min_precision = DBL_MANT_DIG;                            // The minimum precision that we will ever use.
const unsigned int double_precision = DBL_MANT_DIG;                         // The assumed precision of a double.


#if defined(__sparc) || defined(__CYGWIN__)  || defined(__FreeBSD__)
// On sparc solaris long double is bad/broken/different, etc.  E.g.,
// LDBL_MANT_DIG is 113 rather than 106, which causes all kinds of trouble.
// So we only use double_precision.
const unsigned int long_double_precision = double_precision;
#else
const unsigned int long_double_precision = (LDBL_MANT_DIG == 106) ? double_precision : LDBL_MANT_DIG;
#endif
                                                                            // The assumed precision of a long double.
                                                                            // Note: On many systems double_precision = long_double_precision. This is OK, as
                                                                            // the long double stage of the computation will just be skipped.
                                                                            //
                                                                            // NOTE: If long_double_precision > dd_precision, then, again, long doubles
                                                                            //      will not ever be used. It would be nice if this were fixed.


// Second, some constants that control the precision at which we compute.

// When we compute the terms of the Rademacher series, we start by computing with
// mpfr_t variables. But for efficiency we switch to special purpose code when
// we don't need as much precision. These constants control the various
// precision levels at which we switch to different special
// purpose functions.

// const unsigned int level_one_precision = infinity        // We don't actually use this, but if we did it would be infinite.
const unsigned int level_two_precision = long_double_precision;
const unsigned int level_three_precision = long_double_precision;
const unsigned int level_four_precision = long_double_precision;
const unsigned int level_five_precision = double_precision;


const long double ld_pi = 3.141592653589793238462643L;
const double d_pi = ld_pi;

// Third, the rounding mode for mpfr.
const mpfr_rnd_t round_mode = MPFR_RNDN;

/*****************************************************************************
 *
 *  We are finished declaring constants, and next declare semi-constant
 *  variables. These are all set just once per call to part(n).
 *
 *  All of these are set in initialize_constants(), and the mpfr_t variables
 *  are cleared in clear_constants().
 *
 ****************************************************************************/

mpfr_t mp_one_over_12, mp_one_over_24, mp_sqrt2, mp_sqrt3, mp_pi, half, fourth;         // These need to be set at run time
                                                                                        // because we don't know how much precision we will need
                                                                                        // until we know for what n we are computing p(n).

mpfr_t mp_A, mp_B, mp_C, mp_D;                                                          // These "constants" all depend on n, and
double d_A, d_B, d_C, d_D;                                                              //
long double ld_A, ld_B, ld_C, ld_D;                                                     //


/*****************************************************************************
 *
 * We next declare the variables that actually vary. It is cumbersome to have
 * so many global variables, but it is somewhat necessary since we want to call
 * the functions mpz_init(), mpq_init(), mpfr_init(), and the corresponding
 * clear() functions as little as possible.
 *
 ****************************************************************************/

mpz_t ztemp1;                                                                   // These are initialized in the function initialize_mpz_and_mpq_variables(), and
mpq_t qtemps, qtempa, qtempa2;                                                  // cleared in the function clear_mpz_and_mpq_variables().
                                                                                // qtemps is used by q_s() and qtempa and qtempa2() are used by mp_a().
                                                                                //
                                                                                // ztemp1 is only used in one function, and the code where is it used is actually unlikely
                                                                                // to be used for any input, so it should be erased.

mpfr_t tempa1, tempa2, tempf1, tempf2, temps1, temps2;                          // These are all initialized by initialize_mpfr_variables(), and cleared by
                                                                                // clear_mpfr_variables().
                                                                                // NAMING CONVENTIONS:
                                                                                //
                                                                                // -tempa1 and tempa2 are two variables available for use
                                                                                //      in the function mp_a()
                                                                                //
                                                                                // -tempf1 and tempf2 are two variables available for use
                                                                                //      in the function mp_f()
                                                                                //
                                                                                // -etc...

mpfr_t tempc1, tempc2;                                                          // temp variables used by cospi()


/*****************************************************************************
 *
 * End of Global Variables, beginning of function declarations
 *
 ****************************************************************************/

//
// A listing of the main functions, in "conceptual" order.
//

int part(mpz_t answer, unsigned int n);

static void mp_t(mpfr_t result, unsigned int n);                                // Compute t(n,N) for an N large enough, and with a high enough precision, so that
                                                                                // |p(n) - t(n,N)| < .5

static unsigned int compute_initial_precision(unsigned int n);                  // computes the precision required to accurately compute p(n)

static void initialize_constants(unsigned int prec, unsigned int n);            // Once we know the precision that we will need, we precompute
                                                                                // some commonly used constants.

static unsigned int compute_current_precision(unsigned int n, unsigned int N, unsigned int extra);  // Computed the precision required to
                                                                                // accurately compute the tail of the rademacher series
                                                                                // assuming that N terms have already been computed.
                                                                                // This is called after computing each summand, unless
                                                                                // we have already reached the minimum precision.
                                                                                //
                                                                                // Reducing the precision as fast as possible is key to
                                                                                // fast performance.

static double compute_remainder(unsigned int n, unsigned int N);                // Gives an upper bound on the error that occurs
                                                                                // when only N terms of the Rademacher series have been
                                                                                // computed.
                                                                                //
                                                                                // This should only be called when we know that compute_current_precision
                                                                                // will return min_precision. Otherwise the error will be too large for
                                                                                // it to compute.


//
// The following functions are all called in exactly one place, so
// we could ask the compiler to inline everything. Actually making this
// work properly requires changing some compiler options, and doesn't
// provide much improvement, however.
//

static void mp_f(mpfr_t result, unsigned int k);                                // See introduction for an explanation of these functions.
static void q_s(mpq_t result, unsigned int h, unsigned int k);                  //
static void mp_a(mpfr_t result, unsigned int n, unsigned int k);                //

template <class T> static inline T a(unsigned int n, unsigned int k);           // Template versions of the above functions for computing with
template <class T> static inline T f(unsigned int k);                           // low precision. Currently used for computing with qd_real, dd_real,
template <class T> static inline T s(unsigned int h, unsigned int k);           // long double, and double.


//
// The following are a bunch of "fancy macros" designed so that in the templated code
// for a, for example, when we need to use pi, we just call pi<T>() to get pi to
// the proper precision. The compiler should inline these, so using one of these
// functions is just as good as using a constant, and since these functions are static
// they shouldn't even appear in the object code generated.

template <class T> static inline T pi(){ return ld_pi; }

template <class T> static inline T sqrt2() {return sqrt(2.0L);}

template <class T> static inline T sqrt3() {return sqrt(3.0L);}

template <class T> static inline T A() {return 1;}
template <> inline double A() {return d_A;}
template <> inline long double A() {return ld_A;}

template <class T> static inline T B() {return 1;}
template <> inline double B() {return d_B;}
template <> inline long double B() {return ld_B;}

template <class T> static inline T C() {return 1;}
template <> inline double C() {return d_C;}
template <> inline long double C() {return ld_C;}

template <class T> static inline T D() {return 1;}
template <> inline double D() {return d_D;}
template <> inline long double D() {return ld_D;}

template <class T> static inline T pi_sqrt2() {return 1;}
template <> inline double pi_sqrt2() {return d_pi * sqrt(2.0);}
template <> inline long double pi_sqrt2() {return ld_pi * sqrt(2.0L);}

template <class T> static inline T one_over_12() {return 1;}
template <> inline double one_over_12() {return 1.0/12.0;}
template <> inline long double one_over_12() {return 1.0L/12.0L;}

// A few utility functions...

static inline long GCD(long a, long b);
int test(bool longtest = false, bool forever = false);                          // Runs a bunch of tests to make sure
                                                                                // that we are getting the right answers.
                                                                                // Tests are based on a few "known" values that
                                                                                // have been verified by other programs, and
                                                                                // also on known congruences for p(n)

static void cospi (mpfr_t res, mpfr_t x);                                       // Puts cos(pi x) in result. This is not currently
                                                                                // used, but should be faster than using mpfr_cos,
                                                                                // according to comments in the file part.c
                                                                                // mentioned in the introduction.
                                                                                // test

static int grab_last_digits(char * output, int n, mpfr_t x);                    // Might be useful for debugging, but
                                                                                // don't use it for anything else.
                                                                                // (See function definition for more information.)


int main(int argc, char *argv[]);                                               // This program is mainly meant for inclusion
                                                                                // in SAGE (or any other program, if anyone
                                                                                // feels like it). We include a main() function
                                                                                // anyway, because it is useful to compile
                                                                                // this as a standalone program
                                                                                // for testing purposes.

/***********************************************************************
 *
 * That should be the end of both function and variable definitions.
 *
 **********************************************************************/

// The following function can be useful for debugging in come circumstances, but should not be used for anything else
// unless it is rewritten.
int grab_last_digits(char * output, int n, mpfr_t x) {
    // fill output with the n digits of x that occur
    // just before the decimal point
    // Note: this assumes that x has enough digits and enough
    // precision -- otherwise bad things can happen
    //
    // returns: the number of digits to the right of the decimal point

    char * temp;
    mp_exp_t e;

    temp = mpfr_get_str(NULL, &e, 10, 0, x, MPFR_RNDN);

    int retval;

    if(e > 0) {
        strncpy(output, temp + e - n, n);
        retval =  strlen(temp + e);
    }
    else {
        for(int i = 0; i < n; i++)
            output[i] = '0';
        retval = strlen(temp);
    }
    output[n] = '\0';

    mpfr_free_str(temp);

    return retval;
}



unsigned int compute_initial_precision(unsigned int n) {
    // We just want to know how many bits we will need to
    // compute to get an accurate answer.

    // We know that
    //
    //          p(n) ~ exp(pi * sqrt(2n/3))/(4n sqrt(3)),
    //
    // so for now we are assuming that p(n) < exp(pi * sqrt(2n/3))/n,
    // so we need pi*sqrt(2n/3)/log(2) - log(n)/log(2) + EXTRA bits.
    //
    // EXTRA should depend on n, and should be something that ensures
    // that the TOTAL ERROR in all computations is < (something small).
    // This needs to be worked out carefully. EXTRA = log(n)/log(2) + 3
    // is probably good enough, and is convenient...
    //
    // but we really need:
    //
    //                  p(n) < something
    //
    // to be sure that we compute the correct answer

    unsigned int result = (unsigned int)(ceil(3.1415926535897931 * sqrt(2.0 * double(n)/ 3.0) / log(2))) + 3;
    if(debug) cout << "Using initial precision of " << result << " bits." << endl;

    if(result > min_precision) {
        return result;
    }

    else return min_precision;

}

unsigned int compute_current_precision(unsigned int n, unsigned int N, unsigned int extra = 0) {
    // Roughly, we compute
    //
    //      log(A/sqrt(N) + B*sqrt(N/(n-1))*sinh(C * sqrt(n) / N) / log(2)
    //
    // where A, B, and C are the constants listed below. These error bounds
    // are given in the paper by Rademacher listed at the top of this file.
    // We then return this + extra, if extra != 0. If extra == 0, return with
    // what is probably way more extra precision than is needed.
    //
    // extra should probably have been set by a call to compute_extra_precision()
    // before this function was called.


    // n = the number for which we are computing p(n)
    // N = the number of terms that have been computed so far

    // if N is 0, then we can't use the above formula (because we would be
    // dividing by 0).
    if(N == 0) return compute_initial_precision(n) + extra;

    mpfr_t A, B, C;
    mpfr_init2(A, 32);
    mpfr_init2(B, 32);
    mpfr_init2(C, 32);

    mpfr_set_d(A,1.11431833485164,round_mode);
    mpfr_set_d(B,0.059238439175445,round_mode);
    mpfr_set_d(C,2.5650996603238,round_mode);

    mpfr_t error, t1, t2;
    mpfr_init2(error, 32);                                // we shouldn't need much precision here since we just need the most significant bit
    mpfr_init2(t1, 32);
    mpfr_init2(t2, 32);

    mpfr_set(error, A, round_mode);                       // error = A
    mpfr_sqrt_ui(t1, N, round_mode);                      // t1 = sqrt(N)
    mpfr_div(error, error, t1, round_mode);               // error = A/sqrt(N)


    mpfr_sqrt_ui(t1, n, round_mode);                      // t1 = sqrt(n)
    mpfr_mul(t1, t1, C, round_mode);                      // t1 = C * sqrt(n)
    mpfr_div_ui(t1, t1, N, round_mode);                   // t1 = C * sqrt(n) / N
    mpfr_sinh(t1, t1, round_mode);                        // t1 = sinh( ditto )
    mpfr_mul(t1, t1, B, round_mode);                      // t1 = B * sinh( ditto )

    mpfr_set_ui(t2, N, round_mode);                       // t2 = N
    mpfr_div_ui(t2, t2, n-1, round_mode);                 // t2 = N/(n-1)
    mpfr_sqrt(t2, t2, round_mode);                        // t2 = sqrt( ditto )

    mpfr_mul(t1, t1, t2, round_mode);                     // t1 = B * sqrt(N/(n-1)) * sinh(C * sqrt(n)/N)

    mpfr_add(error, error, t1, round_mode);               // error = (ERROR ESTIMATE)

    unsigned int p = mpfr_get_exp(error);                 // I am not 100% certain that this always does the right thing.
                                                          // (It should be the case that p is now the number of bits
                                                          // required to hold the integer part of the error.)


    if(extra == 0) {
        p = p + (unsigned int)ceil(log(n)/log(2));        // This is a stupid case to fall back on,
                                                          // left over from earlier versions of this code.
                                                          // It really should never be used.
    }
    else {
        p = p + extra;                                    // Recall that the extra precision should be
                                                          // large enough so that the accumulated errors
                                                          // in all of the computations that we make
                                                          // are not big enough to matter.
    }

    if(debug) {
        cout << "Error seems to be: ";
        mpfr_out_str(stdout, 10, 0, error, round_mode);
        cout << endl;
        cout.flush();
        cout << "Switching to precision of " << p << " bits. " << endl;
    }

    mpfr_clear(error);
    mpfr_clear(t1);
    mpfr_clear(t2);
    mpfr_clear(A);
    mpfr_clear(B);
    mpfr_clear(C);


    if(p > min_precision) {
        return p;                                         // We don't want to return < min_precision.
                                                          // Note that when the code that calls this
                                                          // function finds that the returned result
                                                          // is min_precision, it should stop calling
                                                          // this function, since the result won't change
                                                          // after that. Also, it should probably switch
                                                          // to computations with doubles, and should
                                                          // start calling compute_remainder().
    }
    return min_precision;
}

int compute_extra_precision(unsigned int n, double error = .25) {
    // Return the number of terms of the Rachemacher series that
    // we will need to compute to get a remainder of less than error
    // in absolute value, and then return the extra precision
    // that will guarantee that the accumulated error after computing
    // that number of steps will be less than .5 - error.
    //
    // How this works:
    // We first need to figure out how many terms of the series we are going to
    // need to compute. That is, we need to know how large k needs to be
    // for compute_remainder(n,k) to return something smaller than error. There
    // might be a clever way to do this, but instead we just keep calling
    // compute_current_precision() until we know that the error will be
    // small enough to call compute_remainder(). Then we just call compute_remainder()
    // until the error is small enough.
    //
    // Now that we know how many terms we will need to compute, k, we compute
    // the number of bits required to accurately store (.5 - error)/k. This ensures
    // that the total error introduced when we add up all k terms of the sum
    // will be less than (.5 - error). This way, if everything else works correctly,
    // then the sum will be within .5 of the correct (integer) answer, and we
    // can correctly find it by rounding.
    unsigned int k = 1;
    for( ; compute_current_precision(n,k,0) > double_precision; k += 100) {
    }

    for( ; compute_remainder(n, k) > error ; k += 100)  {
    }
    if(debug_precision) {
        cout << "To compute p(" << n << ") we will add up approximately " << k << " terms from the Rachemacher series." << endl;
    }
    int bits = (int)((log(k/(.5 - error)))/log(2)) + 5;   // NOTE: reducing the number of bits by 3 here is known to cause errors
                                                          // Why the extra 5 bits? Anytime we call a function, eg mp_a(),
                                                          // we end up doing a bunch of arithmetic operations, and if
                                                          // we want the result of those operations to be accurate
                                                          // within (.5 - error)/k, then we need that function to use
                                                          // a slightly higher working precision, which should be
                                                          // independent of n.
                                                          // TODO:
                                                          // Extensive trial and error has found 3 to be the smallest value
                                                          // that doesn't seem to produce any wrong answers. Thus, to
                                                          // be safe, we use 5 extra bits.
                                                          // (Extensive trial and error means compiling this file to get
                                                          // a.out and then running './a.out testforever' for a few hours.)




    return bits;
}

double compute_remainder(unsigned int n, unsigned int N) {
    // This computes the remainer left after N terms have been computed.
    // The formula is exactly the same as the one used to compute the required
    // precision, but once we know the necessary precision is small, we can
    // call this function to determine the actual error (rather than the precision).
    //
    // Generally, this is only called once we know that the necessary
    // precision is <= min_precision, because then the error is small
    // enough to fit into a double, and also, we know that we are
    // getting close to the correct answer.

    double A = 1.11431833485164;
    double B = 0.059238439175445;
    double C = 2.5650996603238;
    double result;
    result = A/sqrt(N) + B * sqrt(double(N)/double(n-1))*sinh(C * sqrt(double(n))/double(N));
    return result;
}


void initialize_mpfr_variables(unsigned int prec) {
    //
    // Clear and initialize some "temp" variables that are used in the computation of various functions.
    //
    // We do this in this auxilliary function (and endure the pain of
    // the extra global variables) so that we only have to initialize/clear
    // these variables once every time the precision changes.
    //
    // NAMING CONVENTIONS:
    //
    // -tempa1 and tempa2 are two variables available for use
    //      in the function mp_a()
    //
    // -tempf1 and tempf2 are two variables available for use
    //      in the function mp_f()
    //
    // -etc...
    //
    // NOTE: Calls to this function must be paired with calls to clear_mpfr_variables()

    mpfr_init2(tempa1, prec);
    mpfr_init2(tempa2, prec);
    mpfr_init2(tempf1, prec);
    mpfr_init2(tempf2, prec);
    mpfr_init2(temps1, prec);
    mpfr_init2(temps2, prec);
    mpfr_init2(tempc1, prec);
    mpfr_init2(tempc2, prec);
}

void clear_mpfr_variables() {
    mpfr_clear(tempa1);
    mpfr_clear(tempa2);
    mpfr_clear(tempf1);
    mpfr_clear(tempf2);
    mpfr_clear(temps1);
    mpfr_clear(temps2);

    mpfr_clear(tempc1);
    mpfr_clear(tempc2);
}


void initialize_constants(unsigned int prec, unsigned int n) {
    // The variables mp_A, mp_B, mp_C, and mp_D are used for
    // A_n, B_n, C_n, and D_n listed at the top of this file.
    //
    // They depend only on n, so we compute them just once in this function,
    // and then use them many times elsewhere.
    //
    // Also, we precompute some extra constants that we use a lot, such as
    // sqrt2, sqrt3, pi, 1/24, 1/12, etc.
    //
    // NOTE: Calls to this function must be paired with calls to clear_constants()
    static bool init = false;
    mpfr_prec_t p = prec;

    mpfr_init2(mp_one_over_12,p); mpfr_init2(mp_one_over_24,p); mpfr_init2(mp_sqrt2,p); mpfr_init2(mp_sqrt3,p); mpfr_init2(mp_pi,p);
    mpfr_init2(mp_A,p); mpfr_init2(mp_B,p); mpfr_init2(mp_C,p); mpfr_init2(mp_D,p); mpfr_init2(fourth, p); mpfr_init2(half, p);

    init = true;

    mpfr_set_ui(mp_one_over_12, 1, round_mode);                             // mp_one_over_12 = 1/12
    mpfr_div_ui(mp_one_over_12, mp_one_over_12, 12, round_mode);            //

    mpfr_set_ui(mp_one_over_24, 1, round_mode);                             // mp_one_over_24 = 1/24
    mpfr_div_ui(mp_one_over_24, mp_one_over_24, 24, round_mode);            //

    mpfr_set_ui(half, 1, round_mode);                                       //
    mpfr_div_ui(half, half, 2, round_mode);                                 // half = 1/2
    mpfr_div_ui(fourth, half, 2, round_mode);                               // fourth = 1/4

    mpfr_t n_minus;                                                         //
    mpfr_init2(n_minus, p);                                                 //
    mpfr_set_ui(n_minus, n, round_mode);                                    // n_minus = n
    mpfr_sub(n_minus, n_minus, mp_one_over_24,round_mode);                  // n_minus = n - 1/24

    mpfr_t sqrt_n_minus;                                                    //
    mpfr_init2(sqrt_n_minus, p);                                            //
    mpfr_sqrt(sqrt_n_minus, n_minus, round_mode);                           // n_minus = sqrt(n-1/24)


    mpfr_sqrt_ui(mp_sqrt2, 2, round_mode);                                  // mp_sqrt2 = sqrt(2)
    mpfr_sqrt_ui(mp_sqrt3, 3, round_mode);                                  // mp_sqrt3 = sqrt(3)
    mpfr_const_pi(mp_pi, round_mode);                                       // mp_pi = pi

    //mp_A = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0);---------------
    mpfr_set(mp_A, mp_sqrt2, round_mode);                                   // mp_A = sqrt(2)
    mpfr_mul(mp_A, mp_A, mp_pi, round_mode);                                // mp_A = sqrt(2) * pi
    mpfr_mul(mp_A, mp_A, sqrt_n_minus, round_mode);                         // mp_A = sqrt(2) * pi * sqrt(n - 1/24)
    //------------------------------------------------------------------------


    //mp_B = 2.0 * sqrt(3) * (n - 1.0/24.0);----------------------------------
    mpfr_set_ui(mp_B, 2, round_mode);                                       // mp_A = 2
    mpfr_mul(mp_B, mp_B, mp_sqrt3, round_mode);                             // mp_A = 2*sqrt(3)
    mpfr_mul(mp_B, mp_B, n_minus, round_mode);                              // mp_A = 2*sqrt(3)*(n-1/24)
    //------------------------------------------------------------------------


    //mp_C = sqrt(2) * pi * sqrt(n - 1.0/24.0) / sqrt(3);---------------------
    mpfr_set(mp_C, mp_sqrt2, round_mode);                                   // mp_C = sqrt(2)
    mpfr_mul(mp_C, mp_C, mp_pi, round_mode);                                // mp_C = sqrt(2) * pi
    mpfr_mul(mp_C, mp_C, sqrt_n_minus, round_mode);                         // mp_C = sqrt(2) * pi * sqrt(n - 1/24)
    mpfr_div(mp_C, mp_C, mp_sqrt3, round_mode);                             // mp_C = sqrt(2) * pi * sqrt(n - 1/24) / sqrt3
    //------------------------------------------------------------------------


    //mp_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);-----------------------
    mpfr_set_ui(mp_D, 2, round_mode);                                       // mp_D = 2
    mpfr_mul(mp_D, mp_D, n_minus, round_mode);                              // mp_D = 2 * (n - 1/24)
    mpfr_mul(mp_D, mp_D, sqrt_n_minus, round_mode);                         // mp_D = 2 * (n - 1/24) * sqrt(n - 1/24)
    //------------------------------------------------------------------------

    mpfr_clear(n_minus);
    mpfr_clear(sqrt_n_minus);


    // some double precision versions of the above values

    d_A = sqrt(2) * d_pi * sqrt(n - 1.0/24.0);
    d_B = 2.0 * sqrt(3) * (n - 1.0/24.0);
    d_C = sqrt(2) * d_pi * sqrt(n - 1.0/24.0) / sqrt(3);
    d_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);

    ld_A = sqrtl(2) * ld_pi * sqrtl(n - 1.0L/24.0L);
    ld_B = 2.0L * sqrtl(3) * (n - 1.0/24.0);
    ld_C = sqrtl(2) * ld_pi * sqrtl(n - 1.0L/24.0L) / sqrtl(3);
    ld_D = 2.0L * (n - 1.0L/24.0L) * sqrtl(n - 1.0L/24.0L);

}

void clear_constants() {
        mpfr_clear(mp_one_over_12); mpfr_clear(mp_one_over_24); mpfr_clear(mp_sqrt2); mpfr_clear(mp_sqrt3); mpfr_clear(mp_pi);
        mpfr_clear(mp_A); mpfr_clear(mp_B); mpfr_clear(mp_C); mpfr_clear(mp_D); mpfr_clear(half); mpfr_clear(fourth);
}

void initialize_mpz_and_mpq_variables() {
    /*
     * We use a few mpz_t and mpq_t variables which need to be initialized
     * before they can be used. Initialization and clearing take some
     * time, so we initialize just once in this function, and clear in another.
     */
    mpz_init(ztemp1);

    mpq_init(qtemps);
    mpq_init(qtempa);
    mpq_init(qtempa2);
}

void clear_mpz_and_mpq_variables() {
    mpz_clear(ztemp1);

    mpq_clear(qtemps);
    mpq_clear(qtempa);
    mpq_clear(qtempa2);
}

void mp_f(mpfr_t result, unsigned int k) {
    // compute f_n(k) as described in the introduction
    //
    // notice that this doesn't use n - the "constants"
    // A, B, C, and D depend on n, but they are precomputed
    // once before this function is called.

    //result =  pi * sqrt(2) * cosh(A/(sqrt(3)*k))/(B*k) - sinh(C/k)/D;

                                                                    //
    mpfr_set(result, mp_pi, round_mode);                            // result = pi

    mpfr_mul(result, result, mp_sqrt2, round_mode);                 // result = sqrt(2) * pi
                                                                    //
    mpfr_div(tempf1, mp_A, mp_sqrt3, round_mode);                   // temp1 = mp_A/sqrt(3)
    mpfr_div_ui(tempf1, tempf1, k, round_mode);                     // temp1 = mp_A/(sqrt(3) * k)
    mpfr_cosh(tempf1, tempf1, round_mode);                          // temp1 = cosh(mp_A/(sqrt(3) * k))
    mpfr_div(tempf1, tempf1, mp_B, round_mode);                     // temp1 = cosh(mp_A/(sqrt(3) * k))/mp_B
    mpfr_div_ui(tempf1, tempf1, k, round_mode);                     // temp1 = cosh(mp_A/(sqrt(3) * k))/(mp_B*k)
                                                                    //
    mpfr_mul(result, result, tempf1, round_mode);                   // result = sqrt(2) * pi * cosh(mp_A/(sqrt(3) * k))/(mp_B*k)
                                                                    //
    mpfr_div_ui(tempf1, mp_C, k, round_mode);                       // temp1 = mp_C/k
    mpfr_sinh(tempf1, tempf1, round_mode);                          // temp1 = sinh(mp_C/k)
    mpfr_div(tempf1, tempf1, mp_D, round_mode);                     // temp1 = sinh(mp_C/k)/D
                                                                    //
    mpfr_sub(result, result, tempf1, round_mode);                   // result = RESULT!
                                                                    //
    //return pi * sqrt(2) * cosh(A/(sqrt(3)*k))/(B*k) - sinh(C/k)/D;
}

// call the following when 4k < sqrt(UINT_MAX)
//
// TODO: maybe a faster version of this can be written without using mpq_t,
//       or maybe this version can be smarter.
//
//       The actual value of the cosine that we compute using s(h,k)
//       only depends on {s(h,k)/2}, that is, the fractional
//       part of s(h,k)/2. It may be possible to make use of this somehow.
//
void q_s(mpq_t result, unsigned int h, unsigned int k) {

    if(k < 3) {
        mpq_set_ui(result, 0, 1);
        return;
    }

    if (h == 1) {
        unsigned int d = GCD( (k-1)*(k-2), 12*k);
        if(d > 1) {
            mpq_set_ui(result, ((k-1)*(k-2))/d, (12*k)/d);
        }
        else {
            mpq_set_ui(result, (k-1)*(k-2), 12*k);
        }
        return;
    }
    // TODO:
    // It may be advantageous to implement some of the special forms in the comments below,
    // and also some more listed in some of the links mentioned in the introduction, but
    // it seems like there might not be much speed benefit to this, and putting in too
    // many seems to slow things down a little.
    //

    // if h = 2 and k is odd, we have
    // s(h,k) = (k-1)*(k-5)/24k
    //if(h == 2 && k > 5 && k % 2 == 1) {
    //    unsigned int d = GCD( (k-1)*(k-5), 24*k);
    //    if(d > 1) {
    //        mpq_set_ui(result, ((k-1)*(k-5))/d, (24*k)/d);
    //    }
    //    else {
    //        mpq_set_ui(result, (k-1)*(k-5), 24*k);
    //    }
    //    return;
    //}

/*

    // if k % h == 1, then
    //
    //      s(h,k) = (k-1)(k - h^2 - 1)/(12hk)
    //

    // We need to be a little careful here because k - h^2 - 1 can be negative.
    if(k % h == 1) {
        int num = (k-1)*(k - h*h - 1);
        int den = 12*k*h;
        int d = GCD(num, den);
        if(d > 1) {
            mpq_set_si(result, num/d, den/d);
        }
        else {
            mpq_set_si(result, num, den);
        }
        return;
    }

    // if k % h == 2, then
    //
    //      s(h,k) = (k-2)[k - .5(h^2 + 1)]/(12hk)
    //

    //if(k % h == 2) {
    //}
*/




    mpq_set_ui(result, 0, 1);                             // result = 0

    int r1 = k;
    int r2 = h;

    int n = 0;
    int temp3;
    while(r1 > 0 && r2 > 0) {
        unsigned int d = GCD(r1 * r1 + r2 * r2 + 1, r1 * r2);
        if(d > 1) {
            mpq_set_ui(qtemps, (r1 * r1 + r2 * r2 + 1)/d, (r1 * r2)/d);
        }
        else{
            mpq_set_ui(qtemps, r1 * r1 + r2 * r2 + 1, r1 * r2);
        }
        //mpq_canonicalize(qtemps);

        if(n % 2 == 0){                                             //
            mpq_add(result, result, qtemps);                        // result += temp1;
        }                                                           //
        else {                                                      //
            mpq_sub(result, result, qtemps);                        // result -= temp1;
        }                                                           //
        temp3 = r1 % r2;                                            //
        r1 = r2;                                                    //
        r2 = temp3;                                                 //
        n++;                                                        //
    }

    mpq_set_ui(qtemps, 1, 12);
    mpq_mul(result, result, qtemps);                                // result = result * 1.0/12.0;


    if(n % 2 == 1) {
        mpq_set_ui(qtemps, 1, 4);
        mpq_sub(result, result, qtemps);                            // result = result - .25;
    }

}


void mp_a(mpfr_t result, unsigned int n, unsigned int k) {
    // compute a(n,k)

    if (k == 1) {
        mpfr_set_ui(result, 1, round_mode);                         //result = 1
        return;
    }

    mpfr_set_ui(result, 0, round_mode);

    unsigned int h = 0;
    for(h = 1; h < k+1; h++) {
        if(GCD(h,k) == 1) {

            // Note that we compute each term of the summand as
            //      result += cos(pi * ( s(h,k) - (2.0 * h * n)/k) );
            //
            // This is the real part of the exponential that was written
            // down in the introduction, and we don't need to compute
            // the imaginary part because we know that, in the end, the
            // imaginary part will be 0, as we are computing an integer.

            q_s(qtempa, h, k);

            //mpfr_mul_q(tempa1, mp_pi, qtempa, round_mode);
            //mpfr_mul_ui(tempa1, tempa1, k * k, round_mode);

            //mpfr_set_q(tempa1, qtempa, round_mode);
            unsigned int r = n % k;                                     // here we make use of the fact that the
            unsigned int d = GCD(r,k);                                  // cos() term written above only depends
            unsigned int K;                                             // on {hn/k}.
            if(d > 1) {
                r = r/d;
                K = k/d;
            }
            else {
                K = k;
            }
            if(K % 2 == 0) {
                K = K/2;
            }
            else {
                r = r * 2;
            }
            mpq_set_ui(qtempa2, h*r, K);
            mpq_sub(qtempa, qtempa, qtempa2);

            //mpfr_set_q(tempa2, qtempa, round_mode);                   // This might be faster, according to
            //cospi(tempa1, tempa2);                                    // the comments in Ralf Stephan's part.c, but
                                                                        // I haven't noticed a significant speed up.
                                                                        // (Perhaps a different version that takes an mpq_t
                                                                        // as an input might be faster.)

            mpfr_mul_q(tempa1, mp_pi, qtempa, round_mode);
            mpfr_cos(tempa1, tempa1, round_mode);
            mpfr_add(result, result, tempa1, round_mode);

        }

    }

}

template <class T>
inline T partial_sum_of_t(unsigned int n, unsigned int &k, unsigned int exit_precision, unsigned int extra_precision, double error = 0) {
    unsigned int current_precision = compute_current_precision(n, k - 1, extra_precision);
    T result = 0;
    if(error == 0) {
        for(; current_precision > exit_precision; k++) {                        // (don't change k -- it is already the right value)
            result += sqrt(T(int(k))) * a<T>(n,k) * f<T>(k);                    //
            current_precision = compute_current_precision(n,k,extra_precision); // The only reason that we compute the new precision
                                                                                // now is so that we know when we can change to using just doubles.
                                                                                // (There should be a 'long double' version of the compute_current_precision function.
        }
    }
    else {
        double remainder = 1;
        for(; remainder > error; k++) {                                     // (don't change k -- it is already the right value)
            result += sqrt(T(int(k))) * a<T>(n,k) * f<T>(k);                //
            remainder = compute_remainder(n,k);
        }
    }
    return result;
}


void mp_t(mpfr_t result, unsigned int n) {
    // This is the function that actually computes p(n).
    //
    // More specifically, it computes t(n,N) to within 1/2 of p(n), and then
    // we can find p(n) by rounding.
    //
    //
    // NOTE: result should NOT have been initialized when this is called,
    // as we initialize it to the proper precision in this function.

    double error = .25;
    int extra = compute_extra_precision(n, error);

    unsigned int initial_precision = compute_initial_precision(n);  // We begin by computing the precision necessary to hold the final answer.
                                                                    // and then initialize both the result and some temporary variables to that
                                                                    // precision.
    mpfr_t t1, t2;                                                  //
    mpfr_init2(t1, initial_precision);                              //
    mpfr_init2(t2, initial_precision);                              //
    mpfr_init2(result, initial_precision);                          //
    mpfr_set_ui(result, 0, round_mode);                             //

    initialize_mpz_and_mpq_variables();
    initialize_constants(initial_precision, n);                     // Now that we have the precision information, we initialize some constants
                                                                    // that will be used throughout, and also
                                                                    //
    initialize_mpfr_variables(initial_precision);                   // set the precision of the "temp" variables that are used in individual functions.

    unsigned int current_precision = initial_precision;
    unsigned int new_precision;

    // We start by computing with high precision arithmetic, until
    // we are sure enough that we don't need that much precision
    // anymore. Once current_precision == min_precision, we drop
    // out of this loop and switch to a computation
    // that only involves doubles.



    unsigned int k = 1;                                             // (k holds the index of the summand that we are computing.)
    for(k = 1; current_precision > level_two_precision; k++) {      //

        mpfr_sqrt_ui(t1, k, round_mode);                            // t1 = sqrt(k)
                                                                    //
        mp_a(t2, n, k);                                             // t2 = A_k(n)

        if(debuga) {
            cout << "a(" << k <<  ") = ";
            mpfr_out_str(stdout, 10, 10, t2, round_mode);
            cout << endl;
        }

        mpfr_mul(t1, t1, t2, round_mode);                           // t1 = sqrt(k)*A_k(n)
                                                                    //
        mp_f(t2, k);                                                // t2 = f_k(n)

        if(debugf) {
            cout << "f(" << n << "," << k <<  ") = ";
            mpfr_out_str(stdout, 10, 0, t2, round_mode);
            cout << endl;
        }

        mpfr_mul(t1, t1, t2, round_mode);                           // t1 = sqrt(k)*A_k(n)*f_k(n)
                                                                    //
        mpfr_add(result, result, t1, round_mode);                   // result += summand

        if(debugt) {
            int num_digits = 20;
            int num_extra_digits;
            char digits[num_digits + 1];
            num_extra_digits = grab_last_digits(digits, 5, t1);
            grab_last_digits(digits, num_digits, result);

            mpfr_out_str(stdout, 10, 10, t1, round_mode);
            cout << endl;

            cout << k << ": current precision:"  << current_precision << ". 20 last digits of partial result: " << digits << ". Number of extra digits: " << num_extra_digits << endl;
            cout.flush();

        }

        new_precision = compute_current_precision(n,k,extra);       // After computing one summand, check what the new precision should be.
        if(new_precision != current_precision) {                    // If the precision changes, we need to clear
            current_precision = new_precision;                      // and reinitialize all "temp" variables to
            clear_mpfr_variables();                                 // use lower precision.
            initialize_mpfr_variables(current_precision);           //
            mpfr_clear(t1); mpfr_clear(t2);                         //
            mpfr_init2(t1, current_precision);                      //
            mpfr_init2(t2, current_precision);                      //
        }
    }

    mpfr_clear(t1); mpfr_clear(t2);

    mpfr_init2(t1, 200);
    mpfr_init2(t2, 200);

    long double ld_partial_sum = partial_sum_of_t<long double>(n,k,level_five_precision, extra, 0);
    mpfr_set_ld(t1, ld_partial_sum, round_mode);
    mpfr_add(result, result, t1, round_mode);

    double d_partial_sum = partial_sum_of_t<double>(n,k,0,extra,error);


    mpfr_set_d(t1, d_partial_sum, round_mode);                      //
    mpfr_add(result, result, t1, round_mode);                       // We add together the main result and the tail ends'


    mpfr_div(result, result, mp_pi, round_mode);                    // The actual result is the sum that we have computed
    mpfr_div(result, result, mp_sqrt2, round_mode);                 // divided by pi*sqrt(2).

    clear_constants();
    clear_mpz_and_mpq_variables();
    clear_mpfr_variables();
    mpfr_clear(t1);
    mpfr_clear(t2);
}

template <class T>
T f(unsigned int k) {
    return pi_sqrt2<T>() * cosh(A<T>()/(sqrt3<T>()*k))/(B<T>() * k) - sinh(C<T>()/k)/D<T>();
}

template <class T>
T a(unsigned int n, unsigned int k) {
    if(k == 1) {
        return 1;
    }
    T result = 0;

    unsigned int h = 0;
    for(h = 1; h < k+1; h++) {
        if(GCD(h,k) == 1) {
            result += cos( pi<T>() * ( s<T>(h,k) - T(2.0 * double(h) * n)/T(int(k))) ); // be careful to promote 2 and h to Ts
                                                                                        // because the result 2 * h * n could
                                                                                        // be too large.
        }
    }
    return result;
}

template <class T>
T s(unsigned int h, unsigned int k) {
    //mpq_set_ui(result, 1, 1);
    //return;

    //return T(0);  // Uncommenting this line saves 25% or more time in the computation of p(10^9) in the current version (.51)
                    // but, of course, gives wrong results. (This is just to give an idea of how much time could
                    // possibly be saved by optimizing this function.

    if(k < 3) {
        return T(0);
    }

    if (h == 1) {
        return T( int((k-1)*(k-2)) )/T(int(12 * k));
    }
    // TODO: In the function mp_s() there are a few special cases for special forms of h and k.
    // (And there are more special cases listed in one of the references listed in the introduction.)
    //
    // It may be advantageous to implement some here, but I'm not sure
    // if there is any real speed benefit to this.
    //
    // In the mpfr_t version of this function, the speedups didn't seem to help too much, but
    // they might make more of a difference here.
    //
    // Update to the above comments:
    // Actually, a few tests seem to indicate that putting in too many special
    // cases slows things down a little bit.

    // if h = 2 and k is odd, we have
    // s(h,k) = (k-1)*(k-5)/24k
    if(h == 2 && k > 5 && k % 2 == 1) {
        return T( int((k-1)*(k-5)) )/ T(int(24*k));
    }


    // if k % h == 1, then
    //
    //      s(h,k) = (k-1)(k - h^2 - 1)/(12hk)
    //

    // We might need to be a little careful here because k - h^2 - 1 can be negative.
    // THIS CODE DOESN'T WORK.
    //if(k % h == 1) {
    //    int num = (k-1)*(k - h*h - 1);
    //    int den = 12*k*h;
    //    return T(num)/T(den);
    //}

    // if k % h == 2, then
    //
    //      s(h,k) = (k-2)[k - .5(h^2 + 1)]/(12hk)
    //

    //if(k % h == 2) {
    //}

    int r1 = k;
    int r2 = h;

    int n = 1;
    int temp3;

    T result = T(0);
    T temp;
    while(r2 > 0)   // Note that we maintain the invariant r1 >= r2, so
                    // we only need to make sure that r2 > 0
    {
        temp = T(int(r1 * r1 + r2 * r2 + 1))/(n * r1 * r2);
        temp3 = r1 % r2;
        r1 = r2;
        r2 = temp3;

        result += temp;                        // result += temp1;

        n = -n;
    }

    result *= one_over_12<T>();

    if(n < 0) {
        result -= T(.25);
    }
    return result;
}


long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0) {
      a = -a;
   }

   if (b < 0) {
      b = -b;
   }


   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v;
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}



// The following function was copied from Ralf Stephan's code,
// mentioned in the introduction, and then some variable names
// were changed. The function is not currently used, however.
void cospi (mpfr_t res,
        mpfr_t x)
{
//    mpfr_t t, tt, half, fourth;

//    mpfr_init2 (t, prec);
//    mpfr_init2 (tt, prec);
//    mpfr_init2 (half, prec);
//    mpfr_init2 (fourth, prec);

//    mpfr_set_ui (half, 1, r);
//    mpfr_div_2ui (half, half, 1, r);
//    mpfr_div_2ui (fourth, half, 1, r);

        // NOTE: switched t to tempc2
        //       and tt to tempc1


        mpfr_rnd_t r = round_mode;


    mpfr_div_2ui (tempc1, x, 1, r);
    mpfr_floor (tempc1, tempc1);
    mpfr_mul_2ui (tempc1, tempc1, 1, r);
    mpfr_sub (tempc2, x, tempc1, r);
    if (mpfr_cmp_ui (tempc2, 1) > 0)
        mpfr_sub_ui (tempc2, tempc2, 2, r);
    mpfr_abs (tempc1, tempc2, r);
    if (mpfr_cmp (tempc1, half) > 0)
    {
        mpfr_ui_sub (tempc2, 1, tempc1, r);
        mpfr_abs (tempc1, tempc2, r);
        if (mpfr_cmp (tempc1, fourth) > 0)
        {
            if (mpfr_sgn (tempc2) > 0)
                mpfr_sub (tempc2, half, tempc2, r);
            else
                mpfr_add (tempc2, tempc2, half, r);
            mpfr_mul (tempc2, tempc2, mp_pi, r);
            mpfr_sin (tempc2, tempc2, r);
        }
        else
        {
            mpfr_mul (tempc2, tempc2, mp_pi, r);
            mpfr_cos (tempc2, tempc2, r);
        }
        mpfr_neg (res, tempc2, r);
    }
    else
    {
        mpfr_abs (tempc1, tempc2, r);
        if (mpfr_cmp (tempc1, fourth) > 0)
        {
            if (mpfr_sgn (tempc2) > 0)
                mpfr_sub (tempc2, half, tempc2, r);
            else
                mpfr_add (tempc2, tempc2, half, r);
            mpfr_mul (tempc2, tempc2, mp_pi, r);
            mpfr_sin (res, tempc2, r);
        }
        else
        {
            mpfr_mul (tempc2, tempc2, mp_pi, r);
            mpfr_cos (res, tempc2, r);
        }
    }

//    mpfr_clear (half);
//    mpfr_clear (fourth);
//    mpfr_clear (t);
//    mpfr_clear (tt);
}

/* answer must have already been mpz_init'd. */
int part(mpz_t answer, unsigned int n){
    if(n == 1) {
        mpz_set_str(answer, "1", 10);
        return 0;
    }
    mpfr_t result;

    mp_t(result, n);

    mpfr_get_z(answer, result, round_mode);

    mpfr_clear(result);
    return 0;
}


int main(int argc, char *argv[]){
    //init();

    unsigned int n = 10;

    if(argc > 1)
        if(strcmp(argv[1], "test") == 0) {
            n = test(true);
            if(n == 0) {
                cout << "All Tests Passed" << endl;
            }
            else {
                cout << "Error computing p(" << n << ")" << endl;
            }
            return 0;
        }
        else if(strcmp(argv[1], "testforever") == 0) {
            n = test(false, true);
            if(n == 0) {
                cout << "All Tests Passed" << endl;
            }
            else {
                cout << "Error computing p(" << n << ")" << endl;
            }
            return 0;
        }
        else {
            n = atoi(argv[1]);
        }
    else {
        n = test(false);
        if(n == 0) {
            cout << "All short tests passed. Run '" << argv[0] << " test' to run all tests. (This may take some time, but it gives updates as it progresses, and can be interrupted.)" << endl;
            cout << "Run with the argument 'testforever' to run tests until a failure is found (or, hopefully, to run tests forever.)" << endl;
        }
        else {
            cout << "Error computing p(" << n << ")" << endl;
        }
        return 0;
    }
    //mpfr_t result;

    //mp_t(result, n);

    mpz_t answer;
    mpz_init(answer);
    part(answer, n);

    //mpfr_get_z(answer, result, round_mode);

    mpz_out_str (stdout, 10, answer);

    cout << endl;

    return 0;
}


int test(bool longtest, bool forever) {
    // The exact values given below are confirmed by multiple sources, so are probably correct.
    // Other tests rely on known congruences.
    // If longtest is true, then we run some tests that probably take on the order
    // of 10 minutes.
    // If forever is true, then we just do randomized tests until a failure
    // is found. Ideally, this should mean that we test forever, of course.

    int n;

    mpz_t expected_value;
    mpz_t actual_value;

    mpz_init(expected_value);
    mpz_init(actual_value);

    n= 1;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "1", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;

    n = 10;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "42", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;

    n = 100;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "190569292", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;

    n = 1000;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "24061467864032622473692149727991", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;

    n = 10000;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "36167251325636293988820471890953695495016030339315650422081868605887952568754066420592310556052906916435144", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;

    n = 100000;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;

    n = 1000000;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    mpz_set_str(expected_value, "1471684986358223398631004760609895943484030484439142125334612747351666117418918618276330148873983597555842015374130600288095929387347128232270327849578001932784396072064228659048713020170971840761025676479860846908142829356706929785991290519899445490672219997823452874982974022288229850136767566294781887494687879003824699988197729200632068668735996662273816798266213482417208446631027428001918132198177180646511234542595026728424452592296781193448139994664730105742564359154794989181485285351370551399476719981691459022015599101959601417474075715430750022184895815209339012481734469448319323280150665384042994054179587751761294916248142479998802936507195257074485047571662771763903391442495113823298195263008336489826045837712202455304996382144601028531832004519046591968302787537418118486000612016852593542741980215046267245473237321845833427512524227465399130174076941280847400831542217999286071108336303316298289102444649696805395416791875480010852636774022023128467646919775022348562520747741843343657801534130704761975530375169707999287040285677841619347472368171772154046664303121315630003467104673818", 10);
    part(actual_value, n);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return n;

    cout << " OK." << endl;


    // We now run some tests based on the fact that if n = 369 (mod 385) then p(n) = 0 (mod 385).

    n = 369 + 10*385;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    part(actual_value, n);
    if(mpz_divisible_ui_p(actual_value, 385) == 0)
        return n;

    cout << " OK." << endl;

    n = 369 + 10*385;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    part(actual_value, n);
    if(mpz_divisible_ui_p(actual_value, 385) == 0)
        return n;

    cout << " OK." << endl;

    n = 369 + 100*385;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    part(actual_value, n);
    if(mpz_divisible_ui_p(actual_value, 385) == 0)
        return n;

    cout << " OK." << endl;

    n = 369 + 110*385;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    part(actual_value, n);
    if(mpz_divisible_ui_p(actual_value, 385) == 0)
        return n;

    cout << " OK." << endl;

    n = 369 + 120*385;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    part(actual_value, n);
    if(mpz_divisible_ui_p(actual_value, 385) == 0)
        return n;

    cout << " OK." << endl;

    n = 369 + 130*385;
    cout << "Computing p(" << n << ")...";
    cout.flush();
    part(actual_value, n);
    if(mpz_divisible_ui_p(actual_value, 385) == 0)
        return n;

    cout << " OK." << endl;

    // Randomized testing

    srand( time(NULL) );

    for(int i = 0; i < 100; i++) {
        n = int(100000 * double(rand())/double(RAND_MAX) + 1);
        n = n - (n % 385) + 369;
        cout << "Computing p(" << n << ")...";
        cout.flush();
        part(actual_value, n);
        if(mpz_divisible_ui_p(actual_value, 385) == 0) {
            return n;
        }
        cout << " OK." << endl;

    }

    if(longtest) {
        n = 369 + 1000*385;
        cout << "Computing p(" << n << ")...";
        cout.flush();
        part(actual_value, n);
        if(mpz_divisible_ui_p(actual_value, 385) == 0)
            return n;

        cout << " OK." << endl;

        n = 369 + 10000*385;
        cout << "Computing p(" << n << ")...";
        cout.flush();
        part(actual_value, n);
        if(mpz_divisible_ui_p(actual_value, 385) == 0)
            return n;

        cout << " OK." << endl;

        n = 369 + 100000*385;
        cout << "Computing p(" << n << ")...";
        cout.flush();
        part(actual_value, n);
        if(mpz_divisible_ui_p(actual_value, 385) == 0)
            return n;

        cout << " OK." << endl;

        for(int i = 0; i < 20; i++) {
            n = int(100000 * double(rand())/double(RAND_MAX) + 1) + 100000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        for(int i = 0; i < 20; i++) {
            n = int(100000 * double(rand())/double(RAND_MAX) + 1) + 500000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        for(int i = 0; i < 20; i++) {
            n = int(100000 * double(rand())/double(RAND_MAX) + 1) + 1000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        for(int i = 0; i < 10; i++) {
            n = int(100000 * double(rand())/double(RAND_MAX) + 1) + 10000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        n = 369 + 1000000*385;
        cout << "Computing p(" << n << ")...";
        cout.flush();
        part(actual_value, n);
        if(mpz_divisible_ui_p(actual_value, 385) == 0)
            return n;

        cout << " OK." << endl;

        for(int i = 0; i < 10; i++) {
            n = int(100000000 * double(rand())/double(RAND_MAX) + 1) + 100000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        n = 1000000000 + 139;
        cout << "Computing p(" << n << ")...";
        cout.flush();
        part(actual_value, n);
        if(mpz_divisible_ui_p(actual_value, 385) == 0)
            return n;

        cout << " OK." << endl;


        for(int i = 0; i < 10; i++) {
            n = int(100000000 * double(rand())/double(RAND_MAX) + 1) + 1000000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

    }

    if(forever) {
        while(1) {

        for(int i = 0; i < 100; i++) {
            n = int(900000 * double(rand())/double(RAND_MAX) + 1) + 100000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        for(int i = 0; i < 50; i++) {
            n = int(9000000 * double(rand())/double(RAND_MAX) + 1) + 1000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        for(int i = 0; i < 50; i++) {
            n = int(90000000 * double(rand())/double(RAND_MAX) + 1) + 10000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }

        for(int i = 0; i < 10; i++) {
            n = int(900000000 * double(rand())/double(RAND_MAX) + 1) + 100000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }
        for(int i = 0; i < 5; i++) {
            n = int(100000000 * double(rand())/double(RAND_MAX) + 1) + 1000000000;
            n = n - (n % 385) + 369;
            cout << "Computing p(" << n << ")...";
            cout.flush();
            part(actual_value, n);
            if(mpz_divisible_ui_p(actual_value, 385) == 0) {
                return n;
            }
            cout << " OK." << endl;
        }
        }


    }

    mpz_clear(expected_value);
    mpz_clear(actual_value);
    return 0;
}


