#include <cstdint>
#include <x86intrin.h>
#include <array>
#include "config.h"

#ifdef GCC_VECT_CMP
#warning "Using GCC vectors extension syntax"
#else
#warning "Using GCC intrinsics syntax"
#endif

template < unsigned _N, typename Vals=uint8_t >
struct SymGroup {
  static const constexpr unsigned N = _N;
  static const constexpr unsigned vect_len = _N;
  using type = std::array<Vals, N>;
  static type one() { type res; for (unsigned i=0; i<N; i++) res[i]=i; return res; };
  static void mult(type &res, const type &a, const type &b) {
    for (unsigned i=0; i < N; i++) res[i] = a[b[i]];
  };

  static int less_partial(const type &a, const type &b, int k) {
    for (unsigned i=0; i<k; i++) {
      if (a[i] != b[i]) return a[i] > b[i] ? 1 : -1;
    }
    return 0;
  }
};



using vect16 = std::array<uint8_t, 16> ;
using epi8 = uint8_t __attribute__ ((vector_size (16)));

struct VectPerm
{
  union {
    vect16 p;
    __m128i v;
    epi8 v8;
  };
  #ifdef GCC_VECT_CMP
  bool operator==(const VectPerm &vp) const { return _mm_movemask_epi8(v==vp.v)==0xffff; }
  #else
  bool operator==(const VectPerm &vp) const { return _mm_movemask_epi8(_mm_cmpeq_epi8(v,vp.v))==0xffff; }
  #endif
  uint8_t operator[](unsigned long i) const { return p[i]; }
  uint8_t &operator[](unsigned long i) { return p[i]; }
  bool operator<(const VectPerm &b) const {
    const char mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
    uint64_t diff = _mm_cmpestri (v, 16, b.v, 16, mode);
    #ifdef GCC_VECT_CMP
    int diffs = _mm_movemask_epi8(v < b.v);
    #else
    int diffs = _mm_movemask_epi8(_mm_cmplt_epi8(v,b.v));
    #endif
    return (diffs >> diff) & 0x1;;
  }
  bool first_diff(const VectPerm &b, int k) const {
    const char mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
    return _mm_cmpestri (v, k, b.v, k, mode);
  }

  bool less_partial(const VectPerm &b, int diff) const {
    #ifdef GCC_VECT_CMP
    int diffs = _mm_movemask_epi8(v8 < b.v8);
    #else
    int diffs = _mm_movemask_epi8(_mm_cmplt_epi8(__m128i(v8),__m128i(b.v8)));
    #endif
    return (diffs >> diff) & 0x1;
  }
} __attribute__ ((aligned (16)));


template < unsigned _N >
struct SymGroupVect {
  static_assert ( _N <= 16, "Size must be <= 16 because AVX instruction");
  static const constexpr unsigned N = _N;
  static const constexpr unsigned vect_len = 16;
  using type = VectPerm;
    #ifdef GCC_VECT_CMP
    static type one() { return {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; }
    #else
    static type one() { type res; for (int i=0; i<16; i++) res.p[i]=i; return res; }
    #endif
  static void mult(type &res, const type &a, const type &b) {
    res.v = _mm_shuffle_epi8(a.v, b.v);
  };
};


namespace std {

  template <class C, size_t s>
  struct hash< std::array<C, s> > {
    size_t operator () (const std::array<C, s> &ar) const {
      size_t h = 0;
      for (auto e : ar) h = hash<C>()(e) + (h << 6) + (h << 16) - h;
      return h;
    }
  };

  template<>
  struct hash<VectPerm> {
    size_t operator () (const VectPerm &ar) const {
    #ifdef GCC_VECT_CMP
    return ar.v[1] ^ ar.v[0]; // hash<vect16>()(ar.p);
    #else
    return hash<vect16>()(ar.p);
    #endif
    }
  };

  template<>
  struct less<vect16> {
    size_t operator () (const vect16 &v1, const vect16 &v2) const {
      return v1[0] == v2[0] ? v1[1] < v2[1] : v1[0] < v2[0];
    }
  };
}

#include <vector>
#include <set>
#include <unordered_set>

const unsigned N = 8;
// using SGroup = SymGroup<N>;
using SGroup = SymGroupVect<N>;

bool is_canonical(const std::vector< std::vector<SGroup::type> > & sgs,
		  const SGroup::type &v) {
  std::set<SGroup::type> to_analyse({v}), new_to_analyse;
  // std::unordered_set<SGroup::type> to_analyse({v}), new_to_analyse;
  SGroup::type child;

  for (unsigned int i=0; i < SGroup::N-1; i++) {
    new_to_analyse.clear();
    auto &transversal = sgs[i];
    for (SGroup::type list_test : to_analyse) {
      for (SGroup::type x : transversal) {
	SGroup::mult(child, list_test, x);
	int diff = v.first_diff(child, i+1);
	if (diff == 16) new_to_analyse.insert(child);
	if (v.less_partial(child, diff)) return false;
      }
    }
    to_analyse.swap(new_to_analyse);
  }
  return true;
}

