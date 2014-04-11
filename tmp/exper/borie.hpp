#include <cstdint>
#include <x86intrin.h>
#include <iostream>

#include <cilk/cilk_api.h>

#include <array>
#include <vector>
#include <set>
#include <unordered_set>
#include <utility>

#include "config.h"


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
  int less_partial(const VectPerm &b, int k) const {
    const char mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
    int diff = _mm_cmpestri (v, k, b.v, k, mode);
    return (diff == 16) ? 0 : char(p[diff]) - char(b.p[diff]);
  }
} __attribute__ ((aligned (16)));

static std::ostream & operator<<(std::ostream & stream, const VectPerm &term) {
  stream << "[" << unsigned(term[0]);
  for (unsigned i=1; i < 16; i++) stream << "," << unsigned(term[i]);
  stream << "]";
  return stream;
};


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
      for (int i=0; i<s; i++) h = hash<C>()(ar[i]) + (h << 6) + (h << 16) - h;
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


const unsigned N = 16;

// using SGroup = SymGroup<N>;
using SGroup = SymGroupVect<N>;
using StrongGeneratingSet = const std::vector< std::vector<SGroup::type> >;

#ifdef USE_TBB
#include "tbb/scalable_allocator.h"
using allocator = tbb::scalable_allocator< SGroup::type >;
using set = std::set<SGroup::type, std::less<SGroup::type>, allocator>;
#else
using set = std::set<SGroup::type>;
#endif

bool is_canonical(const StrongGeneratingSet & sgs, const SGroup::type &v) {
  set to_analyse({v}), new_to_analyse;
  SGroup::type child;
  for (unsigned int i=0; i < SGroup::N-1; i++) {
    new_to_analyse.clear();
    auto &transversal = sgs[i];
    for (SGroup::type list_test : to_analyse) {
      for (SGroup::type x : transversal) {
	SGroup::mult(child, list_test, x);
	int comp = v.less_partial(child, i+1);
	if (comp == 0) new_to_analyse.insert(child);
	else if (comp < 0) return false;
      }
    }
    to_analyse.swap(new_to_analyse);
  }
  return true;
}


const char LastNonZeroIndex = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_ANY |
  _SIDD_MASKED_NEGATIVE_POLARITY | _SIDD_MOST_SIGNIFICANT;

__m128i zero {0,0};

class children_iterator {

  const SGroup::type &v;
  uint64_t ind;
public:
  children_iterator(const SGroup::type &_v) : v(_v) {
    ind = _mm_cmpestri (zero, 1, _v.v, SGroup::N, LastNonZeroIndex);
    if (ind == 16) ind = 0;
  };

  void operator++() {  ind++;}
  bool is_not_end() { return ind < SGroup::N; }
  SGroup::type operator *() {
    SGroup::type res = v;
    res.p[ind]++;
    return res;
  }
};


// #include <list>
// using list_vect = std::list<SGroup::type>;
// void walk_tree(const SGroup::type &v, list_vect &res, int remaining_depth,
// 	       StrongGeneratingSet & sgs) {
//   if (remaining_depth == 0) res.push_back(v);
//   else for (children_iterator ch(v); ch.is_not_end(); ++ch) {
//       SGroup::type child = *ch;
//       if (is_canonical(sgs, child)) walk_tree(child, res, remaining_depth-1, sgs);
//     }
// }
//
// std::list<SGroup::type> elements_of_depth(int depth, StrongGeneratingSet & sgs) {
//   SGroup::type zero_vect;
//   list_vect list_res;
//   zero_vect.v = zero;
//   walk_tree(zero_vect, list_res, depth, sgs);
//   return list_res;
// }

bool cilk_start(char *nproc) {
  return (__cilkrts_set_param("nworkers", nproc) == __CILKRTS_SET_PARAM_SUCCESS);
}

#include <cilk/reducer_list.h>
using list_vect = cilk::reducer_list_append<SGroup::type>;

void walk_tree(const SGroup::type v, list_vect &res, int remaining_depth,
	       StrongGeneratingSet & sgs) {
  if (remaining_depth == 0) res.push_back(v);
  else for (children_iterator ch(v); ch.is_not_end(); ++ch) {
      SGroup::type child = *ch;
      if (is_canonical(sgs, child))
	cilk_spawn walk_tree(child, res, remaining_depth-1, sgs);
    }
}

std::list<SGroup::type> elements_of_depth(int depth, StrongGeneratingSet & sgs) {
  SGroup::type zero_vect;
  list_vect list_res;
  zero_vect.v = zero;
  walk_tree(zero_vect, list_res, depth, sgs);
  return list_res.get_value();
}

