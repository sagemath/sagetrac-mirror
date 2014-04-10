#include <cstdint>
#include <x86intrin.h>
#include <array>
#include <iostream>
#include <bitset>

using vect16 = std::array<uint8_t, 16> ;

using epi8 = uint8_t __attribute__ ((vector_size (16)));

struct VectPerm
{
  union {
    vect16 p;
    __m128i v;
    epi8 v8;
  };
  bool operator==(const VectPerm &vp) const { return _mm_movemask_epi8(v==vp.v)==0xffff; }
  uint8_t operator[](unsigned long i) const { return p[i]; }
  uint8_t &operator[](unsigned long i) { return p[i]; }
  bool operator<(const VectPerm &b) const { // todo : vectorize
    const char mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
    uint64_t diff = _mm_cmpestri (v, 16, b.v, 16, mode);
    int diffs = _mm_movemask_epi8(v < b.v);
    return (diffs >> diff) & 0x1;;
  }

  int less_partial(const VectPerm &b, int k) {  // vectorize
    const char mode = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
    uint64_t diff = _mm_cmpestri (v, k, b.v, k, mode);
    int diffs = _mm_movemask_epi8(v8 < b.v8);
    std::cout << "diff " << diff << "   diffs " << diffs << ", "
	      << std::bitset<16>(diffs) << std::endl;
    return diff == 16 ? 0 : (1 - ((diffs >> (diff-1)) & 0x2));
  }
} __attribute__ ((aligned (16)));




int main() {
  VectPerm x {0,1,1,3,3,2,0,0 ,0,0,0,0,0,0,0,0};
  VectPerm y {0,3,1,3,5,2,0,0 ,0,0,0,0,0,0,0,1};
  VectPerm z {0,1,0,3,5,2,0,0 ,0,0,0,0,0,0,0,1};
  std::cout << x.less_partial(y, 1) << std::endl; // 0
  std::cout << x.less_partial(y, 2) << std::endl; // -1
  std::cout << y.less_partial(x, 2) << std::endl; // 1
  std::cout << x.less_partial(z, 2) << std::endl; // 0
  std::cout << x.less_partial(z, 4) << std::endl; // 1
}
