/*
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

/*
   C++ codes by Kohei Takeda (k-tak@letter.or.jp)
   Redistribution terms are the same as the original ones above.
*/

#ifndef ___MERSENNE_TWISTER_RNG___
#define ___MERSENNE_TWISTER_RNG___

#include <cassert>
#include <cstdlib>
#include <ctime>

struct Mt64Traits
{
  typedef unsigned long long UINTTYPE;

  static const int INTTYPE_BITS = 64;
  static const unsigned long long MAXDOUBLEVAL = 9007199254740991ULL; // 2^53-1
  static const size_t NN = 312;
  static const size_t MM = 156;
  static const unsigned long long INITVAL = 6364136223846793005ULL;
  static const unsigned long long ARRAYINITVAL_0 = 19650218ULL;
  static const unsigned long long ARRAYINITVAL_1 = 3935559000370003845ULL;
  static const unsigned long long ARRAYINITVAL_2 = 2862933555777941757ULL;

  static unsigned long long twist(const unsigned long long& u, const unsigned long long& v)
  {
    static unsigned long long mag01[2] = { 0ULL, 0xB5026F5AA96619E9ULL };
    return ((((u & 0xFFFFFFFF80000000ULL) | (v & 0x7FFFFFFFULL)) >> 1) ^ mag01[v & 1]);
  }

  static unsigned long long temper(unsigned long long y)
  {
    y ^= (y >> 29) & 0x5555555555555555ULL;
    y ^= (y << 17) & 0x71D67FFFEDA60000ULL;
    y ^= (y << 37) & 0xFFF7EEE000000000ULL;
    y ^= (y >> 43);

    return y;
  }
};

class MtRng64
{
public:
  typedef typename Mt64Traits::UINTTYPE UINTTYPE;

  // member variables
  UINTTYPE state_[Mt64Traits::NN];
  size_t left_;
  UINTTYPE* next_;

protected:
  void nextState()
  {
    UINTTYPE* p = state_;
    size_t j;

    left_ = Mt64Traits::NN - 1;
    next_ = state_;

    for (j = Mt64Traits::NN - Mt64Traits::MM + 1; --j; p++)
      *p = p[Mt64Traits::MM] ^ Mt64Traits::twist(p[0], p[1]);

    for (j = Mt64Traits::MM; --j; p++)
      *p = p[Mt64Traits::MM - Mt64Traits::NN] ^ Mt64Traits::twist(p[0], p[1]);

    *p = p[Mt64Traits::MM - Mt64Traits::NN] ^ Mt64Traits::twist(p[0], state_[0]);
  }

public:
  MtRng64() { init((UINTTYPE)time(NULL)); }

  MtRng64(UINTTYPE seed) { init(seed); }

  void init(UINTTYPE seed)
  {
    assert(sizeof(UINTTYPE) * 8 == (size_t)Mt64Traits::INTTYPE_BITS);

    state_[0] = seed;

    for (size_t j = 1; j < Mt64Traits::NN; j++) {
      state_[j] =
        (Mt64Traits::INITVAL * (state_[j - 1] ^ (state_[j - 1] >> (Mt64Traits::INTTYPE_BITS - 2))) + (UINTTYPE)j);
    }
    left_ = 1;

    next_ = state_;
  }

  /* generates a random number on [0,1)-real-interval */
  double getReal2()
  {
    if (left_-- == 0)
      nextState();
    if (Mt64Traits::INTTYPE_BITS > 53) {
      return ((double)(Mt64Traits::temper(*next_++) >> (Mt64Traits::INTTYPE_BITS - 53)) * (1.0 / 9007199254740992.0));
    } else {
      return ((double)Mt64Traits::temper(*next_++) * (1.0 / ((double)Mt64Traits::MAXDOUBLEVAL + 1.0)));
    }
  }

  // We don't currently need the following constructors, but perhaps we do
  // in the future, and the shallow copying of the state array may create problems.

  MtRng64(const MtRng64& obj)
  {
    for (int i = 0; i < 312; i++) {
      state_[i] = obj.state_[i];
    }

    left_ = obj.left_;
    next_ = state_ + (obj.next_ - obj.state_);
  }

  MtRng64& operator=(const MtRng64& obj)
  {
    for (int i = 0; i < 312; i++) {
      state_[i] = obj.state_[i];
    }

    left_ = obj.left_;
    next_ = state_ + (obj.next_ - obj.state_);
    return *this;
  }
};

#endif
