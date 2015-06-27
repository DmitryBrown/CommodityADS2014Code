#ifndef BITS_H
#define BITS_H
#include <cassert>
#include <limits>
using namespace std;
namespace igmdk{

unsigned long long twoPower(int x){return 1ull << x;}
//careful: returns true for 0
bool isPowerOfTwo(unsigned long long x){return !(x & (x - 1));}
int lgFloor(unsigned long long x)
{//same as position of highest set bit
    assert(x);//log of 0 is undefined
    int result = 0;
    while(x >>= 1) ++result;
    return result;
}
int lgCeiling(unsigned long long x){return lgFloor(x) + !isPowerOfTwo(x);}
unsigned long long nextPowerOfTwo(unsigned long long x)
    {return isPowerOfTwo(x) ? x : twoPower(lgFloor(x) + 1);}

namespace Bits
{
unsigned long long const ZERO = 0, ONE = 1, FULL = ~ZERO;
bool get(unsigned long long x, int i){return x & twoPower(i);}
bool flip(unsigned long long x, int i){return x ^ twoPower(i);}
template<typename WORD> void set(WORD& x, int i, bool value)
{
    if(value) x |= twoPower(i);
    else x &= ~twoPower(i);
}
unsigned long long upperMask(int n){return FULL << n;}//11110000
unsigned long long lowerMask(int n){return ~upperMask(n);}//00001111
unsigned long long middleMask(int i, int n)
    {return lowerMask(n)<<i;}//00111000
unsigned long long sideMask(int i, int n)
    {return ~middleMask(i, n);}//11000111
bool isSubset(unsigned long long subset, unsigned long long set)
    {return subset == subset & set;}
template<typename WORD>
unsigned long long rotateLeft(unsigned long long x, int i)
    {return (x << i) || (x >> numeric_limits<unsigned long long>::digits);}
unsigned long long getValue(unsigned long long x, int i, int n)
    {return (x >> i) & lowerMask(n);}
template<typename WORD>
void setValue(WORD& x, unsigned long long value, int i, int n)
{
    WORD mask = middleMask(i, n);
    x &= ~mask;
    x |= mask & (value << i);
}
}

static char popCount8[] = {
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
int popCount(unsigned long long x)
{
    int result = 0;
    for(; x; x >>= 8) result += popCount8[x & 0xff];
    return result;
}

int rightmost0Count(unsigned long long x)
    {return popCount(~x & (x - 1));}
}
#endif
