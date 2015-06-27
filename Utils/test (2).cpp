#include "Bitset.h"
using namespace igmdk;

int main()
{
    for(int i = 0; i < 5; ++i) DEBUG(nextPowerOfTwo(i));

	Bitset<unsigned char> bs(16);
	bs.output();
	for(int i = 0; i < 16; i+=3)
	{
		bs.set(i, true);
		bs.output();
	}
	bs.set(3, false);
	bs.output();
	bs <<= 9;
	bs.output();
	bs >>= 9;
	bs.output();
	bs.flip();
	bs.output();

    NBitVector<5> x;
    int N = 100000000;
    for(int i = 0; i < N; ++i) x.append(i);
    for(int i = 0; i < N; ++i) x.set(i, i);
    for(int i = 0; i < N; ++i) assert(x[i] == i % 32);
	return 0;
}
