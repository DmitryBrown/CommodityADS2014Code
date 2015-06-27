#ifndef MISC_H
#define MISC_H
#include "../Utils/Bitset.h"
#include "../Utils/Sort.h"
#include "../Heaps/Heap.h"
#include <cmath>

namespace igmdk{

class CRC32
{
    unsigned int polynomial, constant[256];
public:
    CRC32(unsigned int thePolynomial = 0xFA567D89u):
        polynomial(thePolynomial)
    {
        for(int i = 0; i < 256; ++i)
        {
            constant[i] = i << 24;
            for(int j = 0; j < 8; ++j) constant[i] =
                (constant[i] << 1) ^ (constant[i] >> 31 ? polynomial : 0);
        }
    }
    unsigned int hash(unsigned char* array, int size, unsigned int crc = 0)
    {
        for(int i = 0; i < size; ++i)
            crc = (crc << 8) ^ constant[(crc >> 24) ^ array[i]];
        return crc;
    }
};

class LRUCacheIndex
{
    int capacity, accessTime;
    IndexedHeap<int> ids;
public:
    LRUCacheIndex(int theCapacity): capacity(theCapacity),
        accessTime(0) {assert(capacity > 0);}
    pair<bool, int> access(int id)
    {
        int evicted = -1, *index = ids.find(id);
        if(!index && ids.getSize() == capacity) evicted = ids.deleteMin();
        ids.changeKey(accessTime--, id);
        return make_pair(bool(index), evicted);
    }
};

class PrimeTable
{
    long long maxN;
    Bitset<> table;//marks odd numbers starting from 3
    long long nToI(long long n){return (n - 3)/2;}
public:
    PrimeTable(long long primesUpto): maxN(primesUpto - 1),
        table(nToI(maxN) + 1)
    {
        assert(primesUpto > 1);
        table.setAll(true);
        for(long long i = 3; i <= sqrt(maxN); i += 2)
            if(isPrime(i))//set every odd multiple i <= k <= maxN/i to false
                for(long long k = i; i * k <= maxN; k += 2)
                    table.set(nToI(i * k), false);
    }
    bool isPrime(long long n)
    {
        assert(n <= maxN);
        return n == 2 || (n > 2 && n % 2 && table[nToI(n)]);
    }
};

struct Permutator
{
    Vector<int> p;
    Permutator(int size){for(int i = 0; i < size; ++i) p.append(i);}
    bool next()
    {//find largest i such that p[i] < p[i + 1]
        int j = p.getSize() - 1, i = j - 1;
        while(i >= 0 && p[i] >= p[i + 1]) --i;
        bool backToIdentity = i == -1;
        if(!backToIdentity)
        {//find j such that p[j] is next largest element after p[i]
            while(i < j && p[i] >= p[j]) --j;
            swap(p[i], p[j]);
        }
        p.reverse(i + 1, p.getSize() - 1);
        return backToIdentity;//true if returned to smallest permutation
    }
    bool advance(int i)
    {
        assert(i >= 0 && i < p.getSize());
        quickSort(p.getArray(), i + 1, p.getSize() - 1,
            ReverseComparator<int>());
        return next();
    }
};

struct Combinator
{
    int n;
    Vector<int> c;
    Combinator(int m, int theN): n(theN), c(m, m, -1)
    {
        assert(m <= n && m > 0);
        skipAfter(0);
    }
    void skipAfter(int i)
    {//increment c[i] and reset all c[j] for j > i
        assert(i >= 0 && i < c.getSize());
        ++c[i];
        for(int j = i + 1; j < c.getSize(); ++j) c[j] = c[j - 1] + 1;
    }
    bool next()
    {//find rightmost c[i] which can be increased
        int i = c.getSize() - 1;
        while(i >= 0 && c[i] == n - c.getSize() + i) --i;
        bool finished = i == -1;
        if(!finished) skipAfter(i);
        return finished;
    }
};

struct Partitioner
{
    Vector<int> p;
    Partitioner(int n): p(n, n, 0){assert(n > 0);}
    bool skipAfter(int k)
    {//set trailing elements to maximum values and call next
        assert(k >= 0 && k < p.getSize());
        for(int i = k; i < p.getSize(); ++i) p[i] = i;
        return next();
    }
    bool next()
    {//find rightmost p[j] which can be increased
        int m = 0, j = -1;
        for(int i = 0; i < p.getSize(); ++i)
        {
            if(p[i] < m) j = i;
            m = max(m, p[i] + 1);
        }
        bool finished = j == -1;
        if(!finished)
        {//increase it and reset the tail
            ++p[j];
            for(int i = j + 1; i < p.getSize(); ++i) p[i] = 0;
        }
        return finished;
    }
};

}
#endif
