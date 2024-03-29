#ifndef UTILS_H
#define UTILS_H
#include <new>
#include <cassert>
#include <limits>
#include <utility>
#include "Debug.h"
using namespace std;
namespace igmdk{

template<typename ITEM> ITEM* rawMemory(int n)
    {return (ITEM*)::operator new(sizeof(ITEM) * n);}
void rawDelete(void* array){::operator delete(array);}
template<typename ITEM> void rawDestruct(ITEM* array, int size)
{
    for(int i = 0; i < size; ++i) array[i].~ITEM();
    rawDelete(array);
}
long long ceiling(unsigned long long n, long long divisor)
    {return n/divisor + bool(n % divisor);}
template<typename TYPE> TYPE& genericAssign(TYPE& to, TYPE const& rhs)
{
    if(&to != &rhs)
    {
        to.~TYPE();
        new(&to)TYPE(rhs);
    }
    return to;
}
template<typename KEY, typename VALUE> struct KVPair
{
    KEY key;
    VALUE value;
    KVPair(KEY const& theKey = KEY(), VALUE const& theValue = VALUE()):
        key(theKey), value(theValue) {}
};

template<typename ITEM>bool operator<=(ITEM const& lhs,
    ITEM const& rhs){return !(rhs < lhs);}
template<typename ITEM>bool operator>(ITEM const& lhs,
    ITEM const& rhs){return rhs < lhs;}
template<typename ITEM>bool operator>=(ITEM const& lhs,
    ITEM const& rhs){return !(lhs < rhs);}
template<typename ITEM>bool operator==(ITEM const& lhs,
    ITEM const& rhs){return lhs <= rhs && lhs >= rhs;}
template<typename ITEM>bool operator!=(ITEM const& lhs,
    ITEM const& rhs){return !(lhs == rhs);}
template<typename ITEM> struct DefaultComparator
{
    bool isLess(ITEM const& lhs, ITEM const& rhs)const{return lhs < rhs;}
    bool isEqual(ITEM const& lhs, ITEM const& rhs)const{return lhs == rhs;}
};
template<typename ITEM> struct ReverseComparator
{
    bool isLess(ITEM const& lhs, ITEM const& rhs)const{return rhs < lhs;}
    bool isEqual(ITEM const& lhs, ITEM const& rhs)const{return lhs == rhs;}
};
template<typename ITEM> struct PointerComparator
{
    bool isLess(ITEM const& lhs, ITEM const& rhs)const{return *lhs < *rhs;}
    bool isEqual(ITEM const& lhs, ITEM const& rhs)const{return *lhs == *rhs;}
};
template<typename ITEM> struct IndexComparator
{
    ITEM* array;
    IndexComparator(ITEM* theArray): array(theArray){}
    bool isLess(int lhs, int rhs)const{return array[lhs] < array[rhs];}
    bool isEqual(int lhs, int rhs)const{return array[lhs] == array[rhs];}
};
template<typename KEY, typename VALUE, typename COMPARATOR =
    DefaultComparator<KEY> > struct KVComparator
{
    COMPARATOR comparator;
    KVComparator(COMPARATOR const& theComparator = COMPARATOR()):
        comparator(theComparator) {}
    bool isLess(KVPair<KEY, VALUE> const& lhs, KVPair<KEY, VALUE>const& rhs)
        const{return comparator.isLess(lhs.key, rhs.key);}
    bool isEqual(KVPair<KEY, VALUE> const& lhs, KVPair<KEY, VALUE>const& rhs)
        const{return comparator.isEqual(lhs.key, rhs.key);}
};

template<typename VECTOR> struct LexicographicComparator
{
    bool isLess(VECTOR const& lhs, VECTOR const& rhs, int i)const
    {
        return i < lhs.getSize() ?
            i < rhs.getSize() && lhs[i] < rhs[i] : i < rhs.getSize();
    }
    bool isEqual(VECTOR const& lhs, VECTOR const& rhs, int i)const
    {
        return i < lhs.getSize() ?
            i < rhs.getSize() && lhs[i] == rhs[i] : i >= rhs.getSize();
    }
    bool isEqual(VECTOR const& lhs, VECTOR const& rhs)const
    {
        for(int i = 0; i < min(lhs.getSize(), rhs.getSize()); ++i)
            if(lhs[i] != rhs[i]) return false;
        return lhs.getSize() == rhs.getSize();
    }
    bool isLess(VECTOR const& lhs, VECTOR const& rhs)const
    {
        for(int i = 0; i < min(lhs.getSize(), rhs.getSize()); ++i)
        {
            if(lhs[i] < rhs[i]) return true;
            if(rhs[i] < lhs[i]) return false;
        }
        return lhs.getSize() < rhs.getSize();
    }
    int getSize(VECTOR const& value){return value.getSize();}
};

}
#endif
