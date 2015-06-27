#include "Treap.h"
#include "LcpTreap.h"
#include "SkipList.h"
#include "Trie.h"
#include "../Utils/Debug.h"
using namespace igmdk;

struct Fat
{
	enum{SIZE = 10};
	int array[SIZE];
	Fat(int last)
	{
		for(int i = 0; i < SIZE-1; ++i)
		{
			array[i] = i;
		}
		array[SIZE-1] = last;
	}
	bool operator==(Fat const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

struct Fat2
{
	enum{SIZE = 10};
	int array[SIZE];
	Fat2(int last)
	{
		for(int i = 1; i < SIZE; ++i)
		{
			array[i] = i;
		}
		array[0] = last;
	}
	bool operator==(Fat2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

struct Fat4
{
	enum{SIZE = 40};
	unsigned char array[SIZE];
	Fat4(unsigned int last)
	{
		for(int i = 0; i < SIZE; ++i)
		{
			array[i] = i;
		}
		for(int i = 0; i < 4; ++i)
		{
            array[SIZE-1 - i] = last % 256;
            last /= 256;
		}
	}
	bool operator==(Fat4 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat4 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

struct Fat5
{
	enum{SIZE = 40};
	unsigned char array[SIZE];
	Fat5(int last)
	{
		for(int i = 0; i < SIZE; ++i)
		{
			array[i] = i;
		}
		for(int i = 0; i < 4; ++i)
		{
            array[i] = last % 256;
            last /= 256;
		}
	}
	bool operator==(Fat5 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat5 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

void testLCPTreap()
{
    typedef LcpTreap<Vector<int>, int > T;
    T t;
    Vector<int> a, b, c, d;
    a.append(0);
    a.append(1);
    a.append(0);

    b.append(0);
    b.append(0);
    b.append(1);
    b.append(1);

    c.append(0);
    c.append(1);
    c.append(0);
    c.append(0);
    c.append(0);

    t.insert(a,0);
    t.insert(b,1);
    t.insert(c,2);

    d.append(0);
    d.append(1);
    d.append(0);
    d.append(0);
    for(T::Iterator iter = t.buildIterator(t.inclusiveSuccessor(a)),
        end = t.prefixSuccessor(a); iter != end; ++iter)
    {
        DEBUG(iter->value);
    }
}

template<typename KEY> struct DefaultRank
{
    unsigned char* array;
    int size;
	DefaultRank(KEY const& key)
	{//works with pod types only
		array = (unsigned char*)&key;
		size = sizeof(key);
	}
};

template<typename KEY, typename ITEM, typename RANK = DefaultRank<KEY> >
struct TrieWrapper
{
	TernaryTreapTrie<ITEM> trie;
	ITEM* find(KEY const& key)
	{
		RANK rank(key);
		return trie.find(rank.array, rank.size);
	}
	void insert(KEY const& key, ITEM const& item)
	{
		RANK rank(key);
		trie.insert(rank.array, rank.size, item);
	}
	void remove(KEY const& key)
	{
		RANK rank(key);
		trie.remove(rank.array, rank.size);
	}
};

void timeRT()
{
	//LcpTreap<Fat5, int> trie;
	//Treap<int, int> trie;
	SkipList<int, int> trie;
    //TrieWrapper<int, int> trie;
	int N = 1500000;
	for(int i = 0; i < N; ++i)
	{
		trie.insert(i, i);
	}
	for(int j = 0; j < 100; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(trie.find(i));
			assert(*trie.find(i) == i);
			//trie.remove(i);
		}
	}
}

int main()
{
	clock_t start = clock();
    testLCPTreap();
	timeRT();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	return 0;
}
