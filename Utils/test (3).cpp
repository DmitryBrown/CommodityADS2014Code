#include "Stack.h"
#include "Queue.h"
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

int main()
{
    int N = 100000;
	Stack<int> a;

	for(int i=0; i < N; ++i)
	{
		a.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		a.pop();
	}
	for(int i=0; i < N; ++i)
	{
		a.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		a.pop();
	}
	Queue<int> b;
	for(int i=0; i < N; ++i)
	{
		b.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		b.pop();
	}
	for(int i=0; i < N; ++i)
	{
		b.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		b.pop();
	}
    return 0;
}
