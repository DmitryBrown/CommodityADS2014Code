#include "Heap.h"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace igmdk;

void timeSRT()
{
	//Heap<int> heap;
	IndexedHeap<int> heap;
	int N = 1500000;
	//IndexedArrayHeap<int> heap;//(N);
	for(int i = 0; i < N; ++i)
	{
		heap.insert(rand()%10, i);
	}
	for(int i = 0; i < N; ++i)
	{
		heap.deleteMin();
	}
}
int main()
{
	clock_t start = clock();
	timeSRT();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	return 0;
}
