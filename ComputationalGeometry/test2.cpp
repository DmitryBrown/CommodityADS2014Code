#include <iostream>
#include <cmath>
#include "KDTree.h"
#include "../ComputationalGeometry/Point.h"
#include "../RandomNumberGeneration/Random.h"
#include "../NumericalMethods/NumericalMethods.h"
using namespace igmdk;

void testKDTree()
{
    KDTree<Point<int>, int, 2> kdtree;
    int N = 1500000;
    for(int i = 0; i < N; ++i)
    {
        int p1 = (GlobalRNG.next() % 1000), p2 = (GlobalRNG.next() % 1000);
        kdtree.insert(Point<int>(min(p1, p2), max(p1, p2)), i);
    }
    bool dimensions[2];
    dimensions[0] = true;
    dimensions[1] = true;
    for(int k = 0; k < 1; ++k)
    {
        Vector<KDTree<Point<int>, int, 2>::NodeType*> result;
        int point = 0;
        kdtree.rangeQuery(Point<int>(-999999999, point), Point<int>(point, 999999999), dimensions, result);
        DEBUG(result.getSize());
    }
}

void testKDTree2()
{
    KDTree<Point<int>, int, 2> kdtree;
    int N = 1500000;
    for(int i = 0; i < N; ++i)
    {
        kdtree.insert(Point<int>(GlobalRNG.next(), GlobalRNG.next()), i);
    }
    int M = 1500;
    for(int i = 0; i < M; ++i)
    {
        assert(kdtree.nearestNeighbor(Point<int>(GlobalRNG.next(), GlobalRNG.next()), EuclideanDistance<Point<int> >::DistanceIncremental()));
        int k = 2;
        assert(kdtree.kNN(Point<int>(GlobalRNG.next(), GlobalRNG.next()), k, EuclideanDistance<Point<int> >::DistanceIncremental()).getSize() == k);
    }
}

void testVpTree()
{
    VpTree<Point<int>,int, EuclideanDistance<Point<int> >::DistanceIncremental> tree;
    int N = 5;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Point<int>(i, i), i);
    }
    for(int i = 0; i < N; ++i)
    {
        int* result = tree.find(Point<int>(i, i));
        if(result) DEBUG(*result);
        assert(result && *result == i);
    }

    Vector<VpTree<Point<int>,int, EuclideanDistance<Point<int> >::DistanceIncremental>::NodeType*> result2 = tree.distanceQuery(Point<int>(0, 4), sqrt(10));

    for(int i = 0; i < result2.getSize(); ++i)
    {
        DEBUG(result2[i]->key[0]);
        DEBUG(result2[i]->key[1]);
    }

}

void testVpTree2()
{
    VpTree<Point<int>, int, EuclideanDistance<Point<int> >::DistanceIncremental> tree;
    int N = 1500000;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Point<int>(GlobalRNG.next(), GlobalRNG.next()), i);
    }
    int M = 15000;
    for(int i = 0; i < M; ++i)
    {
        assert(tree.nearestNeighbor(Point<int>(GlobalRNG.next(), GlobalRNG.next())));
        int k = 100;
        assert(tree.kNN(Point<int>(GlobalRNG.next(), GlobalRNG.next()), k).getSize() == k);
    }
}

struct FunctionTester
{
    void operator()()const
    {
        testVpTree2();
        testVpTree();
        testKDTree();
        testKDTree2();
    }
};

int main()
{
    NormalSummary result = MonteCarlo::simulate(SpeedTester<FunctionTester>(), 1);
    DEBUG(result.minimum);
    DEBUG(result.maximum);
    DEBUG(result.mean);
    DEBUG(result.variance);
    DEBUG(result.confidence9973());
	return 0;
}
