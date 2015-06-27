#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include "NumericalMethods.h"
#include "../Utils/DEBUG.h"
using namespace std;
using namespace igmdk;


struct SqrtFunction
{
    double value;
    SqrtFunction(double theValue):value(theValue){}
    double operator()(double x)const{return x * x - value;}
};

struct SqrtDerivative
{
    double operator()(double x)const{return 2 * x;}
};

double findSqrt(double x)
{
    return solveFor0(SqrtFunction(x), 1, x);
}

struct Stupid
{
    double x;
    double operator()()const
    {
        double d = (x - 100 + GlobalRNG.uniform01() * 10);
        return d * d;
    }
    Stupid(double theX):x(theX){}
};

struct Stupid2
{
    double operator()(Point<double, 2>const& p)const
    {
        return (p[0] - 100) * (p[0] - 100) +
            (p[1] - 100) * (p[1] - 100) + 200;
    }
};

void testNM()
{
    NelderMead<2, Stupid2> nm;
    NelderMead<2, Stupid2>::P result = nm.minimize(Point<double, 2>(0, 0));
    DEBUG(result.first[0]);
    DEBUG(result.first[1]);
    DEBUG(result.second);
}

void testILSNM()
{
    ILSNelderMead<2, Stupid2> nm;
    ILSNelderMead<2, Stupid2>::P result = nm.minimize(Point<double, 2>(0, 0));
    DEBUG(result.first[0]);
    DEBUG(result.first[1]);
    DEBUG(result.second);
}

struct FunctionOne
{
    template<typename POINT> double operator()(POINT const& point)const{return 1;}
};

void testIntegrate()
{
    double result = Trapezoid<SqrtFunction>::integrate(SqrtFunction(0), 0, 1);
    DEBUG(result);
}

struct UnitSquareTest
{
    template<typename POINT> double operator()(POINT const& point)const
    {return point[0] * point[0] + point[1] * point[1] <= 1;}
};

void testMCI()
{
    Point<double, 2> a(-1, -1), b(1, 1);
    int n = 10000000;
    //sobol 3.1415 va 2.7-e6
    //random 3.14463 va 2.7-e6
    double value = MonteCarloIntegrate<Point<double, 2>, UnitSquareTest, FunctionOne>(make_pair(a, b), n);
    DEBUG(value);
}

int main()
{
    DEBUG(findSqrt(3));
    DEBUG(StochasticMinimizer<Stupid>::min(10000, 0, 200).first);
    testNM();
    testILSNM();
    testMCI();
    testIntegrate();
    return 0;
}
