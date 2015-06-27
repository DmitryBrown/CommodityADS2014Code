#ifndef STATISTICS_H
#define STATISTICS_H
#include "Random.h"
#include "../Utils/Sort.h"
#include "../Utils/Bits.h"
namespace igmdk{

struct NormalSummary
{
    double minimum, maximum, mean, variance;
    double confidence9973()const{return 3 * sqrt(variance);}
    friend NormalSummary operator-(NormalSummary const& a,
        NormalSummary const& b)
    {
        NormalSummary result = {a.maximum - b.minimum, a.minimum - b.maximum,
            a.mean - b.mean, a.variance + b.variance};
        return result;
    }
    double zeroDistance9973()
        {return max(0.0, abs(mean) - confidence9973());}
};

struct IncrementalStatistics
{
    double sum, squaredSum, minimum, maximum;
    long long n;
    IncrementalStatistics(): n(0), sum(0), squaredSum(0),
        minimum(numeric_limits<double>::max()), maximum(-minimum){}
    void addValue(double x)
    {
        ++n;
        maximum = max(maximum, x);
        minimum = min(minimum, x);
        sum += x;
        squaredSum += x * x;
    }
    NormalSummary getSummary()
    {
        NormalSummary result;
        result.minimum = minimum;
        result.maximum = maximum;
        result.mean = sum / n;
        double sampleVariance = max(0.0,
            (squaredSum - sum * result.mean) / (n - 1.0));
        result.variance = sampleVariance / n;
        return result;
    }
};

struct MonteCarlo
{
    template<typename FUNCTION>
    static NormalSummary simulate(FUNCTION const& f, int n)
    {
        IncrementalStatistics s;
        for(int i = 0; i < n; ++i) s.addValue(f());
        return s.getSummary();
    }
    template<typename FUNCTION, typename ITEM>
    static Vector<ITEM> simulateFull(FUNCTION const& f, int n)
    {
        Vector<ITEM> result(n);
        for(double i = 0; i < n; ++i) result.append(f());
        return result;
    }
    template<typename FUNCTION1, typename FUNCTION2> static NormalSummary
        simulateDifference(FUNCTION1 const& f1, FUNCTION2 const& f2, int n)
        {return simulate(f1, n) - simulate(f2, n);}
};

struct BootstrapResult{double mean, minus, plus;};
template<typename FUNCTION, typename DATA> BootstrapResult bootstrap(
    Vector<DATA>const& data, int n, FUNCTION const& f,
    double confidence = 0.95)
{
    assert(confidence > 0 && confidence < 1);
    Vector<double> draws(n, n, 0);
    double sum = 0;
    for(int i = 0; i < n; ++i)
    {
        Vector<DATA> resample;
        for(int j = 0; j < data.getSize(); ++j)
            resample.append(data[GlobalRNG.mod(data.getSize())]);
        sum += draws[i] = f(resample);
    }
    BootstrapResult result;
    result.mean = sum/n;
    double leftIndex = (1 - confidence)/2,
        rightIndex = leftIndex + confidence;
    quickSort(draws.getArray(), 0, draws.getSize() - 1);
    result.minus = result.mean - draws[n * leftIndex];
    result.plus = draws[n * rightIndex] - result.mean;
    return result;
}

template<typename FUNCTION> struct SpeedTester
{
    FUNCTION f;
    SpeedTester(FUNCTION const& theFunction = FUNCTION()): f(theFunction){}
    int operator()()const
    {
        int now = clock();
        f();
        return clock() - now;
    }
};

unsigned char const SobolPolys[] = {0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,19,22,
25,1,4,7,8,14,19,21,28,31,32,37,41,42,50,55,56,59,62,14,21,22,38,47,49,50,52,
56,67,70,84,97,103,115,122};
unsigned char const SobolDegs[] = {1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,
7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
class Sobol
{//SobolPolys do not represent highest and lowest 1s
    enum{B = numeric_limits<double>::digits};
    unsigned long long k;
    Vector<unsigned long long> x, v;
    double factor;
    int index(int d, int b){return d * B + b;}
public:
    bool reachedLimit(){return k >= twoPower(B);}
    static int maxD(){return sizeof(SobolDegs);}
    Sobol(int d): factor(1.0/twoPower(B)), v(d * B, d * B, 0), x(d, d, 0),
        k(1)
    {
        assert(d <= maxD());
        for(int i = 0; i < d; ++i)
            for(int j = 0; j < B; ++j)
            {
                unsigned long long value;
                int l = j - SobolDegs[i];
                if(l < 0) value = (2 * j + 1) * twoPower(B - j - 1);
                else
                {
                    value = v[index(i, l)];
                    value ^= value/twoPower(SobolDegs[i]);
                    for(int k = 1; k < SobolDegs[i]; ++k)
                        if(Bits::get(SobolPolys[i], k - 1))
                            value ^= v[index(i, l + k)];
                }
                v[index(i, j)] = value;
            }
        next();
    }
    void next()
    {
        assert(!reachedLimit());
        for(int i = 0, c = rightmost0Count(k++); i < x.getSize(); ++i)
            x[i] ^= v[index(i, c)];
    }
    double getValue(int i){return x[i] * factor;}
};

}
#endif
