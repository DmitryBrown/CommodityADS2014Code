#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Utils/Sort.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../RandomTreap/Treap.h"
#include "../ComputationalGeometry/Point.h"
#include "../Optimization/Metaheuristics.h"
namespace igmdk{

bool haveDifferentSign(double a, double b){return (a < 0) != (b < 0);}
template<typename FUNCTION> double solveFor0(FUNCTION const& f,
    double xLeft, double xRight)
{
    double yLeft = f(xLeft), xMiddle;
    assert(xRight >= xLeft && haveDifferentSign(yLeft, f(xRight)));
    for(;;)
    {
        xMiddle = (xLeft + xRight) / 2;
        double yMiddle = f(xMiddle), prevDiff = xRight - xLeft;
        if(haveDifferentSign(yLeft, yMiddle)) xRight = xMiddle;
        else
        {
            xLeft = xMiddle;
            yLeft = yMiddle;
        }
        if(xRight - xLeft >= prevDiff) break;
    }
    return xMiddle;
}

template<typename FUNCTION, typename DERIVATIVE> double solveFor0Newton(
    FUNCTION const& f, DERIVATIVE const& derivative, double x)
{
    for(double prevDiff = numeric_limits<double>::max();;)
    {
        double y = f(x), oldX = x;
        x -= y/derivative(x);
        double diff = abs(x - oldX);
        if(diff >= prevDiff) break;
        prevDiff = diff;
    }
    return x;
}

template<typename FUNCTION> pair<double, double> minimize(FUNCTION const& f,
    double xLeft, double xRight)
{
    assert(xLeft <= xRight);
    double GR = 0.618, xMiddle = xLeft * GR + xRight * (1 - GR),
        yMiddle = f(xMiddle);
    for(;;)
    {
        bool chooseR = xRight - xMiddle > xMiddle - xLeft;
        double prevDiff = xRight - xLeft, xNew = GR * xMiddle + (1 - GR) *
            (chooseR ? xRight : xLeft), yNew = f(xNew);
        if(yNew < yMiddle)
        {
            (chooseR ? xLeft : xRight) = xMiddle;
            xMiddle = xNew;
            yMiddle = yNew;
        }
        else (chooseR ? xRight : xLeft) = xNew;
        if(xRight - xLeft >= prevDiff) break;
    }
    return make_pair(xMiddle, yMiddle);
}

template<typename FUNCTION> double findIntervalBound(FUNCTION const& f,
    double guess, double d)
{//run with d < 0 for the left bound and d > 0 for the right bound
    for(double yBest = f(guess);;)
    {
        double xNext = guess + d, yNext = f(xNext);
        if(yNext >= yBest) return xNext;
        yBest = yNext;
    }
}

template<typename FUNCTION> struct StochasticMinimizer
{
    struct Function
    {
        int n;
        Function(int theN): n(theN){}
        double operator()(double x)const
            {return MonteCarlo::simulate(FUNCTION(x), n).mean;}
    };
    static pair<double, double> min(int nSimulations, double xLeft,
        double xRight)
        {return minimize(Function(nSimulations), xLeft, xRight);}
};

template<int D, typename FUNCTION> struct DimentionMinimizer
{
    typedef Point<double, D> P;
    struct Function
    {
        FUNCTION f;
        P point;
        int variable;
        double operator()(double x)
        {
            point[variable] = x;
            return f(point);
        }
    };
    P minimize(P const&initialGuess, double d,
        FUNCTION const& f = FUNCTION())
    {
        assert(d > 0);
        P x = initialGuess;
        for(int i = 0; i < D; ++i)
        {
            Function f1D = {f, x, i};
            x[i] = minimize(f1D, findIntervalBound(f1D, x[i], -d),
                findIntervalBound(f1D, x[i], d)).first;
        }
        return x;
    }
};

template<int D, typename FUNCTION> struct NelderMead
{
    typedef Point<double, D> Vertex;
    Vertex vertexSum;//incremental centroid
    typedef pair<Vertex, double> P;
    P simplex[D + 1];
    FUNCTION f;
    double scale(P& high, double factor)
    {
        P result;
        //convex combination of the high point and the
        //centroid of the remaining vertices
        //centroid = (vertexSum - high)/D and
        //result = centroid * (1 - factor) + high * factor
        double centroidFactor = (1 - factor)/D;
        result.first = vertexSum * centroidFactor +
            high.first * (factor - centroidFactor);
        result.second = f(result.first);
        if(result.second < high.second)
        {//accept scaling if improving
            vertexSum += result.first + high.first * -1;
            high = result;
        }
        return result.second;
    }
    NelderMead(FUNCTION const& theFunction = FUNCTION()): f(theFunction){}

    P minimize(Vertex const& initialGuess, int maxIterations = 1000,
        double precision = numeric_limits<double>::epsilon())
    {
        for(int i = 0; i < D + 1; ++i)
        {
            simplex[i].first = initialGuess;
            if(i > 0) simplex[i].first[i - 1] += GlobalRNG.uniform01();
            simplex[i].second = f(simplex[i].first);
            vertexSum += simplex[i].first;
        }
        for(;;)
        {//calculate high, low, and nextHigh, which must be all different
            int high = 0, nextHigh = 1, low = 2;
            if(simplex[high].second < simplex[nextHigh].second)
                swap(high, nextHigh);
            if(simplex[nextHigh].second < simplex[low].second)
            {
                swap(low, nextHigh);
                if(simplex[high].second < simplex[nextHigh].second)
                    swap(high, nextHigh);
            }
            for(int i = 3; i < D + 1; ++i)
            {
                if(simplex[i].second < simplex[low].second) low = i;
                else if(simplex[i].second > simplex[high].second)
                {
                    high = i;
                    nextHigh = high;
                }
                else if(simplex[i].second > simplex[nextHigh].second)
                    nextHigh = i;
            }
            if(!maxIterations-- || simplex[high].second - simplex[low].second
               < precision) return simplex[low];
            //try to reflect
            double value = scale(simplex[high], -1);
            //try to double if better then low
            if(value <= simplex[low].second) scale(simplex[high], 2);
            else if(value >= simplex[nextHigh].second)
            {//try to halve if worse then next high
                double yHi = simplex[high].second;
                if(scale(simplex[high], 0.5) >= yHi)
                {//contract all to get rid of the high point
                    vertexSum = simplex[low].first;
                    for(int i = 0; i < D + 1; ++i) if(i != low)
                    {
                        vertexSum += simplex[i].first = (simplex[i].first +
                            simplex[low].first) * 0.5;
                        simplex[i].second = f(simplex[i].first);
                    }
                }
            }
        }
    }
};

template<int D, typename FUNCTION> struct ILSNelderMead
{
    typedef NelderMead<D, FUNCTION> NM;
    typedef typename NM::P P;
    struct Move
    {
        NM &nm;
        P current, best;
        int maxIterations;
        double precision;
        void localSearchBest()
            {current = nm.minimize(current.first, maxIterations, precision);}
        void bigMove()
        {
            for(int i = 0; i < D; ++i)
                current.first[i] = GlobalRNG.cauchy01();
        }
        void updateBest()
            {if(best.second > current.second) best = current;}
    };
    NM nelderMead;
    ILSNelderMead(FUNCTION const& theFunction = FUNCTION()):
        nelderMead(theFunction){}
    P minimize(typename NM::Vertex const& initialGuess, int maxJumps =
        1000, int maxIterations = 1000, double precision = 0.001)
    {
        P initial(initialGuess, numeric_limits<double>::max());
        Move move = {nelderMead, initial, initial, maxIterations, precision};
        iteratedLocalSearch(move, maxJumps);
        return move.best;
    }
};

template<typename POINT> double boxVolume(pair<POINT, POINT> const& box)
{
    double result = 1;
    for(int i = 0; i < box.first.getSize(); ++i)
        result *= box.second[i] - box.first[i];
    return result;
}
template<typename POINT, typename TEST, typename FUNCTION>
double MonteCarloIntegrate(pair<POINT, POINT> const& box, int n,
    TEST const& isInside = TEST(), FUNCTION const& f = FUNCTION())
{
    double sum = 0;
    for(int i = 0; i < n; ++i)
    {
        POINT point;
        for(int j = 0; j < point.getSize(); ++j)
            point[j] = GlobalRNG.uniform(box.first[j], box.second[j]);
        if(isInside(point)) sum += f(point);
    }
    return boxVolume(box) * sum/n;
}

template<typename FUNCTION> class Trapezoid
{
    double sum;
    int nIntervals;
public:
    Trapezoid(): nIntervals(1){}
    double addLevel(FUNCTION const& f, double xLeft, double xRight)
    {
        double dX = (xRight - xLeft)/nIntervals;
        if(nIntervals == 1) sum = (f(xLeft) + f(xRight))/2;
        else
        {//nIntervals/2 is the number of added points spaced 2*dX apart
            double x = xLeft + dX;
            for(int i = 0; i < nIntervals/2; ++i, x += 2 * dX) sum += f(x);
        }
        nIntervals *= 2;
        return dX * sum;
    }
    static double integrate(FUNCTION const& f, double xLeft, double xRight,
        double maxXError = 0, int maxIterations = 10, int minIterations = 5)
    {//need minIterations for functions that are same in several places
        Trapezoid t;
        double result = t.addLevel(f, xLeft, xRight), oldResult = 0;
        while(--maxIterations > 0 && (--minIterations > 0 ||
            abs(result - oldResult) >= maxXError))
        {
            oldResult = result;
            result = t.addLevel(f, xLeft, xRight);
        }
        return result;
    }
};

class DynamicLinearInterpolation
{
    Treap<double, double> values;
public:
    double findMin()
    {
        assert(!values.isEmpty());
        return values.findMin()->key;
    }
    double findMax()
    {
        assert(!values.isEmpty());
        return values.findMax()->key;
    }
    double evaluate(double x)
    {
        assert(x >= findMin() && x <= findMax());
        double *y = values.find(x);
        if(y) return *y;
        Treap<double, double>::NodeType* left = values.predecessor(x),
            *right = values.successor(x);
        return left->value + (right->value - left->value) *
            (x - left->key)/(right->key - left->key);
    }
    void remove(double x){values.remove(x);}
    bool contains(double x){return values.find(x);}
    void insert(double x, double y){values.insert(x, y);}
};

}
#endif
