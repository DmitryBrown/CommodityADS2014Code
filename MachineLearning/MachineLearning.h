#ifndef MACHINELEARNING_H
#define MACHINELEARNING_H
#include "../Utils/Utils.h"
#include "../MiscAlgs/Misc.h"
#include "../ComputationalGeometry/KDTree.h"
#include "../ComputationalGeometry/Point.h"
#include "../NumericalMethods/DenseMatrix.h"
#include "../RandomTreap/LCPTreap.h"
#include <cmath>
namespace igmdk{

class NeuralNetwork
{
    Vector<double> inputs;
    struct Neuron
    {
        Vector<int> sources;
        Vector<double> weights;
        double output, error;
    };
    Vector<Vector<Neuron> > layers;
    double activation(double x){return 1/(1 + exp(-x));}
    void propagateInputs(Vector<double> const& theInputs)
    {
        inputs = theInputs;
        for(int i = 0; i < layers.getSize(); ++i)
        {
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& neuron = layers[i][j];
                double sum = neuron.error = 0;
                for(int k = 0; k < neuron.sources.getSize(); ++k)
                    sum += neuron.weights[k] *
                        getInput(i, neuron.sources[k]);
                neuron.output = activation(sum);
            }
        }
    }
    double getInput(int layer, int source)
    {
        return source == -1 ? 1 : layer == 0 ?
            inputs[source] : layers[layer - 1][source].output;
    }
public:
    NeuralNetwork(int nInputs): inputs(nInputs, nInputs){}
    void addLayer(int nNeurons)
        {layers.append(Vector<Neuron>(nNeurons, nNeurons));}
    void addConnection(int layer, int neuron, int to, double weight)
    {
        layers[layer][neuron].sources.append(to);
        layers[layer][neuron].weights.append(weight);
    }
    Vector<double> evaluate(Vector<double> const& theInputs)
    {
        propagateInputs(theInputs);
        Vector<double> result;
        for(int i = 0; i < layers.lastItem().getSize(); ++i)
            result.append(layers.lastItem()[i].output);
        return result;
    }
    void learn(Vector<double> const& theInputs,
        Vector<double> const& results, double learningRate)
    {
        assert(results.getSize() == layers.lastItem().getSize());
        propagateInputs(theInputs);
        for(int j = 0; j < layers.lastItem().getSize(); ++j)
            layers.lastItem()[j].error = learningRate *
                (results[j] - layers.lastItem()[j].output);
        for(int i = layers.getSize() - 1; i >= 0; --i)
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& neuron = layers[i][j];
                double temp = neuron.error * neuron.output *
                    (1 - neuron.output);
                for(int k = 0; k < neuron.sources.getSize(); ++k)
                {
                    int source = neuron.sources[k];
                    if(i > 0 && source != -1) layers[i - 1][source].error +=
                        neuron.weights[k] * temp;
                    neuron.weights[k] += getInput(i, source) * temp;
                }
            }
    }
};

class NaiveBayes
{
    int nClasses, nFeatures;
    struct Feature
    {
        int count;
        Vector<int> valueCounts;
        Feature(int nValues): count(nValues),
            valueCounts(nValues, nValues, 1){}
    };
    Vector<Feature> frequencies;
    int index(int classNumber, int feature)
        {return classNumber * nFeatures + feature;}
public:
    NaiveBayes(int theNClasses, int theNFeatures, int nValues):
        nClasses(theNClasses), nFeatures(theNFeatures), frequencies(
        nClasses * nFeatures, nClasses * nFeatures, Feature(nValues)){}
    void learn(int classNumber, Vector<int> const& featureValues)
    {
        for(int i = 0; i < featureValues.getSize(); ++i)
        {
            Feature& frequency = frequencies[index(classNumber, i)];
            ++frequency.valueCounts[featureValues[i]];
            ++frequency.count;
        }
    }
    int classify(Vector<int> const& featureValues)
    {
        double maxLL;
        int bestClass;
        for(int i = 0; i < nClasses; ++i)
        {
            double ll = 0;
            for(int j = 0; j < featureValues.getSize(); j++)
            {
                Feature& frequency = frequencies[index(i, j)];
                ll += log(frequency.valueCounts[featureValues[j]] * 1.0/
                    frequency.count);
            }
            if(i == 0 || maxLL < ll)
            {
                maxLL = ll;
                bestClass = i;
            }
        }
        return bestClass;
    }
};

template<typename OBJECT, typename INDEX = VpTree<OBJECT, int,
typename EuclideanDistance<OBJECT>::Distance> > class NNClassifier
{
    INDEX instances;
public:
    void learn(int classNumber, OBJECT const& instance)
        {instances.insert(instance, classNumber);}
    int classify(OBJECT const& instance)
        {return instances.nearestNeighbor(instance)->value;}
};

Vector<double> regression(DenseMatrix<double> const& X,
    Vector<double> const& Y)
{
    assert(Y.getSize() == X.rows);
    DenseMatrix<double> Xtranspose = X.transpose(),
        inv = DenseLUP<double>(Xtranspose * X).inverse();
    return inv * Xtranspose * Y;
}

DenseMatrix<double> randomProjection(int fromD, int toD)
{
    DenseMatrix<double> result(toD, fromD);
    for(int i = 0; i < result.rows; ++i)
        for(int j = 0; j < result.columns; ++j)
            result(i, j) = GlobalRNG.bernoulli(0.5) ? -1 : 1;
    return result;
}

template<typename POINT, typename TREE = VpTree<POINT, int,
    typename EuclideanDistance<POINT>::Distance> > struct KMeans
{
    static Vector<int> findClusters(Vector<POINT>& points, int k,
        int maxIterations = 1000)
    {
        assert(k > 0 && k <= points.getSize());
        //generate initial assignment
        Vector<int> assignments(points.getSize());
        //each cluster has at least 1 point, rest are random
        for(int i = 0; i < points.getSize(); ++i)
            assignments.append(i < k ? i : GlobalRNG.next() % k);
        bool converged = false;
        for(int m = 0; !converged && m < maxIterations; ++m)
        {//calculate centroids
            Vector<int> counts(k, k, 0);
            Vector<POINT> centroids(k, k);
            for(int i = 0; i < points.getSize(); ++i)
            {
                ++counts[assignments[i]];
                centroids[assignments[i]] += points[i];
            }
            for(int i = 0; i < k; ++i) centroids[i] *= 1.0/counts[i];
            TREE t;
            for(int i = 0; i < k; ++i) t.insert(centroids[i], i);
            //assign each point to the closest centroid
            converged = true;
            for(int i = 0; i < points.getSize(); ++i)
            {
                int best = t.nearestNeighbor(points[i])->value;
                if(best != assignments[i])
                {
                    converged = false;
                    assignments[i] = best;
                }
            }
        }
        return assignments;
    }
};

double UCB1(double averageValue, int nTries, int totalTries)
    {return averageValue + sqrt(2 * log(totalTries)/nTries);}

template<typename PROBLEM> void TDLearning(PROBLEM& p)
{
    while(p.hasMoreEpisodes())
    {
        double valueCurrent = p.startEpisode();
        while(!p.isInFinalState())
        {
            double valueNext = p.pickNextState();
            p.updateCurrentStateValue(p.learningRate() * (p.reward() +
                p.discountRate() * valueNext - valueCurrent));
            p.goToNextState();
            valueCurrent = valueNext;
        }
        p.updateCurrentStateValue(p.learningRate() *
            (p.reward() - valueCurrent));
    }
}

struct DiscreteValueFunction
{
    Vector<pair<double, int> > values;
    double learningRate(int state){return 1.0/values[state].second;}
    void updateValue(int state, double delta)
    {
        ++values[state].second;
        values[state].first += delta;
    }
    DiscreteValueFunction(int n): values(n, n, make_pair(0.0, 1)){}
};

struct LinearCombinationValueFunction
{
    Vector<double> weights;
    int n;
    double learningRate(){return 1.0/n;}
    void updateWeights(Vector<double> const& stateFeatures, double delta)
    {//set one of the state features to 1 to have a bias weight
        assert(stateFeatures.getSize() == weights.getSize());
        for(int i = 0; i < weights.getSize(); ++i)
            weights[i] += delta * stateFeatures[i];
        ++n;
    }
    LinearCombinationValueFunction(int theN): weights(theN, theN, 0), n(1){}
};

struct APriori
{
    LcpTreap<Vector<int>, int> counts;
    int processBasket(Vector<int> const& basket, int round,
        int rPrevMinCount = 0, int r1MinCount = 0)
    {
        int addedCount = 0;
        if(basket.getSize() > round)
        {
            Combinator c(round, basket.getSize());
            do//prepare the current combination of ids, needn't sort if each
            {//basket is already sorted
                Vector<int> key, single;
                for(int i = 0; i < round; ++i) key.append(basket[c.c[i]]);
                quickSort(key.getArray(), 0, key.getSize() - 1);
                int* count = counts.find(key);
                if(count) ++*count;//combination is frequent if already
                else if(round == 1)//frequent or round is 1
                {
                    counts.insert(key, 1);
                    ++addedCount;
                }
                else//combination is frequent if the last item and
                {//combination without the last item are both frequent
                    single.append(key.lastItem());
                    if(*counts.find(single) >= r1MinCount)
                    {
                        key.removeLast();
                        if(*counts.find(key) >= rPrevMinCount)
                        {
                            key.append(single[0]);
                            counts.insert(key, 1);
                            ++addedCount;
                        }
                    }
                }
            }while(!c.next());
        }
        return addedCount;
    }
    void noCutProcess(Vector<Vector<int> >const& baskets, int nRounds)
    {
        for(int k = 1; k <= nRounds; ++k)
            for(int i = 0; i < baskets.getSize(); ++i)
                processBasket(baskets[i], k);
    }
};

}//end namespace
#endif

