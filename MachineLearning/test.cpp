#include "MachineLearning.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
using namespace igmdk;

//holds the digit data
struct Digit
{
    vector< vector < int > > color;
    int actualDigit;
    Digit(){}
    Digit(vector< vector < int > > color_, int actualDigit_)
    :color(color_),actualDigit(actualDigit_)
    {}
};

//reads in the digits
vector< Digit > readDigits( string filename)
{
    //the data format is 65 comma-separated integer values per line
    //the last value is the actual digit, the 64 are the colors

    vector< Digit > result;
    ifstream fin(filename.c_str());
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        Digit readDigit;
        for(int i = 0; i < 8; i++)
        {
            vector< int > row;
            for(int j = 0; j < 8; j++)
            {
                int currentDigit;
                data >> currentDigit;
                char comma;
                data >> comma;
                row.push_back( currentDigit );
            }
            readDigit.color.push_back(row);
        }
        data >> readDigit.actualDigit;
        result.push_back(readDigit);
    }


    return result;
}

class BayesianDecider
{
    NaiveBayes nb;
    //maps colors fomr 0...16 to 0...5
    int mapColor(int color){return color/3;}
    //learns all probabilities from a digit
    void learnProbColorGivenN(Digit digit)
    {
        Vector<int> values(64);
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                values.append(mapColor(digit.color[i][j]));
        nb.learn(digit.actualDigit, values);
    }
public:
    BayesianDecider( vector< Digit > trainingSet )
    :   nb(10, 64, 6)
    {
        for(int i = 0; i < trainingSet.size(); i++)
        {
            learnProbColorGivenN(trainingSet[i]);
        }
    }
    int mostLikely(Digit const & digit )
    {
        Vector<int> values(64);
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                values.append(mapColor(digit.color[i][j]));
        return nb.classify(values);
    }
};

class NNDecider
{
    NeuralNetwork nn;
    //learns all probabilities from a digit
    void learnProbColorGivenN(Digit digit)
    {
        Vector<double> values(64), results(10, 10, 0);
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                values.append(digit.color[i][j]);
        results[digit.actualDigit] = 1;
        nn.learn(values, results, 0.003);
    }
public:
    NNDecider( vector< Digit > trainingSet )
    :   nn(64)
    {
        nn.addLayer(10);
        for(int i = 0; i < 10; i++)
        {
            for(int j = -1; j < 64; j++)
            {
                nn.addConnection(0, i, j, 0);
            }
        }
        int N = 10;
        while(N--)
        {
            for(int i = 0; i < trainingSet.size(); i++)
            {
                learnProbColorGivenN(trainingSet[i]);
            }
        }
    }
    int mostLikely(Digit const & digit )
    {
        Vector<double> values(64);
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                values.append(digit.color[i][j]);
        Vector<double> result = nn.evaluate(values);
        double maxV = result[0];
        int maxI = 0;
        for(int i = 1; i < result.getSize(); ++i)
        {
            if(maxV < result[i])
            {
                maxV = result[i];
                maxI = i;
            }
        }
        return maxI;
    }
};

class FullStorageDecider
{
    NNClassifier<Point<double, 64> > nn;
    void learnDigitsGivenN(Digit digit)
    {
        Point<double, 64> instance;
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                instance[i * 8 + j] = digit.color[i][j];
        nn.learn(digit.actualDigit, instance);
    }
public:
    //inits and learns digitsGivenN from the trainingSet
    FullStorageDecider(int i, vector< Digit > trainingSet, bool cut = false)//:nn(i)
    {
        DEBUG(trainingSet.size());
        if(cut)//using only /10 of the training data picked at random
        {//cuts success rate to about 94%
            srand(time(0));

            for( int i = 0; i < trainingSet.size()/10; i++ )
            {
                int randomIndex = (int)floor((1.0*rand()/RAND_MAX)*(trainingSet.size()-1));
                learnDigitsGivenN(trainingSet[randomIndex]);
            }
        }
        else
        {
            for( int i = 0; i < trainingSet.size( ); i++ )
            {
                learnDigitsGivenN(trainingSet[i]);
            }
        }

    }
    //does the inference for passed digit image,
    //returning the most likely digit
    int mostLikely(Digit const & digit )
    {
        //changing point to use 2.8-distance and
        //using online learning givens 98.66 success rate
        //but slows down the algorithm
        Point<double, 64> instance;
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                instance[i * 8 + j] = digit.color[i][j];
        return nn.classify(instance);
    }
};

void testKMeans()
{
    Vector<Point2> points;
    int N = 20000;
    for(int i = 0; i < N; ++i)
    {
        points.append(Point2(GlobalRNG.uniform01(), GlobalRNG.uniform01()));
    }
    Vector<int> result = KMeans<Point2>::findClusters(points, 10);
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i]);
    }
}

void testKMeans2()
{
    vector< Digit > testSet = readDigits( "optdigits.tra" );
    Vector<Point<double, 64> > points;
    for(int k = 0; k < testSet.size(); k++)//
    {
        Point<double, 64> instance;
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                instance[i * 8 + j] = testSet[k].color[i][j];
        points.append(instance);
    }
    Vector<int> result = KMeans<Point<double, 64> >::findClusters(points, 10);
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i]);
    }
}



void testNN()
{
    NeuralNetwork nn(2);
    nn.addLayer(3);
    nn.addLayer(1);
    for(int i = 0; i < 3; ++i)
    {
        nn.addConnection(0, i, 0, 0);
        nn.addConnection(0, i, 1, 0);
        nn.addConnection(0, i, -1, 0);
    }
    nn.addConnection(1, 0, 0, 0);
    nn.addConnection(1, 0, 1, 0);
    nn.addConnection(1, 0, -1, 0);
    Vector<double> inputs, result;
    inputs.append(0);
    inputs.append(0);
    result.append(0);
    Vector<double> inputs2, result2;
    inputs2.append(1);
    inputs2.append(0);
    result2.append(1);
    Vector<double> inputs3, result3;
    inputs3.append(0);
    inputs3.append(1);
    result3.append(1);
    Vector<double> inputs4, result4;
    inputs4.append(1);
    inputs4.append(1);
    result4.append(1);
    DEBUG(nn.evaluate(inputs)[0]);
    DEBUG(nn.evaluate(inputs2)[0]);
    DEBUG(nn.evaluate(inputs3)[0]);
    DEBUG(nn.evaluate(inputs4)[0]);
    for(int i = 0; i < 100; ++i)
    {
        nn.learn(inputs, result, 2);
        nn.learn(inputs2, result2, 2);
        nn.learn(inputs3, result3, 2);
        nn.learn(inputs4, result4, 2);
    }
    DEBUG(nn.evaluate(inputs)[0]);
    DEBUG(nn.evaluate(inputs2)[0]);
    DEBUG(nn.evaluate(inputs3)[0]);
    DEBUG(nn.evaluate(inputs4)[0]);
}

void test(int k = 400)
{
    //testNN();
    //testKMeans2();
    //NNDecider decider( readDigits( "optdigits.tra" ) );
    //BayesianDecider decider( readDigits( "optdigits.tra" ) );
    FullStorageDecider decider(k, readDigits( "optdigits.tra" ), false);
    vector< Digit > testSet = readDigits( "optdigits.tes" );
    int successCount = 0, failureCount = 0;
    DEBUG(testSet.size());
    for(int i = 0; i < testSet.size(); i++)//
    {
        if(decider.mostLikely(testSet[i]) == testSet[i].actualDigit)
        {

            successCount++;
        }
        else
        {
            failureCount++;
        }
    }
    cout << successCount << endl;
    cout << failureCount << endl;
    cout << successCount * 1.0 / (successCount + failureCount) << endl;
}

struct FunctionTester
{
    void operator()()const
    {
        test(100);
    }
};

void testAPriori()
{
    Vector<Vector<int> > baskets;
    Vector<int> b1, b2, b3, b4;
    b1.append(0);
    b1.append(1);
    b1.append(2);
    b1.append(3);
    baskets.append(b1);
    b2.append(5);
    b2.append(1);
    b2.append(2);
    b2.append(4);
    baskets.append(b2);
    b3.append(7);
    b3.append(1);
    b3.append(2);
    b3.append(6);
    baskets.append(b3);
    b4.append(1);
    b4.append(0);
    b4.append(4);
    b4.append(6);
    baskets.append(b4);
    APriori ap;
    ap.noCutProcess(baskets, 3);
    for(LcpTreap<Vector<int>, int>::Iterator i(ap.counts.begin()); i != ap.counts.end(); ++i)
    {
        for(int j = 0; j < i->key.getSize(); ++j)
        {
            DEBUG(i->key[j]);
        }
        DEBUG(i->value);
    }
}


struct GridWorld
{
    DiscreteValueFunction u;
    int state, nEpisodes, nextState;
    double reward()
    {
        if(state == 3) return 1;
        if(state == 7) return -1;
        return -0.04;
    }
    double discountRate(){return 1;}
    double goToNextState(){state = nextState;}
    double pickNextState()
    {
        int row = state % 4, column = state / 4;
        int rows[4] = {row+1,row,row,row-1};
        int columns[4] = {column,column+1,column-1,column};
        bool set = false;
        for(int i = 0; i < 4; ++i)
        {
            if(rows[i] >= 0 && rows[i] <= 3 && columns[i] >= 0 && columns[i] <= 2 && !(rows[i] == 1 && columns[i] == 1))
            {
                int newState = rows[i] + columns[i] * 4;
                assert(state !=newState);
                if(!set || u.values[newState].first > u.values[nextState].first) {nextState = newState; set = true;}
            }
        }
        assert(set);
        return u.values[nextState].first;
    }
    bool isInFinalState(){return state == 3 || state == 7;}
    double learningRate(){return u.learningRate(state);}
    bool hasMoreEpisodes(){return nEpisodes;}
    double startEpisode()
    {
        do{state = GlobalRNG.mod(12);} while(state == 5);
        --nEpisodes;
        return u.values[state].first;}
    void updateCurrentStateValue(double delta){u.updateValue(state, delta);}
    GridWorld():nEpisodes(100), u(12){}
    void debug()
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                cout << " " << u.values[j + i * 4].first;
            }
            cout << endl;
        }
    }
};

void testReinforcement()
{
    GridWorld g;
    TDLearning(g);
    g.debug();
}

int main(int argc, char *argv[])
{
    //testReinforcement();
    //testAPriori();
    NormalSummary result = MonteCarlo::simulate(SpeedTester<FunctionTester>(), 1);
    DEBUG(result.minimum);
    DEBUG(result.maximum);
    DEBUG(result.mean);
    DEBUG(result.variance);
    DEBUG(result.confidence9973());
	return 0;
}


