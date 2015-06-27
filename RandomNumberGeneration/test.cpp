#include "Random.h"
#include "Statistics.h"
#include "../Utils/Debug.h"
using namespace igmdk;

struct maxer
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            result += observations[i];
        }
        return result / observations.getSize();
    }
};

struct maxer1
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            if(observations[i] > result) result = observations[i];
        }
        return result;
    }
};

void testSumHeap()
{
    int N = 1000000;
    SumHeap<double> sh;
    sh.add(0);
	sh.add(1.0/36);
	sh.add(2.0/36);
	sh.add(3.0/36);
	sh.add(4.0/36);
	sh.add(5.0/36);
	sh.add(6.0/36);
	sh.add(5.0/36);
	sh.add(4.0/36);
	sh.add(3.0/36);
	sh.add(2.0/36);
	sh.add(1.0/36);
    int sum = 0;
    for(int i = 0 ; i < N; ++i) sum += sh.next();
    DEBUG(sum*1.0/N);
}

int main(int argc, char *argv[])
{
    testSumHeap();
    return 0;

    Vector<double> data;
    //data.append(1000);
    //data.append(2000);
    for(int i = 0; i < 100; ++i) data.append(100 + GlobalRNG.mod(100));
    BootstrapResult r2 = bootstrap(data, 10000, maxer1(), 0.9973);
    DEBUG(r2.mean);
    DEBUG(r2.minus);
    DEBUG(r2.plus);


    DEBUG(Sobol::maxD());
    DEBUG(sizeof(SobolPolys));
    DEBUG(sizeof(SobolDegs));
    Sobol so(1);
    for(int i = 0; i < 10; ++i)
    {
        for(int j = 0; j < 1; ++j)
            DEBUG(so.getValue(j));
        so.next();
    }
    //ARC4 x;
    MRG32k3a x;
    x.jumpAhead();
    unsigned long long N = 1 << 3;
    unsigned long long dummy = 0;
    while(N--) dummy += x.next();
    DEBUG(dummy);

    SumHeap<double> st;
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    DEBUG(st.total());
    DEBUG(st.find(0.1));
    DEBUG(st.find(0.3));
    DEBUG(st.find(0.6));
    DEBUG(st.find(0.9));

    DEBUG(st.find(0));
    DEBUG(st.find(0));
    DEBUG(st.find(0.5));
    DEBUG(st.find(1));

    DEBUG(st.cumulative(0));
    DEBUG(st.cumulative(1));
    DEBUG(st.cumulative(2));
    DEBUG(st.cumulative(3));
    DEBUG(st.cumulative(4));
    DEBUG(st.find(st.cumulative(0)));
    DEBUG(st.find(st.cumulative(1)));
    DEBUG(st.find(st.cumulative(2)));
    DEBUG(st.find(st.cumulative(3)));
    DEBUG(st.find(st.cumulative(4)));
    for(int i = 0; i < 100; ++i)
    {
        DEBUG(st.next());
    }

    clock_t start = clock();
    //Xorshift random;
    //Xorshift64 random;
    //ImprovedXorshift random;
    QualityXorshift64 random;
    //ImprovedXorshift64 random;
    unsigned long long sum = 0;
    for(int i = 0; i < 1000000; ++i)
    {
		sum += random.next();
	}
	DEBUG(sum);
	clock_t end = clock();
	int time = (end - start);
    cout << "IX: " << time << endl;

    if(true)
    {
        DEBUG(GlobalRNG.uniform01());
        DEBUG(GlobalRNG.uniform(10, 20));
        DEBUG(GlobalRNG.normal01());
        DEBUG(GlobalRNG.normal(10, 20));
        DEBUG(GlobalRNG.exponential01());
        DEBUG(GlobalRNG.gamma1(0.5));
        DEBUG(GlobalRNG.gamma1(1.5));
        DEBUG(GlobalRNG.weibull1(20));
        DEBUG(GlobalRNG.erlang(10, 2));
        DEBUG(GlobalRNG.chiSquared(10));
        DEBUG(GlobalRNG.t(10));
        DEBUG(GlobalRNG.logNormal(10, 20));
        DEBUG(GlobalRNG.beta(0.5, 0.5));
        DEBUG(GlobalRNG.F(10 ,20));
        DEBUG(GlobalRNG.cauchy01());
        DEBUG(GlobalRNG.binomial(0.7, 20));
        DEBUG(GlobalRNG.geometric(0.7));
        DEBUG(GlobalRNG.poisson(0.7));
		//system("PAUSE");
	}
	int M = 100000;
	double average = 0;
	for(int i = 0; i < M; ++i)
	{
        average += GlobalRNG.beta(0.5, 0.5);
	}
	DEBUG(average/M);
    return 0;
}
