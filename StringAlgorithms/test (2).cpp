#include "SuffixArray.h"
#include <iostream>
#include <string>
using namespace std;
using namespace igmdk;

void timeSRT()
{
    string w = "awe.yfdlyiwhrkiyufsduyrwoihsiyufkshrseyfdshfriwuyriosyufhkshrweuirfysuhrfkahrfs";
    Vector<int> sa = suffixArray(w.c_str(), w.length(), SARank());
    for(int i = 0; i < w.length() ; ++i)
    {
        cout << sa[i] << " " <<w[sa[i]] << endl;
    }
}

void timeSRT2()
{
    int n = 100000;
    Vector<char> w(n, n, 0);
    for(int i = 0; i < n; ++i){
        w[i] = GlobalRNG.next();
    }
    Vector<int> sa = suffixArray(w.getArray(), n, BWTRank());
}

void timeSRT3()
{
    string s = "aaa";
    int n = s.length();
    Vector<char> w(n, n, 0);
    for(int i = 0; i < n; ++i){
        w[i] = s[i];
    }
    SuffixIndex<char> index(w);
    for(int i = 0; i < n; ++i)
    {
        DEBUG(index.sa[i]);
        DEBUG(index.lcpa[i]);
    }
    string p = "a";
    pair<int, int> lr = index.interval((char*)p.c_str(), p.length());
    DEBUG(lr.first);
    DEBUG(lr.second);
}

int main()
{
	clock_t start = clock();
	timeSRT();
	timeSRT2();
    timeSRT3();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	return 0;
}
