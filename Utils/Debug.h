#ifndef DEBUG_H
#define DEBUG_H
#include <iostream>
using namespace std;
namespace igmdk{

#define DEBUG(var) cout << #var " "<< (var) << endl;

}//end namespace
#endif
