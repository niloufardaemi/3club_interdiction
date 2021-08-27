#include "gurobi_c++.h"
#include "GRBInterface.h"
#include "LazyConstraints.h"
#include "KGraph.h"
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace std::chrono;


long s = 3;    //the parameter s in s-club

// global variables to count the number of callbacks and lazy cuts
long num_callbacks_interdiction;
long num_Lazycuts_interdiction_1;   // lazy cuts for s-clubs that form a star
long num_Lazycuts_interdiction_2;	// lazy cuts for s-clubs with leavs
long num_Lazycuts_interdiction_3;	// regular lazy cuts

double SclubTime = 0;             // the total time to solve the maximum s-club problem in all the iterations
double Finding_LCDS_time = 0;     // the total time to solve the LCDS in all the iterations

// count the number of iterations the solution is found by heuristics
long HS_Counter = 0;
long ICUT_Counter = 0;
