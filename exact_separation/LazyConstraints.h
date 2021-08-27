#ifndef LAZYCONSTRAINTS_H
#define LAZYCONSTRAINTS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "GRBInterface.h"
#include "gurobi_c++.h"
#include <ctime>
#include <memory>
#include <sstream>
#include "KGraph.h"

extern double copy_time;
using namespace std;

class LazyConstraints : public GRBCallback //Lazy constraints for solving minimum CDS problem with no latency constraints (i.e., latency-s CDS with s=n-1).
{
public:
	LazyConstraints(GRBVar *xvars, KGraph &g, long &s);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;

};

class LazyConstraints2 : public GRBCallback // Lazy constraints for r-robust latency-s CDS (r=2 only)
{
public:
	LazyConstraints2(GRBVar *xvars, KGraph &g, long &s);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;

};

class LazyConstraints3 : public GRBCallback // Lazy constraints for latency-s CDS with node-weighted distances
{
public:
	LazyConstraints3(GRBVar *xvars, KGraph &g, long &s, vector<long> W);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;

};

class LazyConstraints4 : public GRBCallback // Lazy constraints for solving reinforced minimum latency-s CDS
{
public:
	LazyConstraints4(GRBVar *xvars, KGraph &g, long &s);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;

};

class LazyConstraints5 : public GRBCallback // Lazy constraints for solving robust minimum latency-s CDS with r=2 and s=2
{
public:
	LazyConstraints5(GRBVar *xvars, KGraph &g, long &s);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;
	static long numCallbacks2;
	static double TotalCallbackTime2;
	static long numLazyCuts2;

};

#endif
