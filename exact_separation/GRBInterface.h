#ifndef GRBINTERFACE_H
#define GRBINTERFACE_H
#include "gurobi_c++.h"
#include "KGraph.h"
#include "LazyConstraints.h"
#include <stack>

string itos(int i);

vector<bool> boolify(vector<long> &S, long n);

/* Find a large k-club and then do preprocessing, as follows:
1. Find a large distance k-clique using k-th power graph G^k and clique heuristic in G^k
2. Use DROP to find a large k-club within the large k-clique. This gives a lower bound LB.
3. Use the LB to do preprocessing in a k-core "peeling" approach.
Here, k = LB-1 (and not the k from k-club).
*/
vector<long> HeuristicAndPreprocess(KGraph &g, long k);

// DROP heuristic for max k-club, due to Bourjolly et al. (2000).
vector<long> DROP(KGraph &g, long k);
vector<long> DROP(KGraph &g, vector<bool> &W, long k);


// returns the size of the k-hop neighborhood of vertex v in G (or in G[W]).
long KNeighborhoodSize(KGraph &g, long k, long v);
long KNeighborhoodSize(KGraph &g, vector<bool> &W, long k, long v);

// Our cut-like formulation for max k-club. Uses the routine MINIMALIZE to strengthen the cuts.
vector<long> solveMaxKClub_CutLike(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt);
vector<long> MINIMALIZE(KGraph &g, long a, long b, vector<long> length_k_ab_separator, long k);


vector<long>  ICUT(KGraph &g, long k, vector<long> &BestKClub);
vector<long> ICUT_Subproblem(KGraph &g, long k, long v, long LowerBound, long &SubProblemLB, long &SubProblemUB);


// callback functions for cut-like formulation
class Kclub_callback : public GRBCallback
{
public:
	GRBVar *vars;
	KGraph g1;
	long k1;

	Kclub_callback(GRBVar *xvars, KGraph &g, long k)
	{
		vars = xvars;
		g1.Duplicate(g);
		k1 = k;
	}
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;
};





// Functions to find minimum LCDS

string itos(int i);
string statusNumtoString(int num);
// random generator function:
int myrandom(int i);

long smallestFeasibleLatency2Robust(KGraph &g);

// overall latency-s CDS function, which calls others depending on values of (r,s). Assumes hop-based distances
vector<long> solveLatencyCDS(KGraph &g, long r, long s);

// latency-s CDS funtions for specific values of (r,s) and hop-based distances
vector<long> solveMCDS(KGraph &g, long s);			 // r=1, s=arbitrary


// functions for verifying a latency-s CDS.
bool IsLatencyConstrainedCDS(KGraph &g, vector<bool> &D, long s);
bool IsLatencyConstrainedRCDS(KGraph &g, vector<bool> D, long s);
long EccentricitySubroutine(KGraph &g, vector<bool> &D, long v);
vector<long> ComputeSSSPinGBv(KGraph &g, vector<bool> &B, long v);

// functions used to minimalize length-s vertex cuts (to strengthen the formulation)
vector<long> Minimalize(KGraph &g, vector<bool> &B, long s);	// B is the infeasible CDS "solution"
vector<long> Minimalize(KGraph &g, vector<long> &CUT, long s);  // CUT are the vertices not in the CDS solution
vector< vector<long> > EnumerateFarPairs(KGraph &g, vector<bool> &B, long s);
vector< vector<long> > EnumerateFarPairsRob(KGraph &g, vector<bool> &B, long s);
vector<long> ComplementVector(vector<bool> &B, long n);

// slower versions of minimalize functions
vector<long> MinimalizeBasic(KGraph &g, vector<bool> &B, long s);
vector<long> MinimalizeBasic(KGraph &g, vector<long> &CUT, long s);
vector<long> MinimalizeNone(KGraph &g, vector<bool> &b, long s);
vector<long> MinimalizeNone(KGraph &g, vector<long> &CUT, long s);

// heuristics that we use
vector<bool> HeuristicLCDS(KGraph &g, long s);
vector<bool> HeuristicLCDS(KGraph &g, vector<bool> &heuristicSoln, long s);
vector<bool> HeuristicLCDSBestIn(KGraph &g, long s);


// heuristics that we use for r-robust latency-s CDS (case r=2)
vector<bool> HeuristicRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s);
vector<bool> MinimalizeRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s);


// Tarjan's linear-time algorithm for finding biconnected components
vector< vector<long> > FindBiconnectedComponents(KGraph &g, vector<bool> &AV);
void Bico_Sub(long v, long u, long &i, KGraph &g, vector<long> &number, vector<long> &lowopt, stack<long> &le, stack<long> &re, vector< vector<long> > &BC);

#endif
