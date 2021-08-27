#include "GRBInterface.h"
#include "LazyConstraints.h"
#include "KGraph.h"
#include <sstream>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <queue>
#include <limits.h>
#include <string.h>
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

string itos(int i) { stringstream s; s << i; return s.str(); }

vector<bool> boolify(vector<long> &S, long n)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++) Sbool[S[i]] = true;
	return Sbool;
}

vector<long> HeuristicAndPreprocess(KGraph &g, long k)
{
	KGraph gk = g.CreatePowerGraph(k);
	vector<long> rd;  // right-degree of degeneracy ordering.
	vector<long> degeneracyordering = gk.FindDegeneracyOrdering(rd);
	vector<long> kclique = gk.FindHeuristicClique(degeneracyordering, rd);
	cerr << "k-clique (heuristic) size is = " << kclique.size() << endl;

	// perform DROP heuristic, with kclique as starting point.
	vector<bool> HeuristicSolution = boolify(kclique, g.n); // make the heuristic kclique into a boolean form
	vector<long> DROP_Solution = DROP(g, HeuristicSolution, k);
	long lb = DROP_Solution.size();
	cerr << "Drop heuristic gives LB = " << lb << endl;

	// perform preprocessing
	vector<long> kcorevertices = gk.FindVerticesOfKCore(degeneracyordering, rd, lb - 1);
	cerr << "Preprocessed instances has this many vertices = " << kcorevertices.size() << endl;
	g.FindInducedGraph(kcorevertices);

	return DROP_Solution;
}

vector<long> DROP(KGraph &g, long k)
{
	vector<bool> W(g.n, true);	// our k-club (eventually)
	return DROP(g, W, k);
}
vector<long> DROP(KGraph &g, vector<bool> &W, long k)
{
	vector<long> near(g.n, 0);	// the number of vertices "nearby" to a vertex (within k-hops)

	long Wsize = count(W.begin(), W.end(), true);
	// while W is not an k-club, delete some "bad" vertex. The number of vertices in W is size.
	for (long size = Wsize; size >= 0; size--)
	{
		// compute how many vertices are "nearby" to each vertex.
		for (long i = 0; i < g.n; i++)
		{
			if (!W[i]) continue;
			near[i] = KNeighborhoodSize(g, W, k, i);
		}

		// pick vertex w (in W) with smallest "nearby" vertices
		long w;	// vertex with smallest number of nearby vertices
		long smallestNearby = size; // number of vertices nearby to w.
		for (long i = 0; i < g.n; i++)
		{
			if (!W[i]) continue;
			if (near[i] < smallestNearby)
			{
				w = i;
				smallestNearby = near[i];
			}
		}

		// check if k-club. If not, then remove vertex w.
		if (smallestNearby == size) break;
		W[w] = false;
	}

	// convert W to vector<long> form, and return
	vector<long> Wvector;
	for (long i = 0; i < g.n; i++)
	{
		if (W[i])
		{
			Wvector.push_back(i);
		}
	}
	return Wvector;
}

long Kclub_callback::numCallbacks = 0;
double Kclub_callback::TotalCallbackTime = 0;
long Kclub_callback::numLazyCuts = 0;



// calback function for our main method - using MINIMALIZE 
// to obtain minimal subset of length-k a,b-separator
void Kclub_callback::callback()
{
	try {
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks++;
			time_t start = clock();

			// get the ''solution (S)'' from Gurobi
			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			// make it boolean and call it D
			vector<bool> D(g1.n, false);

			for (long i = 0; i < g1.n; i++)
			{
				if (x[i] > 0.5) D[i] = true;
			}

			if (count(D.begin(), D.end(), false) != g1.n)
			{
				vector<long> SelectedVertices;
				vector<long> C_Prime;

				for (long i = 0; i < g1.n; i++)
				{
					if (D[i]) SelectedVertices.push_back(i);
					else C_Prime.push_back(i);
				}

				// create G[D], which we call g2
				vector<long> Rmap;
				KGraph g2 = g1.CreateInducedGraph(D, Rmap);
				for (long i = 0; i < g2.n; i++)
				{
					vector<long> dist_from_i_to = g2.ShortestPathsUnweighted(i);
					for (long j = i + 1; j < g2.n; j++)
					{
						if (dist_from_i_to[j] > k1)
						{
							long a = SelectedVertices[i];
							long b = SelectedVertices[j];

							// C_Prime is a length-s a,b-separator. now make it minimal
							vector<long> minimal_length_k_ab_separator = MINIMALIZE(g1, a, b, C_Prime, k1);
							GRBLinExpr expr = vars[a] + vars[b];
							for (long q = 0; q < minimal_length_k_ab_separator.size(); q++)
							{
								long v = minimal_length_k_ab_separator[q];
								expr -= vars[v];
							}

							addLazy(expr <= 1);
							numLazyCuts++;
						}
					}
				}
			}
			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] x;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}


vector<long> MINIMALIZE(KGraph &g, long a, long b, vector<long> ab_Separator, long k)
{
	vector<long> MinimalCut; // what we return at end of function
	vector<bool> Cut(g.n, false); // a boolean representation of the cut. 

								  // first do some linear-time preprocessing to remove vertices that are not on length-k a,b-path
								  //		for example, this removes vertices that belong to a different component in G.
	vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
	vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
	for (long i = 0; i < ab_Separator.size(); i++)
	{
		long v = ab_Separator[i];
		if (dist_from_a[v] + dist_from_b[v] <= k)
		{
			Cut[v] = true;  // vertex v was in ab_Separator AND belongs to a length-k a,b-path in graph G
		}
	}

	// initialize VerticesInGc = V \ Cut
	vector<bool> VerticesInGc(g.n, true);
	for (long i = 0; i < g.n; i++) if (Cut[i]) VerticesInGc[i] = false;

	// now run the normal minimalize algorithm.
	for (long c = 0; c < g.n; c++)
	{
		if (!Cut[c]) continue; // only test for cut vertices
		if (dist_from_a[c] == 1 && dist_from_b[c] == 1) continue; // if c \in N(a)\cap N(b), then c belongs to every minimal cut.

																  // put vertex c in G_c
		VerticesInGc[c] = true;

		// compute distances from c in G_c
		vector<long> dist_from_c = g.ShortestPathsUnweighted(c, VerticesInGc);

		// check if vertex c would close the distance from a to b to be at most k.
		if (dist_from_c[a] + dist_from_c[b] <= k) // vertex c remains in cut
		{
			VerticesInGc[c] = false;
		}
		else // vertex c is pulled from the cut.
		{
			Cut[c] = false;
		}
	}
	for (long i = 0; i < g.n; i++) if (Cut[i]) MinimalCut.push_back(i);
	return MinimalCut;
}


long KNeighborhoodSize(KGraph &g, vector<bool> &W, long k, long v)
{
	if (!W[v])
	{
		cerr << "\n ERROR: K-neighborhoodSize. You are calculating distances across W nodes starting from some vertex v which does not belong to W.";
		return 0;
	}
	vector<long> dist = g.ShortestPathsUnweighted(v, W);
	long nsize = 0;		// the number of nodes j with dist(v,j) <= s (including node v).
	for (long i = 0; i < g.n; i++) if (dist[i] <= k) nsize++;
	return nsize;
}
long KNeighborhoodSize(KGraph &g, long k, long v)
{
	vector<bool> W(g.n, true);
	return KNeighborhoodSize(g, W, k, v);
}

vector<long> solveMaxKClub_CutLike(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt)
{
	vector<long> S;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = HeuristicSolution.size();

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
		model.update();


		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();


		//Adding lazy constraints.
		Kclub_callback cb = Kclub_callback(X, g, k);
		model.setCallback(&cb);

		// Providing initial solution
		for (long i = 0; i < g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < HeuristicSolution.size(); i++)
		{
			long v = HeuristicSolution[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		model.optimize();

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << Kclub_callback::numCallbacks << endl;
		cerr << "# lazy cuts = " << Kclub_callback::numLazyCuts << endl;

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}



vector<long> ICUT_Subproblem(KGraph &g, long k, long v_i, long LowerBound, long &SubProblemLB, long &SubProblemUB)
{
	vector<long> S;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		env.set(GRB_DoubleParam_Cutoff, LowerBound);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < LowerBound)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
		model.update();

		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}

		// Fixing X[v_i] to 1.
		model.addConstr(X[v_i] == 1);

		model.update();


		//Adding lazy constraints.
		Kclub_callback cb = Kclub_callback(X, g, k);
		model.setCallback(&cb);


		model.optimize();


		SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
		SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_CUTOFF)
		{
			SubProblemLB = LowerBound;
			SubProblemUB = LowerBound;
			return S;
		}

		else if (status == GRB_OPTIMAL)
		{
			SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
			SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}


vector <long> ICUT(KGraph &g, long k, vector<long> &BestKClub)
{
	time_t start = clock();
	KGraph gk = g.CreatePowerGraph(k);

	long BestLowerBound = BestKClub.size();
	long BestUpperBound = BestKClub.size();

	vector<long> degeneracyordering = gk.FindDegeneracyOrdering();
	vector<bool> T(g.n, false);

	for (long i = g.n - 1; i >= 0; i--)
	{
		vector<long> SubproblemVertices;

		long v = degeneracyordering[i];
		T[v] = true;
		vector<long> dist_from_v = g.ShortestPathsUnweighted(v, T);

		for (long j = 0; j < g.n; j++)
			if (dist_from_v[j] <= k)
				SubproblemVertices.push_back(j);

		if (SubproblemVertices.size() <= BestLowerBound) continue;

		vector<long> Rmap;
		KGraph g_subproblem = g.CreateInducedGraph(SubproblemVertices, Rmap);

		long SubProblemLB;
		long SubProblemUB;
		vector<long> SubproblemSolution = ICUT_Subproblem(g_subproblem, k, Rmap[v], BestLowerBound, SubProblemLB, SubProblemUB);

		BestLowerBound = max(BestLowerBound, SubProblemLB);
		BestUpperBound = max(BestUpperBound, SubProblemUB);

		if (SubproblemSolution.size() >= BestLowerBound)
		{
			BestKClub = SubproblemSolution;
			for (long q = 0; q < BestLowerBound; q++)
			{
				long v = BestKClub[q];
				long w = SubproblemVertices[v];
				BestKClub[q] = w;
			}
		}
	}

	cout << "bestLB = " << BestLowerBound << ", bestUB = " << BestUpperBound << " ";
	return BestKClub;
}






// For finding minimumm LCDS

//string itos(int i) { stringstream s; s << i; return s.str(); }

string statusNumtoString(int num)
{
	if (num == 1) return "Model is loaded, but no solution information is available.";
	else if (num == 2) return "Model was solved to optimality (subject to tolerances), and an optimal solution is available.";
	else if (num == 3) return "Model was proven to be infeasible.";
	else if (num == 4) return "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.";
	else if (num == 5) return "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.";
	else if (num == 6) return "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.";
	else if (num == 7) return "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.";
	else if (num == 8) return "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.";
	else if (num == 9) return "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.";
	else if (num == 10) return "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.";
	else if (num == 11) return "Optimization was terminated by the user.";
	else if (num == 12) return "Optimization was terminated due to unrecoverable numerical difficulties.";
	else if (num == 13) return "Unable to satisfy optimality tolerances; a sub-optimal solution is available.";
	else if (num == 14) return "An asynchronous optimization call was made, but the associated optimization run is not yet complete.";
	else if (num >= 15 || num <= 0) return "No specific error could be recognized.";
}

long smallestFeasibleLatency2Robust(KGraph &g)
{
	long s = g.DiameterUnweighted();
	vector<bool> allNodes(g.n, true);
	vector<long> dist;

	for (long i = 0; i < g.n; i++)
	{
		allNodes[i] = false;

		// find the smallest s such that V-v is a latency-s CDS
		for (long k = 0; k < g.n; k++)
		{
			dist = ComputeSSSPinGBv(g, allNodes, k);
			for (long j = 0; j < g.n; j++) if (dist[j] > s) s = dist[j];
		}
		allNodes[i] = true;
	}
	return s;
}

vector<long> solveLatencyCDS(KGraph &g, long r, long s)
{
	bool subOpt;
	vector<long> D;
	if (r == 1 && s > 2) D = solveMCDS(g, s);
	else cout << "ERROR: Your values of (r,s) are NOT supported." << endl;

	cerr << "Nodes in CDS are : ";
	for (long i = 0; i < D.size(); i++) cerr << D[i] << " ";
	cerr << endl;

	return D;
}


vector<long> SortVerticesByIncreasingDegree(KGraph &g) //Sorting vertices by increasing degree
{
	vector<long> emptyBucket;
	vector< vector<long> > degreeBuckets(g.n, emptyBucket);

	long deg;
	for (long i = 0; i < g.n; i++)
	{
		deg = g.degree[i];
		degreeBuckets[deg].push_back(i);
	}

	vector<long> sortedList(g.n, -1);
	long count = 0;
	long v;
	for (long i = 0; i < degreeBuckets.size(); i++)
	{
		for (long j = 0; j < degreeBuckets[i].size(); j++)
		{
			v = degreeBuckets[i][j];
			sortedList[count] = v;
			count++;
		}
	}

	return sortedList;
}

bool IsLatencyConstrainedCDS(KGraph &g, vector<bool> &D, long s)
{
	long ev;
	for (long v = 0; v < g.n; v++)
	{
		ev = EccentricitySubroutine(g, D, v);
		if (ev > s) return false;
	}
	return true;
}


vector<long> ComplementVector(vector<bool> &B, long n)
{
	vector<long> Q;					// complement of B (i.e., V\B), where V = {0, 1, ..., n-1 }
	for (long i = 0; i < n; i++)
	{
		if (!B[i]) Q.push_back(i);
	}
	return Q;
}
vector<long> Minimalize(KGraph &g, vector<long> &CUT, long s)
{
	vector<bool> B(g.n, true);
	long v;
	for (long i = 0; i < CUT.size(); i++)
	{
		v = CUT[i];
		B[v] = false;
	}
	return Minimalize(g, B, s);
}
vector<long> MinimalizeBasic(KGraph &g, vector<long> &CUT, long s)
{
	vector<bool> B(g.n, true);
	long v;
	for (long i = 0; i < CUT.size(); i++)
	{
		v = CUT[i];
		B[v] = false;
	}
	return MinimalizeBasic(g, B, s);
}
vector< vector<long> > EnumerateFarPairs(KGraph &g, vector<bool> &B, long s)
{
	vector< vector<long> > FarPairs;
	vector<long> dist;
	vector<long> Pair;

	for (long i = 0; i < g.n; i++)
	{
		dist = ComputeSSSPinGBv(g, B, i);
		for (long j = i + 1; j < g.n; j++)
		{
			if (dist[j] > s)
			{
				Pair.push_back(i);
				Pair.push_back(j);
				FarPairs.push_back(Pair);
				Pair.clear();
			}
		}
	}

	return FarPairs;
}

vector< vector<long> > EnumerateFarPairsRob(KGraph &g, vector<bool> &B, long s)
{
	vector< vector<long> > FarPairs;
	vector<long> dist;
	vector<long> Pair;

	for (long i = 0; i < g.n; i++)
	{
		dist = ComputeSSSPinGBv(g, B, i);
		for (long j = i + 1; j < g.n; j++)
		{
			if (dist[j] > s)
			{
				Pair.push_back(i);
				Pair.push_back(j);
				FarPairs.push_back(Pair);
				Pair.clear();
			}
			else
			{
				for (long k = 0; k < B.size(); k++)
				{
					if (!B[k]) continue;
					B[k] = false;
					dist = ComputeSSSPinGBv(g, B, i);
					if (dist[j] > s)
					{
						Pair.push_back(i);
						Pair.push_back(j);
						FarPairs.push_back(Pair);
						Pair.clear();
						break;
					}
				}
			}
		}

	}

	return FarPairs;
}


vector<long> Minimalize(KGraph &g, vector<bool> &B, long s)
{
	vector<long> Q = ComplementVector(B, g.n);					// complement of initial B (i.e., V\B)
	long q;
	long a, b;

	vector< vector<long> > FarPairs = EnumerateFarPairs(g, B, s);
	vector<long> dist;

	// make B maximal but not a latency-s CDS
	for (long i = 0; i < Q.size(); i++)
	{
		q = Q[i];

		// compute SSSP in graph G^B_q
		dist = ComputeSSSPinGBv(g, B, q);

		// check if FarPairs will become empty if we add q
		bool StillTooFar = false;
		for (long j = 0; j < FarPairs.size() && !StillTooFar; j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];

			if (dist[a] + dist[b] > s)
			{
				StillTooFar = true;
			}
		}

		if (StillTooFar)
		{
			// move q to set B
			B[q] = true;

			// create new set of FarPairs
			vector< vector<long> > newFarPairs;
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];

				if (dist[a] + dist[b] > s)
				{
					newFarPairs.push_back(FarPairs[j]);
				}
			}
			FarPairs = newFarPairs;
		}
	}

	// return complement of *final* B
	return ComplementVector(B, g.n);
}
vector<long> MinimalizeNone(KGraph &g, vector<bool> &B, long s)
{
	vector<long> Q = ComplementVector(B, g.n);
	return Q;
}

vector<long> MinimalizeNone(KGraph &g, vector<long> &CUT, long s)
{
	return CUT;
}

vector<long> MinimalizeBasic(KGraph &g, vector<bool> &B, long s)
{
	vector<long> Q = ComplementVector(B, g.n);		// complement of initial B (i.e., V\B)

	long q;
	// make B maximal but not a latency-s CDS
	for (long i = 0; i < Q.size(); i++)
	{
		q = Q[i];
		B[q] = true;
		if (IsLatencyConstrainedCDS(g, B, s))
		{
			B[q] = false;
		}
	}

	// return complement of *final* B
	return ComplementVector(B, g.n);
}

long EccentricitySubroutine(KGraph &g, vector<bool> &D, long v)
{
	vector<long> dist = ComputeSSSPinGBv(g, D, v);
	long ev = 0;
	for (long i = 0; i < g.n; i++)
	{
		if (dist[i] > ev) ev = dist[i];
	}
	return ev;

}
vector<long> ComputeSSSPinGBv(KGraph &g, vector<bool> &B, long v) //Computing single source shortest path
{
	vector<long> dist(g.n, g.n);
	vector<long> parents;
	vector<long> children;

	dist[v] = 0;
	long u, w;

	children.push_back(v);

	while (!children.empty())
	{
		parents = children;
		children.clear();
		for (long i = 0; i < parents.size(); i++)
		{
			w = parents[i];
			for (long j = 0; j < g.degree[w]; j++) {
				u = g.adj[w][j];
				if (dist[u] == g.n) {
					dist[u] = dist[w] + 1;
					if (B[u])
						children.push_back(u);
				}
			}
		}
	}
	return dist;
}

vector<long> ComputeSSSPinGBvWeighted(KGraph &g, vector<bool> &B, long origin, vector<long> W) //Computing weighted single source shortest path
{
	long infty = long(floor((double)0.4*LONG_MAX));
	vector<long> dist(g.n, infty); //shortest distance from origin node to each other node. dist[i] = \infty means i not reachable
	vector<bool> Q(g.n, true); // the set of vertices whose distance labels are not yet permanent
	//W[origin] = 0;   // for convenience, since we are not counting the weight of the first vertex in the path as part of its length
	dist[origin] = 0; //the origin node is distance 0 from itself
	long minDistance;
	long minVertex;

	// do node-weighted version of Dijkstra's shortest path algorithm.
	for (long i = 0; i < g.n; i++)
	{
		// find a vertex u from Q of minimum distance label dist[u]
		minDistance = LONG_MAX;
		minVertex = -1;
		for (long v = 0; v < g.n; v++)
		{
			if (!Q[v]) continue;
			if (dist[v] < minDistance)
			{
				minDistance = dist[v];
				minVertex = v;
			}
		}

		// remove minVertex from Q
		Q[minVertex] = false;

		if (B[minVertex] || minVertex == origin)  // only allow communication paths to emanate from the origin or from nodes in the CDS B
		{
			// update distance labels of neighbors
			for (long j = 0; j < g.degree[minVertex]; j++)
			{
				long v = g.adj[minVertex][j];
				if (Q[v] && dist[minVertex] + W[minVertex] < dist[v])
				{
					dist[v] = dist[minVertex] + W[minVertex];
				}
			}
		}
	}
	return dist;
}

vector<bool> HeuristicLCDS(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, true);
	return HeuristicLCDS(g, heuristicSoln, s);
}
vector<bool> HeuristicLCDS(KGraph &g, vector<bool> &heuristicSoln, long s)
{
	vector<long> sortedList = SortVerticesByIncreasingDegree(g);
	long v;

	for (int i = 0; i < g.n; i++) {
		v = sortedList[i];
		if (!heuristicSoln[v]) continue;
		heuristicSoln[v] = false;
		if (!IsLatencyConstrainedCDS(g, heuristicSoln, s))
		{
			heuristicSoln[v] = true;
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}


vector<bool> HeuristicLCDS2(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, true);
	vector<long> sortedList = SortVerticesByIncreasingDegree(g);
	long v;

	for (int i = 0; i < g.n; i++)
	{
		v = sortedList[i];
		heuristicSoln[v] = false;
		if (!IsLatencyConstrainedCDS(g, heuristicSoln, s))
		{
			heuristicSoln[v] = true;
		}
		else
		{
			for (long j = 0; j < g.n && !heuristicSoln[v]; j++)
			{
				if (j != v && heuristicSoln[j])
				{
					heuristicSoln[j] = false;
					if (!IsLatencyConstrainedCDS(g, heuristicSoln, s)) heuristicSoln[v] = true;
					heuristicSoln[j] = true;
				}
			}
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

bool IsLatencyConstrainedRCDS(KGraph &g, vector<bool> D, long s)
{
	for (long i = 0; i < g.n; i++)
	{
		if (!D[i]) continue;
		D[i] = false;
		if (!IsLatencyConstrainedCDS(g, D, s))
		{
			return false;
		}
		D[i] = true;
	}
	return true;
}

vector<bool> MinimalizeRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s)
{
	for (long i = 0; i < g.n; i++)
	{
		if (!heuristicSoln[i]) continue;
		// now node i belongs to the heuristic solution for r=2. Can we remove it and still be feasibe?
		heuristicSoln[i] = false;
		if (!IsLatencyConstrainedRCDS(g, heuristicSoln, s)) heuristicSoln[i] = true;
		// else we didn't need node i in the solution, and we can keep it set to false.
	}
	return heuristicSoln;
}

vector<bool> HeuristicRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s)
{
	vector<bool> initialHeuristicSoln = heuristicSoln;
	vector< vector<long> > FarPairs;
	vector< vector<long> > newFarPairs;
	vector<long> dist;
	vector<long> newPair;
	vector<long> score;
	long a, b;
	long max;
	long imax;
	vector<long> scoreZero(g.n, 0);

	for (long v = 0; v < g.n; v++)
	{
		if (!initialHeuristicSoln[v]) continue;  //if node was not in initial heuristic soln, then removing it will not make current heuristic soln infeasibe
		heuristicSoln[v] = false;
		FarPairs = EnumerateFarPairs(g, heuristicSoln, s);
		//cerr << "v = " << v << endl;
		while (!FarPairs.empty())
		{
			score = scoreZero;
			score[v] = -1;
			max = -1;
			imax = -1;
			for (long i = 0; i < g.n; i++)
			{
				if (heuristicSoln[i] || i == v) continue;
				dist = ComputeSSSPinGBv(g, heuristicSoln, i);
				for (long j = 0; j < FarPairs.size(); j++)
				{
					a = FarPairs[j][0];
					b = FarPairs[j][1];
					if (dist[a] + dist[b] <= s)
					{
						score[i]++;
					}

				}

				if (score[i] > max) {
					max = score[i];
					imax = i;
				}
			}
			heuristicSoln[imax] = true;
			newFarPairs.clear();

			dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] > s)
				{
					newPair.clear();
					newPair.push_back(a);
					newPair.push_back(b);
					newFarPairs.push_back(newPair);
				}
			}

			FarPairs = newFarPairs;

		}
		heuristicSoln[v] = true;
	}

	return heuristicSoln;
}


vector<bool> HeuristicRLCDS2(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, false);
	vector< vector<long> > FarPairs = EnumerateFarPairs(g, heuristicSoln, 1);
	vector<long> dist;
	vector<long> score;
	vector<long> scoreZero(g.n, 0);
	vector< vector<long> > newFarPairs;
	vector<long> newPair;

	long a, b;
	long max;
	long imax;
	while (!FarPairs.empty())
	{
		score = scoreZero;
		max = -1;
		imax = -1;
		for (long i = 0; i < g.n; i++)
		{
			if (heuristicSoln[i]) continue;
			dist = ComputeSSSPinGBv(g, heuristicSoln, i);
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] <= s)
				{
					score[i]++;
				}

			}

			if (score[i] > max) {
				max = score[i];
				imax = i;
			}
		}
		heuristicSoln[imax] = true;

		dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
		for (long j = 0; j < FarPairs.size(); j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];
			if (dist[a] + dist[b] > s)
			{
				newPair.clear();
				newPair.push_back(a);
				newPair.push_back(b);
				newFarPairs.push_back(newPair);
			}
			else
			{
				for (long k = 0; k < heuristicSoln.size(); k++)
				{
					if (k == imax) continue;
					if (!heuristicSoln[k]) continue;
					heuristicSoln[k] = false;
					dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
					heuristicSoln[k] = true;
					if (dist[a] + dist[b] > s)
					{
						newPair.clear();
						newPair.push_back(a);
						newPair.push_back(b);
						newFarPairs.push_back(newPair);
						break;
					}
				}
			}
		}

		FarPairs = newFarPairs;

	}

	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<bool> HeuristicLCDSBestIn(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, false);
	vector< vector<long> > FarPairs = EnumerateFarPairs(g, heuristicSoln, 1);
	vector<long> dist;
	vector<long> score;
	vector<long> scoreZero(g.n, 0);

	long a, b;
	long max;
	long imax;
	while (!FarPairs.empty())
	{
		score = scoreZero;
		max = -1;
		imax = -1;
		for (long i = 0; i < g.n; i++)
		{
			if (heuristicSoln[i]) continue;
			dist = ComputeSSSPinGBv(g, heuristicSoln, i);
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] <= s)
				{
					score[i]++;
				}

			}

			if (score[i] > max) {
				max = score[i];
				imax = i;
			}
		}
		heuristicSoln[imax] = true;

		vector< vector<long> > newFarPairs;
		vector<long> newPair;
		dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
		for (long j = 0; j < FarPairs.size(); j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];
			if (dist[a] + dist[b] > s)
			{
				newPair.clear();
				newPair.push_back(a);
				newPair.push_back(b);
				newFarPairs.push_back(newPair);
			}
		}

		FarPairs = newFarPairs;

	}

	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<long> BiasedShuffle(KGraph &g)
{
	vector<long> order(g.n, -1);
	long rn;
	long remdeg = 2 * g.m; // remaining degree
	long tempdeg;
	vector<bool> selected(g.n, false);  // whether or not the vertex is in order yet.
	bool numFound;

	for (long i = 0; i < g.n; i++)	// randomly pick a vertex that has not been selected, biased by degree
	{
		rn = rand() % remdeg;

		numFound = false;
		tempdeg = 0;
		for (long j = 0; j < g.n && !numFound; j++)
		{
			if (selected[j]) continue;
			tempdeg += g.degree[j];
			if (rn < tempdeg)
			{
				order[g.n - i - 1] = j;
				selected[j] = true;
				remdeg -= g.degree[j];
				numFound = true;
			}
		}
	}
	return order;
}
vector<bool> RandomHeuristicLCDS(KGraph &g, long s, long iterations)
{
	vector<bool> incumbent;
	long incsize = g.n + 1;
	long newsize;
	vector<long> sortedList;

	for (long iter = 0; iter < iterations; iter++)
	{
		vector<bool> heuristicSoln(g.n, true);
		sortedList = BiasedShuffle(g);
		long v;

		for (int i = 0; i < g.n; i++)
		{
			v = sortedList[i];
			heuristicSoln[v] = false;
			if (!IsLatencyConstrainedCDS(g, heuristicSoln, s))
			{
				heuristicSoln[v] = true;
			}
		}
		newsize = count(heuristicSoln.begin(), heuristicSoln.end(), true);
		if (newsize < incsize)
		{
			incumbent = heuristicSoln;
			incsize = newsize;
		}
	}

	cerr << "Heuristic solution found,  number of nodes = " << count(incumbent.begin(), incumbent.end(), true) << endl;
	return incumbent;
}


vector< vector<long> > FindBiconnectedComponents(KGraph &g, vector<bool> &AV)
{
	/* I tried to use the naming conventions presented in Tarjan's 1972 paper.
	I assume that the graph is connected, so that only one call to the recursion is necessary. */

	// declare vars
	long u = -1, v = 0, i = 0;
	vector<long> number(g.n, (long)-1);
	vector<long> lowopt(g.n, g.n);
	vector< vector<long> > BC;		// biconnected components
	stack<long> le, re;				// used to store a stack of edges. le is "left edge" and re is "right edge". An edge is stored (le[i],re[i]). 

									// perform DFS-based algorithm
	Bico_Sub(v, u, i, g, number, lowopt, le, re, BC);

	vector<long> countComp(g.n, 0);
	for (long p = 0; p < BC.size(); p++) // count how many components each vertex belongs to
		for (long q = 0; q < BC[p].size(); q++)
			countComp[BC[p][q]]++;

	vector<bool> AV_temp(g.n, false);
	AV = AV_temp;
	for (long p = 0; p < g.n; p++) // if a vertex belongs to >1 component, then it is a cut vertex
		if (countComp[p] > 1)
			AV[p] = true;

	return BC;
}

void Bico_Sub(long v, long u, long &i, KGraph &g, vector<long> &number, vector<long> &lowopt, stack<long> &le, stack<long> &re, vector< vector<long> > &BC)
{
	i++;
	number[v] = i;
	lowopt[v] = number[v];
	long w;
	for (long j = 0; j < g.degree[v]; j++)
	{
		w = g.adj[v][j];
		if (number[w] == -1)
		{
			le.push(v);
			re.push(w);
			Bico_Sub(w, v, i, g, number, lowopt, le, re, BC);
			lowopt[v] = (long)min(lowopt[v], lowopt[w]);
			if (lowopt[w] >= number[v])
			{
				vector<long> temp_BC;
				vector<bool> bBC(g.n, false);
				while (!le.empty() && !re.empty() && number[le.top()] >= number[w])
				{
					bBC[le.top()] = true;
					bBC[re.top()] = true;
					le.pop();
					re.pop();
				}
				if (!le.empty() && le.top() == v)
				{
					bBC[le.top()] = true;
					bBC[re.top()] = true;
					le.pop();
					re.pop();
				}
				else
				{
					cerr << "ERROR: edge (v,w) not on top of stack" << endl;
				}
				for (long p = 0; p < g.n; p++)
					if (bBC[p])
						temp_BC.push_back(p);
				BC.push_back(temp_BC);
			}
		}
		else if (number[w] < number[v] && w != u)
		{
			le.push(v);
			re.push(w);
			lowopt[v] = min(lowopt[v], number[w]);
		}
	}
}

int myrandom(int i)
{
	return std::rand() % i;
}

vector<long> solveMCDS(KGraph &g, long s)
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}
	try {
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);

		model.update();

		//cerr << "Adding objective function" << endl;
		for (int i = 0; i < g.n; i++)
		{
			X[i].set(GRB_DoubleAttr_Obj, 1);
		}
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();

		//cerr << "Adding domination constraints" << endl;
		long v;
		vector<long> minimalCut;
		string MinimalizeSetting = "Basic";
		for (int i = 0; i < g.n; i++) // for each vertex i, make sure you choose at least one vertex from N(i). This is a vertex cut!
		{
			if (MinimalizeSetting == "None")
				minimalCut = MinimalizeNone(g, g.adj[i], s);
			else if (MinimalizeSetting == "Basic")
				minimalCut = MinimalizeBasic(g, g.adj[i], s);
			else if (MinimalizeSetting == "Fast")
				minimalCut = Minimalize(g, g.adj[i], s);   // a minimal subset of N(i) that is a length-s cut.
			else cerr << "ERROR: not a supported value for MinimalizeSetting." << endl;
			GRBLinExpr expr = 0;
			for (long j = 0; j < minimalCut.size(); j++)
			{
				v = minimalCut[j];
				expr += X[v];
			}
			model.addConstr(expr >= 1);
		}
		model.update();

		//Initial Solution
		double TotaltimeinHeuristic = 0;
		time_t start1 = clock();
		vector<bool> soln2 = HeuristicLCDSBestIn(g, s);
		//	cerr << "Best-in heuristic soln = " << count(soln2.begin(), soln2.end(), true) << endl;
		soln2 = HeuristicLCDS(g, soln2, s);
		//	cerr << "Best-in heuristic soln (after minimalizing) = " << count(soln2.begin(), soln2.end(), true) << endl;
		vector<bool> initialSolution = soln2;
		TotaltimeinHeuristic = (double)(clock() - start1) / CLOCKS_PER_SEC;

		for (long i = 0; i < g.n; i++)
		{
			if (initialSolution[i])
			{
				X[i].set(GRB_DoubleAttr_Start, 1.0);
			}
			else
			{
				X[i].set(GRB_DoubleAttr_Start, 0.0);
			}
		}

		// fix articulation vertices in solution
		vector<bool> ArticulationVertices; // (g.n, false);
		FindBiconnectedComponents(g, ArticulationVertices);

		for (int i = 0; i < g.n; i++)
		{
			if (ArticulationVertices[i])
			{
				X[i].set(GRB_DoubleAttr_LB, 1.0);
			}
		}

		//	cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		model.update();

		//End of initial solution and articulation vertex cut

	//	cerr << "Adding lazy constraints (the vertex cut inequalities)\n";

		LazyConstraints cb = LazyConstraints(X, g, s);	// tell Gurobi which function generates the lazy cuts.
		model.setCallback(&cb);

		//	cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;

		cerr << endl;
		cerr << "Number of callbacks : " << LazyConstraints::numCallbacks << endl;
		cerr << "Time in callbacks : " << LazyConstraints::TotalCallbackTime << endl;
		cerr << "Number of lazy cuts : " << LazyConstraints::numLazyCuts << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;

		for (int i = 0; i < g.n; i++)
		{
			if (X[i].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				cds.push_back(i);
			}
		}

		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << count(initialSolution.begin(), initialSolution.end(), true) << " " << cds.size() << " " << lb << " " << NumBBNodes << " " << LazyConstraints::numLazyCuts << " " << TotaltimeinHeuristic << " " << LazyConstraints::TotalCallbackTime << " ";

		delete[] X;


	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return cds;
}
