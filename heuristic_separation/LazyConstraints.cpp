#include "LazyConstraints.h"
#include "gurobi_c++.h"
#include "KGraph.h"
#include "GRBInterface.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <queue>
#include <set>
#include <map>
#include <omp.h>
#include <memory>
#include <sstream>
#include <vector>
using namespace std;

long LazyConstraints::numCallbacks = 0;
double LazyConstraints::TotalCallbackTime = 0;
long LazyConstraints::numLazyCuts = 0;

long LazyConstraints2::numCallbacks = 0;
double LazyConstraints2::TotalCallbackTime = 0;
long LazyConstraints2::numLazyCuts = 0;

long LazyConstraints3::numCallbacks = 0;
double LazyConstraints3::TotalCallbackTime = 0;
long LazyConstraints3::numLazyCuts = 0;

long LazyConstraints4::numCallbacks = 0;
double LazyConstraints4::TotalCallbackTime = 0;
long LazyConstraints4::numLazyCuts = 0;

long LazyConstraints5::numCallbacks = 0;
double LazyConstraints5::TotalCallbackTime = 0;
long LazyConstraints5::numLazyCuts = 0;

long LazyConstraints5::numCallbacks2 = 0;
double LazyConstraints5::TotalCallbackTime2 = 0;
long LazyConstraints5::numLazyCuts2 = 0;


GRBVar *vars;
KGraph g1;
vector<long> W1;
long s1;

LazyConstraints3::LazyConstraints3(GRBVar * xvars, KGraph & g, long &s, vector<long> W)
{

	vars = xvars;
	g1.Duplicate(g);
	s1 = s;
	W1 = W;
};

void LazyConstraints3::callback() {
	try {
		if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			numCallbacks++;
			time_t start = clock();

			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			vector<bool> selectedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)			// convert the ``solution'' to a boolean vector
			{
				if (x[i] > 0.5)
				{
					selectedVertices[i] = true;
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


LazyConstraints4::LazyConstraints4(GRBVar * xvars, KGraph & g, long &s)
{
	vars = xvars;
	g1.Duplicate(g);
	s1 = s;
};

void LazyConstraints4::callback() {
	try {
		if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			numCallbacks++;
			time_t start = clock();

			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			vector<bool> selectedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)			// convert the ``solution'' to a boolean vector
			{
				if (x[i] > 0.5)
				{
					selectedVertices[i] = true;
				}
			}

			if (!IsLatencyConstrainedCDS(g1, selectedVertices, s1)) // If the ``solution'' is not really a CDS, then add a new vertex cut constraint to the problem. (Note: it is enough to see if it is connected, since our initial constraints enforce that it is dominating.
			{
				vector<long> minimalCut = Minimalize(g1, selectedVertices, s1);

				long v;
				GRBLinExpr expr = 0;
				for (long i = 0; i < minimalCut.size(); i++)
				{
					v = minimalCut[i];
					expr += vars[v];
				}
				addLazy(expr >= 1);
				numLazyCuts++;
			}
			else
			{
				bool cutAdded = false;
				for (long i = 0; i < g1.n && !cutAdded; i++)
				{
					if (selectedVertices[i]) {
						selectedVertices[i] = false;
						if (!IsLatencyConstrainedCDS(g1, selectedVertices, g1.n - 1)) {
							vector<bool> tempSelectedVertices = selectedVertices;
							vector<long> minimalCut = Minimalize(g1, tempSelectedVertices, g1.n - 1);

							long v;
							GRBLinExpr expr = 0;
							for (long i = 0; i < minimalCut.size(); i++)
							{
								v = minimalCut[i];
								expr += vars[v];
							}
							addLazy(expr >= 2);
							numLazyCuts++;
							cutAdded = true;
						}
						selectedVertices[i] = true;
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



LazyConstraints5::LazyConstraints5(GRBVar * xvars, KGraph & g, long &s)
{
	vars = xvars;
	g1.Duplicate(g);
	s1 = s;
};

void LazyConstraints5::callback() {
	try {
		if (where == GRB_CB_MIPNODE && GRB_CB_MIPNODE_STATUS == GRB_OPTIMAL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			numCallbacks++;
			long MaxNumCuts = 10;
			long NumCutsAdded = 0;
			time_t start = clock();
			double u;
			double epsilon = 0.01;
			double *x = new double[g1.n];
			x = getNodeRel(vars, g1.n);

			for (long i = 0; i < g1.n && NumCutsAdded < MaxNumCuts; i++)
			{
				vector<bool> neighbors(g1.n, false);
				long v;
				for (long j = 0; j < g1.degree[i]; j++)
				{
					v = g1.adj[i][j];
					neighbors[v] = true;
				}

				for (long j = i + 1; j < g1.n && NumCutsAdded < MaxNumCuts; j++)
				{
					if (neighbors[j]) continue;

					vector<long> p = g1.CommonNeighborsList(i, j);

					for (long k = 0; k < p.size(); k++)
					{
						u += x[p[k]];
					}

					if (u < 1 - epsilon)
					{
						GRBLinExpr expr = 0;
						for (long k = 0; k < p.size(); k++)
						{
							expr += vars[p[k]];
						}
						addLazy(expr >= 1);
						numLazyCuts++;
						NumCutsAdded++;
					}

				}
			}
			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] x;
		}
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks2++;
			time_t start2 = clock();

			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			vector<bool> selectedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)			// convert the ``solution'' to a boolean vector
			{
				if (x[i] > 0.5)
				{
					selectedVertices[i] = true;
				}
			}
			vector< vector<long> > B = FindBiconnectedComponents(g1, selectedVertices);
			if (B.empty())
			{
				vector<long> minimalCut = Minimalize(g1, selectedVertices, s1);
				long v;
				GRBLinExpr expr = 0;
				for (long i = 0; i < minimalCut.size(); i++)
				{
					v = minimalCut[i];
					expr += vars[v];
				}
				addLazy(expr >= 2);
				numLazyCuts2++;
			}
			TotalCallbackTime2 += (double)(clock() - start2) / CLOCKS_PER_SEC;
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

LazyConstraints2::LazyConstraints2(GRBVar * xvars, KGraph & g, long &s)
{
	vars = xvars;
	g1.Duplicate(g);
	s1 = s;
};



void LazyConstraints2::callback() {
	try {
		if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			numCallbacks++;
			time_t start = clock();

			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			vector<bool> selectedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)			// convert the ``solution'' to a boolean vector
			{
				if (x[i] > 0.5)
				{
					selectedVertices[i] = true;
				}
			}

			if (!IsLatencyConstrainedCDS(g1, selectedVertices, s1)) // If the ``solution'' is not really a CDS, then add a new vertex cut constraint to the problem. (Note: it is enough to see if it is connected, since our initial constraints enforce that it is dominating.
			{
				vector<long> minimalCut = Minimalize(g1, selectedVertices, s1);

				long v;
				GRBLinExpr expr = 0;
				for (long i = 0; i < minimalCut.size(); i++)
				{
					v = minimalCut[i];
					expr += vars[v];
				}
				addLazy(expr >= 2);
				numLazyCuts++;
			}
			else
			{
				bool cutAdded = false;
				for (long i = 0; i < g1.n && !cutAdded; i++)
				{
					if (selectedVertices[i]) {
						selectedVertices[i] = false;
						if (!IsLatencyConstrainedCDS(g1, selectedVertices, s1)) {
							vector<bool> tempSelectedVertices = selectedVertices;
							vector<long> minimalCut = Minimalize(g1, tempSelectedVertices, s1);

							long v;
							GRBLinExpr expr = 0;
							for (long i = 0; i < minimalCut.size(); i++)
							{
								v = minimalCut[i];
								expr += vars[v];
							}
							addLazy(expr >= 2);
							numLazyCuts++;
							cutAdded = true;
						}
						selectedVertices[i] = true;
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

LazyConstraints::LazyConstraints(GRBVar * xvars, KGraph & g, long &s)
{
	vars = xvars;
	g1.Duplicate(g);
	s1 = s;
};

void LazyConstraints::callback() {
	try {
		if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			numCallbacks++;
			time_t start = clock();

			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			vector<bool> selectedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)			// convert the ``solution'' to a boolean vector
			{
				if (x[i] > 0.5)
				{
					selectedVertices[i] = true;
				}
			}

			if (!IsLatencyConstrainedCDS(g1, selectedVertices, s1)) // If the ``solution'' is not really a CDS, then add a new vertex cut constraint to the problem. (Note: it is enough to see if it is connected, since our initial constraints enforce that it is dominating.
			{
				vector<long> minimalCut;
				string MinimalizeSetting = "Basic";
				if (MinimalizeSetting == "None")
					minimalCut = MinimalizeNone(g1, selectedVertices, s1);
				else if (MinimalizeSetting == "Basic")
					minimalCut = MinimalizeBasic(g1, selectedVertices, s1);
				else if (MinimalizeSetting == "Fast")
					minimalCut = Minimalize(g1, selectedVertices, s1);   // a minimal subset of N(i) that is a length-s cut.
				else cerr << "ERROR: not a supported value for MinimalizeSetting." << endl;
				long v;
				GRBLinExpr expr = 0;
				for (long i = 0; i < minimalCut.size(); i++)
				{
					v = minimalCut[i];
					expr += vars[v];
				}
				addLazy(expr >= 1);
				numLazyCuts++;
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
