#include "gurobi_c++.h"
#include "GRBInterface.h"
#include "KGraph.h"
#include "LazyConstraints.h"
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include<iostream>
#include<algorithm>

using namespace std;
using namespace std::chrono;


long s = 3;     //the parameter s in s-club

// global variables to count the number of callbacks and lazy cuts
long num_callbacks_interdiction;
long num_Lazycuts_interdiction;

// total time to solve the maximum s-club problem and minimum LCDS in all the iterations
double SclubTime = 0;
double Finding_LCDS_time = 0;


class addlazycut_theta : public GRBCallback
{
public:

	GRBVar* Xvar;
	GRBVar Theta;
	long n_Xvar;
	KGraph graph_1;

	addlazycut_theta(GRBVar* xvar, GRBVar Teta, long N_xvar, KGraph &graph_2)
	{
		Xvar = xvar;
		Theta = Teta;
		n_Xvar = N_xvar;
		graph_1 = graph_2;
	}

protected:

	void callback() {
		try
		{
			if (where == GRB_CB_MIPSOL)
			{
				num_callbacks_interdiction++;    // count the number of callbacks

				// get the solution of the master problem
				double *x_master = new double[n_Xvar];
				x_master = getSolution(Xvar, n_Xvar);
				double THETA = getSolution(Theta);

				// if value of variable x is less than or equal to 0.5, it is not interdicted so x = 0
				vector <long> non_interdicted_vertices;
				for (int p2 = 0; p2 < n_Xvar; p2++)
				{
					if (x_master[p2] <= 0.5)
					{
						non_interdicted_vertices.push_back(p2);
					}
				}

				// mapping to find the adjacency list of the graph after interdiction		
				vector<long> ReverseMap;
				KGraph induced_g = graph_1.CreateInducedGraph(non_interdicted_vertices, ReverseMap);

				// solve the separation only if the interdicted graph has at least 2 vertices
				if (induced_g.n >= 2)
				{

					auto start_sclub = chrono::steady_clock::now();  // begin to compute the time for solving the separtion problem
					
					vector <long> sclb_index;
					vector <long> HS;
					HS = HeuristicAndPreprocess(induced_g, s);   // find the heuristic s-club in the interdicted graph
					sclb_index = ICUT(induced_g, s, HS);         // find maximum s-club in the interdicted graph					
					
					chrono::duration <double> duration_sclb = chrono::steady_clock::now() - start_sclub;   // duration of solving the separation
					SclubTime += duration_sclb.count();

					// add a cut if the solution of master problem is not feasible, i.e., theta < |s-club|
					if (THETA < sclb_index.size())
					{
						// converting the index of vertices to the indices in the original graph (before interdiction)
						vector <long> sclb_original_index;
						for (long i = 0; i < sclb_index.size(); i++)
						{
							sclb_original_index.push_back(non_interdicted_vertices[sclb_index[i]]);
						}
						KGraph induced_kclb = graph_1.CreateInducedGraph(sclb_original_index, ReverseMap);

						// define requrired structures to add the lazy cut	
						GRBLinExpr Critical_vertices = 0;
						GRBLinExpr Hereditary_vertices = 0;
						vector<long> Critical_set;
						vector<long> Hereditary_set;
						vector<long> Critical_set_original;
						
						//find critical set
						auto start_lcds = chrono::steady_clock::now();
						long r = 1;
						Critical_set = solveMCDS(induced_kclb, s);
						chrono::duration <double> duration_lcds = chrono::steady_clock::now() - start_lcds;
						Finding_LCDS_time += duration_lcds.count();

						// converting the index of vertices in the critical set to the indices in the original graph (before interdiction)
						for (long i = 0; i < Critical_set.size(); i++)
						{
							Critical_set_original.push_back(sclb_original_index[Critical_set[i]]);
						}

						//find hereditary set
						long key;
						long c = 0;
						for (long i = 0; i < sclb_original_index.size(); i++)
						{
							key = sclb_original_index[i];
							c = std::count(Critical_set_original.begin(), Critical_set_original.end(), key);
							if (c == 0)
							{
								Hereditary_set.push_back(key);
							}
							c = 0;
						}

						// add the lazy cut
						for (long i = 0; i < Critical_set_original.size(); i++)
						{
							Critical_vertices += Xvar[Critical_set_original[i]];
						}

						for (long i = 0; i < Hereditary_set.size(); i++)
						{
							Hereditary_vertices += Xvar[Hereditary_set[i]];
						}
						addLazy(Theta >= ((1 - Critical_vertices) * sclb_index.size()) - Hereditary_vertices);
						num_Lazycuts_interdiction++;
						Critical_vertices = 0;
						Hereditary_vertices = 0;
						Critical_set.clear();
						Critical_set_original.clear();
						Hereditary_set.clear();
					}
				}				
				delete[] x_master;
			}
		}  
		catch (GRBException e)
		{
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...)
		{
			cout << "Error during callback" << endl;
		}
	}     
}; 



int main(int argc, char *argv[])
{
	auto start = chrono::steady_clock::now();
	if (argc < 2)
		cerr << "ERROR: Not enough arguments.";

	// read the input graph, format of the dataset and the value of the penalty (alpha)
	KGraph grph(argv[2], argv[2], argv[1]);
	double alpha = atof(argv[3]);
	long num_interdicted_vertices = 0;
	
	// Master (interdiction) Problem 
	try {
		GRBEnv env_master = GRBEnv();
		GRBModel model_Master = GRBModel(env_master);

		// veriables
		GRBVar* X_Master = model_Master.addVars(grph.n, GRB_BINARY);
		GRBVar theta = model_Master.addVar(0.0, INFINITY, 1.0, GRB_CONTINUOUS);
		model_Master.update();

		//set obj coefficient
		for (long p1 = 0; p1 < grph.n; ++p1)
		{
			X_Master[p1].set(GRB_DoubleAttr_Obj, alpha);
		}
		
		// add constraints for edges and stars
		long counter = 0;
		GRBLinExpr neighbors_u = GRBLinExpr();
		GRBLinExpr neighbors_v = GRBLinExpr();
		GRBLinExpr intersection_of_uv = GRBLinExpr();
		vector<long> common_vertices;
		long v;
		long edge_size;
		vector <bool> star_cut(grph.n, false);
		for (long u = 0; u < grph.n; u++)
		{
			neighbors_u = 0;
			for (long i = 0; i < grph.degree[u]; i++)
			{
				neighbors_u += X_Master[grph.adj[u][i]];
			}
			if (star_cut[u] == false)
			{
				model_Master.addConstr(theta >= ((1 - X_Master[u]) * (grph.degree[u] + 1)) - neighbors_u);
				star_cut[u] = true;
			}
			for (long j = 0; j < grph.degree[u]; j++)
			{
				v = grph.adj[u][j];
				if (u < v)
				{
					for (long i = 0; i < grph.degree[v]; i++)
					{
						neighbors_v += X_Master[grph.adj[v][i]];
					}
					if (star_cut[v] == false)
					{
						model_Master.addConstr(theta >= ((1 - X_Master[v]) * (grph.degree[v] + 1)) - neighbors_v);
						star_cut[v] = true;
					}
					common_vertices = grph.CommonNeighborsList(u, v);
					for (long i = 0; i < common_vertices.size(); i++)
					{
						intersection_of_uv += X_Master[common_vertices[i]];
					}
					edge_size = grph.degree[u] + grph.degree[v] - common_vertices.size();
					model_Master.addConstr(theta >= ((1 - X_Master[u] - X_Master[v]) * edge_size) - neighbors_u - neighbors_v  + intersection_of_uv);
					counter++;
					neighbors_v = 0;
					intersection_of_uv = 0;
					common_vertices.clear();
				}
			}
		}

		//SET GUROBI PARAMETERS

		//Specify the use of lazy constraints
		model_Master.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		//Set feasibility vs optimality balance
		model_Master.getEnv().set(GRB_IntParam_MIPFocus, 0);
		//1-feasible sols quickly;2-prove optimality;3-focus on MIP bound; default is 0

		//Set threads; review guidance on Gurobi.com; 0-default;
		model_Master.getEnv().set(GRB_IntParam_Threads, 0);

		//Set root node LPR solver
		model_Master.getEnv().set(GRB_IntParam_Method, -1);
		//-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent

		//Set BC node LPR solver
		model_Master.getEnv().set(GRB_IntParam_NodeMethod, 1);
		//0=primal simplex, 1=dual simplex, 2=barrier

		//Set global cut aggressiveness; over-ridden by individual cut settings
		model_Master.getEnv().set(GRB_IntParam_Cuts, -1);
		//0=no cuts;1=moderate;2=aggressive;3=very aggressive;-1=default

		//Set maximum time limit
		model_Master.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

		//Set termination gap limit; as needed; default is 1e-4
		model_Master.getEnv().set(GRB_DoubleParam_MIPGap, 1e-4);

		//Set Gurobi screen display flag
		model_Master.getEnv().set(GRB_IntParam_OutputFlag, 1);
		//0=switch off; 1=default

		// Save log file
		model_Master.set(GRB_StringParam_LogFile, "logfile");

		//Set objective to minimize
		model_Master.set(GRB_IntAttr_ModelSense, 1);

		model_Master.update();


		// set callback
		addlazycut_theta cb1 = addlazycut_theta(X_Master, theta, grph.n, grph);
		model_Master.setCallback(&cb1);
		model_Master.optimize();

		if (model_Master.get(GRB_IntAttr_SolCount) == 0)
		{
			cout << "No solution found, Gurobi optimization status = " << model_Master.get(GRB_IntAttr_Status) << endl;
		}
		else
		{
			double obj_master = model_Master.get(GRB_DoubleAttr_ObjVal);   // objective value of the master problem
			cout << endl << "obj = " << obj_master << endl;
			for (long i = 0; i < grph.n; i++)
			{
				if (X_Master[i].get(GRB_DoubleAttr_X) > 0.5)
				{
					num_interdicted_vertices++;
				}
			}		
			cout << "num_interdicted_vertices : " << num_interdicted_vertices << endl;
			cout << "theta : " << theta.get(GRB_DoubleAttr_X) << endl;
			cout << "# B&B nodes in interdiction = " << (long)model_Master.get(GRB_DoubleAttr_NodeCount) << endl;
			cout << "# of callbacks in interdiction  = " << num_callbacks_interdiction << endl;
			cout << "# of lazy cuts in interdiction = " << num_Lazycuts_interdiction << endl;
		}
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

	// Total Time 
	chrono::duration <double> duration = chrono::steady_clock::now() - start;
	printf("Total Time : %.2fs\n", duration.count());

	//print time for max sclub problem
	printf("sclb Time : %.2fs\n", SclubTime);

	//print time for LCDS
	printf("Finding LCDS Time : %.2fs\n", Finding_LCDS_time);

	return 0;
}
