#ifndef KGRAPH_H
#define KGRAPH_H
//#include <ilcplex/cplex.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdio.h>
#include <string.h>


extern double copy_time;
using namespace std;

class KGraph
{
public:
	bool IsConnected(vector<long> S);
	bool IsConnected(vector<bool> S);
	bool IsConnected();

	long DiameterUnweighted();
	long DiameterUnweighted(vector<long> S);
	long DiameterUnweighted(vector<bool> S1);
	long LongestShortestPathUnweighted(long origin);
	vector<long> ShortestPathsUnweighted(long origin);
	vector<long> ShortestPathsUnweighted(long origin, vector<bool> &S);
	vector<long> MultiSourceShortestPaths(vector<long> &D);
	vector<long> MultiSourceShortestPaths(vector<long> &D, vector<bool> &S);

	bool DeleteNode(long i);

	vector<long> *adj;  // stores the adj lists. Each list is maintained as a sorted vector.
	long *degree;       // stores the degree seq
	long n;             // num of nodes
	long m;             // num of edges
	long Delta;         // higest degree. As of now, is instanciated when graph is read or copied, not maintained afterwards
	string name;        // name of the graph. Could be anything.

						/* Constructors */
	KGraph();
	KGraph(string nm);
	KGraph(long n);
	KGraph(string nm, string file, string type);
	KGraph(const KGraph &rhs);

	/* File IO utility functions */
	void ReadDIMACSGraph(string file);
	void ReadDIMACSColorGraph(string file);
	void ReadDIMACSGraphParallel(string file);
	void ReadDATGraph(string file);
	void ReadSNAPGraph(string file);
	bool WriteGVizGraph(string outfile);
	void WriteSNAPGraph(string outfile);
	void WriteDIMACSGraph(string outfile);

	/* General purpose utility functions */
	bool CheckValid();
	void Duplicate(const KGraph &rhs);
	void DuplicateConnected(const KGraph &rhs, map<long, long> &node_map);
	bool AddEdge(long i, long j, bool reverseToo, bool safe);
	bool DeleteEdge(long i, long j, bool reverseToo, bool safe);
	long ConnectedVertices();  // number of non-isolated nodes

	vector<long> CommonNeighborsList(long u, long v);
	vector<long> CommonNeighborsList(long u, vector<long> &v);

	/* functions for induced subgraphs */
	// deletes the edges adjacent to nodes in V\S. Keeps them as isolated nodes
	void FindInducedGraph(vector<bool> &S);
	void FindInducedGraph(vector<long> &S);

	// Creates a brand new graph G[S] with new node count.
	KGraph CreateInducedGraph(vector<long> &S);
	KGraph CreateInducedGraph(vector<long> &S, vector<long> &ReverseMap);
	KGraph CreateInducedGraph(vector<bool> &S);
	KGraph CreateInducedGraph(vector<bool> &S, vector<long> &ReverseMap);

	void FindConnectedComponents(vector< vector< long> > &components, vector<long> &degreeZero);

	KGraph CreatePowerGraph(long s);
	vector<long> FindHeuristicClique(vector<long> &degeneracyorder, vector<long> &rightdegree);
	vector<long> FindDegeneracyOrdering(vector<long> &rightdegree);
	vector<long> FindDegeneracyOrdering();
	vector<long> FindVerticesOfKCore(vector<long> &degeneracyorder, vector<long> &rightdegree, long k);

	bool IsKClub(vector<bool> &S, long k);
	bool IsKClub(vector<long> &S, long k);

	/* Destructors */
	void clear();
	~KGraph();
};

#endif
