#include "KGraph.h"
#include "LazyConstraints.h"
#include "GRBInterface.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <queue>
#include <set>
#include <map>
#include <omp.h>
#include <string>      
#include <iostream>    
#include <sstream>    

using namespace std;

bool KGraph::IsConnected(vector<bool> S)
{
	/* Is the subgraph induced by S connected? */
	long u, v;
	vector<bool> reached(n, false);
	vector<long> children, parents;

	for (long i = 0; i < n; i++) //find a root node
		if (S[i]) {
			children.push_back(i);
			reached[i] = true;
			break;
		}

	for (; !children.empty(); ) { //do BFS
		parents = children;
		children.clear();
		for (long i = 0; i < parents.size(); i++) { //for each parent, examine the children
			u = parents[i];
			for (long j = 0; j < degree[u]; j++) {
				v = adj[u][j];
				if (S[v] && !reached[v]) { //can only use vertices in S
					reached[v] = true;
					children.push_back(v);
				}
			}
		}
	}
	for (long i = 0; i < n; i++) //if a vertex in S hasn't been reached, return false
		if (S[i] && !reached[i])
			return false;

	return true;
}
bool KGraph::IsConnected(vector<long> S1)
{
	vector<bool> S(n, false); //convert S1 to bool
	for (long i = 0; i < S1.size(); i++)
		S[S1[i]] = true;
	return IsConnected(S);
}

bool KGraph::IsConnected()
{
	//Is G connected?
	vector<bool> S(n, true);
	return IsConnected(S);
}

long KGraph::DiameterUnweighted()
{
	/* Returns the (unweighted) diameter of a connected graph.
	If the graph is not connected it returns n.
	Solves a series of shortest paths problems. */
	vector<long> ShortestPaths = ShortestPathsUnweighted(0); //check to see if the graph is connected
	for (long i = 0; i < n; i++)
		if (ShortestPaths[i] == n)
			return n;

	long diameter = 0, temp_longest;
	for (long i = 0; i < n; i++) { //solve shortest paths problem, originating from node i
		temp_longest = LongestShortestPathUnweighted(i);
		if (temp_longest > diameter)
			diameter = temp_longest;
	}
	return diameter;
}

long KGraph::DiameterUnweighted(vector<long> S)
{
	/* returns the diameter of the graph induced by S */
	KGraph g = CreateInducedGraph(S);
	return g.DiameterUnweighted();
}

long KGraph::DiameterUnweighted(vector<bool> S1)
{
	vector<long> S;
	for (long i = 0; i < n; i++)
		if (S1[i])
			S.push_back(i);
	return DiameterUnweighted(S);
}

vector<long> KGraph::ShortestPathsUnweighted(long origin)
{
	vector<bool> S(n, true);
	return ShortestPathsUnweighted(origin, S);
}


vector<long> KGraph::ShortestPathsUnweighted(long origin, vector<bool> &S)
{
	vector<long> D;
	D.push_back(origin);
	return MultiSourceShortestPaths(D, S);
}


vector<long> KGraph::MultiSourceShortestPaths(vector<long> &D)
{
	vector<bool> S(n, true);
	return MultiSourceShortestPaths(D, S);
}


vector<long> KGraph::MultiSourceShortestPaths(vector<long> &D, vector<bool> &S)
{
	/*Finds the shortest paths from node v to set D.*/
	long u, v;
	vector<long> dist(n, n);
	vector<bool> reached(n, false);
	vector<long> children, parents;
	bool status = false;

	for (long i = 0; i < D.size(); i++)
	{
		if (S[D[i]]) status = true;
		else continue;
		children.push_back(D[i]);
		dist[D[i]] = 0;
		reached[D[i]] = true;
	}

	if (!status) return dist; // if none of sources are in S, return infinities.

	for (long d = 1; !children.empty(); d++)
	{
		parents = children;
		children.clear();
		for (long i = 0; i < parents.size(); i++)
		{
			u = parents[i];
			for (long j = 0; j < degree[u]; j++)
			{
				v = adj[u][j];
				if (!reached[v] && S[v])
				{
					reached[v] = true;
					dist[v] = d;
					children.push_back(v);
				}
			}
		}
	}
	return dist;
}


long KGraph::LongestShortestPathUnweighted(long origin)
{
	vector<long> SP = ShortestPathsUnweighted(origin);
	return (long)*max_element(SP.begin(), SP.end());
}

bool KGraph::DeleteNode(long i)
{
	/* Deletes all edges which are incident to node i.
	This does NOT actually remove the node per se.
	The value g.n remains the same. */
	for (long j = 0; j < degree[i]; j++)
	{
		long v = adj[i][j];
		DeleteEdge(v, i, false, false);
	}
	m -= degree[i];
	adj[i].clear();
	degree[i] = 0;
	return true;
}



/* Default constructor. Does nothing but initializing n and m. */
KGraph::KGraph()
{
	n = 0;
	m = 0;
	Delta = 0;
}

/* Default constructor, but names the graph in addition. */
KGraph::KGraph(string nm)
{
	name = nm;
	n = 0;
	m = 0;
	Delta = 0;
}

/* Useful constructor. Names the graph, and reads it from a files. */
KGraph::KGraph(string nm, string file, string type)
{
	n = 0;
	m = 0;
	Delta = 0;
	name = nm;
	if (type == "dimacs")
		ReadDIMACSGraph(file);
	else if (type == "snap_d")
		ReadSNAPGraph(file);
	else if (type == "dimacs_color")
		ReadDIMACSColorGraph(file);
	else if (type == "DAT")
		ReadDATGraph(file);
	else
		cerr << "Format " << type << " not found.\n";
}

/* Default constructor, but initializes the number of nodes and the vectors in addition */
KGraph::KGraph(long nodes)
{
	n = nodes;
	m = 0;
	Delta = 0;
	degree = new long[n];
	memset(degree, 0, sizeof(long)*n);
	adj = new vector<long>[n];
}

/* Copy constructor */
KGraph::KGraph(const KGraph &rhs)
{
	Duplicate(rhs);
}

/* Copying function. Makes the calling graph same as the passed graph. */
void KGraph::Duplicate(const KGraph &rhs)
{
	if (n > 0)
		clear();
	n = rhs.n;
	m = rhs.m;
	name = rhs.name;
	Delta = rhs.Delta;
	degree = new long[n];
	adj = new vector<long>[n];
	long i = 0;
	//#pragma omp parallel for	
	for (i = 0; i < n; i++)
	{
		degree[i] = rhs.degree[i];
		adj[i] = rhs.adj[i];
	}
}

/* Copying function. Makes the calling graph same as the passed graph, but removes isolated vertices
and renumbers the remaining vertices, making n smaller. */
void KGraph::DuplicateConnected(const KGraph &rhs, map<long, long> &node_map)
{
	if (n > 0)
		clear();

	map<long, long> node;
	for (long i = 0; i < rhs.n; i++)
		if (rhs.degree[i] > 0)
		{
			node[i] = n;
			node_map[n] = i;
			n++;
		}
	m = rhs.m;
	name = rhs.name;
	degree = new long[n];
	adj = new vector<long>[n];


	for (long i = 0; i < rhs.n; i++)
	{
		long t = 0;
		if (rhs.degree[i] <= 0)
			continue;
		t = node.find(i)->second;
		degree[t] = rhs.degree[i];
		adj[t].resize(rhs.degree[i]);
		for (long s = 0; s < rhs.degree[i]; s++)
		{
			adj[t][s] = node.find(rhs.adj[i][s])->second;
		}
	}
}

long KGraph::ConnectedVertices()
{
	/* Returns the number of nodes which have positive degree. */
	long connectedVertices = 0;
	for (long i = 0; i < n; i++)
		if (degree[i] > 0)
			connectedVertices++;
	return connectedVertices;
}

/* reverseToo: if false, j will be added to i's adj, but not the other way round.
* safe: if true, a check will be performed to make sure the edge does not already exists. */
bool KGraph::AddEdge(long i, long j, bool reverseToo, bool safe)
{
	vector<long>::iterator it;
	if (degree[i] == 0 || j > adj[i][degree[i] - 1])
		adj[i].push_back(j);
	else
	{
		it = lower_bound(adj[i].begin(), adj[i].end(), j);
		if (!safe)
			adj[i].insert(it, j);
		else
		{
			if (*it != j)
				adj[i].insert(it, j);
			else
				return false;
		}
	}
	degree[i]++;

	if (reverseToo)
	{
		if (degree[j] == 0)
			adj[j].push_back(i);
		else
		{
			it = lower_bound(adj[j].begin(), adj[j].end(), i);
			if (!safe)
				adj[j].insert(it, i);
			else
			{
				if (*it != i)
					adj[j].insert(it, i);
				else return false;
			}
		}
		degree[j]++;
	}
	return true;
}

/* reverseToo: if false, j will be removed from i's adj, but not the other way round.
* safe: if true, a check will be performed to make sure the edge exists. */
bool KGraph::DeleteEdge(long i, long j, bool reverseToo, bool safe)
{
	vector<long>::iterator it = lower_bound(adj[i].begin(), adj[i].end(), j);
	if (!safe)
		adj[i].erase(it);
	else
	{
		if (it != adj[i].end() && *it == j)
			adj[i].erase(it);
		else return false;
	}
	degree[i]--;
	if (reverseToo)
	{
		it = lower_bound(adj[j].begin(), adj[j].end(), i);
		if (!safe)
			adj[j].erase(it);
		else
		{
			if (it != adj[j].end() && *it == i)
				adj[j].erase(it);
			else return false;
		}
		degree[j]--;
	}
	return true;
}

bool KGraph::CheckValid()
{
	long m1 = 0;
	for (long i = 0; i < n; i++)
		m1 += degree[i];
	if (m1 != 2 * m)
	{
		cerr << "ERROR: " << m1 << " != " << 2 * m << endl;
		return false;
	}
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < degree[i] - 1; j++)
		{
			if (adj[i][j] >= adj[i][j + 1])
			{
				cerr << "ERROR: " << "adj[i][j]=" << adj[i][j] << " >= " << adj[i][j + 1] << "adj[i][j+1]" << endl;
				return false;
			}
		}
	}
	return true;
}


/* Returns the sorted list of common neighbors node u has with v */
vector<long> KGraph::CommonNeighborsList(long u, long v)
{
	vector<long> t;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t.push_back(t1);
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the sorted list of common neighbors node u has with v */
vector<long> KGraph::CommonNeighborsList(long u, vector<long> &v)
{
	vector<long> t;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = v.size();
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &v.front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t.push_back(t1);
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	return t;
}

bool KGraph::IsKClub(vector<long> &S, long k)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++)
		Sbool[S[i]] = true;
	return IsKClub(Sbool, k);
}

bool KGraph::IsKClub(vector<bool> &S, long k)
{
	KGraph g1 = CreateInducedGraph(S);
	if (g1.DiameterUnweighted() <= k) return true;
	else return false;
}

/* Finds the subgraph induced by all the nodes in S. S is sorted */
void KGraph::FindInducedGraph(vector<bool> &S)
{
	for (long i = 0; i < n; i++)
		if (!S[i])
			DeleteNode(i);
}
void KGraph::FindInducedGraph(vector<long> &S)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++) Sbool[S[i]] = true;
	return FindInducedGraph(Sbool);
}

KGraph KGraph::CreateInducedGraph(vector<bool> &S)
{
	vector<long> S1;
	for (long i = 0; i < n; i++)
	{
		if (S[i]) S1.push_back(i);
	}
	return CreateInducedGraph(S1);
}

KGraph KGraph::CreateInducedGraph(vector<long> &S)
{
	vector<long> ReverseMap;
	return CreateInducedGraph(S, ReverseMap);
}
KGraph KGraph::CreateInducedGraph(vector<bool> &S, vector<long> &ReverseMap)
{
	vector<long> vecS;
	for (long i = 0; i < n; i++)
	{
		if (S[i])
		{
			vecS.push_back(i);
		}
	}
	return CreateInducedGraph(vecS, ReverseMap);
}
KGraph KGraph::CreateInducedGraph(vector<long> &S, vector<long> &ReverseMap)
{
	/* Finds the subgraph induced by all the nodes in S. S is sorted */
	unsigned long S_size = S.size();
	KGraph g(S_size);
	ReverseMap.resize(n, -1);
	for (long i = 0; i < S_size; i++)
		ReverseMap[S[i]] = i;
	for (unsigned long i = 0; i < S_size; i++)
	{
		g.adj[i] = CommonNeighborsList(S[i], S);
		g.degree[i] = g.adj[i].size();
		for (long j = 0; j < g.degree[i]; j++) //relabel the vertices for the new, smaller graph
			g.adj[i][j] = ReverseMap[g.adj[i][j]];
		g.m += g.degree[i];
	}
	g.m /= 2;
	return g;
}

KGraph KGraph::CreatePowerGraph(long s)
{
	KGraph g(n);
	g.m = 0;

	for (long i = 0; i < n; i++)
	{
		vector<long> dist = ShortestPathsUnweighted(i);
		for (long j = 0; j < n; j++)
		{
			if (i != j && dist[j] <= s)
			{
				g.adj[i].push_back(j);
				g.degree[i]++;
				g.m++;
			}
		}
	}
	g.m /= 2;
	return g;
}

vector<long> KGraph::FindHeuristicClique(vector<long> &degeneracyorder, vector<long> &rightdegree)
{
	vector<long> clique;
	for (long i = 0; i < n && clique.empty(); i++)
	{
		long v = degeneracyorder[i];
		// if v neighbors all vertices after it in the ordering, 
		//		i.e., rightdegree[v] = (n-1) - (i+1) + 1 = n-i-1, 
		//		then v and all vertices after it form a clique.
		if (rightdegree[v] == n - i - 1)
		{
			clique.resize(n - i);
			for (long j = i; j < n; j++)	clique[j - i] = degeneracyorder[j];
			sort(clique.begin(), clique.end());
		}
	}
	return clique;
}
vector<long> KGraph::FindDegeneracyOrdering(vector<long> &rightdegree)
{
	long degeneracy = 0;
	rightdegree.resize(n);

	// initialize deg. Also update max degree Delta just in case.
	Delta = 0;
	for (long i = 0; i < n; i++)
	{
		rightdegree[i] = degree[i];
		Delta = max(Delta, rightdegree[i]);
	}

	// prepare the bins
	vector<long> bin(Delta + 1, (long)0);
	for (long i = 0; i < n; i++)	bin[rightdegree[i]]++;
	long start = 0;
	for (long d = 0; d <= Delta; d++)
	{
		long num = bin[d];
		bin[d] = start;
		start += num;
	}

	// initialize the ordering & position vectors
	vector<long> pos(n);	// pos[v] is position of vertex v in vert
	vector<long> vert(n);	// vert[v] is the v-th vertex in the ordering
	for (long i = 0; i < n; i++)
	{
		pos[i] = bin[rightdegree[i]];
		vert[pos[i]] = i;
		bin[rightdegree[i]]++;
	}

	// reset the bin starting points
	for (long d = Delta; d >= 1; d--) bin[d] = bin[d - 1];
	bin[0] = 0;

	// start peeling away minimum degree nodes
	for (long i = 0; i < n; i++)
	{
		long minv = vert[i];	// this is a min-degree vertex in the remaining graph
		bin[rightdegree[minv]]++;
		degeneracy = max(degeneracy, rightdegree[minv]);

		for (long j = 0; j < degree[minv]; j++) // adjust the degrees of the neighbors of v
		{
			long u = adj[minv][j];
			if (pos[u] > pos[minv])	// this means vertex u is still "in the graph" so we need to update its degree and its bucket
			{
				if (rightdegree[u] == rightdegree[minv])
				{
					long pw = bin[rightdegree[minv]];	// the first position of the bin that contains vertex minv
					long w = vert[pw];					// the vertex in that position
					if (u != w)						// if u is not the first vertex in the bin, swap u and w
					{
						vert[pw] = u;
						vert[pos[u]] = w;
						pos[w] = pos[u];
						pos[u] = pw;
					}
					bin[rightdegree[minv] - 1] = pos[minv] + 1;
					bin[rightdegree[u]]++;
					rightdegree[u]--;
				}
				else
				{
					long pw = bin[rightdegree[u]];
					long w = vert[pw];

					if (u != w)
					{
						vert[pw] = u;
						vert[pos[u]] = w;
						pos[w] = pos[u];
						pos[u] = pw;
					}
					bin[rightdegree[u]]++;
					rightdegree[u]--;
				}
			}
		}
	}
	//cerr << "\n Degeneracy = " << degeneracy << endl;
	return vert;
}

vector<long> KGraph::FindDegeneracyOrdering()
{
	vector<long> rightdegree;
	return FindDegeneracyOrdering(rightdegree);
}

vector<long> KGraph::FindVerticesOfKCore(vector<long> &degeneracyorder, vector<long> &rightdegree, long k)
{
	vector<long> vertices;
	for (long i = 0; i < n && vertices.empty(); i++)
	{
		long v = degeneracyorder[i];
		if (rightdegree[v] >= k)	// k-core is this vertex and all to the right in degeneracy order.
		{
			vertices.resize(n - i);
			for (long j = i; j < n; j++)	vertices[j - i] = degeneracyorder[j];
			sort(vertices.begin(), vertices.end());
		}
	}
	return vertices;
}

//check if B is subset of A
bool Issubset(vector<long> A, vector<long> B)
{
	sort(A.begin(), A.end());
	sort(B.begin(), B.end());
	return includes(A.begin(), A.end(), B.begin(), B.end());
}
/* This function finds all the connected components in a graph and places them in separate clusters.
* Nodes with degree 0 are placed in a vector called degreeZero.
* Arguments:
* clusters: the vector in which the components are stored. Could be non-empty, signifying
previously found components (For consistency, vertices in old clusters should have degree 0).
* degreeZero: array where singleton nodes are placed. Cleared before being filled. */
void KGraph::FindConnectedComponents(vector< vector< long> > &clusters, vector<long> &degreeZero)
{
	//cerr<<m<<" edges to start off in FindConnectedComponents."<<endl;
	long v;
	degreeZero.clear();
	bool* label = new bool[n];
	for (long i = 0; i < n; i++)
		label[i] = false;
	for (long i = 0; i < n; i++)
	{
		if (degree[i] != 0 && label[i] == false)
		{
			vector<long> cluster;
			cluster.push_back(i);
			label[i] = true;
			long c = 0;
			while (c != cluster.size())
			{
				long j = cluster[c];
				if (label[j] == false)
				{
					cluster.push_back(j);
					label[j] = true;
				}
				for (long t = 0; t < degree[j]; t++)
				{
					v = adj[j][t];
					if (label[v] == false)
					{
						cluster.push_back(v);
						label[v] = true;
					}
				}
				c++;
			}
			sort(cluster.begin(), cluster.end());
			clusters.push_back(cluster);
		}
	}
	for (long i = 0; i < n; i++)
		if (label[i] == false)
			degreeZero.push_back(i);
	delete[] label;
}

/* Reads a graph from a file in the DIMACS format.
* The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::ReadDIMACSGraph(string file)
{
	if (n > 0)
		clear();
	cerr << "ReadDIMACSGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	long u, v;
	ifstream input;

	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line[0] != '%')
			lineread = true;
	}
	istringstream nLine = istringstream(line);

	nLine >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i < n; i++)
	{
		lineread = false;
		line.clear();
		while (!lineread)
		{
			getline(input, line);
			if (line[0] != '%')
				lineread = true;
		}
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		if (line == "") continue;
		istringstream iLine = istringstream(line);

		v = -1;
		while (!iLine.eof())
		{
			iLine >> u;
			if (u != v)
			{
				adj[i].push_back(u - 1);
				degree[i]++;
				m++;
				//if(i>(u-1) && !binary_search(adj[u-1].begin(), adj[u-1].end(), i))
				//	cerr<<i<<"-"<<u<<" found, but "<<u<<"-"<<i<<"wasn't\n";
			}
			v = u;
		}
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		if (degree[i] != adj[i].size()) { cerr << " Error on line number " << i << endl; exit(0); }
		Delta = max(degree[i], Delta);
	}
	cerr << endl;
	m = m / 2;
	if (m1 != m)
	{
		cerr << "WRONG DATA!!!!!!!!!!!!!!!!!!!!! " << m << " != " << m1 << endl;
		exit(0);
	}
}

/* Reads a graph from a file in the DIMACS-2 format (clique/coloring challenges).
* The format is described at  */
void KGraph::ReadDIMACSColorGraph(string file)
{
	if (n > 0)
		clear();
	cerr << "ReadDIMACSColoringGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	long u, v;
	ifstream input;

	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line != "" && line[0] != 'c')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	string p;
	nLine >> p >> temp >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m1 / 1000000 << " dots: ";
	for (long i = 0; i < m1; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		getline(input, line);
		istringstream nLine = istringstream(line);

		nLine >> temp;
		if (temp == "e")
		{
			nLine >> u >> v;
			if (u == v)
				continue;
			adj[u - 1].push_back(v - 1);
			adj[v - 1].push_back(u - 1);
		}
		else i--;
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	if (m1 != m)
	{
		cerr << "Possible error in ReadDirectedGraphFromFile: " << m1 << "!=" << m << "\n";
		//	exit(0);
	}
}

/* Reads a graph from a file in the DIMACS format.
* The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::ReadDIMACSGraphParallel(string file)
{
	if (n > 0)
		clear();
	cerr << "ReadDIMACSGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	ifstream input;

	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line[0] != '%')
			lineread = true;
	}
	istringstream nLine = istringstream(line);

	nLine >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << n / 100000 << " dots: ";
	int nthreads = 0;
	int tid = 0;
	long lineNum = 0;

#pragma omp parallel shared(input,lineNum) private(line,tid)
	{
		long u, v;
		long share = 0;
		tid = nthreads++;
		cout << "Number of threads = " << tid << endl;
		string line;
		bool lineread;
		long i = 0;

		while (i < n && lineNum < n)
		{
			lineread = false;
			line.clear();
#pragma omp critical
			if (lineNum < n)
			{
#pragma omp flush (lineNum)
				while (!lineread && lineNum < n)
				{
					getline(input, line);
					if (line[0] != '%')
						lineread = true;
				}
				i = lineNum;
				lineNum++;
			}

			if ((i + 1) % 100000 == 0)
				cerr << ".";

			share++;
			if (i < n)
			{
				//cerr<<tid<<"------Using \t"<<i<<"\t"<<lineNum<<endl;
				if (line == "") continue;
				istringstream iLine = istringstream(line);

				v = -1;

				while (!iLine.eof())
				{
					iLine >> u;
					if (u != v)
					{
						adj[i].push_back(u - 1);
						degree[i]++;
#pragma omp atomic
						m++;
						//if(i>(u-1) && !binary_search(adj[u-1].begin(), adj[u-1].end(), i))
						//	cerr<<i<<"-"<<u<<" found, but "<<u<<"-"<<i<<"wasn't\n";
					}
					v = u;
				}
				sort(adj[i].begin(), adj[i].end());
				adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
				if (degree[i] != adj[i].size()) { cerr << " Error on line number " << i << " " << degree[i] << " " << adj[i].size() << endl; exit(0); }
#pragma omp critical(delta)
				Delta = max(degree[i], Delta);
			}
		}
		cerr << tid << " share = " << (double)share / (double)n << endl;
	}
	cerr << endl;
	m = m / 2;
	if (m1 != m)
	{
		cerr << "WRONG DATA!!!!!!!!!!!!!!!!!!!!! " << m << " != " << m1 << endl;
		exit(0);
	}
}
// reads the .dat graphs from Simonetti et al (2011) The Minimum Connected Dominating Set Problem: Formulation, Valid Inequalities and a Branch-and-Cut Algorithm
void KGraph::ReadDATGraph(string file)
{
	cerr << "ReadDATGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;

	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> m;
	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i < m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		v--; u--;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}
/* Reads a graph from a file in the SNAP format.
* The format is described at http://snap.stanford.edu/data/index.html
* Note: if (i,j) is read from file, both (i,j) and (j,i) are added to the graph to make sure its undirected*/
void KGraph::ReadSNAPGraph(string file)
{
	cerr << "ReadSNAPGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;

	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> temp >> m >> temp;

	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i < m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}

/* Writes the graph to a file in the SNAP format.
* The format is described at http://snap.stanford.edu/data/index.html
* Note: if (i,j) is read from file, both (i,j) and (j,i) are added to the graph to make sure its undirected*/
void KGraph::WriteSNAPGraph(string file)
{
	cerr << "WriteSNAPGraph ";
	ofstream output;
	output.open(file.c_str(), ios::out);
	if (!output.is_open())
	{
		cout << "File " << file << " could not be opened!!!\n";
		return;
	}

	output << n << " nodes, " << m << " edges.\n";

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i < n; i++)
	{
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		for (long j = 0; j < degree[i]; j++)
			output << i << "\t" << adj[i][j] << endl;
	}
}

/* Writes the graph to a file in the DIMACS-10 format.
* The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::WriteDIMACSGraph(string file)
{
	cerr << "WriteDIMACSGraph: " << file << endl;
	ofstream output;
	output.open(file.c_str(), ios::out);
	if (!output.is_open())
	{
		cout << "File " << file << " could not be opened!!!\n";
		return;
	}

	output << n << " " << m << endl;

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i < n; i++)
	{
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		for (long j = 0; j < degree[i]; j++)
			output << adj[i][j] + 1 << " ";
		output << endl;
	}
}

/* Writes the graph in the format such that GraphViz can plot it.
* No position is specified, so graphviz determines whats best */
bool KGraph::WriteGVizGraph(string outfile)
{
	ofstream gviz;
	gviz.open(outfile.c_str(), ios::out);
	gviz << "graph test\n{\nnode [shape=point, pin=true, fontsize=1];\n";
	double scale = 10;
	for (long i = 0; i < n; i++)
		gviz << i + 1 << ";\n";

	long temp;
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < degree[i]; j++)
		{
			temp = adj[i][j];
			if (i < temp)
				gviz << i + 1 << " -- " << temp + 1
				<< "[" << "color = gray55" << "]"
				<< ";\n";
		}
	}
	gviz << "}";
	gviz.close();
	return true;
}

void KGraph::clear()
{
	//cout<<"===Clearing DG "<<name<<"\n";
	if (n > 0)
	{
		for (long i = 0; i < n; i++)
			adj[i].clear();
		delete[]adj;
		delete[]degree;
		n = 0;
		m = 0;
	}
}

KGraph::~KGraph()
{
	//cout<<"===Destructing DG "<<name<<"\n";
	if (n > 0)
	{
		for (long i = 0; i < n; i++)
			adj[i].clear();
		delete[]adj;
		delete[]degree;
		n = 0;
		m = 0;
	}
}
