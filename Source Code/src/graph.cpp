#include "graph.h"

thread_local vector<vector<EdgeInfo>> Graph::thEdgeStatus;
thread_local long long Graph::th_version;

void Graph::readNM()
{
	ifstream cin((folder + "/attribute.txt").c_str());
	string s;
	while (cin >> s)
	{
		if (s.substr(0, 2) == "n=")
		{
			n = atoi(s.substr(2).c_str());
			continue;
		}
		if (s.substr(0, 2) == "m=")
		{
			m = atoi(s.substr(2).c_str());
			continue;
		}
		/*
		if (s.substr(0, 2) == "A=")
		{
			A_ratio = atof(s.substr(2).c_str());
			continue;
		}
		if (s.substr(0, 2) == "B=")
		{
			B_ratio = atof(s.substr(2).c_str());
			continue;
		}
		*/
	}
	cin.close();
	//cout << "[info] n = " << n << ", m = " << m << ", A ratio =" << A_ratio << ", B ration =" << B_ratio << endl;
}

void Graph::add_in_edge(int a, int b, double p)
{
	gT[b].push_back(a);
	probT[b].push_back(p);
	inDeg[b]++;

	EdgeInfo e;
	inEdgeStatus[b].push_back(e);
}

void Graph::add_out_edge(int a, int b, double p)
{
	gO[a].push_back(b);
	probO[a].push_back(p);
	outDeg[a]++;

	EdgeInfo e;
	outEdgeStatus[a].push_back(e);
	
}

void Graph::readGraph_v2()
{
	string fileName = folder + graph_file;
	ifstream myfile(fileName.c_str(), ios::in);
	string delim = " ";

	//data format
	// line1; a b p q
	// line2; b a q p
	if (myfile.is_open()) {
		string line;
		while (getline(myfile, line)) {
			if (line.empty()) continue;

			istringstream iss(line);
			vector<string> tokens;
			string token;


			while (iss >> token) {
				tokens.push_back(token);
			}

			if (tokens.size() < 4) {
				cout << "[warning] Invalid line format: " << line << endl;
				continue;
			}

			int a = std::stoi(tokens[0]);
			int b = std::stoi(tokens[1]);
			double p = std::stod(tokens[2]);
			double q = std::stod(tokens[3]);
			//cout << "a="<<a << endl;
			//cout << "b=" << b << endl;

			add_in_edge(a, b, p);
			add_out_edge(a, b, p);
		}
		myfile.close();

	}
	else {
		cout << "[error] can't open graph file " << fileName << endl;
		exit(1);
	}

	cout << "[info] finish reading graph data" << endl;
}

void Graph::readGraph()
{
	FILE* fin = nullptr;
	fopen_s(&fin, (folder + graph_file).c_str(), "r");
;

	int readCnt = 0;
	for (int i = 0; i < m; i++)
	{
		readCnt++;
		//cout << readCnt << endl;
		int a, b;
		double p;
		int c = fscanf_s(fin, "%d%d%lf", &a, &b, &p);
		//cout << a << ", " << b << ", " << p << endl;

		add_in_edge(a, b, p);
		add_out_edge(a, b, p);
	}

	fclose(fin);
	cout << "[info] finish reading graph data" << endl;
}

void Graph::printGraph(int num_nodes)
{
	if (num_nodes > n)
		num_nodes = n;
	for (int u = 0; u < num_nodes; u++) {
		printf("%d: ", u);
		for (int i = 0; i < outDeg[u]; i++) {
			printf("(%d, %g) ", gO[u][i], probO[u][i]);
		}
		printf("\n");
	}
}


// find top-k highest degree nodes
vector<int> Graph::findHighDegree(int k)
{
	cout << "[info] ranking nodes by degrees." << endl;
	vector<int> ret;
	vector<std::pair<int, int>> degVec;
	for (int i = 0; i < n; i++) {
		std::pair<int, int> p;
		p = std::make_pair(i, outDeg[i]);
		degVec.push_back(p);
	}

	std::sort(degVec.begin(), degVec.end(), [](const std::pair<int, int>& left, const std::pair<int, int>& right) {
		return left.second < right.second;
		});

	for (int i = 0; i < k; i++) {
		std::pair<int, int> p = degVec.at(degVec.size() - 1 - i);
		//cout << "degree = " << p.second << endl;
		ret.push_back(p.first);
	}

	return ret;
}

//PageRank
vector<int> Graph::findPageRank(int k)
{

	double damping = 0.85;
	int maxIter = 10000;
	double tol = 1e-6;
	cout << "[info] ranking nodes by PageRank." << endl;
	vector<int> retVec;
	vector<double> pr(n, 1.0 / n); // Initialize PageRank values for each node (uniform distribution)
	vector<double> new_pr(n, 0.0); // To store updated PageRank values
	for (int iter = 0; iter < maxIter; iter++) {
		// Initialize new PageRank values with the teleportation probability
		for (int v = 0; v < n; v++) {
			new_pr[v] = (1.0 - damping) / n;
		}

		// Update each node's PageRank   v->u->k
		for (int v = 0; v < n; v++) {
			for (int j = 0; j < inDeg[v]; j++) { // iterate over its In-neighbors
				int u = gT[v][j];
				int outDegree = outDeg[u];
				/*
				for (int k = 0; k < outDeg[u]; k++) {
					double coin = (double)rand() / (double)RAND_MAX;
					if (coin > probO[u][k]) {
						outDegree--;
					}
				}*/
				if (outDegree > 0) {
					//cout << "outDegree" << outDegree << endl;
					new_pr[v] += damping * (pr[u] / outDegree);
				}
			}
		}

		// Check for convergence by comparing old and new PageRank values
		double diff = 0.0;
		for (int v = 0; v < n; v++) {
			diff += fabs(new_pr[v] - pr[v]);
		}
		if (diff < tol) {
			cout << "Converged after " << iter + 1 << " iterations." << endl;
			pr = new_pr;
			break;
		}

		// Update the PageRank values for the next iteration
		pr = new_pr;
	}
	vector<pair<int, double>> pageVec;
	for (int v = 0; v < n; v++) 
		pageVec.push_back(make_pair(v, pr[v]));
	
	std::sort(pageVec.begin(), pageVec.end(), [](const pair<int, double>& left, const pair<int, double>& right) {
		return left.second > right.second;
		});

	for (int i = 0; i < k; i++) {

		pair<int, double> p = pageVec[i];
		cout << "node "<<p.first<<", PageRank = " << p.second << endl;
		retVec.push_back(p.first);
	}

	return retVec;
}

