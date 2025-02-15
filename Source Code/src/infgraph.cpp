#include "infgraph.h"
#include <algorithm>
#include <chrono>
using namespace std;
using namespace std::chrono;

thread_local vector<int> InfGraph::visited;
thread_local int InfGraph::time_stamp;

bool hasDuplicateElements(const vector<int>& vec) {
	unordered_set<int> uniqueElements;
	for (int element : vec) {
		if (uniqueElements.count(element) > 0) {
			return true;
		}
		uniqueElements.insert(element);
	}
	return false;
}

void InfGraph::setParametersIG(int kA, int kB, vector<double> qq, bool ignore, string folder, string bseeds)
{

	this->kA = kA;
	this->kB = kB;
	this->qao = qq[0];
	this->qab = qq[1];
	this->qbo = qq[2];
	this->qba = qq[3];
	ignore_B = ignore;
	b_seeds_file_name = bseeds;

	output_file_name = folder + "/output/seeds_tim_" + to_string(k);
	output_file_name += "_" + double_to_string(qao) + "_" + double_to_string(qab) + "_" + double_to_string(qbo) + "_" + double_to_string(qba);
	if (ignore_B)
		output_file_name += "_1";
	else
		output_file_name += "_0";
	output_file_name += ".txt";
}


void InfGraph::readBSeedsIG()
{
	bSeeds.clear();
	if (ignore_B) {
		cout << "[info] config file indicates that B seeds should be ignored!" << endl;
		return;
	}

	ifstream myfile(b_seeds_file_name.c_str(), ios::in);
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			string line;
			getline(myfile, line);
			if (line == "")
				continue;
			int u = atoi(line.c_str());
			bSeeds.push_back(u);
		}
		myfile.close();
	}
	else {
		cout << "[error] unable to open file: " << graph->folder << "/bseeds.txt" << endl;
		exit(1);
	}

	std::sort(bSeeds.begin(), bSeeds.end());//sort the B seeds 

	cout << "[info] B-seeds are: ";
	for (auto i : bSeeds)
		cout << i << " ";
	cout << endl;
}


// compute influence spread
double InfGraph::computeInfluenceHyperGraph()
{
	unordered_set<int> s;
	for (auto t : seedSet) {
		for (auto tt : hyperG[t]) {
			s.insert(tt);
		}
	}
	double inf = (double)n * s.size() / hyperId;
	return inf;
}
double InfGraph::computeInfHyperGraph_Order(int theta) {
	if (theta <= 0) {
		throw std::invalid_argument("theta must be greater than 0");
	}
	unordered_set<int> s;
	for (auto t : aSeeds) {
		for (auto tt : hyperG_A[t]) {
			s.insert(tt); 
		}
	}
	for (auto t : bSeeds) {
		for (auto tt : hyperG_B[t]) {
			s.insert(tt); 
		}
	}
	double inf = static_cast<double>(n) * static_cast<double>(s.size()) / static_cast<double>(theta);
	return inf;
}

// print seed sets
void InfGraph::printSeeds()
{
	int count = 1;
	for (auto s : seedSet) {
		cout << "[seed] rank = " << count << ", id = " << s << endl;
		count += 1;
	}
}




long long RA_cost = 0;
long long RB_cost = 0;
long long total = 0;


int InfGraph::BuildSingleRRSet_q1(int uStart, int rrSetID, bool isrA, bool q1, vector<vector<int>>& hyperGT, vector<vector<int>>& hyperG) {

	auto& probT = graph->probT;
	double coin = randomDouble();
	double qo;  vector<int> temp;
	if (isrA) {
		qo = q1 ? 1.0 : qao;
	}
	else
		qo = q1 ? 1.0 : qbo;

	if (coin > qo) {
		return 0;
	}
	mutex& currentLock = isrA ? lockA : lockB;
	//add to RR set
	hyperGT[rrSetID].push_back(uStart);
	{
		lock_guard<mutex> lock(currentLock);
		hyperG[uStart].push_back(rrSetID);
	}

	queue<int> Q;
	Q.push(uStart);
	visited[uStart] = time_stamp;
	while (!Q.empty()) {
		int u = Q.front();
		Q.pop();
		auto& gTu = graph->gT[u];
		// go through all neighbours
		for (int j = 0; j < gTu.size(); j++) {
			int w = gTu[j]; // edge w --> u

			if (visited[w] == time_stamp)
				continue;
			
			/***************** test edge w->u **********************/
			int prevEdgeStatus = graph->get_thEdgeStatus(u, j);
			if (prevEdgeStatus == BLOCKED) {
				continue;
			}
			else if (prevEdgeStatus == INACTIVE) {
				coin = randomDouble();
				if (coin > probT[u][j]) {
					graph->set_thEdgeStatus(u, j, BLOCKED);
					continue;
				}
				else {
					graph->set_thEdgeStatus(u, j, LIVE);
				}
			}
			/*******************************************************/
			visited[w] = time_stamp;

			coin = randomDouble();
			if (coin <= qo) {
				Q.push(w);
				temp.push_back(w);
				hyperGT[rrSetID].push_back(w);;

			}
		}
	}//end-while
	{
		lock_guard<mutex> lock(currentLock); 
		for (int w : temp)
			hyperG[w].push_back(rrSetID);
	}
	return 0;
}
int InfGraph::BuildSingleRRSet_idp(int uStart, int rrSetID, bool isrA, bool q1, vector<vector<int>>& hyperGT, vector<vector<int>>& hyperG) {
	auto& probT = graph->probT;
	double coin = randomDouble();
	double qo;  vector<int> temp;

	qo = isrA ? qao : qbo;
	if (coin > qo) {
		return 0;
	}
	mutex& currentLock = isrA ? lockA : lockB;
	//add to RR set
	hyperGT[rrSetID].push_back(uStart);
	{
		lock_guard<mutex> lock(currentLock);
		hyperG[uStart].push_back(rrSetID);
	}

	queue<int> Q;
	Q.push(uStart);
	visited[uStart] = time_stamp;
	while (!Q.empty()) {
		int u = Q.front();
		Q.pop();
		auto& gTu = graph->gT[u];
		// go through all neighbours
		for (int j = 0; j < gTu.size(); j++) {
			int w = gTu[j]; // edge w --> u

			if (visited[w] == time_stamp)
				continue;

			/***************** test edge w->u **********************/
			int prevEdgeStatus = graph->get_thEdgeStatus(u, j);
			if (prevEdgeStatus == BLOCKED) {
				continue;
			}
			else if (prevEdgeStatus == INACTIVE) {
				coin = randomDouble();
				if (coin > probT[u][j]) {
					graph->set_thEdgeStatus(u, j, BLOCKED);
					continue;
				}
				else {
					graph->set_thEdgeStatus(u, j, LIVE);
				}
			}
			/*******************************************************/
			visited[w] = time_stamp;

			coin = randomDouble();
			if (coin <= qo) {
				Q.push(w);
				temp.push_back(w);
				hyperGT[rrSetID].push_back(w);;

			}
		}
	}//end-while
	{
		lock_guard<mutex> lock(currentLock);
		for (int w : temp)
			hyperG[w].push_back(rrSetID);
	}
	return 0;
}
void InfGraph::BuildHyperGraphR_q1() {
	reset_hyper();

	////////////////////////////////////////////Parallel////////////////////////////////////////////
	omp_set_num_threads(4);
	long long progress = 0;
#pragma omp parallel
	{
		graph->init_EdgeStatus();
		init_visit();

#pragma omp for
		for (int i = 0; i < theta; i++) {
			int r = randomInt(n);
			//Afirst
			//Bfirst
			BuildSingleRRSet_q1(r, i, true, true, hyperGT_A, hyperG_A);
			//BuildSingleRRSet_q1(r, i, true, false, hyperGT_A, hyperG_A);
			time_stamp++;
			BuildSingleRRSet_q1(r, i, false, false, hyperGT_B, hyperG_B);
			//BuildSingleRRSet_q1(r, i, false, true, hyperGT_B, hyperG_B);
			time_stamp++;

			graph->reset_thEdgeStatus();

#pragma omp atomic
			progress++;

			if (progress % (theta / 10) == 0) {
				cout << "Progress: " << (100 * progress / theta) << "% completed" << endl;
			}
		}
		cout << "Thread ID: " << omp_get_thread_num() << " finished!" << endl;
		///////////////////////////////////////////////////////////////////////////////////////////
		
	}

}
void InfGraph::BuildHyperGraphR_idp() {
	reset_hyper();

	////////////////////////////////////////////Parallel////////////////////////////////////////////
	omp_set_num_threads(4);
	long long progress = 0;
#pragma omp parallel
	{
		graph->init_EdgeStatus();
		init_visit();

#pragma omp for
		for (int i = 0; i < theta; i++) {
			int r = randomInt(n);
			//Afirst
			//Bfirst
			BuildSingleRRSet_idp(r, i, true, true, hyperGT_A, hyperG_A);
			//BuildSingleRRSet_idp(r, i, true, false, hyperGT_A, hyperG_A);
			time_stamp++;
			BuildSingleRRSet_idp(r, i, false, false, hyperGT_B, hyperG_B);
			//BuildSingleRRSet_idp(r, i, false, true, hyperGT_B, hyperG_B);
			time_stamp++;

			graph->reset_thEdgeStatus();

#pragma omp atomic
			progress++;

			if (progress % (theta / 10) == 0) {
				cout << "Progress: " << (100 * progress / theta) << "% completed" << endl;
			}
		}
		cout << "Thread ID: " << omp_get_thread_num() << " finished!" << endl;
		///////////////////////////////////////////////////////////////////////////////////////////

	}

}
void InfGraph::BuildSeedSetIG_q1(bool output)
{
	int a = hyperGT_A.size(), b = hyperGT_B.size();
	vector<int> visit_A(a);vector<int> visit_B(b);
	degree_A.clear(); degree_B.clear();
	aSeeds.clear();	bSeeds.clear();

	for (int i = 0; i < n; i++) {
		degree_A.push_back(hyperG_A[i].size());
		degree_B.push_back(hyperG_B[i].size());
	}
	int cA = 0, cB = 0;
	for (int i = 0; i < kA + kB; i++) {
		auto tA = max_element(degree_A.begin(), degree_A.end());
		auto tB = max_element(degree_B.begin(), degree_B.end());
		int aseed = tA - degree_A.begin(), bseed = tB - degree_B.begin();
		int coverA = *tA, coverB = *tB;
		
		if(output)
			cout << "[info] A seed <" << aseed << "> coverage : " << coverA << ", B seed <" << bseed << "> coverage : " << coverB << endl;

		if (cB < kB && (coverB >= coverA || cA == kA)) {

			cB++;
			if(output)
				cout << "Round: " << i << ", new B seed : " << bseed << endl;
			bSeeds.push_back(bseed);

			for (int t : hyperG_B[bseed]) {
				if (!visit_B[t]) {
					visit_B[t] = true;
					for (int item : hyperGT_B[t]) {
						degree_B[item]--;
					}
				}
				if (!visit_A[t]) {
					visit_A[t] = true;
					for (int item : hyperGT_A[t]) {
						degree_A[item]--;
					}
				}

			}
		}
		//select A seed
		else {

			cA++;
			if (output)
				cout << "Round: " << i << ", new A seed : " << aseed << endl;
			aSeeds.push_back(aseed);

			for (int t : hyperG_A[aseed]) {

				if (!visit_A[t]) {
					visit_A[t] = true;
					for (int item : hyperGT_A[t]) {
						degree_A[item]--;
					}
				}
			}
			for (int t : hyperG_A[aseed]) {
				if (!visit_B[t]) {
					visit_B[t] = true;
					for (int item : hyperGT_B[t]) {
						degree_B[item]--;
					}
				}
			}
		}
		if (output)
			cout << endl;
	}//end for
}
void InfGraph::BuildSeedSetIG_idp() {
	int a = hyperGT_A.size(), b = hyperGT_B.size();
	vector<int> visit_A(a); vector<int> visit_B(b);
	degree_A.clear(); degree_B.clear();
	aSeeds.clear();	bSeeds.clear();
	for (int i = 0; i < n; i++) {
		degree_A.push_back(hyperG_A[i].size());
		degree_B.push_back(hyperG_B[i].size());
	}
	//mine A Seeds
	for (int i = 0; i < kA; i++) {
		auto t = max_element(degree_A.begin(), degree_A.end());
		int id = t - degree_A.begin();
		aSeeds.push_back(id);
		degree_A[id] = 0;
		for (int t : hyperG_A[id]) {
			if (!visit_A[t]) {
				visit_A[t] = true;
				for (int item : hyperGT_A[t]) {
					degree_A[item]--;
				}
			}
		}
	}
	//mine B Seeds
	for (int i = 0; i < kB; i++) {
		auto t = max_element(degree_B.begin(), degree_B.end());
		int id = t - degree_B.begin();
		bSeeds.push_back(id);
		degree_B[id] = 0;
		for (int t : hyperG_B[id]) {
			if (!visit_B[t]) {
				visit_B[t] = true;
				for (int item : hyperGT_B[t]) {
					degree_B[item]--;
				}
			}
		}
	}
}

double  InfGraph::calculateRepeatRate(const vector<int>& aSeeds, const vector<int>& bSeeds) {
	int count = 0;
	for (int i = 0; i < aSeeds.size(); i++) {
		if (find(bSeeds.begin(), bSeeds.end(), aSeeds[i]) != bSeeds.end()) {
			count++;
		}
	}
	return static_cast<double>(count) / aSeeds.size();
}
void InfGraph::printTree(RNode* root) {
	if (root == nullptr) {
		cout << "ERROR: RRTree is empty!" << endl;
		return;
	}

	stack<pair<RNode*, int>> sta;  // Stack stores node and depth level
	sta.push({ root, 0 });

	while (!sta.empty()) {
		auto element = sta.top();
		RNode* cur = element.first;
		int depth = element.second;
		sta.pop();

		// Print current node with depth indentation
		cout << string(depth * 2, ' ') << cur->vex << endl;

		// Push each child to the stack with increased depth
		for (auto it = cur->children.rbegin(); it != cur->children.rend(); ++it) {
			// Print edge relationship
			cout << string((depth + 1) * 2, ' ') << "|-- " << cur->vex << " -> " << (*it)->vex << endl;

			// Push the child into the stack with its depth
			sta.push({ (*it).get(), depth + 1 });
		}
	}
}
unique_ptr<RNode> InfGraph::BuildSingleARRSet_Order(int uStart, int rrSetID, bool Afirst) {
	// Create root node and manage it with unique_ptr
	auto root = std::make_unique<RNode>(uStart);
	auto& probT = graph->probT;

	// Temporary storage for data to be written
	vector<int> temp;

	// Test the root node
	double coin = randomDouble();
	if (coin > qao) {
		if (Afirst) examined[rrSetID].insert(uStart);
		return nullptr; // Return an empty tree
	}

	hyperGT_A[rrSetID].push_back(uStart); // Add to RR-set
	RRP_A[rrSetID].emplace(uStart, root.get()); // Store the raw pointer of root in the index container
	{
		// Critical section to protect write operations on hyperG_A
		lock_guard<mutex> lock(lockA);
		hyperG_A[uStart].push_back(rrSetID);
	}

	// Initialize BFS queue
	queue<RNode*> Q; // Queue stores raw pointers, does not own the nodes
	Q.push(root.get());
	visited[uStart] = time_stamp;

	// BFS traversal
	while (!Q.empty()) {
		RNode* curNode = Q.front();
		Q.pop();
		int u = curNode->vex;
		auto& gTu = graph->gT[u];

		// Traverse neighbor nodes
		for (int j = 0; j < gTu.size(); j++) {
			int w = gTu[j]; // edge w --> u
			if (visited[w] == time_stamp)
				continue;

			/***************** Test edge w->u **********************/
			int prevEdgeStatus = graph->get_thEdgeStatus(u, j);
			if (prevEdgeStatus == BLOCKED) {
				continue;
			}
			else if (prevEdgeStatus == INACTIVE) {
				coin = randomDouble();
				if (coin > probT[u][j]) {
					graph->set_thEdgeStatus(u, j, BLOCKED);
					continue;
				}
				else {
					graph->set_thEdgeStatus(u, j, LIVE);
				}
			}
			/*******************************************************/
			visited[w] = time_stamp;

			// Check node and construct RR set
			coin = randomDouble();
			if (coin <= qao) {
				// Create child node
				auto newNode = std::make_unique<RNode>(w);
				newNode->parent = curNode; // Set parent node
				RNode* rawNewNode = newNode.get(); // Get raw pointer
				curNode->children.push_back(std::move(newNode)); // Add child node to the current node's children list

				// Store content requiring lock in temporary variable
				temp.push_back(w);
				// Write content not requiring lock directly
				hyperGT_A[rrSetID].push_back(w);
				RRP_A[rrSetID].emplace(w, rawNewNode);
				Q.push(rawNewNode);
			}
			else if (Afirst) {
				examined[rrSetID].insert(w);
			}
		}
	}
	{
		lock_guard<mutex> lock(lockA); // Lock to protect write operations
		for (int w : temp)
			hyperG_A[w].push_back(rrSetID);
	}

	return root; // Return the unique_ptr of root
}
unique_ptr<RNode> InfGraph::BuildSingleBRRSet_Order(int uStart, int rrSetID, bool Afirst) {
	//initailization

	auto root = std::make_unique<RNode>(uStart);
	auto& probT = graph->probT;

	vector<int> temp;     

	//test root v
	double coin = randomDouble();
	if (coin > qbo) {
		if (!Afirst) examined[rrSetID].insert(uStart);
		return NULL;
	}
	hyperGT_B[rrSetID].push_back(uStart); 
	RRP_B[rrSetID].emplace(uStart, root.get());
	{
		lock_guard<mutex> lock(lockB);
		hyperG_B[uStart].push_back(rrSetID);
	}

	queue<RNode*> Q; 
	Q.push(root.get());

	visited[uStart] = time_stamp;

	while (!Q.empty()) {
		RNode* curNode = Q.front();
		Q.pop();
		int u = curNode->vex;
		auto& gTu = graph->gT[u];

		for (int j = 0; j < gTu.size(); j++) {
			int w = gTu[j]; // edge w-->u
			if (visited[w] == time_stamp)
				continue;

			/***************** test edge w->u **********************/
			int prevEdgeStatus = graph->get_thEdgeStatus(u, j);
			if (prevEdgeStatus == BLOCKED) {
				continue;
			}
			else if (prevEdgeStatus == INACTIVE) {
				coin = randomDouble();
				if (coin > probT[u][j]) {
					graph->set_thEdgeStatus(u, j, BLOCKED);
					continue;
				}
				else {
					graph->set_thEdgeStatus(u, j, LIVE);
				}
			}
			/*******************************************************/

			visited[w] = time_stamp;

			coin = randomDouble();
			if (coin <= qbo) {

				auto newNode = std::make_unique<RNode>(w);
				newNode->parent = curNode;
				RNode* rawNewNode = newNode.get();
				curNode->children.push_back(std::move(newNode));

				temp.push_back(w);
				hyperGT_B[rrSetID].push_back(w);
				RRP_B[rrSetID].emplace(w,rawNewNode);
				Q.push(rawNewNode);
			}
			else if (!Afirst) {
				examined[rrSetID].insert(w);
			}
		}
	}
	{
		std::lock_guard<std::mutex> lock(lockB); 
		for (int w : temp)
			hyperG_B[w].push_back(rrSetID);
	}
	return root;
}
void InfGraph::BuildHypergraphR_Order(int run, bool Afirst) {

	if (run > 0) {
		reset_order();
		resetRNG();
	}
	else init_order();

	reset_hyper();

	omp_set_num_threads(1);
	int progress = 0;

	////////////////////////////////////////////Parallel////////////////////////////////////////////
	#pragma omp parallel
	{
		graph->init_EdgeStatus();
		init_visit();
		
		#pragma omp for
		for (int i = 0; i < theta; i++) {
			int r = randomInt(n);

			unique_ptr<RNode> RaTree = BuildSingleARRSet_Order(r, i, Afirst);
			RRTree_A[i] = move(RaTree);
			time_stamp++;

			unique_ptr<RNode> RbTree = BuildSingleBRRSet_Order(r, i, Afirst);
			RRTree_B[i] = move(RbTree);
			time_stamp++;

			graph->reset_thEdgeStatus();
			#pragma omp atomic
			progress++;

			if (progress % (theta / 10) == 0) {
				cout << "Progress: " << (100 * progress / theta) << "% completed" << endl;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		cout << "Thread ID: " << omp_get_thread_num() << " finished!" << endl;
	}


}
////////////////////////////////////Afirst//////////////////////////////////////////////////////
void InfGraph::BuildSeedSetIG_Order_Afirst()
{
	int a = hyperGT_A.size(), b = hyperGT_B.size();
	vector<int> visit_A(a);
	vector<int> visit_B(b);
	unordered_set<int> delayDel;
	degree_A.clear(); degree_B.clear();
	aSeeds.clear();	bSeeds.clear();
	for (int i = 0; i < n; i++) {
		degree_A.push_back(hyperG_A[i].size());
		degree_B.push_back(hyperG_B[i].size());

	}

	int cA = 0, cB = 0; double coin;
	for (int i = 0; i < kA + kB; i++) {

		auto tA = max_element(degree_A.begin(), degree_A.end());
		auto tB = max_element(degree_B.begin(), degree_B.end());
		int aseed = tA - degree_A.begin(), bseed = tB - degree_B.begin();
		int coverA = *tA, coverB = *tB;

		cout << "[info] A seed <" << aseed << "> coverage : " << coverA << ", B seed <" << bseed << "> coverage : " << coverB << endl;

		/*
		bool priorityA = false;
		if (i == 0 || tA - degree_A.begin() == tB - degree_B.begin())
			priorityA = true;
		*/

		//select A seed
		delayDel.clear();
		if (cA < kA && (coverA >= coverB || cB == kB)) {
			//	degree_B[aseed] = 0;//won't be B seed
			cA++;
			cout << "Round: " << i << ", new A seed : " << aseed << endl;
			aSeeds.push_back(aseed);

			//delete all rA which contains aseed
			for (int t : hyperG_A[aseed]) {
				if (!visit_A[t]) {
					visit_A[t] = true;
					for (int item : hyperGT_A[t]) {
						degree_A[item]--;
					}
				}
				//trim all pair which rA contains aseed
				//Case1: aseed exist in rA, regardless of rB
				if (!visit_B[t] && RRTree_B[t]) {
					visit_B[t] = pruningRB_01(RRTree_A[t].get(), RRP_A[t][aseed], t);
				}

			}

			for (int t : hyperG_B[aseed]) {
				//RPointer prob = RPointer(t, nullptr);
				//auto it = RRP_A[aseed].find(prob);
				auto it = RRP_A[t].find(aseed);
				if (!visit_B[t] && it == RRP_A[t].end()) {//Case2: aseed exist in rB but not in rA
					pruningRB_02(RRTree_B[t].get(), RRP_B[t][aseed], t, delayDel);

				}
			}
			for (int elem : delayDel) {
				auto& arr = hyperG_B[aseed];
				arr.erase(remove(arr.begin(), arr.end(), elem), arr.end());
			}

		}
		//select B seed
		else {

			// degree_A[bseed] = 0;//won't be A seed
			cB++;
			cout << "Round: " << i << ", new B seed : " << bseed << endl;
			bSeeds.push_back(bseed);
			if (bseed == 641)
				cout << "there";
			//calculate raw discount for A seed
			int discount1 = coverB * (1 - qba / qbo);//raw discount
			int discount2 = 0;                     //precise discount

			for (int t : hyperG_B[bseed]) {
				if (!visit_B[t]) {

					// Mark this `t` as visited
					visit_B[t] = true;

					//calculate precise discount for A seed

					if (tested[t].find(bseed) == tested[t].end()) {
						double coin = randomDouble();
						if (coin > qba / qbo) {
							discount2++;
						}
					}

					// Update the degree of items in `hyperGT_B[t]`
					for (int item : hyperGT_B[t]) {
						degree_B[item]--;
					}
				}
			}
		}

		cout << endl;
	}//end for


	cout << "A seed set: ";
	for (int i = 0; i < kA; i++)
		cout << aSeeds[i] << " ";
	cout << endl;
	cout << "B seed set: ";
	for (int i = 0; i < kB; i++)
		cout << bSeeds[i] << " ";
	cout << endl;


}
bool InfGraph::pruningRB_01(RNode* root, RNode* aseed, int rrSetID) {
	if (root == nullptr)
		return false;
	RNode* cur_a = aseed;
	int v, t; double coin;
	double hold = qba / qbo;
	auto& test = tested[rrSetID];
	// Backtrack from aseed to root in ra
	while (cur_a != root) {
		v = cur_a->vex;

		if (test.find(v) == test.end()) {
			test.insert(v);
			// If v exists in rb 
			// RPointer prob = RPointer(rrSetID, nullptr);
			// auto it = RRP_B[v].find(prob);
			// if (it != RRP_B[v].end()) 
			auto loc = RRP_B[rrSetID].find(v);
			if (loc != RRP_B[rrSetID].end()) {
				coin = randomDouble();
				if (coin > qba) {
					stack<RNode*> sta;
					RNode* cur_b = loc->second;
					sta.push(cur_b);
					vector<unique_ptr<RNode>> nodesToDelete;

					while (!sta.empty()) {
						RNode* s = sta.top();
						sta.pop();

						// Update related statistics
						t = s->vex;
						degree_B[t]--;

						auto& vec_p1 = hyperGT_B[rrSetID];
						vec_p1.erase(remove(vec_p1.begin(), vec_p1.end(), t), vec_p1.end());

						auto& vec_p2 = hyperG_B[t];
						vec_p2.erase(remove(vec_p2.begin(), vec_p2.end(), rrSetID), vec_p2.end());

						RRP_B[rrSetID].erase(t);

						for (auto& child : s->children) {
							sta.push(child.get());
						}

						// Delay deletion
						if (!s->children.empty()) {
							nodesToDelete.push_back(std::move(s->children.back())); // Transfer ownership to nodesToDelete
							s->children.pop_back(); // Remove from children list
						}
					}

					// Remove cur_b from its parent's siblings
					auto& siblings = cur_b->parent->children;  // Get parent's children list
					siblings.erase(remove_if(siblings.begin(), siblings.end(),
						[cur_b](const unique_ptr<RNode>& node) { return node.get() == cur_b; }),
						siblings.end());
				}
			}
		} // End if
		// cout << cur_a->vex << "->";
		cur_a = cur_a->parent;
	} // End while
	if (cur_a == root) {
		// cout << "Root reached" << endl;
		v = cur_a->vex;
		if (test.find(v) == test.end()) {
			test.insert(v);
			double coin = randomDouble();
			if (coin > qba) { // Delete the whole RR set
				for (int t : hyperGT_B[rrSetID]) {
					degree_B[t]--;
					// RPointer temp(rrSetID, nullptr);
					// RRP_B[t].erase(temp);
					RRP_B[rrSetID].erase(t);
				}
				return true;
			}
		}
	}
	else {
		cout << "Error: cur_A is not root!" << endl;
	}
	return false;
}

bool InfGraph::pruningRB_02(RNode* root, RNode* aseed, int rrSetID, unordered_set<int>& delay) {
	if (root == nullptr)
		return true;

	RNode* cur = aseed;
	RNode* parent;
	int v, t;
	double coin1, coin2;
	double hold = qba / qbo;
	auto& test = tested[rrSetID];
	auto& examine = examined[rrSetID];

	// Container for delayed deletion, storing subtrees to be deleted
	std::vector<std::unique_ptr<RNode>> nodesToDelete;

	// Backtrack from `aseed` to the break point
	while (cur != root) {
		if (cur == nullptr) {
			std::cerr << "Error: nullptr cur !" << std::endl;
			break;
		}

		v = cur->vex;
		parent = cur->parent;

		if (examine.find(v) != examine.end())
			break;

		if (test.find(v) == test.end()) {
			test.insert(v);

			coin1 = randomDouble();
			coin2 = randomDouble();
			if (coin1 > qao) {
				break;
			}

			if (coin2 > qba) { // Pruning condition met
				std::stack<RNode*> sta;
				RNode* p = cur;

				sta.push(cur);

				// Delete current node from its parent's children
				auto& siblings = parent->children;
				auto it = std::find_if(siblings.begin(), siblings.end(),
					[p](const std::unique_ptr<RNode>& node) { return node.get() == p; });

				if (it != siblings.end()) {
					// Transfer ownership of the node to be deleted to the delayed deletion container
					nodesToDelete.push_back(std::move(*it));
					siblings.erase(it);
				}

				while (!sta.empty()) {
					RNode* s = sta.top();
					t = s->vex;
					sta.pop();

					if (s == nullptr) {
						std::cerr << "[Warning]: nullptr is in stack!" << std::endl;
						continue;
					}

					if (t >= 0 && t < degree_B.size()) {
						degree_B[t]--;
					}
					else {
						std::cerr << "[Warning]: Index t = " << t << " is out of bounds for degree_B." << std::endl;
						continue;
					}

					auto& vec_p1 = hyperGT_B[rrSetID];
					vec_p1.erase(std::remove(vec_p1.begin(), vec_p1.end(), t), vec_p1.end());
					// Delayed deletion to avoid issues in the calling loop
					if (t == aseed->vex)
						delay.insert(rrSetID);
					else {
						auto& vec_p2 = hyperG_B[t];
						vec_p2.erase(std::remove(vec_p2.begin(), vec_p2.end(), rrSetID), vec_p2.end());
					}

					RRP_B[rrSetID].erase(t);

					for (auto& child : s->children) {
						if (child != nullptr) {
							if (child->vex < 0) {
								std::cerr << "Error: child->vex is negative! Value = " << child->vex << std::endl;
							}
							else {
								sta.push(child.get());
							}
						}
					}

					// Delayed deletion: clear the current node's children without manual deletion, memory managed by `unique_ptr`
					if (!s->children.empty()) {
						nodesToDelete.push_back(std::move(s->children.back())); // Transfer ownership to nodesToDelete
						s->children.pop_back(); // Remove from children list
					}
				} // End while
			}
		}
		cur = parent;
	}

	// Delayed deletion nodes will be automatically released when `nodesToDelete` goes out of scope
	return true;
}
////////////////////////////////////Bfirst//////////////////////////////////////////////////////
void InfGraph::BuildSeedSetIG_Order_Bfirst()
{
	int a = hyperGT_A.size(), b = hyperGT_B.size();
	vector<int> visit_A(a);
	vector<int> visit_B(b);
	unordered_set<int> delayDel;
	degree_A.clear(); degree_B.clear();
	aSeeds.clear();	bSeeds.clear();
	for (int i = 0; i < n; i++) {
		degree_A.push_back(hyperG_A[i].size());
		degree_B.push_back(hyperG_B[i].size());
	}


	int cA = 0, cB = 0; double coin;
	for (int i = 0; i < kA + kB; i++) {

		auto tA = max_element(degree_A.begin(), degree_A.end());
		auto tB = max_element(degree_B.begin(), degree_B.end());
		int aseed = tA - degree_A.begin(), bseed = tB - degree_B.begin();
		int coverA = *tA, coverB = *tB;

		cout << "[info] A seed <" << aseed << "> coverage : " << coverA << ", B seed <" << bseed << "> coverage : " << coverB << endl;

		/*
		bool priorityA = false;
		if (i == 0 || tA - degree_A.begin() == tB - degree_B.begin())
			priorityA = true;
		*/
	
		delayDel.clear();
		if (cB < kB && (coverB >= coverA || cA == kA)) {
			degree_A[bseed] = 0;
			cB++;
			cout << "Round: " << i << ", new B seed : " << bseed << endl;
			bSeeds.push_back(bseed);

			for (int t : hyperG_B[bseed]) {
				if (!visit_B[t]) {
					visit_B[t] = true;
					for (int item : hyperGT_B[t]) {
						degree_B[item]--;
					}
				}

				if (!visit_A[t] && RRTree_A[t]) {
					visit_A[t] = pruningRA_01(RRTree_B[t].get(), RRP_B[t][bseed], t);
				}

			}

			for (int t : hyperG_A[bseed]) {

				auto it = RRP_B[t].find(bseed);
				if (!visit_A[t] && it == RRP_B[t].end()) {
					pruningRA_02(RRTree_A[t].get(), RRP_A[t][bseed], t, delayDel);
				}
			}
			for (int elem : delayDel) {
				auto& arr = hyperG_A[bseed];
				arr.erase(remove(arr.begin(), arr.end(), elem), arr.end());
			}
			


		}
		//select A seed
		else {
			degree_B[aseed] = 0;
			cA++;
			cout << "Round: " << i << ", new A seed : " << aseed << endl;
			aSeeds.push_back(aseed);
			int discount = 0;                  

			for (int t : hyperG_A[aseed]) {
				if (!visit_A[t]) {

					visit_A[t] = true;

					if (tested[t].find(aseed) == tested[t].end()) {
						double coin = randomDouble();
						if (coin > qab / qao) {
							discount++;
						}
					}

					for (int item : hyperGT_A[t]) {
						degree_A[item]--;
					}
				}
			}
			degree_B[aseed] -= discount * qbo;
		}

		cout << endl;
	}//end for


	cout << "A seed set: ";
	for (int i = 0; i < kA; i++)
		cout << aSeeds[i] << " ";
	cout << endl;
	cout << "B seed set: ";
	for (int i = 0; i < kB; i++)
		cout << bSeeds[i] << " ";
	cout << endl;


}
bool InfGraph::pruningRA_01(RNode* root, RNode* bseed, int rrSetID) {
	if (root == nullptr)
		return false;
	RNode* cur_b = bseed;
	int v, t; double coin;
	double hold = qab / qao;
	auto& test = tested[rrSetID];
	//backtrack from aseed to root in ra
	while (cur_b != root) {
		v = cur_b->vex;

		if (test.find(v) == test.end()) {
			test.insert(v);
			//if v exists in ra 
			auto loc = RRP_A[rrSetID].find(v);
			if (loc != RRP_A[rrSetID].end())
			{
				coin = randomDouble();
				if (coin > qab) {
					stack<RNode*> sta;
					RNode* cur_a = loc->second;
					sta.push(cur_a);
					vector<unique_ptr<RNode>> nodesToDelete;

					while (!sta.empty()) {
						RNode* s = sta.top();
						sta.pop();

						t = s->vex;
						degree_A[t]--;

						auto& vec_p1 = hyperGT_A[rrSetID];
						vec_p1.erase(remove(vec_p1.begin(), vec_p1.end(), t), vec_p1.end());
						auto& vec_p2 = hyperG_A[t];
						vec_p2.erase(remove(vec_p2.begin(), vec_p2.end(), rrSetID), vec_p2.end());

						RRP_A[rrSetID].erase(t);

						for (auto& child : s->children) {
							sta.push(child.get());
						}

						if (!s->children.empty()) {
							nodesToDelete.push_back(move(s->children.back())); 
							s->children.pop_back(); 
						}

					}

					auto& siblings = cur_a->parent->children;  
					siblings.erase(remove_if(siblings.begin(), siblings.end(),
						[cur_a](const unique_ptr<RNode>& node) { return node.get() == cur_a; }),
						siblings.end());

				}
			}

		}//end if
		//cout << cur_a->vex << "->";
		cur_b = cur_b->parent;
	}//end while
	if (cur_b == root) {
		//cout << "Root reached" << endl;
		v = cur_b->vex;
		if (test.find(v) == test.end()) {
			test.insert(v);
			double coin = randomDouble();
			if (coin > qab) {//delete the whole RR set
				for (int t : hyperGT_A[rrSetID]) {
					degree_A[t]--;
					RRP_A[rrSetID].erase(t);
				}
				return true;
			}

		}
	}
	else {
		cout << "[Error]: cur_B is not root!" << endl;
	}
	return false;
}
bool InfGraph::pruningRA_02(RNode* root, RNode* bseed, int rrSetID, unordered_set<int>& delay) {

	if (root == nullptr)
		return true;

	RNode* cur = bseed;
	RNode* parent;
	int v, t;
	double coin1, coin2;
	double hold = qab / qao;
	auto& test = tested[rrSetID];
	auto& examine = examined[rrSetID];

	std::vector<std::unique_ptr<RNode>> nodesToDelete;

	while (cur != root) {
		if (cur == nullptr) {
			std::cerr << "Error: nullptr cur !" << std::endl;
			break;
		}

		v = cur->vex;
		parent = cur->parent;

		if (examine.find(v) != examine.end())
			break;

		if (test.find(v) == test.end()) {
			test.insert(v);

			coin1 = randomDouble();
			coin2 = randomDouble();
			if (coin1 > qbo) {
				break;
			}

			if (coin2 > qab) { // Pruning condition met
				std::stack<RNode*> sta;
				RNode* p = cur;

				sta.push(cur);

				// Delete current node from its parent's children
				auto& siblings = parent->children;
				auto it = std::find_if(siblings.begin(), siblings.end(),
					[p](const std::unique_ptr<RNode>& node) { return node.get() == p; });

				if (it != siblings.end()) {
					nodesToDelete.push_back(move(*it));
					siblings.erase(it);
				}

				while (!sta.empty()) {
					RNode* s = sta.top();
					t = s->vex;
					sta.pop();

					if (s == nullptr) {
						std::cerr << "[Warning]: nullptr is in stack!" << std::endl;
						continue;
					}

					if (t >= 0 && t < degree_A.size()) {
						degree_A[t]--;
					}
					else {
						std::cerr << "[Warning]: Index t = " << t << " is out of bounds for degree_A." << std::endl;
						continue;
					}

					auto& vec_p1 = hyperGT_A[rrSetID];
					vec_p1.erase(remove(vec_p1.begin(), vec_p1.end(), t), vec_p1.end());

					if (t == bseed->vex)
						delay.insert(rrSetID);
					else {
						auto& vec_p2 = hyperG_A[t];
						vec_p2.erase(remove(vec_p2.begin(), vec_p2.end(), rrSetID), vec_p2.end());
					}

					RRP_A[rrSetID].erase(t);

					for (auto& child : s->children) {
						if (child != nullptr) {
							if (child->vex < 0) {
								std::cerr << "Error: child->vex is negative! Value = " << child->vex << std::endl;
							}
							else {
								sta.push(child.get());
							}
						}
					}

					if (!s->children.empty()) {
						nodesToDelete.push_back(move(s->children.back())); 
						s->children.pop_back(); 
					}
				}//end-while

			}
		}
		cur = parent;
	}
	return true;
}