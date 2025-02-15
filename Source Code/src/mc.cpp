#include "mc.h"
#include <algorithm>
#include <unordered_set>
#include<random>
void MonteCarlo::setParametersMC(int k_A, vector<double>qq, bool ignore, string folder, string bseeds)
{
	k = k_A;
	qao = qq[0];
	qab = qq[1];
	qbo = qq[2];
	qba = qq[3];

	
	reconsider_a = max((qab - qao), 0) / (1 - qao);
	reconsider_b = max((qba - qbo), 0) / (1 - qbo);
	cout << "reconsider_a: " << reconsider_a << " reconsider_b: " << reconsider_b << endl;
	qa_aware = (qao + qab) / 2;
	qb_aware = (qbo + qba) / 2;
	cout<< "qa_aware: " << qa_aware << " qb_aware: " << qb_aware << endl;
	ignore_B = ignore;
	b_seeds_file_name = bseeds;
	output_file_name = folder + "/output/seeds_mc_" + to_string(k);
	output_file_name += "_" + double_to_string(qao) + "_" + double_to_string(qab) + "_" + double_to_string(qbo) + "_" + double_to_string(qba);
	if (ignore_B)
		output_file_name += "_1";
	else
		output_file_name += "_0";
	output_file_name += ".txt";
}


void MonteCarlo::readBSeedsMC()
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
			if (line == "") continue;
			int u = atoi(line.c_str());
			bSeeds.push_back(u);
		}
		myfile.close();
	}
	else {
		cout << "[error] unable to open file: " << graph->folder << "/bseeds.txt" << endl;
		exit(1);
	}

	// we should not need sorting in any case
	//std::sort(bSeeds.begin(), bSeeds.end());

	cout << "[info] B-seeds (from input) are: ";
	for (auto i : bSeeds)
		cout << i << " ";
	cout << endl;

	//ASSERT((int)aSeeds.size() == 0);
}
void MonteCarlo::readBSeedsMC(string str)
{
	bSeeds.clear();
	if (ignore_B) {
		cout << "[info] config file indicates that B seeds should be ignored!" << endl;
		return;
	}

	ifstream myfile(str.c_str(), ios::in);
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			string line;
			getline(myfile, line);
			if (line == "") continue;
			int u = atoi(line.c_str());
			bSeeds.push_back(u);
		}
		myfile.close();
	}
	else {
		cout << "[error] unable to open file: " << graph->folder << "/bseeds.txt" << endl;
		exit(1);
	}


	cout << "[info] B-seeds (from input) are: ";
	for (auto i : bSeeds)
		cout << i << " ";
	cout << endl;

}

void MonteCarlo::readASeedsMC(string str)
{
	aSeeds.clear();
	ifstream myfile(str.c_str(), ios::in);
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			string line;
			getline(myfile, line);
			if (line == "") continue;
			int u = atoi(line.c_str());
			aSeeds.push_back(u);
		}
		myfile.close();
	}
	else {
		cout << "[error] unable to open file: " << graph->folder << "/" << str << endl;
		exit(1);
	}

	cout << "[info] A-seeds (from input) are: ";
	for (auto i : aSeeds)
		cout << i << " ";
	cout << endl;

}

void adjustArray(vector<int>& arr, const vector<int>& common_order) {
	vector<int> result;
	unordered_set<int> common_set(common_order.begin(), common_order.end());
	unordered_map<int, int> order_map;


	for (int i = 0; i < common_order.size(); ++i) {
		order_map[common_order[i]] = i;
	}


	for (int num : common_order) {
		if (find(arr.begin(), arr.end(), num) != arr.end()) {
			result.push_back(num);
		}
	}

	for (int num : arr) {
		if (common_set.find(num) == common_set.end()) {
			result.push_back(num);
		}
	}

	arr = result;
}
vector<int> computeCommonOrder(const vector<int>& arrA, const vector<int>& arrB) {
	unordered_set<int> setA(arrA.begin(), arrA.end());
	vector<int> common_order;

	for (int num : arrB) {
		if (setA.find(num) != setA.end()) {
			common_order.push_back(num);
		}
	}

	return common_order;
}
void preProcess(deque<int>& listA, deque<int>& listB) {

	vector<int> tempA(listA.begin(), listA.end());
	vector<int> tempB(listB.begin(), listB.end());

	vector<int> common_order = computeCommonOrder(tempA, tempB);

	adjustArray(tempA, common_order);
	adjustArray(tempB, common_order);

	listA.assign(tempA.begin(), tempA.end());
	listB.assign(tempB.begin(), tempB.end());
}
void printList(deque<int>& list) {
	vector<int> temp(list.begin(), list.end());
	cout << "[";
	for (int i = 0; i < temp.size(); i++) {
		cout << temp[i] << " ";
	}
	cout << "]";
	cout << endl;
}

/**
 *  baseSpread: \sigma_A(S_A, \emptyset)
 */



void MonteCarlo::adjustHeap(int* heap, double* improve, int i) {
	int x = 0;
	while (x * 2 + 2 <= n - i) {
		int newx = x * 2 + 1;
		if ((newx + 1 < n - i) && (improve[heap[newx]] < improve[heap[newx + 1]]))
			newx++;
		if (improve[heap[x]] < improve[heap[newx]]) {
			int t = heap[x];
			heap[x] = heap[newx];
			heap[newx] = t;
			x = newx;
		}
		else {
			break;
		}
	} // endwhile
}

void MonteCarlo::mineABseed_Greedy_Afirst(int kA, int kB, string file_SA,string file_SB) {
	double* improveSA = new double[graph->n];
	double* improveSB = new double[graph->n];
	int* last_update_A = new int[graph->n];
	int* last_update_B = new int[graph->n];
	int* heap_A = new int[graph->n];
	int* heap_B = new int[graph->n];
	vector<int> SA(kA,-1); vector<int> SB(kB,-1);
	double a, b, c, d;//empty argument

	for (int i = 0; i < n; i++) {
		heap_A[i] = i;
		heap_B[i] = i;
		last_update_A[i] = -1;
		last_update_B[i] = -1;
		improveSA[i] = (double)(n + 1);
		improveSB[i] = (double)(n + 1);
	}
	double old = 0;
	srand(time(NULL));

	aSeeds.clear();
	bSeeds.clear();
	mg.clear();
	int ca = 0, cb = 0;
	for (int i = 0; i < kA + kB; i++) {
		//construct A heap
		bool doneA = ca == kA; bool doneB = cb == kB;
		cout << "adjusting A heap ..." << endl;
		
		cout << "i=" << i << ", ca=" << ca << ", cb=" << cb << endl;
		while (!doneA && last_update_A[heap_A[0]] != i) {
			int v = heap_A[0];
			last_update_A[v] = i;

			SA[ca] = v;
				improveSA[v] = compute_coverage_seq_Afirst(SA, ca + 1, SB, cb, a, b, c, d) - old;
			SA[ca] = -1;

			adjustHeap(heap_A, improveSA, ca);
		}
		//construct B heap
		cout << "adjusting B heap ..." << endl;
		while (!doneB && last_update_B[heap_B[0]] != i) {
			int v = heap_B[0];
			last_update_B[v] = i;

			SB[cb] = v;
				improveSB[v] = compute_coverage_seq_Afirst(SA, ca, SB, cb + 1, a, b, c, d) - old;
			SB[cb] = -1;


			adjustHeap(heap_B, improveSB, cb);
		}
		//add one aseed
		if (!doneA && (improveSA[heap_A[0]] >= improveSB[heap_B[0]]||doneB)) {
			int chosed = heap_A[0];
			SA[ca] = chosed;
			aSeeds.push_back(chosed);
			old += improveSA[chosed];
			mg.push_back(improveSA[chosed]);
			heap_A[0] = heap_A[n - 1 - ca];
			ca++;

			cout << "\tround " << i + 1 << ": A seed = " << aSeeds[ca - 1] << ", mg = " << mg[i] << ", total = " << old << endl;
			adjustHeap(heap_A, improveSA, ca);
		}
		//add one bseed
		else if(!doneB && (improveSA[heap_A[0]] < improveSB[heap_B[0]]||doneA)) {
			int chosed = heap_B[0];
			SB[cb] = chosed;
			bSeeds.push_back(chosed);
			old += improveSB[chosed];
			mg.push_back(improveSB[chosed]);
			heap_B[0] = heap_B[n - 1 - cb];
			cb++;

			cout << "\tround " << i + 1 << ": B seed = " << bSeeds[cb - 1] << ", mg = " << mg[i] << ", total = " << old << endl;
			adjustHeap(heap_B, improveSB, cb);
		}
		
		
	}//endfor

	//write seeds to file 
	write_seeds_to_file(aSeeds, file_SA);
	write_seeds_to_file(bSeeds, file_SB);
	cout << "Seed set repeat rate: " << graph->calculateRepeatRate(aSeeds, bSeeds) * 100 << '%' << endl;
}
void MonteCarlo::mineABseed_Greedy_Bfirst(int kA, int kB, string file_SA, string file_SB) {
	double* improveSA = new double[graph->n];
	double* improveSB = new double[graph->n];
	int* last_update_A = new int[graph->n];
	int* last_update_B = new int[graph->n];
	int* heap_A = new int[graph->n];
	int* heap_B = new int[graph->n];
	vector<int> SA(kA, -1); vector<int> SB(kB, -1);

	double a, b, c, d;//empty argument

	for (int i = 0; i < n; i++) {
		heap_A[i] = i;
		heap_B[i] = i;
		last_update_A[i] = -1;
		last_update_B[i] = -1;
		improveSA[i] = (double)(n + 1);
		improveSB[i] = (double)(n + 1);
	}
	double old = 0;
	srand(time(NULL));

	aSeeds.clear();
	bSeeds.clear();
	mg.clear();
	int ca = 0, cb = 0;
	for (int i = 0; i < kA + kB; i++) {
		//construct A heap
		bool doneA = ca == kA; bool doneB = cb == kB;
		cout << "adjusting A heap ..." << endl;

		cout << "i=" << i << ", ca=" << ca << ", cb=" << cb << endl;
		while (!doneA && last_update_A[heap_A[0]] != i) {
			int v = heap_A[0];
			last_update_A[v] = i;

			SA[ca] = v;
			improveSA[v] = compute_coverage_seq_Bfirst(SA, ca + 1, SB, cb, a,b,c,d) - old;
			SA[ca] = -1;

			adjustHeap(heap_A, improveSA, ca);
		}
		//construct B heap
		cout << "adjusting B heap ..." << endl;
		while (!doneB && last_update_B[heap_B[0]] != i) {
			int v = heap_B[0];
			last_update_B[v] = i;

			SB[cb] = v;
			improveSB[v] = compute_coverage_seq_Bfirst(SA, ca, SB, cb + 1, a, b, c, d) - old;
			SB[cb] = -1;


			adjustHeap(heap_B, improveSB, cb);
		}
		//add one aseed
		if (!doneA && (improveSA[heap_A[0]] > improveSB[heap_B[0]] || doneB)) {
			int chosed = heap_A[0];
			SA[ca] = chosed;
			aSeeds.push_back(chosed);
			old += improveSA[chosed];
			mg.push_back(improveSA[chosed]);
			heap_A[0] = heap_A[n - 1 - ca];
			ca++;

			cout << "\tround " << i + 1 << ": A seed = " << aSeeds[ca - 1] << ", mg = " << mg[i] << ", total = " << old << endl;
			adjustHeap(heap_A, improveSA, ca);
		}
		//add one bseed
		else if (!doneB && (improveSA[heap_A[0]] <= improveSB[heap_B[0]] || doneA)) {
			int chosed = heap_B[0];
			SB[cb] = chosed;
			bSeeds.push_back(chosed);
			old += improveSB[chosed];
			mg.push_back(improveSB[chosed]);
			heap_B[0] = heap_B[n - 1 - cb];
			cb++;

			cout << "\tround " << i + 1 << ": B seed = " << bSeeds[cb - 1] << ", mg = " << mg[i] << ", total = " << old << endl;
			adjustHeap(heap_B, improveSB, cb);
		}

		
		
	}//endfor

	//write seeds to file 
	write_seeds_to_file(aSeeds, file_SA);
	write_seeds_to_file(bSeeds, file_SB);
	cout << "Seed set repeat rate: " << graph->calculateRepeatRate(aSeeds, bSeeds) * 100 << '%' << endl;
}



int MonteCarlo::oneRound_A_Spread(deque<int>& list_A, int* status_A, int* status_B, double* alpha_A, int* next_A, double& covA, int& cnt_qab) {
	bool finish = 1;
	if (!list_A.empty())
	{
		finish = 0;
		int v = list_A.front();
		if (status_A[v] != ADOPTED && status_A[v] != REJECTED) {

			if (status_B[v] != ADOPTED) {
				if (alpha_A[v] <= qao) {
					status_A[v] = ADOPTED;  
					covA++;
				}
				else status_A[v] = REJECTED;
			}
			else {
				cnt_qab++;
				
				if (alpha_A[v] <= qab) {
					status_A[v] = ADOPTED;
					covA++;
				}
				else status_A[v] = REJECTED;
			}
			
			if (status_A[v] == ADOPTED) {
				examine_out_neighbors(v, &list_A, next_A, status_A);
			}
		}//END-IF

	}//END-IF
	return finish;
}
int MonteCarlo::oneRound_B_Spread(deque<int>& list_B, int* status_A, int* status_B, double* alpha_B, int* next_B, double& covB, int& cnt_qba) {
	bool finish = 1;
	if (!list_B.empty())
	{
		finish = 0;
		int v = list_B.front();
		if (status_B[v] != ADOPTED && status_B[v] != REJECTED) {

			if (status_A[v] != ADOPTED) {
				if (alpha_B[v] <= qbo) {
					status_B[v] = ADOPTED;
					covB++;

				}
				else status_B[v] = REJECTED;
			}
			else {
				cnt_qba++;

				if (alpha_B[v] <= qba) {
					status_B[v] = ADOPTED;
					covB++;
				}
				else status_B[v] = REJECTED;
			}

			if (status_B[v] == ADOPTED) {
				examine_out_neighbors(v, &list_B, next_B, status_B);
			}
		}//END-IF

	}//END-IF
	return finish;
}

double MonteCarlo::compute_coverage_simu_v1(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, double& Acover, double& Bcover, double& num_qab, double& num_qba)
{


	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int* status_B = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	bool SeedTestFlag = 1;

	deque<int> q1, q2;


	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();


		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		for (int i = 0; i < k1; i++) {
			int u = aSeedSet.at(i);
			q1.push_back(u);
		}
		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);
			q2.push_back(u);
		}

		preProcess(q1, q2);
		
		////////////////////////////////////Tie breaking in seed test/////////////////////////////////////////
		double coin;
		while (!q1.empty() && !q2.empty()) {
			int v= q1.front(), u = q2.front();
			if (v == u) {
				q1.pop_front(); q2.pop_front();
				coin = randomDouble();
				if (coin > 0.5) {//A seed first

					if (SeedTestFlag && alpha_A[v] > qao) 
						status_A[v] = REJECTED;
					else {
						
						/////////////////////////////////////////////////////////////////////////////////////
						status_A[v] = ADOPTED;
						covA++;
						for (int j = 0; j < graph->outDeg[v]; j++) { // iterate over its out-neighbors
							int w = graph->gO[v][j];
							int preEdgeStatus = graph->get_outedge_status(v, j); 

							if (preEdgeStatus == INACTIVE) { 
								coin = randomDouble();
								if (coin <= graph->probO[v][j]) {
									graph->set_outedge_status(v, j, LIVE); 
									if (status_A[w] != ADOPTED) {
										list_A.push_back(w);
									}
								}
								else {
									graph->set_outedge_status(v, j, BLOCKED); 
								}
							}
							else if (preEdgeStatus == LIVE && status_A[w] != ADOPTED) {
								list_A.push_back(w);
							}
						}
					}

					if (SeedTestFlag && ((alpha_B[u] > qbo) || (status_A[u] == ADOPTED && alpha_B[u] > qba))) {
						status_B[u] = REJECTED;
						continue;
					}
					/////////////////////////////////////////////////////////////////////////////////////
					status_B[u] = ADOPTED;
					covB++;
					for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
						int w = graph->gO[u][j];


						int prvEdgeStatus = graph->get_outedge_status(u, j);

						if (prvEdgeStatus == INACTIVE) { 
							double coin = randomDouble();
							if (coin <= graph->probO[u][j]) {
				
								graph->set_outedge_status(u, j, LIVE);

								if (status_B[w] != ADOPTED) {
									list_B.push_back(w);
								}
							}
							else {
					
								graph->set_outedge_status(u, j, BLOCKED);
							}
						}
						else if (prvEdgeStatus == LIVE && status_B[w] != ADOPTED) {
							list_B.push_back(w);
						}
					}
				}
				else {//B seed first
					if (SeedTestFlag && alpha_B[u] > qbo)
						status_B[u] = REJECTED;
					else {
						/////////////////////////////////////////////////////////////////////////////////////
						status_B[u] = ADOPTED;
						covB++;
						for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
							int w = graph->gO[u][j];

				
							int prvEdgeStatus = graph->get_outedge_status(u, j);

							if (prvEdgeStatus == INACTIVE) { 
								double coin = randomDouble();
								if (coin <= graph->probO[u][j]) {
								
									graph->set_outedge_status(u, j, LIVE);

									if (status_B[w] != ADOPTED) {
										list_B.push_back(w);
									}
								}
								else {
							
									graph->set_outedge_status(u, j, BLOCKED);
								}
							}
							else if (prvEdgeStatus == LIVE && status_B[w] != ADOPTED) {
								list_B.push_back(w);
							}
						}
					}
					if (SeedTestFlag && ((alpha_A[v] > qao) || (status_B[v] == ADOPTED && alpha_A[v] > qab))) {
						status_A[v] = REJECTED;
						continue;
					}
					/////////////////////////////////////////////////////////////////////////////////////
					status_A[v] = ADOPTED;
					covA++;
					for (int j = 0; j < graph->outDeg[v]; j++) { // iterate over its out-neighbors
						int w = graph->gO[v][j];
						int preEdgeStatus = graph->get_outedge_status(v, j); 

						if (preEdgeStatus == INACTIVE) { 
							coin = randomDouble();
							if (coin <= graph->probO[v][j]) {
								graph->set_outedge_status(v, j, LIVE); 
								if (status_A[w] != ADOPTED) {
									list_A.push_back(w);
								}
							}
							else {
								graph->set_outedge_status(v, j, BLOCKED); 
							}
						}
						else if (preEdgeStatus == LIVE && status_A[w] != ADOPTED) {
							list_A.push_back(w);
						}
					}
				}
			}
			else break;
		}//end while

		while (!q1.empty()) {
			int v = q1.front(); q1.pop_front();
			if (SeedTestFlag && alpha_A[v] > qao) {
				status_A[v] = REJECTED;
				continue;
			}
			/////////////////////////////////////////////////////////////////////////////////////
			status_A[v] = ADOPTED;
			covA++;
			for (int j = 0; j < graph->outDeg[v]; j++) { // iterate over its out-neighbors
				int w = graph->gO[v][j];
				int preEdgeStatus = graph->get_outedge_status(v, j); 

				if (preEdgeStatus == INACTIVE) { 
					coin = randomDouble();
					if (coin <= graph->probO[v][j]) {
						graph->set_outedge_status(v, j, LIVE); 
						if (status_A[w] != ADOPTED) {
							list_A.push_back(w);
						}
					}
					else {
						graph->set_outedge_status(v, j, BLOCKED); 
					}
				}
				else if (preEdgeStatus == LIVE && status_A[w] != ADOPTED) {
					list_A.push_back(w);
				}
			}
				
		}
		while (!q2.empty()) {
			int u = q2.front(); q2.pop_front();
			if (SeedTestFlag && alpha_B[u] > qbo) {
				status_B[u] = REJECTED;
				continue;
			}
			/////////////////////////////////////////////////////////////////////////////////////
			status_B[u] = ADOPTED;
			covB++;
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int w = graph->gO[u][j];


				int prvEdgeStatus = graph->get_outedge_status(u, j);

				if (prvEdgeStatus == INACTIVE) { 
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
					
						graph->set_outedge_status(u, j, LIVE);

						if (status_B[w] != ADOPTED) {
							list_B.push_back(w);
						}
					}
					else {

						graph->set_outedge_status(u, j, BLOCKED);
					}
				}
				else if (prvEdgeStatus == LIVE && status_B[w] != ADOPTED) {
					list_B.push_back(w);
				}
			}
		}

		//cout << endl;
		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;

		while (curr_A > 0 || curr_B > 0) {

			
			preProcess(list_A, list_B);

			for (int i = 0; i < max(curr_A, curr_B); i++) {
				// A-adoption test
				double coin = randomDouble();
				//cout << "coin=" << coin << endl;
				if (coin <= 0.5) {

					oneRound_A_Spread(list_A, status_A, status_B, alpha_A, &next_A, covA, cnt_qab);
					oneRound_B_Spread(list_B, status_A, status_B, alpha_B, &next_B, covB, cnt_qba);
				}//End branch 1

				else {
					oneRound_B_Spread(list_B, status_A, status_B, alpha_B, &next_B, covB, cnt_qba);
					oneRound_A_Spread(list_A, status_A, status_B, alpha_A, &next_A, covA, cnt_qab);
				}//End branch2

				//reach the next element
				if (!list_A.empty())
					list_A.pop_front();
				if (!list_B.empty())
					list_B.pop_front();


			} // ENDFOR
			curr_A = next_A;
			curr_B = next_B;
			next_A = next_B = 0;

		} // END-WHILE
	}

	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;


	num_qab += cnt_qab / (double)MC_RUNS;
	num_qba += cnt_qba / (double)MC_RUNS;
	Acover += covA / (double)MC_RUNS;
	Bcover += covB / (double)MC_RUNS;
	return (covA + covB) / (double)MC_RUNS;
}

double MonteCarlo::compute_coverage_seq_Afirst(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, double& Acover, double& Bcover, double& num_qab, double& num_qba) {
	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int* status_B = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	double covB_qbo=0, covB_qba=0; double covA_qao=0,covA_reconsider = 0;

	bool SeedTestFlag = 1;
;
	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		// scan all a-seeds
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && alpha_A[u] > qao) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			covA++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);
					if (status_A[v] != ADOPTED)
						list_A.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED);
				}
			}
		}
		/*adpption test in sequence*/
		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;
		
		while (curr_A > 0) {
			for (int i = 0; i < curr_A; i++) {
				
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == ADOPTED || status_A[v] == REJECTED)
					continue;


				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						covA++; covA_qao++;
					}
					else {status_A[v] = REJECTED; }// A-rejected
				}//end qao adoption

				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) {examine_out_neighbors(v, &list_A, &next_A, status_A);} 
			} // ENDFOR
			curr_A = next_A; next_A = 0;
		}//END-WHILE

		// scan all B-seeds

		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);

			if (SeedTestFlag && ((alpha_B[u] > qbo) || (status_A[u] == ADOPTED && alpha_B[u] > qba) )) {
				status_B[u] = REJECTED;
				continue;
			}

			status_B[u] = ADOPTED;
			covB++;
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				int preEdgeStatus = graph->get_outedge_status(u, j);
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);  // edge is live
						if (status_B[v] != ADOPTED)
							list_B.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
					}
				}
				else if (preEdgeStatus == LIVE && status_B[v] != ADOPTED) {
					list_B.push_back(v);
				}

			}
		}

		curr_B = list_B.size(); next_B = 0;

		while (curr_B > 0) {
		
			for (int i = 0; i < curr_B; i++) {
				int v = list_B.front();
				list_B.pop_front();
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
				
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						covB++; covB_qbo++;

					}
					else { status_B[v] = REJECTED; }
				}//end qbo adoption


				else {
				
					cnt_qba++;
					if (alpha_B[v] <= qba) {
						status_B[v] = ADOPTED;
						covB++;covB_qba++;
					}
					else {status_B[v] = REJECTED;}
				}//end qba adoption


				if (status_B[v] == ADOPTED) {examine_out_neighbors(v, &list_B, &next_B, status_B);} 
			} // END-FOR
			curr_B = next_B;next_B = 0;
		} // END-WHILE
	}
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;


	num_qab += cnt_qab / (double)MC_RUNS;
	num_qba += cnt_qba / (double)MC_RUNS;
	Acover += covA / (double)MC_RUNS;
	Bcover += covB / (double)MC_RUNS;
	return (covA + covB) / (double)MC_RUNS;
}
double MonteCarlo::compute_coverage_seq_Afirst_rhythm(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, int timing, double& Acover, double& Bcover, double& num_qab, double& num_qba) {
	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int* status_B = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	double covB_qbo = 0, covB_qba = 0; double covA_qao = 0, covA_qab = 0;

	bool SeedTestFlag = 1;
	;
	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);

		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		///////////////////////// scan all a-seeds///////////////////////////////////
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && alpha_A[u] > qao) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			covA++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);
					if (status_A[v] != ADOPTED)
						//cout << "node=" << v << " statusA=" << status_A[v] << endl;
						list_A.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED);
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;
		int ticks = 0, curr;
		bool Bready = false, scaned = false;

		while (curr_A > 0 || curr_B > 0 || !Bready) {

			Bready = ticks >= timing;
			if (Bready && !scaned) {
				///////////////////////// scan all b-seeds///////////////////////////////////
				for (int i = 0; i < k2; i++) {
					int u = bSeedSet.at(i);

					if (SeedTestFlag && status_A[u] == ADOPTED && alpha_B[u] > qba) {
						//cout << "A";
						status_B[u] = REJECTED;
						continue;
					}

					status_B[u] = ADOPTED;
					covB++;
					for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
						int v = graph->gO[u][j];
						int preEdgeStatus = graph->get_outedge_status(u, j);
						if (preEdgeStatus == INACTIVE) {
							double coin = randomDouble();
							if (coin <= graph->probO[u][j]) {
								graph->set_outedge_status(u, j, LIVE);  // edge is live
								if (status_B[v] != ADOPTED)
									list_B.push_back(v);
							}
							else {
								graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
							}
						}
						else if (preEdgeStatus == LIVE && status_B[v] != ADOPTED) {
							list_B.push_back(v);
						}

					}
				}
				scaned = true;
				curr_B = list_B.size();
				if (r % 1000 == 0)
					cout << "after scan all B-seeds,curr_B: " << list_B.size() << endl;
				////////////////////////////////////////////////////////////////////////////////
			}

	
			curr = !Bready ? curr_A : max(curr_A, curr_B);
			for (int i = 0; i < curr; i++) {

				double coin = randomDouble();
				if (coin <= 0.5) {
					oneRound_A_Spread(list_A, status_A, status_B, alpha_A, &next_A, covA, cnt_qab);
					if (Bready)
						oneRound_B_Spread(list_B, status_A, status_B, alpha_B, &next_B, covB, cnt_qba);

				}
				else {
					if (Bready)
						oneRound_B_Spread(list_B, status_A, status_B, alpha_B, &next_B, covB, cnt_qba);
					oneRound_A_Spread(list_A, status_A, status_B, alpha_A, &next_A, covA, cnt_qab);
				}
				if (!list_A.empty())
					list_A.pop_front();
				if (!list_B.empty())
					list_B.pop_front();

			}//EDN-FOR

			curr_A = next_A;
			curr_B = next_B;
			next_A = next_B = 0;
			ticks++;

		}//END-WHILE	


	}//END-RUNS
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;
	//

	num_qab += cnt_qab / (double)MC_RUNS;
	num_qba += cnt_qba / (double)MC_RUNS;
	Acover += covA / (double)MC_RUNS;
	Bcover += covB / (double)MC_RUNS;
	return (covA + covB) / (double)MC_RUNS;
}
double MonteCarlo::compute_coverage_seq_Afirst_nodeCount(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, bool q1) {
	double qao = this->qao;
	if (q1)qao = 1;
	double  nodeCnt = 0;

	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int* status_B = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B


	bool SeedTestFlag = 1;
	;
	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);

		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		// scan all a-seeds
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && alpha_A[u] > qao) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			nodeCnt++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);
					if (status_A[v] != ADOPTED)
						//cout << "node=" << v << " statusA=" << status_A[v] << endl;
						list_A.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED);
				}
			}
		}
		/*adpption test in sequence*/
		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;

		while (curr_A > 0) {

			for (int i = 0; i < curr_A; i++) {

				int v = list_A.front();
				list_A.pop_front();

				if (status_A[v] == ADOPTED || status_A[v] == REJECTED)
					continue;

				if (alpha_A[v] <= qao) {
					status_A[v] = ADOPTED;  // A-adopted
					nodeCnt++;
				}
				else { status_A[v] = REJECTED; }// A-rejected

				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) { examine_out_neighbors(v, &list_A, &next_A, status_A); }
			} // ENDFOR
			curr_A = next_A; next_A = 0;
		}//END-WHILE

		// scan all B-seeds
		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);

			if (SeedTestFlag && ((alpha_B[u] > qbo) || (status_A[u] == ADOPTED && alpha_B[u] > qba))) {
				status_B[u] = REJECTED;
				continue;
			}
			if (status_A[u] != ADOPTED) nodeCnt++;

			status_B[u] = ADOPTED;

			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				int preEdgeStatus = graph->get_outedge_status(u, j);
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);  // edge is live
						if (status_B[v] != ADOPTED)
							list_B.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
					}
				}
				else if (preEdgeStatus == LIVE && status_B[v] != ADOPTED) {
					list_B.push_back(v);
				}

			}
		}
		
		curr_B = list_B.size(); next_B = 0;
		
		while (curr_B > 0) {

			for (int i = 0; i < curr_B; i++) {
				int v = list_B.front();
				list_B.pop_front();
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						nodeCnt++;
					}
					else { status_B[v] = REJECTED; }
				}//end qbo adoption


				else {
					if (alpha_B[v] <= qba) {
						status_B[v] = ADOPTED;
					}
					else { status_B[v] = REJECTED; }
				}//end qba adoption


				if (status_B[v] == ADOPTED) { examine_out_neighbors(v, &list_B, &next_B, status_B); }
			} // END-FOR
			curr_B = next_B; next_B = 0;
		} // END-WHILE
	}
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;

	return nodeCnt / (double)MC_RUNS;
}
double MonteCarlo::compute_coverage_seq_Afirst_upper(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2) {
	double qba = this->qbo;
	cout << "qba changed to " << qba << endl;
	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int* status_B = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	double covB_qbo = 0, covB_qba = 0; double covA_qao = 0, covA_reconsider = 0;

	bool SeedTestFlag = 1;
	;
	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);

		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		// scan all a-seeds
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && alpha_A[u] > qao) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			covA++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);
					if (status_A[v] != ADOPTED)
						list_A.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED);
				}
			}
		}
		/*adpption test in sequence*/
		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;

		while (curr_A > 0) {
			
			for (int i = 0; i < curr_A; i++) {

				int v = list_A.front();
				
				list_A.pop_front();
				if (status_A[v] == ADOPTED || status_A[v] == REJECTED)
					continue;


				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0

					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						covA++; covA_qao++;
					}
					else { status_A[v] = REJECTED; }// A-rejected
				}//end qao adoption

				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) { examine_out_neighbors(v, &list_A, &next_A, status_A); }
				
			} // ENDFOR
			curr_A = next_A; next_A = 0;
		}//END-WHILE

		// scan all B-seeds
		//for (auto it = bSeedSet.begin(); it != bSeedSet.end(); ++it) 
		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);

			if (SeedTestFlag && ((alpha_B[u] > qbo) || (status_A[u] == ADOPTED && alpha_B[u] > qba))) {
				status_B[u] = REJECTED;
				continue;
			}

			status_B[u] = ADOPTED;
			covB++;
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				int preEdgeStatus = graph->get_outedge_status(u, j);
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);  // edge is live
						if (status_B[v] != ADOPTED)
							list_B.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
					}
				}
				else if (preEdgeStatus == LIVE && status_B[v] != ADOPTED) {
					list_B.push_back(v);
				}

			}
		}
	
		curr_B = list_B.size(); next_B = 0;

		while (curr_B > 0) {
			
			for (int i = 0; i < curr_B; i++) {
				int v = list_B.front();
				list_B.pop_front();
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					//cout << v << "test with qbo=" << qbo << ", alphaB=" << alpha_B[v] << endl;
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						covB++; covB_qbo++;

					}
					else { status_B[v] = REJECTED; }
				}//end qbo adoption


				else {
					//cout << v << "test with qba=" << qba << ", alphaB=" << alpha_B[v] << endl;
					cnt_qba++;
					if (alpha_B[v] <= qba) {
						status_B[v] = ADOPTED;
						covB++; covB_qba++;
					}
					else { status_B[v] = REJECTED; }
				}//end qba adoption


				if (status_B[v] == ADOPTED) { examine_out_neighbors(v, &list_B, &next_B, status_B); }
			} // END-FOR
			curr_B = next_B; next_B = 0;
		} // END-WHILE
	}
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;
	return (covA + covB) / (double)MC_RUNS;
}


double MonteCarlo::compute_coverage_seq_Bfirst(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, double& Acover, double& Bcover, double& num_qab, double& num_qba) {
	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; 
	int* status_B = new int[graph->n]; 

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	double covB_qbo = 0, covB_qba = 0, covB_reconsider = 0;
	double covA_qao = 0, covA_qab = 0, covA_reconsider = 0;

	bool SeedTestFlag = 1;

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();


		// scan all B-seeds
		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);

			if (SeedTestFlag && alpha_B[u] > qbo) {
				status_B[u] = REJECTED;
				continue;
			}

			status_B[u] = ADOPTED;
			covB++;
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);  // edge is live
					if (status_B[v] != ADOPTED)
						list_B.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
				}


			}
		}
		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;
		
		while (curr_B > 0) {

			for (int i = 0; i < curr_B; i++) {

				int v = list_B.front();  list_B.pop_front();
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						covB++;covB_qbo++;
					}
					else {
						status_B[v] = REJECTED;
					}
				}

				if (status_B[v] == ADOPTED) {examine_out_neighbors(v, &list_B, &next_B, status_B);} // END-IF
			} // END-FOR
			curr_B = next_B; next_B = 0;
		} // END-WHILE

		// scan all A-seeds
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && ((alpha_A[u] > qao) || (status_B[u] == ADOPTED && alpha_A[u] > qab))) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			covA++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				int preEdgeStatus = graph->get_outedge_status(u, j);
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);
						if (status_A[v] != ADOPTED)
							list_A.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED);
					}
				}
				else if (preEdgeStatus == LIVE && status_A[v] != ADOPTED) {
					list_A.push_back(v);
				}
			}
		}

		curr_A = list_A.size();
		next_A = 0;
		while (curr_A > 0) {
			
			for (int i = 0; i < curr_A; i++) {
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == ADOPTED || status_A[v] == REJECTED)
					continue;


				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
					
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						covA++; covA_qao++;

					}
					else {status_A[v] = REJECTED;}
				}
				

				else {
					
					cnt_qab++;
					// v is already B-adopted, test with q_A|B
					if (alpha_A[v] <= qab) {
						status_A[v] = ADOPTED;
						covA++; covA_qab++;
					}
					else {
						status_A[v] = REJECTED;
					}
				}
				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) {examine_out_neighbors(v, &list_A, &next_A, status_A);} 
			} // ENDFOR
			curr_A = next_A; next_A = 0;
		}//END-WHILE
	

	}
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;

	num_qab += cnt_qab / (double)MC_RUNS;
	num_qba += cnt_qba / (double)MC_RUNS;
	Acover += covA / (double)MC_RUNS;
	Bcover += covB / (double)MC_RUNS;

	return (covA + covB) / (double)MC_RUNS;
}
double MonteCarlo::compute_coverage_seq_Bfirst_rhythm(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, int timing, double& Acover, double& Bcover, double& num_qab, double& num_qba) {
	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n];
	int* status_B = new int[graph->n];

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	double covB_qbo = 0, covB_qba = 0, covB_reconsider = 0;
	double covA_qao = 0, covA_qab = 0, covA_reconsider = 0;

	bool SeedTestFlag = 1;

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();


		/////////////////////////////////// scan all B-seeds////////////////////////////////////
		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);
			status_B[u] = ADOPTED;
			covB++;
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);  // edge is live
					if (status_B[v] != ADOPTED)
						list_B.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
				}

			}
		}

		///////////////////////////////////////////////////////////////////////////////////////

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;
		int ticks = 0, curr;
		bool Aready = false, scaned = false;

		while (curr_B > 0 || curr_A > 0 || !Aready) {

			Aready = ticks >= timing;
			if (Aready && !scaned) {
				///////////////////////// scan all A-seeds///////////////////////////////////
				for (int i = 0; i < k1; ++i) {
					int u = aSeedSet.at(i);

					if (SeedTestFlag && status_B[u] == ADOPTED && alpha_A[u] > qab) {
						status_A[u] = REJECTED;
						continue;
					}

					status_A[u] = ADOPTED;
					covA++;
					// iterate over its out-neighbors
					for (int j = 0; j < graph->outDeg[u]; j++) {
						int v = graph->gO[u][j];
						int preEdgeStatus = graph->get_outedge_status(u, j);
						if (preEdgeStatus == INACTIVE) {
							double coin = randomDouble();
							if (coin <= graph->probO[u][j]) {
								graph->set_outedge_status(u, j, LIVE);
								if (status_A[v] != ADOPTED)
									list_A.push_back(v);
							}
							else {
								graph->set_outedge_status(u, j, BLOCKED);
							}
						}
						else if (preEdgeStatus == LIVE && status_A[v] != ADOPTED) {
							list_A.push_back(v);
						}
					}
				}
				scaned = true;
				curr_A = list_A.size();

				////////////////////////////////////////////////////////////////////////////////
			}
			curr = !Aready ? curr_B : max(curr_A, curr_B);
			for (int i = 0; i < curr; i++) {

				double coin = randomDouble();
				if (coin <= 0.5) {
					oneRound_B_Spread(list_B, status_A, status_B, alpha_B, &next_B, covB, cnt_qba);
					if (Aready)
						oneRound_A_Spread(list_A, status_A, status_B, alpha_A, &next_A, covA, cnt_qab);

				}
				else {
					if (Aready)
						oneRound_A_Spread(list_A, status_A, status_B, alpha_A, &next_A, covA, cnt_qab);
					oneRound_B_Spread(list_B, status_A, status_B, alpha_B, &next_B, covB, cnt_qba);
				}
				if (!list_A.empty())
					list_A.pop_front();
				if (!list_B.empty())
					list_B.pop_front();

			}//EDN-FOR
			curr_A = next_A;
			curr_B = next_B;
			next_A = next_B = 0;
			ticks++;
		} // END-WHILE



	}//END-RUNS
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;

	cout << endl;

	num_qab += cnt_qab / (double)MC_RUNS;
	num_qba += cnt_qba / (double)MC_RUNS;
	Acover += covA / (double)MC_RUNS;
	Bcover += covB / (double)MC_RUNS;

	return (covA + covB) / (double)MC_RUNS;
}
double MonteCarlo::compute_coverage_seq_Bfirst_nodeCount(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, bool q1) {
	double qbo = this->qbo;
	if (q1) qbo = 1;
	double nodeCnt = 0;

	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n];
	int* status_B = new int[graph->n];

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B


	bool SeedTestFlag = 1;

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		//for (auto it = bSeedSet.begin(); it != bSeedSet.end(); ++it) 
		// scan all B-seeds
		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);

			if (SeedTestFlag && alpha_B[u] > qbo) {
				status_B[u] = REJECTED;
				continue;
			}

			status_B[u] = ADOPTED;
			nodeCnt++;

			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);  // edge is live
					if (status_B[v] != ADOPTED)
						list_B.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
				}

			}
		}

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;
		
		while (curr_B > 0) {
			
			for (int i = 0; i < curr_B; i++) {

				int v = list_B.front();  list_B.pop_front();
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;

				if (alpha_B[v] <= qbo) {
					status_B[v] = ADOPTED;
					nodeCnt++;
				}
				else {
					status_B[v] = REJECTED;
				}

				if (status_B[v] == ADOPTED) { examine_out_neighbors(v, &list_B, &next_B, status_B); } // END-IF
			} // END-FOR
			curr_B = next_B; next_B = 0;
		} // END-WHILE

		// scan all A-seeds
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && ((alpha_A[u] > qao) || (status_B[u] == ADOPTED && alpha_A[u] > qab))) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			if (status_B[u] != ADOPTED) nodeCnt++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				int preEdgeStatus = graph->get_outedge_status(u, j);
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);
						if (status_A[v] != ADOPTED)
							list_A.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED);
					}
				}
				else if (preEdgeStatus == LIVE && status_A[v] != ADOPTED) {
					list_A.push_back(v);
				}
			}
		}

		curr_A = list_A.size();
		next_A = 0;
		while (curr_A > 0) {
			
			for (int i = 0; i < curr_A; i++) {
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == ADOPTED || status_A[v] == REJECTED)
					continue;


				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
			
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						nodeCnt++;

					}
					else { status_A[v] = REJECTED; }
				}


				else {

					// v is already B-adopted, test with q_A|B
					if (alpha_A[v] <= qab) {
						status_A[v] = ADOPTED;
					}
					else {
						status_A[v] = REJECTED;
					}
				}
				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) { examine_out_neighbors(v, &list_A, &next_A, status_A); }
			} // ENDFOR
			curr_A = next_A; next_A = 0;
		}//END-WHILE

	}
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;



	return nodeCnt / (double)MC_RUNS;
}
double MonteCarlo::compute_coverage_seq_Bfirst_upper(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2) {
	int qab = this->qao;
	double  covA = 0, covB = 0;
	int cnt_qab = 0, cnt_qba = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n];
	int* status_B = new int[graph->n];

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	double covB_qbo = 0, covB_qba = 0, covB_reconsider = 0;
	double covA_qao = 0, covA_qab = 0, covA_reconsider = 0;

	bool SeedTestFlag = 1;

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();


		for (int i = 0; i < k2; i++) {
			int u = bSeedSet.at(i);

			if (SeedTestFlag && alpha_B[u] > qbo) {
				status_B[u] = REJECTED;
				continue;
			}

			status_B[u] = ADOPTED;
			covB++;
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);  // edge is live
					if (status_B[v] != ADOPTED)
						list_B.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
				}


			}
		}

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;
		//second adoption test
	
		while (curr_B > 0) {
		
			for (int i = 0; i < curr_B; i++) {

				int v = list_B.front();  list_B.pop_front();
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						covB++; covB_qbo++;
					}
					else {
						status_B[v] = REJECTED;
					}
				}

				if (status_B[v] == ADOPTED) { examine_out_neighbors(v, &list_B, &next_B, status_B); } // END-IF
			} // END-FOR
			curr_B = next_B; next_B = 0;
		} // END-WHILE

		// scan all A-seeds
		for (int i = 0; i < k1; ++i) {
			int u = aSeedSet.at(i);

			if (SeedTestFlag && ((alpha_A[u] > qao) || (status_B[u] == ADOPTED && alpha_A[u] > qab))) {
				status_A[u] = REJECTED;
				continue;
			}

			status_A[u] = ADOPTED;
			covA++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				int preEdgeStatus = graph->get_outedge_status(u, j);
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);
						if (status_A[v] != ADOPTED)
							list_A.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED);
					}
				}
				else if (preEdgeStatus == LIVE && status_A[v] != ADOPTED) {
					list_A.push_back(v);
				}
			}
		}

		curr_A = list_A.size();
		next_A = 0;
		while (curr_A > 0) {
			
			for (int i = 0; i < curr_A; i++) {
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == ADOPTED || status_A[v] == REJECTED)
					continue;

				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						covA++; covA_qao++;

					}
					else { status_A[v] = REJECTED; }
				}
				else {
					cnt_qab++;
					// v is already B-adopted, test with q_A|B
					if (alpha_A[v] <= qab) {
						status_A[v] = ADOPTED;
						covA++; covA_qab++;
					}
					else {
						status_A[v] = REJECTED;
					}
				}
				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) { examine_out_neighbors(v, &list_A, &next_A, status_A); }
			} // ENDFOR
			curr_A = next_A; next_A = 0;
		}//END-WHILE
	}
	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;
	return (covA + covB) / (double)MC_RUNS;
}

double MonteCarlo::compute_coverage(int* set_A, int size_A)
{
	double  cov = 0;
	double* alpha_A = new double[graph->n];
	double* alpha_B = new double[graph->n];
	int* status_A = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int* status_B = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++) {
			alpha_A[i] = randomDouble();
			alpha_B[i] = randomDouble();
		}
		graph->reset_outedge_status();

		// scan all A-seeds
		for (int i = 0; i < size_A; ++i) {
			//cout<<"u="<<set_A[i]<<endl;
			int u = set_A[i];
			status_A[u] = ADOPTED;
			cov++;
			// iterate over its out-neighbors
			for (int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = randomDouble();
				if (coin <= graph->probO[u][j]) {
					graph->set_outedge_status(u, j, LIVE);
					if (status_A[v] != ADOPTED)
						list_A.push_back(v);
				}
				else {
					graph->set_outedge_status(u, j, BLOCKED);
				}
			}
		}

		// scan all B-seeds  (skip if bSeeds is empty)
		for (auto it = bSeeds.begin(); it != bSeeds.end(); ++it) {
			int u = *it;
			status_B[u] = ADOPTED;
			
			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int preEdgeStatus = graph->get_outedge_status(u, j);
				int v = graph->gO[u][j];
				if (preEdgeStatus == INACTIVE) {
					double coin = randomDouble();
					if (coin <= graph->probO[u][j]) {
						graph->set_outedge_status(u, j, LIVE);  // edge is live
						if (status_B[v] != ADOPTED)
							list_B.push_back(v);
					}
					else {
						graph->set_outedge_status(u, j, BLOCKED); // edge is blocked
					}
				}
				else if (preEdgeStatus == LIVE && status_B[v] != ADOPTED) {
					list_B.push_back(v);
				}

			}
		}

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;

		while (curr_A > 0 || curr_B > 0) {
			// A-adoption test
			for (int i = 0; i < curr_A; i++) {
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == SUSPENDED || status_A[v] == ADOPTED)
					continue;

				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						cov++;

						if (status_B[v] == SUSPENDED && alpha_B[v] <= qba) { //reconsider only happen once for the same node.
							status_B[v] = ADOPTED; // reconsider to adopt B
							examine_out_neighbors(v, &list_B, &next_B, status_B);
						}
					}
					else {
						status_A[v] = SUSPENDED; // A-suspended
					}

				}
				else {
					// v is already B-adopted, test with q_A|B
					if (alpha_A[v] <= qab) {
						status_A[v] = ADOPTED;
						cov++;
					}
					else {
						status_A[v] = SUSPENDED;
					}
				}

				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) {
					examine_out_neighbors(v, &list_A, &next_A, status_A);
				} // END-IF
			} // ENDFOR

			// B adoption test
			for (int i = 0; i < curr_B; i++) {
				int v = list_B.front();
				list_B.pop_front();
				if (status_B[v] == SUSPENDED || status_B[v] == ADOPTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						if (status_A[v] == SUSPENDED && alpha_A[v] <= qab) {
							status_A[v] = ADOPTED; // reconsideration for A!
							cov++;
							examine_out_neighbors(v, &list_A, &next_A, status_A);
						}
					}
					else {
						status_B[v] = SUSPENDED;
					}

				}
				else {
					status_B[v] = (alpha_B[v] <= qba) ? ADOPTED : SUSPENDED; // already A-adopted
				}

				if (status_B[v] == ADOPTED) {
					examine_out_neighbors(v, &list_B, &next_B, status_B);
				} // END-IF
			} // END-FOR

			curr_A = next_A;
			curr_B = next_B;
			next_A = next_B = 0;

		} // END-WHILE
	}

	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;

	return cov / (double)MC_RUNS;
}

void MonteCarlo::examine_out_neighbors(int v, deque<int>* p_list, int* p_next, int* status)
{
	for (int j = 0; j < graph->outDeg[v]; j++) {
		int w = graph->gO[v][j];
		int preEdgeStatus = graph->get_outedge_status(v, j);
		int status_w = status[w];
		/*potential repeadted calculation*/
		if (preEdgeStatus == LIVE && status_w != ADOPTED && status_w != REJECTED && status_w != SUSPENDED) {
			p_list->push_back(w);
			(*p_next)++;
		}
		else if (preEdgeStatus == INACTIVE) {
			double coin = randomDouble();

			if (coin <= graph->probO[v][j]) {
				graph->set_outedge_status(v, j, LIVE);
				if (status_w != ADOPTED && status_w != REJECTED && status_w != SUSPENDED) {
					p_list->push_back(w);
					(*p_next)++;
				}
			}
			else {
				graph->set_outedge_status(v, j, BLOCKED);
			}
		} // ENDIF
	} // ENDFOR
}
void MonteCarlo::estimate_spread_from_file(string seed_file_name,int k)
{
	string aseeds = "dataset/A-SEED/aseeds_" + to_string(k) + ".txt";
	int* seedArr = new int[k];
	ifstream afile(aseeds.c_str(), ios::in);
	if (afile.is_open()) {
		int idx = 0;
		while (!afile.eof()) {
			string line;
			getline(afile, line);
			if (line == "")
				continue;
			int u = atoi(line.c_str());
			seedArr[idx++] = u;
		}
		afile.close();
		cout << "[info] A-seeds read from file: " << endl;
		for (int i = 0; i < k; i++)
			cout << seedArr[i] << " ";
		cout << endl;
	}
	else {
		cout << "[error] unable to open file: " << aseeds << endl;
	}

	int found = (int)seed_file_name.find_last_of(".");
	string ofile_name = seed_file_name.substr(0, found) + "_TRUE.txt";
	ofstream myfile;
	myfile.open(ofile_name.c_str());
	if (myfile.is_open()) {
		for (int i = 0; i < k; i++) {
			if (i == 0 || (i + 1) % 10 == 0) {
				double spread = compute_coverage(seedArr, i + 1);
				myfile << (i + 1) << "\t" << spread << endl;
				cout << (i + 1) << "\t" << spread << endl;
			}
		}
		myfile.close();
	}
	else {
		cerr << "Unable to open output file; no output was written in: " << ofile_name << endl;
		exit(1);
	}

	delete[] seedArr;
}


void MonteCarlo::write_seeds_to_file(bool isB)
{
	ofstream myfile;
	myfile.open(output_file_name.c_str());
	double spread = 0;
	if (myfile.is_open()) {
		for (int i = 0; i < aSeeds.size(); ++i) {
			spread += mg.at(i);
			if (!isB)
				myfile << i + 1 << "\t" << aSeeds.at(i) << "\t" << mg.at(i) << "\t" << spread << endl;
			else
				myfile << i + 1 << "\t" << bSeeds.at(i) << "\t" << mg.at(i) << "\t" << spread << endl;
		}
		myfile.close();
	}
	else {
		cerr << "[error] unable to open output file; no output was written " << output_file_name << endl;
		exit(1);
	}
}
void MonteCarlo::write_seeds_to_file(vector<int> seeds,string file_name) {
	ofstream myfile;
	myfile.open(file_name.c_str());
	if (myfile.is_open()) {
		for (int i = 0; i < seeds.size(); ++i) {
			myfile << seeds.at(i) << endl;

		}
		myfile.close();
	}
	else {
		cerr << "[error] unable to open output file; no output was written " << file_name << endl;
		exit(1);
	}
}


