#pragma once
#ifndef COMIC_INFGRAPH_H
#define COMIC_INFGRAPH_H

#include "graph.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <deque>
#include <cstring>
#include <stack>
#include "random_utils.h"
using namespace std;

//Node for RR tree
struct RNode {
public:
	int vex;
	RNode* parent;  // 弱引用，指向父节点
	vector<unique_ptr<RNode>> children;  // 子节点由 unique_ptr 管理生命周期

	RNode(int v) : vex(v), parent(nullptr) {}
	RNode() : vex(-1), parent(nullptr) {}
};

class InfGraph //: public Graph
{
public:
	vector<vector<int>> hyperG;
	vector<vector<int>> hyperGT;  // all RR sets
	int64 hyperId;
	vector<int> seedSet; // seed set to be mined (A-seeds for Self-Inf-Max and B-seeds for Comp-Inf-Max)
	vector<int> aSeeds;  // S_A from input files (for Comp-Inf-Max)
	vector<int> bSeeds;  // S_B from input files (for Self-Inf-Max)

	double qao;
	double qab;
	double qbo;
	double qba;
	int k;
	int kA;
	int kB;
	int64 theta;
	//double* alpha_A;
	//double* alpha_B;
	int* status_A;
	int* status_B;
	deque<int> Q;
	//bool *drawn; // TRUE if its alpha values are drawn in RR-set generation process
	//bool* visited;  // TRUE if being visited in the backward BFS
	string b_seeds_file_name;
	string output_file_name;
	string dataset;
	bool ignore_B;
	bool* used; // used[v] is true if v has been added to RR-set
	bool* discovered;

	Graph* graph;
	int n, m;
	
	///////////////////////////////////////////////////////////////并行区///////////////////////////////////////////////////
	mutex lockA,lockB;

	thread_local static vector<int> visited;
	thread_local static int time_stamp;

	void init_visit() {
		visited.clear();
		time_stamp = 1;
		visited.resize(n, 0);
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
	
	
	vector<vector<int>> hyperG_A, hyperG_B;
	vector<vector<int>> hyperGT_A, hyperGT_B;
	//tree form of RR set
	vector<unique_ptr<RNode>>RRTree_A, RRTree_B;// size = n,store the root of each RR tree
	//pointer for each node in it's covered RR sets
	//vector<unordered_set<RPointer, RPointerHash>> RRP_A, RRP_B;
	vector<unordered_map<int, RNode*>> RRP_A, RRP_B;

	//tested mark
	vector<unordered_set<int>> tested;
	vector<unordered_set<int>> examined;




	vector<int> degree_A, degree_B;

	InfGraph(string folder, string graph_file, int64 theta) //:Graph(folder, graph_file) 
	{
		dataset = folder;
		graph = new Graph(folder, graph_file);
		this->n = graph->n;
		this->m = graph->m;
		this->theta = theta;

		init_infgraph();

		/*
		hyperG.clear();
		for (int i = 0; i < n; i++)
			hyperG.push_back(vector<int>());
		*/

		//alpha_A = new double[n];
		//alpha_B = new double[n];
		//status_A = new int[n];
		//status_B = new int[n];
		//visited = new bool[n];
		//used = new bool[n];
		//discovered = new bool[n];
	}
	void init_infgraph() {
		hyperG_A.reserve(n);
		hyperG_B.reserve(n);
		degree_A.reserve(n);
		degree_B.reserve(n);
		for (int i = 0; i < n; i++) {
			hyperG_A.push_back(vector<int>());
			hyperG_B.push_back(vector<int>());
		}
		hyperGT_A.reserve(theta);
		hyperGT_B.reserve(theta);
		for (int i = 0; i < theta; i++) {
			hyperGT_A.push_back(vector<int>());
			hyperGT_B.push_back(vector<int>());
		}
	}
	void init_order() {
		hyperId = theta;
		RRTree_A.resize(theta);
		RRTree_B.resize(theta);
		RRP_A.reserve(theta);
		RRP_B.reserve(theta);
		tested.reserve(theta);
		examined.reserve(theta);
		for (int i = 0; i < theta; i++) {

			RRP_A.push_back(unordered_map<int, RNode*>());
			RRP_B.push_back(unordered_map<int, RNode*>());
			tested.push_back(unordered_set<int>());
			examined.push_back(unordered_set<int>());
		}
	}
	void reset_hyper() {
		// 清空 hyperG_A 和 hyperG_B（array form）
		for (auto& vec : hyperG_A) {
			vec.clear(); // 清除每个节点的 RR-set 列表
		}
		for (auto& vec : hyperG_B) {
			vec.clear(); // 清除每个节点的 RR-set 列表
		}
		// 清空 hyperGT_A 和 hyperGT_B（reverse RR sets）
		for (auto& vec : hyperGT_A) {
			vec.clear(); // 清除每个 RR-set 的节点列表
		}
		for (auto& vec : hyperGT_B) {
			vec.clear(); // 清除每个 RR-set 的节点列表
		}
		hyperGT_A.resize(theta);
		hyperGT_B.resize(theta);
	}
	void reset_order() {


		// 重置 RR 树 (tree form)
		for (auto& root : RRTree_A) {
			root.reset(); // 释放 unique_ptr 管理的内存
		}
		for (auto& root : RRTree_B) {
			root.reset(); // 释放 unique_ptr 管理的内存
		}
		

		// 清空 tested 和 examined 集合
		for (auto& set : tested) {
			set.clear(); // 清除每个 RR-set 的测试标记
		}
		for (auto& set : examined) {
			set.clear(); // 清除每个 RR-set 的检查标记
		}

		// 清空 RRP_A 和 RRP_B（指针容器）
		for (auto& map : RRP_A) {
			map.clear(); // 清除每个节点的指针映射
		}
		for (auto& map : RRP_B) {
			map.clear(); // 清除每个节点的指针映射
		}
	}
	~InfGraph()
	{
		//delete[] alpha_A;
		//delete[] alpha_B;
		//delete[] status_A;
		//delete[] status_B;
		//delete[] visited;
		//delete[] used;
		//delete[] discovered;
		delete graph;

	}

	void setParametersIG(int kA, int kB, vector<double> qq, bool ignore, string folder, string bseeds);
	double computeInfluenceHyperGraph();  // compute influence
	void readBSeedsIG();
	void printSeeds();
	
	double calculateRepeatRate(const vector<int>& aSeeds, const vector<int>& bSeeds);
;
	//for theory
	void BuildHyperGraphR_idp();
	void BuildHyperGraphR_q1();
	int BuildSingleRRSet_q1(int uStart, int rrSetID, bool isrA, bool q1, vector<vector<int>>& hyperGT, vector<vector<int>>& hyperG);
	int BuildSingleRRSet_idp(int uStart, int rrSetID, bool isrA, bool q1, vector<vector<int>>& hyperGT, vector<vector<int>>& hyperG);
	void BuildSeedSetIG_q1(bool output);
	void BuildSeedSetIG_idp();
	double computeInfHyperGraph_Order(int theta);
	//for Order
	void BuildHypergraphR_Order(int run,bool Afirst);
	void BuildSeedSetIG_Order_Afirst();
	void BuildSeedSetIG_Order_Bfirst();
	unique_ptr<RNode> BuildSingleARRSet_Order(int uStart, int rrSetID, bool Afirst);
	unique_ptr<RNode> BuildSingleBRRSet_Order(int uStart, int rrSetID, bool Afirst);
	void printTree(RNode* root);
	bool pruningRB_01(RNode* root, RNode* aseed, int rrsetID);
	bool pruningRB_02(RNode* root, RNode* aseed, int rrsetID, unordered_set<int>& delay);
	bool pruningRA_01(RNode* root, RNode* bseed, int rrsetID);
	bool pruningRA_02(RNode* root, RNode* bseed, int rrsetID, unordered_set<int>& delay);

	void printQQ();
};
#endif
