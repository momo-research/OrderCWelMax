#pragma once
#ifndef __MC_H
#define __MC_H

#include "graph.h"
#include <fstream>
#include <sstream>
#include <algorithm>

class MonteCarlo // : public Graph
{
public:
	vector<int>     aSeeds; // seed set to be found (for item A)
	vector<int>     bSeeds; // input: the other company's seed set
	vector<double>  mg;  	// marginal gain of each seed added (in greedy sequence)
	/*
	vector<double>  qao;
	vector<double>  qab;
	vector<double>  qbo;
	vector<double>  qba;
	*/
	double qao;
	double qab;
	double qbo;
	double qba;;
	
	double          reconsider_a;
	double          reconsider_b;
	double			qa_aware;
	double			qb_aware;
	int             k;
	bool            ignore_B;
	string          output_file_name;
	string 			b_seeds_file_name;
	string			a_seeds_file_name;

	Graph* graph;
	int n;
	int m;

	MonteCarlo(string folder, string graph_file)
	{
		ignore_B = false;
		k = 0;
		output_file_name = folder + "/output/default_output.txt";

		graph = new Graph(folder, graph_file);
		this->n = graph->n;
		this->m = graph->n;
	}

	~MonteCarlo()
	{
		delete graph;
	}

	void setParametersMC(int k_A, vector<double> qq, bool ignore, string folder, string bseeds);
	void readBSeedsMC();
	void readBSeedsMC(string str);
	void readASeedsMC(string str);

	double mineSeedsMC();
	double compute_coverage(int* set_A, int size_A);
	/* if v becomes X-adopted (X = A or X = B), we examine v's out-neighbours */
	void examine_out_neighbors(int v, deque<int>* p_list, int* p_next, int* status);

	void estimate_spread_from_file(string seed_file_name, int k);
	void write_seeds_to_file(bool isB);
	void write_seeds_to_file(vector<int> seeds, string file_name);
	
	//for sequence


	void mineABseed_Greedy_Afirst(int kA, int kB, string file_SA, string file_SB);
	void mineABseed_Greedy_Bfirst(int kA, int kB, string file_SA, string file_SB);
	void adjustHeap(int* heap, double* improve, int i);

	int oneRound_A_Spread(deque<int>& list_A, int* status_A, int* status_B, double* alpha_A, int* next_A, double& covA, int& cnt_qab);
	int oneRound_B_Spread(deque<int>& list_B, int* status_A, int* status_B, double* alpha_B, int* next_B, double& covB, int& cnt_qba);

	
	double compute_coverage_simu_v1(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, double& Acover, double& Bcover, double& num_qab, double& num_qba2);
	

	double compute_coverage_seq_Afirst(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, double& Acover, double& Bcover, double& num_qab, double& num_qba);
	double compute_coverage_seq_Afirst_rhythm(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, int timing, double& Acover, double& Bcover, double& num_qab, double& num_qba);
	double compute_coverage_seq_Afirst_nodeCount(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, bool q1);
	double compute_coverage_seq_Afirst_upper(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2);

	double compute_coverage_seq_Bfirst(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, double& Acover, double& Bcover, double& num_qab, double& num_qba);
	double compute_coverage_seq_Bfirst_rhythm(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, int timing, double& Acover, double& Bcover, double& num_qab, double& num_qba);
	double compute_coverage_seq_Bfirst_nodeCount(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2, bool q1);
	double compute_coverage_seq_Bfirst_upper(vector<int> aSeedSet, int k1, vector<int>bSeedSet, int k2);
};


#endif
