#ifndef COMIC_TIMGRAPH_H
#define COMIC_TIMGRAPH_H

#include "infgraph.h"
#include <random>

class TimGraph : public InfGraph {
public:
	string output_file_name;

	TimGraph(string folder, string graph_file, int64 theta) :InfGraph(folder, graph_file,theta) {}
	~TimGraph() {}

	//for Order
	vector<vector<vector<int>>> mineSeeds_Order_Afirst(int runs, double epsilon, string file_SA, string file_SB);
	vector<vector<vector<int>>> mineSeeds_Order_Bfirst(int runs, double epsilon, string file_SA, string file_SB);
	void minseeds_Order_q1(double epsilon, string file_SA, string file_SB);
	void minseeds_idp(double epsilon, string file_SA, string file_SB);
	double LambdaPrime(double epsprime, double ell, int n);
	double LambdaStar(double eps, double ell);
	void EstimateTheta(double eps, double ell);

	inline double logcnk(int n, int k) {
		double ans = 0;
		for (int i = n - k + 1; i <= n; i++) {
			ans += log(i);
		}
		for (int i = 1; i <= k; i++) {
			ans -= log(i);
		}
		return ans;
	}

};

#endif //COMIC_TIMGRAPH_H
#pragma once
