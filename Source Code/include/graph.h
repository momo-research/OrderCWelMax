#pragma once
#ifndef COMIC_GRAPH_H
#define COMIC_GRAPH_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include "head.h"
#include "random_utils.h"
using namespace std;

typedef double (*pf)(int, int);

struct EdgeInfo {
    int status;    // Current status (0: INACTIVE, 1: LIVE, 2: BLOCKED)
    int version;   // Current version number
    struct EdgeInfo() :status(INACTIVE), version(0) {}
};

class Graph
{
public:
    int n, m, k;
    //int A_size, B_size, AB_size, NoneAB_size;
    //double A_ratio, B_ratio;
    vector<int> inDeg;
    vector<int> outDeg;
    vector<vector<int>> gT; // Transpose graph, maintaining in-neighbors
    vector<vector<double>> probT;
    vector<vector<int>> gO; // Non-transpose graph, maintaining out-neighbors
    vector<vector<double>> probO;
    vector<vector<EdgeInfo>>  outEdgeStatus;
    //unordered_map<int, int> edgeStatusMap;
    vector<vector<EdgeInfo>> inEdgeStatus;
    long long in_version;
    long long out_version;

    /////////////////////////////////////////////////////////////// Parallel Section ///////////////////////////////////////////////////
    thread_local static vector<vector<EdgeInfo>> thEdgeStatus;
    thread_local static long long th_version;

    void init_EdgeStatus() {
        thEdgeStatus.clear();
        th_version = in_version + 1;
        thEdgeStatus.reserve(inEdgeStatus.size());
        for (const auto& edges : inEdgeStatus) {
            thEdgeStatus.emplace_back(edges); // Directly copy each row
        }
    }
    inline void reset_thEdgeStatus() {
        th_version++;
    }
    inline int get_thEdgeStatus(int u, int j) {
        auto& edge = thEdgeStatus[u][j];
        if (edge.version < th_version) {
            edge.status = INACTIVE;
            edge.version = th_version;
        }
        return edge.status;
    }
    inline void set_thEdgeStatus(int u, int j, int status) {
        auto& edge = thEdgeStatus[u][j];
        edge.status = status;
        edge.version = th_version;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    string folder;
    string graph_file;

    Graph(string folder, string graph_file) :folder(folder), graph_file(graph_file) {
        // Read # nodes and # edges
        readNM();
        // Initialize vectors
        for (int i = 0; i < n; i++) {
            gT.push_back(vector<int>());  // Transpose graph
            probT.push_back(vector<double>()); // Incoming edge weight
            inEdgeStatus.push_back(vector<EdgeInfo>());
            inDeg.push_back(0); // In-degree of each node

            gO.push_back(vector<int>()); // Normal graph
            probO.push_back(vector<double>()); // Outgoing edge weight
            outEdgeStatus.push_back(vector<EdgeInfo>());
            outDeg.push_back(0); // Out-degree of each node

            //group.push_back(-1);
        }
        in_version = 0;
        out_version = 0;
        // Read edge list from file
        cout << "Begin to read graph_v2..." << endl;
        readGraph_v2();
    }

    ~Graph() {}

    void readNM();
    void add_in_edge(int a, int b, double p);  // Add edge (a,b): a as b's in-neighbour
    void add_out_edge(int a, int b, double p); // Add edge (a,b): b as a's out-neighbour
    void readGraph();
    void readGraph_v2();
    void printGraph(int num_nodes);
    vector<int> findHighDegree(int k);
    vector<int> findPageRank(int k);

    inline int cantor(int x, int y) {
        return ((x + y) * (x + y + 1)) / 2 + y;
    }

    inline vector<int> inv_cantor(int z) {
        vector<int> ret;
        int w = (int)floor((sqrt(8 * z + 1) - 1) / 2);
        int t = (w * w + w) / 2;
        int y = z - t;
        int x = w - y;
        ret.push_back(x);
        ret.push_back(y);
        return ret;
    }

    inline void reset_inedge_status() {
        in_version++; // Increment global version number
    }

    inline int get_inedge_status(int u, int j) {
        auto& edge = inEdgeStatus[u][j];
        if (edge.version < in_version) {
            edge.status = INACTIVE;       // Default status is INACTIVE
            edge.version = in_version;
        }
        return edge.status;
    }

    inline void set_inedge_status(int u, int j, int status) {
        auto& edge = inEdgeStatus[u][j];
        edge.status = status;
        edge.version = in_version;
    }

    inline void reset_outedge_status() {
        out_version++; // Increment global version number
    }

    inline int get_outedge_status(int u, int j) {
        auto& edge = outEdgeStatus[u][j];
        if (edge.version < out_version) {
            edge.status = INACTIVE;       // Default status is INACTIVE
            edge.version = out_version;
        }
        return edge.status;
    }

    inline void set_outedge_status(int u, int j, int status) {
        auto& edge = outEdgeStatus[u][j];
        edge.status = status;
        edge.version = out_version;
    }

    /*
    inline void reset_out_edge_status() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < outDeg[i]; j++) {
                outEdgeStatus[i][j] = INACTIVE;
            }
        }
    }*/

    double calculateRepeatRate(const vector<int>& aSeeds, const vector<int>& bSeeds) {
        int count = 0;
        for (int i = 0; i < aSeeds.size(); i++) {
            if (find(bSeeds.begin(), bSeeds.end(), aSeeds[i]) != bSeeds.end()) {
                count++;
            }
        }
        return static_cast<double>(count) / aSeeds.size();
    }

};
#endif