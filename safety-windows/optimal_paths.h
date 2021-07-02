#pragma once
#include <vector>
#include <map>
#include <string>
#include <array>

#include <gmpxx.h>

// Dag here is supposed to be the subgraph constructed with the gen_dag function
struct Dag {
	std::vector<std::vector<int>> adj;
	int src, sink;
	std::map<std::pair<int, int>, std::array<int, 3>> trans;
	std::map<int, std::pair<int, int>> transr;
};

struct Node {
	int N_index, M_index, type;
	int cost;
	Node (int N_index, int M_index, int type, int cost) : N_index(N_index),
			M_index(M_index), type(type), cost(cost)
	{}
};

// Just plain dijkstra
std::vector<std::vector<std::vector<int>>>
dijkstra(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj,
		const int sn, const int sm);

// Construct alignment paths
std::vector<std::vector<std::vector<std::vector<Node>>>> build_dp_matrix(const std::string &a,
		const std::string &b, const int GAP_COST, const int START_GAP);

// Construct optimal alignment score matrix
std::vector<std::vector<std::vector<int>>>
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, int sn, int sm);

// Find the sub-graph of the alignment paths with (sub-)optimal paths
Dag gen_dag(const std::string &a, const std::string &b, const mpq_class TH = 1.0,
		const int GAP_COST = 1, const int START_GAP = 11); // BLOSUM62 costs
