#pragma once
#include <vector>
#include <map>
#include <string>
#include <array>
#include <assert.h>

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
		const std::string &b, const int GAP_COST, const int START_GAP, const int cost_matrix[21][21], int sign);

// Construct optimal alignment score matrix
std::vector<std::vector<std::vector<int>>>
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, bool dir);

// Find the sub-graph of the alignment paths with (sub-)optimal paths
template<class K>
Dag gen_dag(const std::string &a, const std::string &b, const int cost_matrix[21][21],
		const K TH, const int GAP_COST, const int START_GAP) {
	std::vector<std::vector<std::vector<std::vector<Node>>>> e = build_dp_matrix(a, b, GAP_COST, START_GAP, cost_matrix, 1);
	int n = (int) e.size() - 1;
	assert(n > 0);
	int m = (int) e[0].size() - 1;
	assert(m > 0);

	// find smallest distances for each node
	std::vector<std::vector<std::vector<std::vector<Node>>>> er(n + 1, std::vector<std::vector<std::vector<Node>>>(m + 1, std::vector<std::vector<Node>>(3)));
	for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int k = 0; k <= 2; k++) {
		for (Node nxt: e[i][j][k])
			er[nxt.N_index][nxt.M_index][nxt.type].push_back(Node(i, j, k, nxt.cost));
	}

	std::vector<std::vector<std::vector<int>>> dp = opt_alignment(e, 0);
	std::vector<std::vector<std::vector<int>>> dpr = opt_alignment(er, 1);
	assert(dp[n][m][0] == dpr[0][0][0]);

	const int OPT = dp[n][m][0];
	int WORST;
	if (TH != 0) {
		// find the largest distance
		std::vector<std::vector<std::vector<std::vector<Node>>>> el = build_dp_matrix(a, b, GAP_COST, START_GAP, cost_matrix, -1);
		std::vector<std::vector<std::vector<int>>> dpl = opt_alignment(el, 0);
		WORST = (-1) * dpl[n][m][0];
	} else {
		WORST = OPT;
	}

	assert(OPT <= WORST);
	int current = 0;
	std::vector<std::vector<int>> adj(1);
	std::map<std::pair<int, int>, std::array<int, 3>> trans; // translate to index
	std::map<int, std::pair<int, int>> transr; // translate index to pair
	trans[std::make_pair(0, 0)] = {0, -1, -1};
	transr[0] = std::make_pair(0, 0);

	auto add_node = [&](const Node &node) {
		if (trans.find(std::make_pair(node.N_index, node.M_index)) == trans.end()) {
			trans[std::make_pair(node.N_index, node.M_index)] = {-1, -1, -1};
		}
		if (trans[std::make_pair(node.N_index, node.M_index)][node.type] == -1) {
			trans[std::make_pair(node.N_index, node.M_index)][node.type] = ++current;
			transr[current] = std::make_pair(node.N_index, node.M_index);
			adj.push_back(std::vector<int>());
		}
	};

	// TH: v 0%      v 50%     v 100%
	//    [OPT, .........., WORST]
	auto check_th = [&](const int k, const int OPT, const int WORST, K TH) {
		//return k == OPT;
		return k <= OPT + TH * (WORST - OPT);
	};

	for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int k = 0; k <= 2; k++) {
		if (!check_th(dp[i][j][k] + dpr[i][j][k], OPT, WORST, TH)) continue;
		add_node(Node(i, j, k, 0));
		for (const Node &nxt: e[i][j][k]) {
			if (!check_th(dp[nxt.N_index][nxt.M_index][nxt.type] + dpr[nxt.N_index][nxt.M_index][nxt.type], OPT, WORST, TH)) continue;
			if (check_th(dpr[nxt.N_index][nxt.M_index][nxt.type] + dp[i][j][k] + nxt.cost, OPT, WORST, TH)) {
				add_node(nxt);
				adj[trans[std::make_pair(i, j)][k]].push_back(trans[std::make_pair(nxt.N_index, nxt.M_index)][nxt.type]);
			}
		}
	}

	return { adj, 0, trans[std::make_pair(n, m)][0], trans, transr };
}

