#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <string>
#include <assert.h>
#include <array>

#include <iostream>

#include "optimal_paths.h"

std::vector<std::vector<std::vector<int>>>
dijkstra(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj,
		const int sn, const int sm) {
	int n = (int) adj.size() - 1;
	assert(n > 0);
	int m = (int) adj[0].size() - 1;
	assert(m > 0);

	auto cmp = [&](const Node &a, const Node &b) {
		if (a.cost > b.cost) return true;
		if (a.cost < b.cost) return false;

		// compare the rest (TODO: can this be done more elegantly?)
		std::tuple<int, int, int> A = { a.N_index, a.M_index, a.type };
		std::tuple<int, int, int> B = { b.N_index, b.M_index, b.type };
		std::vector<std::tuple<int, int, int>> tmp = { A, B };
		sort(tmp.begin(), tmp.end(), std::greater<std::tuple<int, int, int>>());
		return tmp[0] == A;
	};

	std::priority_queue<Node, std::vector<Node>, decltype(cmp)> q(cmp);
	q.push(Node(sn, sm, 0, 0));

	std::vector<std::vector<std::vector<int>>> dist(n + 1,
		std::vector<std::vector<int>>(m + 1, std::vector<int>(3, (1 << 30))));
	dist[sn][sm][0] = 0;

	while (!q.empty()) {
		Node cur = q.top();
		q.pop();

		if (cur.cost > dist[cur.N_index][cur.M_index][cur.type])
			continue;

		for (Node nxt: adj[cur.N_index][cur.M_index][cur.type]) {
			nxt.cost += cur.cost;
			if (nxt.cost < dist[nxt.N_index][nxt.M_index][nxt.type]) {
				dist[nxt.N_index][nxt.M_index][nxt.type] = nxt.cost;
				q.push(nxt);
			}
		}
	}
	return dist;
}

std::vector<std::vector<std::vector<std::vector<Node>>>> build_dp_matrix(const std::string &a,
		const std::string &b) {
	int n = (int) a.size();
	int m = (int) b.size();

	const int SAME_COST = 0;
	const int DIFF_COST = 1;
	const int GAP_COST = 1;
	const int START_GAP = 0;

	// create affine linear gap cost graph (see README for illustration)
	std::vector<std::vector<std::vector<std::vector<Node>>>> adj(n + 1,
			std::vector<std::vector<std::vector<Node>>>(m + 1,
			std::vector<std::vector<Node>>(3)));
	for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
		if (i + 1 <= n && j + 1 <= m)
			adj[i][j][0].push_back(Node(i + 1, j + 1, 0,
					(a[i] == b[j] ? SAME_COST : DIFF_COST)));

		if (i + 1 <= n) {
			adj[i][j][0].push_back(Node(i + 1, j, 1, START_GAP + GAP_COST));
			adj[i][j][1].push_back(Node(i + 1, j, 1, GAP_COST));
		}

		if (j + 1 <= m) {
			adj[i][j][0].push_back(Node(i, j + 1, 2, START_GAP + GAP_COST));
			adj[i][j][2].push_back(Node(i, j + 1, 2, GAP_COST));
		}

		adj[i][j][1].push_back(Node(i, j, 0, 0));
		adj[i][j][2].push_back(Node(i, j, 0, 0));
	}

	return adj;
}

std::vector<std::vector<std::vector<int>>>
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, int sn, int sm) {
	return dijkstra(adj, sn, sm);
}

Dag gen_dag(const std::string &a, const std::string &b) {
	std::vector<std::vector<std::vector<std::vector<Node>>>> e = build_dp_matrix(a, b);
	int n = (int) e.size() - 1;
	assert(n > 0);
	int m = (int) e[0].size() - 1;
	assert(m > 0);
	std::vector<std::vector<std::vector<std::vector<Node>>>> er(n + 1, std::vector<std::vector<std::vector<Node>>>(m + 1, std::vector<std::vector<Node>>(3)));
	for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int k = 0; k <= 2; k++) {
		for (Node nxt: e[i][j][k])
			er[nxt.N_index][nxt.M_index][nxt.type].push_back(Node(i, j, k, nxt.cost));
	}

	std::vector<std::vector<std::vector<int>>> dp = opt_alignment(e, 0, 0);
	std::vector<std::vector<std::vector<int>>> dpr = opt_alignment(er, n, m);
	assert(dp[n][m][0] == dpr[0][0][0]);

	const int OPT = dp[n][m][0];
	std::vector<std::vector<int>> res;
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

	for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int k = 0; k <= 2; k++) {
		if (dp[i][j][k] + dpr[i][j][k] != OPT) continue;
		add_node(Node(i, j, k, 0));
		for (const Node &nxt: e[i][j][k]) {
			if (dp[nxt.N_index][nxt.M_index][nxt.type] + dpr[nxt.N_index][nxt.M_index][nxt.type] != OPT) continue;
			if (dp[nxt.N_index][nxt.M_index][nxt.type] == dp[i][j][k] + nxt.cost) {
				add_node(nxt);
				adj[trans[std::make_pair(i, j)][k]].push_back(trans[std::make_pair(nxt.N_index, nxt.M_index)][nxt.type]);
			}
		}
	}

	return { adj, 0, trans[std::make_pair(n, m)][0], trans, transr };
}
