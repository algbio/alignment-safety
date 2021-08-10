#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <string>
#include <assert.h>
#include <array>

#include <iostream>

#include "optimal_paths.h"

// translate fasta file letters to amino acid symbols (see http://www.math.utep.edu/Faculty/mleung/bioinformatics/aacodon.html)
//                    A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z
const int LTA[26] = { 0, 20, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, 20, 18, 20 };

std::vector<std::vector<std::vector<int>>>
dijkstra(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj,
		const int sn, const int sm) {
	int n = (int) adj.size() - 1;
	assert(n > 0);
	int m = (int) adj[0].size() - 1;
	assert(m > 0);

	auto cmp = [&](const Node &a, const Node &b) {
		if (a.cost > b.cost) return true;
		return false;
	};

	std::priority_queue<Node, std::vector<Node>, decltype(cmp)> q(cmp);
	q.push(Node(sn, sm, 0, 0));

	std::vector<std::vector<std::vector<int>>> dist(n + 1,
		std::vector<std::vector<int>>(m + 1, std::vector<int>(3, (1 << 30))));
	dist[sn][sm][0] = 0;

	while (!q.empty()) {
		Node cur = q.top();
		q.pop();

		if (cur.cost > (1 << 30) || cur.cost < - (1 << 29)) {
			std::cerr << "ERROR. COST TOO HIGH" << std::endl;
		}

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
		const std::string &b, const int GAP_COST, const int START_GAP, const int cost_matrix[21][21],
		int sign) {
	int n = (int) a.size();
	int m = (int) b.size();

	// create affine linear gap cost graph (see README for illustration)
	std::vector<std::vector<std::vector<std::vector<Node>>>> adj(n + 1,
			std::vector<std::vector<std::vector<Node>>>(m + 1,
			std::vector<std::vector<Node>>(3)));
	for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
		if (i < n) {
			if (a[i] < 'A' || a[i] > 'Z') {
				std::cerr << "ERROR: INVALID CHARACTER in the reference string: " << a[i] << std::endl;
			}
			if (LTA[a[i] - 'A'] == -1) {
				std::cerr << "ERROR: WRONG CHARACTER in the reference string: " << a[i] << std::endl;
			}
		}
		if (j < m) {
			if (b[j] < 'A' || b[j] > 'Z') {
				std::cerr << "ERROR: INVALID CHARACTER in the comparing string: " << b[j] << std::endl;
			}
			if (LTA[b[j] - 'A'] == -1) {
				std::cerr << "ERROR: WRONG CHARACTER in the comparing string: " << b[j] << std::endl;
			}
		}
		if (i + 1 <= n && j + 1 <= m)
			adj[i][j][0].push_back(Node(i + 1, j + 1, 0,
					sign * cost_matrix[LTA[a[i] - 'A']][LTA[b[j] - 'A']]));
					//a[i] != b[j]));

		if (i + 1 <= n) {
			adj[i][j][0].push_back(Node(i + 1, j, 1, sign * (START_GAP + GAP_COST)));
			adj[i][j][1].push_back(Node(i + 1, j, 1, sign * GAP_COST));
		}

		if (j + 1 <= m) {
			adj[i][j][0].push_back(Node(i, j + 1, 2, sign * (START_GAP + GAP_COST)));
			adj[i][j][2].push_back(Node(i, j + 1, 2, sign * GAP_COST));
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

Dag gen_dag(const std::string &a, const std::string &b, const int cost_matrix[21][21],
		const mpq_class TH, const int GAP_COST, const int START_GAP) {
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

	std::vector<std::vector<std::vector<int>>> dp = opt_alignment(e, 0, 0);
	std::vector<std::vector<std::vector<int>>> dpr = opt_alignment(er, n, m);
	assert(dp[n][m][0] == dpr[0][0][0]);

	const int OPT = dp[n][m][0];
	int WORST;
	if (TH != 0) {
		// find the largest distance
		std::vector<std::vector<std::vector<std::vector<Node>>>> el = build_dp_matrix(a, b, GAP_COST, START_GAP, cost_matrix, -1);
		std::vector<std::vector<std::vector<int>>> dpl = opt_alignment(el, 0, 0);
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
	auto check_th = [&](const int k, const int OPT, const int WORST, mpq_class TH) {
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
