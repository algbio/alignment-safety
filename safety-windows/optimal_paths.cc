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
const int LTA[26] = { 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20 };

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
		// Those errors shouldn't appear anymore
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
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, bool dir) {
	int n = (int) adj.size() - 1;
	assert(n > 0);
	int m = (int) adj[0].size() - 1;
	assert(m > 0);
//	return dijkstra(adj, (dir ? n : 0), (dir ? m : 0));

	std::vector<std::vector<std::vector<int>>> d(n + 1, std::vector<std::vector<int>>(m + 1, std::vector<int>(3, (1 << 30))));
	(dir ? d[n][m][0] : d[0][0][0]) = 0;
	auto update_dist = [&](int i, int j, int k) {
		for (const Node &node: adj[i][j][k]) {
			d[node.N_index][node.M_index][node.type] = std::min(d[node.N_index][node.M_index][node.type], d[i][j][k] + node.cost);
		}
	};
	if (!dir) {
		for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int k = 2; k >= 0; k--)
			update_dist(i, j, k);
	} else {
		for (int i = n; i >= 0; i--) for (int j = m; j >= 0; j--) for (int k = 0; k <= 2; k++)
			update_dist(i, j, k);
	}
	return d;
}
