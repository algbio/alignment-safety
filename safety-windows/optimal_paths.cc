#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <string>
#include <assert.h>
#include <array>

#include <fstream>
#include <iostream>

#include "optimal_paths.h"

// translate fasta file letters to amino acid symbols (see http://www.math.utep.edu/Faculty/mleung/bioinformatics/aacodon.html)
//                        A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z
const int64_t LTA[26] = { 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20 };

void alignments_into_fasta(int64_t number_of_paths, Dag &d, const std::string &a, const std::string &fasta_file, const std::string &id) {
	std::string alignment = "";
	std::ofstream fasta_stream;
	fasta_stream.open(fasta_file, std::ofstream::app);
	int64_t al_idx = 0;
	std::function<void(std::string, int64_t)> print_alignment = [&](std::string s, int64_t c) {
		fasta_stream << ">Alignment_" << id << '_' << al_idx++ << '_' << c << '\n';
		fasta_stream << s << '\n';
	};

	// calculate set of paths P as the top suboptimal alignments
	int64_t n = (int64_t) d.adj.size();
	std::vector<std::vector<int64_t>> P;
	std::vector<int64_t> costs_of_final_paths;
	std::vector<int64_t> count(n, 0LL);

	std::vector<std::pair<int64_t, std::vector<int64_t>>> B; // pair (cost, path)
	std::vector<int64_t> just_src = { d.src };
	B.push_back(std::make_pair(0LL, just_src));
	while (!B.empty() && count[d.sink] < number_of_paths) {
		std::pop_heap(B.begin(), B.end());
		auto [cost, Pp] = B.back();
		B.pop_back();
		count[Pp.back()]++;
		if (Pp.back() == d.sink)
			P.push_back(Pp), costs_of_final_paths.push_back(cost);
		if (count[Pp.back()] <= number_of_paths) {
			for (auto [v, nxtc]: d.adj_costs[Pp.back()]) {
				std::vector<int64_t> Pv = Pp;
				Pv.push_back(v);
				B.push_back(std::make_pair(cost + nxtc, Pv));
				std::push_heap(B.begin(), B.end());
			}
		}
	}

	int64_t ci = 0;
	for (const std::vector<int64_t> &path: P) {
		std::string alignment = "";
		for (int64_t i = 1; i < (int64_t) path.size(); i++) {
			std::pair<int64_t, int64_t> p = d.transr[path[i - 1]];
			std::pair<int64_t, int64_t> n = d.transr[path[i]];
			alignment += (n.first == p.first + 1 ? a[p.first] : '_');
		}
		print_alignment(alignment, costs_of_final_paths[ci++]);
	}

	fasta_stream.close();
}

std::vector<std::vector<std::vector<std::vector<Node>>>> build_dp_matrix(const std::string &a,
		const std::string &b, const int64_t GAP_COST, const int64_t START_GAP, const int64_t cost_matrix[21][21],
		int64_t sign) {
	int64_t n = (int64_t) a.size();
	int64_t m = (int64_t) b.size();

	// create affine linear gap cost graph (see README for illustration)
	std::vector<std::vector<std::vector<std::vector<Node>>>> adj(n + 1,
			std::vector<std::vector<std::vector<Node>>>(m + 1,
			std::vector<std::vector<Node>>(3)));
	for (int64_t i = 0; i <= n; i++) for (int64_t j = 0; j <= m; j++) {
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

std::vector<std::vector<std::vector<int64_t>>>
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, bool dir) {
	int64_t n = (int64_t) adj.size() - 1;
	assert(n > 0);
	int64_t m = (int64_t) adj[0].size() - 1;
	assert(m > 0);
//	return dijkstra(adj, (dir ? n : 0), (dir ? m : 0));

	std::vector<std::vector<std::vector<int64_t>>> d(n + 1, std::vector<std::vector<int64_t>>(m + 1, std::vector<int64_t>(3, -(1 << 30))));
	(dir ? d[n][m][0] : d[0][0][0]) = 0;
	auto update_dist = [&](int64_t i, int64_t j, int64_t k) {
		for (const Node &node: adj[i][j][k]) {
			d[node.N_index][node.M_index][node.type] = std::max(d[node.N_index][node.M_index][node.type], d[i][j][k] + node.cost);
		}
	};
	if (!dir) {
		for (int64_t i = 0; i <= n; i++) for (int64_t j = 0; j <= m; j++) for (int64_t k = 2; k >= 0; k--)
			update_dist(i, j, k);
	} else {
		for (int64_t i = n; i >= 0; i--) for (int64_t j = m; j >= 0; j--) for (int64_t k = 0; k <= 2; k++)
			update_dist(i, j, k);
	}
	return d;
}
