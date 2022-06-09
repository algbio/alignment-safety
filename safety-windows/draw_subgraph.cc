#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "draw_subgraph.h"

std::string draw_subgraph(const int64_t IDX, const int64_t n, const int64_t m, const Dag &d, const std::vector<int64_t> &alpha_path, const std::vector<std::pair<int64_t, int64_t>> &windows, const std::string &a, const std::string &b) {
	std::unordered_map<int64_t, int64_t> order;
	for (int64_t i = 0; i < (int) alpha_path.size(); i++) order[alpha_path[i]] = i;

	auto has = [&](int64_t i, int64_t j) {
		// check if pair (i, j) is part of the suboptimal dag
		return d.trans.find(std::make_pair(i, j)) != d.trans.end();
	};
	auto is_safe = [&](int64_t i1, int64_t j1, int64_t i2, int64_t j2) {
		if (!has(i1, j1) || !has(i2, j2)) return 0;
		int64_t a1 = -1, a2 = -1;
		for (int64_t a: alpha_path) {
			if (d.transr.at(a) == std::make_pair(i1, j1)) a1 = a;
			if (d.transr.at(a) == std::make_pair(i2, j2)) a2 = a;
		}
		if (a1 < 0 || a2 < 0) return 0;
		if (order[a1] > order[a2]) {
			std::cerr << "Error. order of left index larger than order of right index.\n";
			std::cerr << "DBG: " << i1 << ' ' << j1 << ' ' << i2 << ' ' << j2 << ' ' << a1 << ' ' << a2 << std::endl;
			assert(false);
		}
		int idx = 1;
		for (auto [L, R]: windows) {
			if (order[a1] >= order[L] && order[a1] <= order[R] && order[a2] >= order[L] && order[a2] <= order[R]) return idx;
			idx++;
		}
		return 0;
	};
	auto outside = [&](int64_t i, int64_t j) {
		return i < 0 || i >= n || j < 0 || j >= m;
	};
	auto in_alpha_path = [&](int64_t i, int64_t j) {
		if (!has(i, j)) return false;
		std::array<int64_t, 3> as = d.trans.at(std::make_pair(i, j));
		for (int64_t a: as) if (std::find(alpha_path.begin(), alpha_path.end(), a) != alpha_path.end()) return true;
		return false;
	};
	std::string output = "graph fasta_" + std::to_string(IDX) + " {\nnode [shape=none, width=.2, height=.05, fixedsize=true];\n";//rankdir=LR;\n";
	std::vector<std::string> nodes, edges;
	for (int64_t i = 0; i < n; i++) {
		for (int64_t j = 0; j < m; j++) {
			std::string node_label;
			if (i == 0 && j > 0) node_label = std::string(1, b[j - 1]);
			else if (i > 0 && j == 0) node_label = std::string(1, a[i - 1]);
			else node_label = std::to_string(i) + "_" + std::to_string(j);
			std::string node = "\"" + std::to_string(i) + "_" + std::to_string(j) + "\" [label=\"" + node_label + "\"";
			if ((i > 0 && j > 0) || (i == 0 && j == 0)) node += ", style=invis]"; else node += "]";
			nodes.push_back(node);

			std::vector<std::pair<int64_t, int64_t>> from = { {i - 1, j - 1}, {i - 1, j}, {i, j - 1} };
			for (auto [l, k]: from) if (!outside(l, k)) {
				std::string edge = "\"" + std::to_string(l) + "_" + std::to_string(k) + "\" -- \"" + std::to_string(i) + "_" + std::to_string(j) + "\"";
				if (!has(i, j) || !has(l, k)) {
					edge += " [style=invis]";
					assert(!in_alpha_path(i, j) || !in_alpha_path(l, k));
				} else {
					int s;
					if ((s = is_safe(l, k, i, j))) edge += ((s % 2) ? " [color=red, penwidth=8]" : " [color=red, penwidth=8]"); // If wished, we can differentiate between two seperate SW by color here. Though, some SW are intersecting each other.
					//else if (in_alpha_path(i, j) && in_alpha_path(l, k)) edge += " [color=blue, penwidth=3]";
					else edge += " [penwidth=3]";
				}
				edges.push_back(edge);
			}
		}
	}

	for (int i = 0, k = 0; i < n; i++) {
		output += "subgraph row_" + std::to_string(i) + " {\n";
		output += "rank=same;\n";
		for (int j = 0; j < m; j++) {
			output += nodes[k++] + ";\n";
		}
		output += "}\n";
	}
	for (std::string edge: edges) output += edge + ";\n";
	output += "}\n";

	return output;
}

