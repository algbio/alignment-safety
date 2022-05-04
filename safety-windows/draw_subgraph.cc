#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "draw_subgraph.h"

std::string draw_subgraph(const int64_t IDX, const int64_t n, const int64_t m, const Dag &d, const std::vector<int64_t> &alpha_path, const std::vector<std::pair<int64_t, int64_t>> &windows) {
	auto has = [&](int64_t i, int64_t j) {
		// check if pair (i, j) is part of the suboptimal dag
		return d.trans.find(std::make_pair(i, j)) != d.trans.end();
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
	std::string output = "graph fasta_" + std::to_string(IDX) + " {\nnode [shape=none];\n";//rankdir=LR;\n";
	std::vector<std::string> nodes, edges;
	for (int64_t i = 0; i < n; i++) {
		for (int64_t j = 0; j < m; j++) {
			std::string node = "\"" + std::to_string(i) + "_" + std::to_string(j) + "\" [label=\"" + std::to_string(i) + "_" + std::to_string(j) + "\", style=invis]";
			nodes.push_back(node);

			std::vector<std::pair<int64_t, int64_t>> from = { {i - 1, j - 1}, {i - 1, j}, {i, j - 1} };
			for (auto [l, k]: from) if (!outside(l, k)) {
				std::string edge = "\"" + std::to_string(l) + "_" + std::to_string(k) + "\" -- \"" + std::to_string(i) + "_" + std::to_string(j) + "\"";
				if (!has(i, j) || !has(l, k)) {
					edge += " [style=invis]";
					assert(!in_alpha_path(i, j) || !in_alpha_path(l, k));
				}
				else if (in_alpha_path(i, j) && in_alpha_path(l, k)) edge += " [color=blue]"; // this can be made faster easily if needed
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

