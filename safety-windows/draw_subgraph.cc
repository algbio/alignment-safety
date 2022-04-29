#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "draw_subgraph.h"

// TODO:
// Draw *all* nodes of G(A,B) in a grid, have all nodes invisible
// Only make arcs visible if they are part of G_Delta(A,B)
// Highlight OPT arcs and alpha safe arcs (with different alpha values)

void draw_subgraph(const int64_t IDX, const Dag &d, const std::vector<int64_t> &alpha_path, const std::vector<std::pair<int64_t, int64_t>> &windows) {
	std::cout << "digraph fasta_" << IDX << " {\n";
	std::cout << "node [shape=none];" << '\n';
	std::cout << "rankdir=LR;" << '\n';
	std::vector<std::string> edges;
	std::map<std::pair<int64_t, int64_t>, int64_t> where;
	for (int64_t i = 0; i < (int64_t) d.adj.size(); i++) {
		std::pair<int64_t, int64_t> p = d.transr.at(i);
		for (int64_t k: d.adj[i]) {
			std::pair<int64_t, int64_t> pp = d.transr.at(k);
			if (p == pp) continue;
			std::string edge = "\"" + std::to_string(p.first) + "_" + std::to_string(p.second) + "\" -> \"" + std::to_string(pp.first) + "_" + std::to_string(pp.second) + "\"";
			//std::string edge = "\"" + std::to_string(i) + "\" -> \"" + std::to_string(k) + "\"";
			where[std::make_pair(i, k)] = (int64_t) edges.size();
			edges.push_back(edge);
		}
	}

	int64_t cw = 0;
	bool inside = false;
	for (int64_t i = 0; i < (int64_t) alpha_path.size() - 1; i++) {
		int64_t u = alpha_path[i], v = alpha_path[i + 1];
		if (where.find(std::make_pair(u, v)) == where.end()) continue;
		if (windows[cw].second == u) {
			cw++;
			inside = false;
		}
		if (cw < (int64_t) windows.size() && windows[cw].first == u) {
			inside = true;
		}
		if (inside) {
			edges[where[std::make_pair(u, v)]] += (cw % 2 ? " [color=red]" : " [color=blue]");
		} else {
			edges[where[std::make_pair(u, v)]] += " [color=pink]";
		}
	}

	for (std::string edge: edges) {
		std::cout << edge << '\n';
	}
	std::cout << "}\n";
}
