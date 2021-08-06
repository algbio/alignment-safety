#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "draw_subgraph.h"

void draw_subgraph(const int IDX, const Dag &d, const std::vector<int> &alpha_path, const std::vector<std::pair<int, int>> &windows) {
	std::cout << "digraph fasta_" << IDX << " {\n";
	std::vector<std::string> edges;
	std::map<std::pair<int, int>, int> where;
	for (int i = 0; i < (int) d.adj.size(); i++) {
		std::pair<int, int> p = d.transr.at(i);
		for (int k: d.adj[i]) {
			std::pair<int, int> pp = d.transr.at(k);
			if (p == pp) continue;
			std::string edge = "\"" + std::to_string(p.first) + "_" + std::to_string(p.second) + "\" -> \"" + std::to_string(pp.first) + "_" + std::to_string(pp.second) + "\"";
			//std::string edge = "\"" + std::to_string(i) + "\" -> \"" + std::to_string(k) + "\"";
			where[std::make_pair(i, k)] = (int) edges.size();
			edges.push_back(edge);
		}
	}

	int cw = 0;
	bool inside = false;
	for (int i = 0; i < (int) alpha_path.size() - 1; i++) {
		int u = alpha_path[i], v = alpha_path[i + 1];
		if (where.find(std::make_pair(u, v)) == where.end()) continue;
		if (windows[cw].second == u) {
			cw++;
			inside = false;
		}
		if (cw < (int) windows.size() && windows[cw].first == u) {
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
