#include <iostream>
#include <vector>
#include <set>
#include <stack>
#include <unordered_map>
#include <assert.h>

#include <gmpxx.h>

#include "alpha_safe_paths.h"
#include "optimal_paths.h"


// given a dag of optimal paths, find a path with almost safe (> alpha) paths

std::vector<int> topsort(std::vector<std::vector<int>> &dag) {
	int n = (int) dag.size();
	std::vector<int> indeg(n, 0);
	for (int i = 0; i < n; i++) {
		for (int v: dag[i]) indeg[v]++;
	}
	std::stack<int> nxt;
	for (int i = 0; i < n; i++) if (indeg[i] == 0) nxt.push(i);

	std::vector<int> sorted;
	while (!nxt.empty()) {
		int v = nxt.top();
		nxt.pop();
		sorted.push_back(v);

		for(int u: dag[v]) if (--indeg[u] == 0) nxt.push(u);
	}

	assert((int) sorted.size() == n); // true iff input is a dag
	return sorted;
}

void find_path(int src, int dest, std::vector<int> &path, std::vector<std::vector<int>> &dag,
		std::vector<int> &order) {
	//std::cerr << "FIND PATH DBG: " << src << ' ' << dest << std::endl;
	std::unordered_map<int, bool> vis;
	std::function<bool(int)> dfs = [&](int current) {
		if (order[current] > order[dest]) return false;
		if (vis[current]) return false;
		path.push_back(current);
		vis[current] = true;
		if (current == dest) return true;
		for (int nxt: dag[current]) {
			if (dfs(nxt)) return true;
		}
		path.pop_back();
		return false;
	};
	assert(dfs(src));
}
