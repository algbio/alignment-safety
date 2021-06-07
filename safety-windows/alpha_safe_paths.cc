#include <iostream>
#include <vector>
#include <set>
#include <stack>
#include <unordered_map>
#include <gmp.h>

#include "alpha_safe_paths.h"

#define SRC 0
#define DEST (n-1)

const double eps = 1e-6;

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

// TODO: consider returning vector<double>, or use a bigint library
std::vector<double> amount_paths(std::vector<std::vector<int>> &dag, int sink) {
	int n = (int) dag.size();
	std::vector<int> sorted = topsort(dag);

	std::vector<double> am(n, 0);
	am[sink] = 1;
	for (int i = n - 1; i >= 0; i--) {
		for (int v: dag[sorted[i]]) am[sorted[i]] += am[v];
	}
	return am;
}

std::vector<std::vector<double>> path_ratios(std::vector<std::vector<int>> &dag) {
	int n = (int) dag.size();

	std::vector<std::vector<int>> rdag(n);
	for (int i = 0; i < n; i++) for (int v: dag[i]) rdag[v].push_back(i);

	// TODO: perhaps create a dag class with source/sink variables
	std::vector<double> am = amount_paths(dag, DEST);
	std::vector<double> ram = amount_paths(rdag, SRC);

	std::vector<std::vector<double>> ratios(n);
	for (int i = 0; i < n; i++) for (int v: dag[i]) {
		ratios[i].push_back((double) am[v] * ram[i] / am[SRC]);
	}
	return ratios;
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

std::vector<int> find_alpha_path(std::vector<std::vector<int>> &dag,
		std::vector<std::vector<double>> &ratios, double alpha) {
	int n = (int) ratios.size();
	std::vector<int> sorted = topsort(dag);

	std::vector<std::pair<int, int>> needed;
	for (int i = 0; i < n; i++) {
		int u = sorted[i];
		int am = 0;
		for (int j = 0; j < (int) dag[u].size(); j++) {
			int v = dag[u][j];
			double d = ratios[u][j];
			if (d > alpha) needed.emplace_back(u, v), am++;
		}
		assert(am <= 1);
	}

	std::vector<int> order(n);
	for (int i = 0; i < n; i++) order[sorted[i]] = i;

	std::vector<int> path;
	int last = SRC;
	for (auto [u,v]: needed) {
		find_path(last, u, path, dag, order);
		last = v;
	}
	if (last != DEST)
		find_path(last, DEST, path, dag, order);

	return path;
}
