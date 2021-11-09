#pragma once
#include <vector>
#include <gmpxx.h>

#include "optimal_paths.h"

// given a dag of optimal paths, find a path with almost safe (> alpha) paths

std::vector<int> topsort(std::vector<std::vector<int>> &dag);

// Find path from src to dest. Assertion error in case if such path does not exist
void find_path(int src, int dest, std::vector<int> &path,
		std::vector<std::vector<int>> &dag, std::vector<int> &order);

// For each vertex, save the amount of paths starting from the vertex to the sink
template<class T>
std::vector<T> amount_paths(std::vector<std::vector<int>> &dag) {
	int n = (int) dag.size();
	std::vector<int> sorted = topsort(dag);
	int sink = sorted.back();

	//std::vector<mpz_class> am(n, mpz_class(0));
	std::vector<T> am(n, T(0));
	am[sink] = 1;
	for (int i = n - 1; i >= 0; i--) {
		for (int v: dag[sorted[i]]) am[sorted[i]] += am[v];
	}
	return am;
}

// For each edge, calculate the % of s-t paths they are part in
template<class T, class K>
std::vector<std::vector<K>> path_ratios(Dag &d) {
	std::vector<std::vector<int>> &dag = d.adj;
	int n = (int) dag.size();

	std::vector<std::vector<int>> rdag(n);
	for (int i = 0; i < n; i++) for (int v: dag[i]) rdag[v].push_back(i);

	// TODO: perhaps create a dag class with source/sink variables
	//std::vector<mpz_class> am = amount_paths(dag);
	//std::vector<mpz_class> ram = amount_paths(rdag);
	std::vector<T> am = amount_paths<T>(dag);
	std::vector<T> ram = amount_paths<T>(rdag);

	for (int i = 0; i < n; i++) assert(am[i] <= am[d.src]);

	std::vector<std::vector<K>> ratios(n);

	int SRC = d.src;

	for (int i = 0; i < n; i++) for (int v: dag[i]) {
		// nxt = am[v] * ram[i] / am[SRC]
		//mpq_class nxt = am[v] * ram[i];
		K nxt = am[v] * ram[i];
		nxt /= am[SRC];
		ratios[i].push_back(nxt);
	}
	return ratios;
}

// Find s--t path that contains all edges with occurence ratio > alpha. Might fail if alpha < 0.5
template<class K>
std::vector<int> find_alpha_path(Dag &d,
		std::vector<std::vector<K>> &ratios, K alpha) {
	std::vector<std::vector<int>> &dag = d.adj;
	int n = (int) ratios.size();
	std::vector<int> sorted = topsort(dag);

	std::vector<std::pair<int, int>> needed;
	for (int i = 0; i < n; i++) {
		int u = sorted[i];
		int am = 0;
		for (int j = 0; j < (int) dag[u].size(); j++) {
			int v = dag[u][j];
			K d = ratios[u][j];
			assert(d <= 1);
			if (d > alpha) needed.emplace_back(u, v), am++;
		}
		assert(am <= 1);
	}

	/*std::cerr << "NEEDED DBG" << std::endl;
	for (auto [x,y]: needed) std::cerr << x << ' ' << y << std::endl;
	std::cerr << "NEEDED DBG END" << std::endl;*/

	std::vector<int> order(n);
	for (int i = 0; i < n; i++) order[sorted[i]] = i;

	std::vector<int> path;
	int last = d.src;
	for (auto [u,v]: needed) {
		find_path(last, u, path, dag, order);
		last = v;
	}
	find_path(last, d.sink, path, dag, order);

	return path;
}

