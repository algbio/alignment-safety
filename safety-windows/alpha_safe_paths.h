#pragma once
#include <vector>
#include <gmpxx.h>

#include "optimal_paths.h"

// given a dag of optimal paths, find a path with almost safe (> alpha) paths

std::vector<int> topsort(std::vector<std::vector<int>> &dag);

// For each vertex, save the amount of paths starting from the vertex to the sink
std::vector<mpz_class> amount_paths(std::vector<std::vector<int>> &dag);

// For each edge, calculate the % of s-t paths they are part in
std::vector<std::vector<mpq_class>> path_ratios(Dag &d);

// Find path from src to dest. Assertion error in case if such path does not exist
void find_path(int src, int dest, std::vector<int> &path,
		std::vector<std::vector<int>> &dag, std::vector<int> &order);

// Find s--t path that contains all edges with occurence ratio > alpha. Might fail if alpha < 0.5
std::vector<int> find_alpha_path(Dag &d,
		std::vector<std::vector<mpq_class>> &ratios, mpq_class alpha);
