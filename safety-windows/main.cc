#include <iostream>
#include <iomanip>
#include <vector>

#include "alpha_safe_paths.h"

int main() {
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output

	const double alpha = 0.5;
	
	// test implementations
	int n = 6;
	std::vector<std::vector<int>> dag(n);
	dag[0] = { 1, 2 };
	dag[1] = { 3, 4 };
	dag[2] = { 5 };
	dag[3] = { 4, 5 };
	dag[4] = { 5 };

	std::vector<int> paths = amount_paths(dag, 5);
	for (int i = 0; i < n; i++) std::cout << paths[i] << " \n"[i + 1 == n];

	std::vector<std::vector<double>> ratios = path_ratios(dag);
	for (int i = 0; i < n; i++) {
		std::cout << "\nout from " << i << ":\n";
		for (int j = 0; j < (int) dag[i].size(); j++) {
			std::cout << dag[i][j] << ' ' << ratios[i][j] << '\n';
		}
	}

	// test find_path function
	std::vector<int> path;
	/*std::vector<int> sorted = topsort(dag);
	std::vector<int> order(n);
	for (int i = 0; i < n; i++) order[sorted[i]] = i;
	for (int i = 0; i < n; i++) std::cerr << order[i] << " \n"[i + 1 == n];
	find_path(2, 5, path, dag, order);
	std::cout << '\n';
	for (int v: path) std::cout << v << ' ';
	std::cout << '\n';*/

	// test find_alpha_path function
	std::cout << "\nFinal alpha safe path:\n";
	path = find_alpha_path(dag, ratios, alpha);
	for (int v: path) std::cout << v << ' ';
	std::cout << '\n';
	return 0;
}
