#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"

void print_usage(char **argv) {
	// TODO: Use cerr if -h argument is not used
	std::cout << "How to run: " << argv[0] << " -f <clusterfile>\n";
	std::cout << "Run \"" << argv[0] << " -h\" to show this help message.\n";
}

struct Protein {
	std::string descriptor;
	std::string sequence;
	Protein(std::string descriptor) : descriptor(descriptor) {}
};

int main(int argc, char **argv) {
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output

	const double alpha = 0.5;

	if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-f") != 0) {
		print_usage(argv);
		return 0;
	}

	if (argc <= 2) {
		std::cerr << "Error: File missing\n";
		print_usage(argv);
		return 0;
	}

	std::ifstream input(argv[2]);
	std::vector<Protein> proteins;
	for (std::string line; std::getline(input, line); ) {
		if ((int) line.size() <= 0) continue;
		if (line[0] == '>') {
			proteins.push_back(Protein(line));
		} else {
			proteins.back().sequence += line;
		}
	}
	// proteins[0] will be the reference

	int PS = (int) proteins.size();

	std::cout << 0 << '\n' << PS << '\n';
	for (int i = 1; i < PS; i++) {
		std::cout << i << ' ';
		const std::string &a = proteins[0].sequence;
		const std::string &b = proteins[i].sequence;
		std::vector<std::vector<int>> dp = opt_alignment(a, b);

		Dag d = gen_dag(dp, a, b);
		std::vector<std::vector<int>> adj = d.adj;
		int k = (int) adj.size();

		std::vector<std::vector<double>> ratios = path_ratios(adj);

		std::vector<int> path = find_alpha_path(adj, ratios, alpha);

		std::vector<double> r = find_ratios(path, adj, ratios);
		std::vector<std::pair<int, int>> windows = safety_windows(adj, path, r, alpha);

		/*std::map<std::pair<int, int>, int> trans = d.trans;
		for (auto [a, b]: trans) {
			std::cout << a.first << ' ' << a.second << ' ' << b << '\n';
		}*/

		std::map<int, std::pair<int, int>> transr = d.transr;
		std::cout << windows.size();
		for (auto [x, y]: windows) {
			int a = transr[x].first, b = transr[y].first;
			std::cout << ' ' << x << ' ' << y;
		}
		std::cout << '\n';
	}


	
	
	// test implementations
	/*int n = 6;
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
	std::vector<int> sorted = topsort(dag);
	std::vector<int> order(n);
	for (int i = 0; i < n; i++) order[sorted[i]] = i;
	for (int i = 0; i < n; i++) std::cerr << order[i] << " \n"[i + 1 == n];
	find_path(2, 5, path, dag, order);
	std::cout << '\n';
	for (int v: path) std::cout << v << ' ';
	std::cout << '\n';

	// test find_alpha_path function
	std::cout << "\nFinal alpha safe path:\n";
	path = find_alpha_path(dag, ratios, alpha);
	for (int v: path) std::cout << v << ' ';
	std::cout << '\n';

	// Calculate the safety windows
	std::cout << "Safety intervals:\n";
	std::vector<double> rat_path = find_ratios(path, dag, ratios);
	std::vector<std::pair<int, int>> windows = safety_windows(dag, path, rat_path, alpha);
	for (auto [x,y]: windows) {
		std::cout << x << ' ' << y << '\n';
	}*/
}
