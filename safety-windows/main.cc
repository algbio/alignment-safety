#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <unordered_map>

#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"

int print_usage(char **argv, int help) {
	std::cout << "How to run: " << argv[0] << " -f <clusterfile>\n";
	std::cout << "Run \"" << argv[0] << " -h\" to show this help message.\n";
	return help;
}

struct Protein {
	std::string descriptor;
	std::string sequence;
	Protein(std::string descriptor) : descriptor(descriptor) {}
};

int main(int argc, char **argv) {
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output

	const mpq_class alpha = 0.5;

	/*
	std::string a = "EBCDE";
	std::string b = "ABCDFE";
	std::vector<std::vector<std::vector<int>>> dp = opt_alignment(build_dp_matrix(a, b), 0, 0);
	int n = a.size(), m = b.size();
	std::cout << dp[n][m][0] << std::endl;
	*/

	if (argc >= 2 && strcmp(argv[1], "-h") == 0) {
		return print_usage(argv, 0);
	}

	if (argc <= 1 || argc > 3 || strcmp(argv[1], "-f") != 0) {
		return print_usage(argv, 1);
	}

	if (argc <= 2) {
		std::cerr << "Error: File missing\n";
		return print_usage(argv, 1);
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

	// reference protein and amount of proteins in the cluster
	std::cout << 0 << ' ' << proteins[0].sequence << '\n' << PS << '\n';
	for (int i = 1; i < PS; i++) {
		const std::string &a = proteins[0].sequence;
		const std::string &b = proteins[i].sequence;
		std::cout << i << ' ' << b << ' ';

		Dag d = gen_dag(a, b);
		std::vector<std::vector<int>> adj = d.adj;
		int k = (int) adj.size();

		std::vector<std::vector<mpq_class>> ratios = path_ratios(d);

		std::vector<int> path = find_alpha_path(d, ratios, alpha);
		std::unordered_map<int, int> cnt;
		for (int v: path) {
			assert(cnt[v] == 0);
			cnt[v]++;
		}

		std::vector<mpq_class> r = find_ratios(path, adj, ratios);
		std::vector<std::pair<int, int>> windows_tmp = safety_windows(r, path, alpha);

		std::vector<std::pair<int, int>> windows;
		auto outside = [&](const int &L, const int &R) {
			if (windows.empty()) return false;
			auto [bL, bR] = windows.back();
			return bL >= L && bR <= R;
		};
		auto inside = [&](const int &L, const int &R) {
			if (windows.empty()) return false;
			auto [bL, bR] = windows.back();
			return L >= bL && R <= bR;
		};
		
		std::map<int, std::pair<int, int>> transr = d.transr;
		for (int i = 0; i < (int) windows_tmp.size(); i++) {
			auto [LT, RT] = windows_tmp[i];
			int L = transr[LT].first, R = transr[RT].first;
			while (outside(L, R)) windows.pop_back();
			if (!inside(L, R)) windows.emplace_back(L, R);
		}


		std::cout << windows.size();
		for (auto [x, y]: windows) {
			std::cout << ' ' << x << ' ' << y;
		}
		std::cout << '\n';
	}


	
	/*
	// test implementations
	int n = 6;
	std::vector<std::vector<int>> dag(n);
	dag[0] = { 4 };
	dag[1] = { 5 };
	dag[2] = { 1 };
	dag[3] = { 2 };
	dag[4] = { 3 };

	std::vector<int> sorted = topsort(dag);
	for (int v: sorted) std::cerr << v << ' ';
	std::cerr << std::endl;

	std::vector<int> paths = amount_paths(dag, 5);

	std::vector<std::vector<double>> ratios = path_ratios(dag);

	// test find_path function
	std::vector<int> path;
	std::vector<int> order(n);
	for (int i = 0; i < n; i++) order[sorted[i]] = i;
	find_path(2, 5, path, dag, order);

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
