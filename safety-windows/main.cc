#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <string>


#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"

int print_usage(char **argv, int help) {
	std::cout << "How to run: " << argv[0] << " -f <clusterfile> [OPTION...]\n\n";
	std::cout << "\t-a, --alpha\tFloating value, choose edges that appear in (alpha*100)% of all\n\t            \t(sub-)optimal paths in the alpha-safe path. (Default: 0.75)\n";
	std::cout << "\t-b, --samecost\tInteger, set the cost of aligning two equal characters. (Default: 0)\n";
	std::cout << "\t-c, --diffcost\tInteger, set the cost of aligning two unequal characters. (Default: 1)\n";
	std::cout << "\t-d, --gapcost\tInteger, set the cost of aligning a character to a gap. (Default: 1)\n";
	std::cout << "\t-e, --startgap\tInteger, set the cost of starting a gap alignment. (Default: 0)\n";
	std::cout << "\t-g, --threshold\tFloating value, set the treshold for suboptimality. (Default: 1.0, Range: [1.0,infinity))\n";
	std::cout << "\t-h, --help\tShows this help message.\n";
	return help;
}

struct Protein {
	std::string descriptor;
	std::string sequence;
	Protein(std::string descriptor) : descriptor(descriptor) {}
};

static int verbose_flag; // this does nothing for now

int main(int argc, char **argv) {
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output

	mpq_class alpha = 0.75, TH = 1;

	/*
	std::string a = "EBCDE";
	std::string b = "ABCDFE";
	std::vector<std::vector<std::vector<int>>> dp = opt_alignment(build_dp_matrix(a, b), 0, 0);
	int n = a.size(), m = b.size();
	std::cout << dp[n][m][0] << std::endl;
	*/

    int SAME_COST = 0;
	int DIFF_COST = 1;
	int GAP_COST = 1;
	int START_GAP = 0;

	std::string file;
	bool help_flag = false;
	bool read_file = false;

	int c;
	while (1) {
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 }, // this does nothing for now
			{ "alpha", required_argument, 0, 'a' },
			{ "samecost", required_argument, 0, 'b' },
			{ "diffcost", required_argument, 0, 'c' },
			{ "gapcost", required_argument, 0, 'd' },
			{ "startgap", required_argument, 0, 'e' },
			{ "threshold", required_argument, 0, 'g' },
			{ "file", required_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};
	
		int option_index = 0;
		c = getopt_long(argc, argv, "a:b:c:d:e:g:f:h", long_options, &option_index);
		if (c == -1) break;

		switch (c) {
			case 'a':
				alpha = std::stof(optarg);
				break;
			case 'b':
				SAME_COST = atoi(optarg);
				break;
			case 'c':
				DIFF_COST = atoi(optarg);
				break;
			case 'd':
				GAP_COST = atoi(optarg);
				break;
			case 'e':
				START_GAP = atoi(optarg);
				break;
			case 'g':
				TH = std::stof(optarg);
				break;
			case 'f':
				read_file = true;
				file = optarg;
				break;
			case 'h':
				help_flag = true;
				break;
			case '?': break;
			default: abort();
		}
	}
	
	if (help_flag) {
		return print_usage(argv, 0);
	}
	if (!read_file) {
		return print_usage(argv, 1);
	}

	if (alpha >= 1.0) {
		std::cerr << "Warning: for alpha values >= 1.0, the program will not return any safety windows.\n";
	} else if (alpha < 0.5) {
		std::cerr << "Warning: for alpha values < 0.5, the program will not behave well defined and might crash.\n";
	}


	std::ifstream input(file);
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

		Dag d = gen_dag(a, b, TH, SAME_COST, DIFF_COST, GAP_COST, START_GAP); // fix custom threshold bug!
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
