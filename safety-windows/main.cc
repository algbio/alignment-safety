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
	std::cout << "\t-c, --costmat\tReads the aligning cost of two symbols from a text file. The text file is a lower triangular matrix with 20 lines. (Default: BLOSUM62)\n";
	std::cout << "\t-d, --gapcost\tInteger, set the cost of aligning a character to a gap. (Default: 1)\n";
	std::cout << "\t-e, --startgap\tInteger, set the cost of starting a gap alignment. (Default: 11)\n";
	//std::cout << "\t-g, --threshold\tFloating value, set the treshold for suboptimality. (Default: 0.0, Range: [0.0, 1.0])\n";
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

	mpq_class alpha = 0.75, TH = 0;

	/*
	std::string a = "EBCDE";
	std::string b = "ABCDFE";
	std::vector<std::vector<std::vector<int>>> dp = opt_alignment(build_dp_matrix(a, b), 0, 0);
	int n = a.size(), m = b.size();
	std::cout << dp[n][m][0] << std::endl;
	*/

	int GAP_COST = 1;
	int START_GAP = 11;
	std::string file, cost_matrix_file;
	bool help_flag = false;
	bool read_file = false;
	bool read_cost_matrix = false;

	int c;
	while (1) {
		// TODO: Read cost_matrix file
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 }, // this does nothing for now
			{ "alpha", required_argument, 0, 'a' },
			{ "costmat", required_argument, 0, 'c' },
			{ "gapcost", required_argument, 0, 'd' },
			{ "startgap", required_argument, 0, 'e' },
			//{ "threshold", required_argument, 0, 'g' },
			{ "file", required_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};
	
		int option_index = 0;
		c = getopt_long(argc, argv, "a:c:d:e:f:h", long_options, &option_index);
		if (c == -1) break;

		switch (c) {
			case 'a':
				alpha = std::stof(optarg);
				break;
			case 'c':
				read_cost_matrix = true;
				cost_matrix_file = optarg;
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

	// BLOSUM62 matrix
	int cost_matrix[21][21] = {
		// Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val  Def
		{  -4,   1,   2,   2,   0,   1,   1,   0,   2,   1,   1,   1,   1,   2,   1,  -1,   0,   3,   2,   0,   1 }, // Ala
		{   1,  -5,   0,   2,   3,  -1,   0,   2,   0,   3,   2,  -2,   1,   3,   2,   1,   1,   3,   2,   3,   1 }, // Arg
		{   2,   0,  -6,  -1,   3,   0,   0,   0,  -1,   3,   3,   0,   2,   3,   2,  -1,   0,   4,   2,   3,   1 }, // Asn
		{   2,   2,  -1,  -6,   3,   0,  -2,   1,   1,   3,   4,   1,   3,   3,   1,   0,   1,   4,   3,   3,   1 }, // Asp
		{   0,   3,   3,   3,  -9,   3,   4,   3,   3,   1,   1,   3,   1,   2,   3,   1,   1,   2,   2,   1,   1 }, // Cys
		{   1,  -1,   0,   0,   3,  -5,  -2,   2,   0,   3,   2,  -1,   0,   3,   1,   0,   1,   2,   1,   2,   1 }, // Gln
		{   1,   0,   0,  -2,   4,  -2,  -5,   2,   0,   3,   3,  -1,   2,   3,   1,   0,   1,   3,   2,   2,   1 }, // Glu
		{   0,   2,   0,   1,   3,   2,   2,  -6,   2,   4,   4,   2,   3,   3,   2,   0,   2,   2,   3,   3,   1 }, // Gly
		{   2,   0,  -1,   1,   3,   0,   0,   2,  -8,   3,   3,   1,   2,   1,   2,   1,   2,   2,  -2,   3,   1 }, // His
		{   1,   3,   3,   3,   1,   3,   3,   4,   3,  -4,  -2,   3,  -1,   0,   3,   2,   1,   3,   1,  -3,   1 }, // Ile
		{   1,   2,   3,   4,   1,   2,   3,   4,   3,  -2,  -4,   2,  -2,   0,   3,   2,   1,   2,   1,  -1,   1 }, // Leu
		{   1,  -2,   0,   1,   3,  -1,  -1,   2,   1,   3,   2,  -5,   1,   3,   1,   0,   1,   3,   2,   2,   1 }, // Lys
		{   1,   1,   2,   3,   1,   0,   2,   3,   2,  -1,  -2,   1,  -5,   0,   2,   1,   1,   1,   1,  -1,   1 }, // Met
		{   2,   3,   3,   3,   2,   3,   3,   3,   1,   0,   0,   3,   0,  -6,   4,   2,   2,  -1,  -3,   1,   1 }, // Phe
		{   1,   2,   2,   1,   3,   1,   1,   2,   2,   3,   3,   1,   2,   4,  -7,   1,   1,   4,   3,   2,   1 }, // Pro
		{  -1,   1,  -1,   0,   1,   0,   0,   0,   1,   2,   2,   0,   1,   2,   1,  -4,  -1,   3,   2,   2,   1 }, // Ser
		{   0,   1,   0,   1,   1,   1,   1,   2,   2,   1,   1,   1,   1,   2,   1,  -1,  -5,   2,   2,   0,   1 }, // Thr
		{   3,   3,   4,   4,   2,   2,   3,   2,   2,   3,   2,   3,   1,  -1,   4,   3,   2,  -11, -2,   3,   1 }, // Trp
		{   2,   2,   2,   3,   2,   1,   2,   3,  -2,   1,   1,   2,   1,  -3,   3,   2,   2,  -2,  -7,   1,   1 }, // Tyr
		{   0,   3,   3,   3,   1,   2,   2,   3,   3,  -3,  -1,   2,  -1,   1,   2,   2,   0,   3,   1,  -4,   1 }, // Val
		{   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 },
	};

	if (read_cost_matrix) {
		std::ifstream costmat(cost_matrix_file);
		for (int i = 0; i < 20; i++) for (int j = 0; j <= i; j++) {
			costmat >> cost_matrix[i][j];
			cost_matrix[j][i] = cost_matrix[i][j];
		}
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

		Dag d = gen_dag(a, b, cost_matrix, TH, GAP_COST, START_GAP);
		std::vector<std::vector<int>> adj = d.adj;

		std::vector<std::vector<mpq_class>> ratios = path_ratios(d);

		std::vector<int> path = find_alpha_path(d, ratios, alpha);
		std::unordered_map<int, int> cnt;
		for (int v: path) {
			assert(cnt[v] == 0);
			cnt[v]++;
		}

		std::vector<mpq_class> r = find_ratios(path, adj, ratios);
		std::vector<std::pair<int, int>> windows_tmp = safety_windows(r, path, alpha);

		std::vector<std::pair<int, int>> windows, windowsp;
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
			int Lp = transr[LT].second, Rp = transr[RT].second;
			/*while (outside(L, R)) windows.pop_back();
			if (!inside(L, R)) windows.emplace_back(L, R);*/
			windows.emplace_back(L, R);
			windowsp.emplace_back(Lp, Rp);
		}


		std::cout << windows.size() << '\n';;
		for (int i = 0; i < (int) windows.size(); i++) {
			auto [x, y] = windows[i];
			auto [xp, yp] = windowsp[i];
			std::cout << x << ' ' << y << ' ' << xp << ' ' << yp << '\n';;
		}
	}
}
