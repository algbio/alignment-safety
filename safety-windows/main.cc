#include <assert.h>
#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <string>

#include <omp.h>

#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"
#include "draw_subgraph.h"

int64_t print64_t_usage(char **argv, int64_t help) {
	std::cout << "How to run: " << argv[0] << " -f <clusterfile> [OPTION...]\n\n";
	std::cout << "\t-a, --alpha\tFloating value, choose edges that appear in (alpha*100)% of all\n\t            \t(sub-)optimal paths in the alpha-safe path. (Default: 0.75)\n";
	std::cout << "\t-d, --delta\tInteger value, defines suboptimal paths to be in the delta neighborhood of the optimal. (Default: 0, Range: [0.0, inf))\n";
	std::cout << "\t-c, --costmat\tReads the aligning cost of two symbols from a text file.\n\t            \tThe text file is a lower triangular matrix with 20 lines. (Default: BLOSUM62)\n";
	std::cout << "\t-g, --gapcost\tInteger, set the cost of aligning a character to a gap. (Default: 1)\n";
	std::cout << "\t-e, --startgap\tInteger, set the cost of starting a gap alignment. (Default: 11)\n";
	std::cout << "\t-s, --special\tInteger, sets the cost of aligning symbols with special characters.\n\t            \tINF value ignores these charachters. (Default: 1)\n";
	std::cout << "\t-i, --threads\tInteger, specifies the number of threads (Default: 1).\n";
	std::cout << "\t-r, --reference\tProtein identity, selects reference protein. By default, this is the first protein.\n";
	std::cout << "\t-h, --help\tShows this help message.\n";
	return help;
}

struct Protein {
	std::string descriptor;
	std::string sequence;
	Protein(std::string descriptor) : descriptor(descriptor) {}
};

static int verbose_flag; // this does nothing for now
bool use_approx = false;

float alpha = 0.75, TH = 0;
int64_t delta = 0;

int64_t GAP_COST = 1;
int64_t START_GAP = 11;
int64_t SP = 1;
bool ignore_special = false;
std::string file, cost_matrix_file;
bool help_flag = false;
bool read_file = false;
bool read_cost_matrix = false;
std::string reference = "-";
int64_t threads = 1;

// BLOSUM62 matrix
int64_t cost_matrix[21][21] = {
	// Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val  Def
	{  -4,   1,   2,   2,   0,   1,   1,   0,   2,   1,   1,   1,   1,   2,   1,  -1,   0,   3,   2,   0,  SP }, // Ala
	{   1,  -5,   0,   2,   3,  -1,   0,   2,   0,   3,   2,  -2,   1,   3,   2,   1,   1,   3,   2,   3,  SP }, // Arg
	{   2,   0,  -6,  -1,   3,   0,   0,   0,  -1,   3,   3,   0,   2,   3,   2,  -1,   0,   4,   2,   3,  SP }, // Asn
	{   2,   2,  -1,  -6,   3,   0,  -2,   1,   1,   3,   4,   1,   3,   3,   1,   0,   1,   4,   3,   3,  SP }, // Asp
	{   0,   3,   3,   3,  -9,   3,   4,   3,   3,   1,   1,   3,   1,   2,   3,   1,   1,   2,   2,   1,  SP }, // Cys
	{   1,  -1,   0,   0,   3,  -5,  -2,   2,   0,   3,   2,  -1,   0,   3,   1,   0,   1,   2,   1,   2,  SP }, // Gln
	{   1,   0,   0,  -2,   4,  -2,  -5,   2,   0,   3,   3,  -1,   2,   3,   1,   0,   1,   3,   2,   2,  SP }, // Glu
	{   0,   2,   0,   1,   3,   2,   2,  -6,   2,   4,   4,   2,   3,   3,   2,   0,   2,   2,   3,   3,  SP }, // Gly
	{   2,   0,  -1,   1,   3,   0,   0,   2,  -8,   3,   3,   1,   2,   1,   2,   1,   2,   2,  -2,   3,  SP }, // His
	{   1,   3,   3,   3,   1,   3,   3,   4,   3,  -4,  -2,   3,  -1,   0,   3,   2,   1,   3,   1,  -3,  SP }, // Ile
	{   1,   2,   3,   4,   1,   2,   3,   4,   3,  -2,  -4,   2,  -2,   0,   3,   2,   1,   2,   1,  -1,  SP }, // Leu
	{   1,  -2,   0,   1,   3,  -1,  -1,   2,   1,   3,   2,  -5,   1,   3,   1,   0,   1,   3,   2,   2,  SP }, // Lys
	{   1,   1,   2,   3,   1,   0,   2,   3,   2,  -1,  -2,   1,  -5,   0,   2,   1,   1,   1,   1,  -1,  SP }, // Met
	{   2,   3,   3,   3,   2,   3,   3,   3,   1,   0,   0,   3,   0,  -6,   4,   2,   2,  -1,  -3,   1,  SP }, // Phe
	{   1,   2,   2,   1,   3,   1,   1,   2,   2,   3,   3,   1,   2,   4,  -7,   1,   1,   4,   3,   2,  SP }, // Pro
	{  -1,   1,  -1,   0,   1,   0,   0,   0,   1,   2,   2,   0,   1,   2,   1,  -4,  -1,   3,   2,   2,  SP }, // Ser
	{   0,   1,   0,   1,   1,   1,   1,   2,   2,   1,   1,   1,   1,   2,   1,  -1,  -5,   2,   2,   0,  SP }, // Thr
	{   3,   3,   4,   4,   2,   2,   3,   2,   2,   3,   2,   3,   1,  -1,   4,   3,   2,  -11, -2,   3,  SP }, // Trp
	{   2,   2,   2,   3,   2,   1,   2,   3,  -2,   1,   1,   2,   1,  -3,   3,   2,   2,  -2,  -7,   1,  SP }, // Tyr
	{   0,   3,   3,   3,   1,   2,   2,   3,   3,  -3,  -1,   2,  -1,   1,   2,   2,   0,   3,   1,  -4,  SP }, // Val
	{  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP }, // Def
};

std::vector<Protein> proteins;
int64_t ref = 0; // reference protein
std::vector<int64_t> random_order;

template<class T, class K>
void run_case(const int64_t j, std::vector<std::string> &output) {
	int64_t i = random_order[j];
//		std::cerr << "handling id " << i << std::endl;
//		std::cerr << omp_get_thread_num() << std::endl;
	const std::string &a = proteins[ref].sequence;
	const std::string &b = proteins[i].sequence;
	output[i] += std::to_string(i) + ' ' + b + ' ';
//		std::cout << i << ' ' << b << ' ' << std::flush;

	Dag d = gen_dag<K>(a, b, cost_matrix, delta, GAP_COST, START_GAP, verbose_flag);
	std::vector<std::vector<int64_t>> adj = d.adj;

	std::vector<std::vector<K>> ratios = path_ratios<T, K>(d);

	std::vector<int64_t> path = find_alpha_path<K>(d, ratios, alpha);
	std::unordered_map<int64_t, int64_t> cnt;
	for (int64_t v: path) {
		assert(cnt[v] == 0);
		cnt[v]++;
	}

	std::vector<K> r = find_ratios<K>(path, adj, ratios);
	std::vector<std::pair<int64_t, int64_t>> windows_tmp = safety_windows<K>(r, path, alpha);

	std::vector<std::pair<int64_t, int64_t>> windows, windowsp;
	/*auto outside = [&](const int64_t &L, const int64_t &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return bL >= L && bR <= R;
	};
	auto inside = [&](const int64_t &L, const int64_t &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return L >= bL && R <= bR;
	};*/
	
	std::map<int64_t, std::pair<int64_t, int64_t>> transr = d.transr;
	for (int64_t i = 0; i < (int64_t) windows_tmp.size(); i++) {
		auto [LT, RT] = windows_tmp[i];
		int64_t L = transr[LT].first, R = transr[RT].first;
		int64_t Lp = transr[LT].second, Rp = transr[RT].second;
		
		// This removes safety windows, that are a subset of another safety window.
		// Here, this is only the case, if gaps are being used.
		// As we print64_t the safety windows wrt. both strings, we don't want to remove these kind
		// of subsets, though, as we'd lose to one-to-one correspondence between the safety-windows
		// of both strings.
		/*while (outside(L, R)) windows.pop_back();
		if (!inside(L, R)) windows.emplace_back(L, R);*/

		windows.emplace_back(L, R);
		windowsp.emplace_back(Lp, Rp);
	}


	output[i] += std::to_string(windows.size()) + '\n' + proteins[i].descriptor + '\n';
	//std::cout << windows.size() << std::endl;
	for (int64_t k = 0; k < (int64_t) windows.size(); k++) {
		auto [x, y] = windows[k];
		auto [xp, yp] = windowsp[k];
		//std::cout << x << ' ' << y << ' ' << xp << ' ' << yp << '\n';
		output[i] += std::to_string(x) + ' ' + std::to_string(y) + ' ' + std::to_string(xp) + ' ' + std::to_string(yp) + '\n';
	}
	//std::cout << std::flush;
}

signed main(int argc, char **argv) {
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output


	int64_t c;
	while (1) {
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 }, // this does nothing for now
			{ "alpha", required_argument, 0, 'a' },
			{ "delta", required_argument, 0, 'd' },
			{ "costmat", required_argument, 0, 'c' },
			{ "gapcost", required_argument, 0, 'g' },
			{ "startgap", required_argument, 0, 'e' },
			{ "special", required_argument, 0, 's' },
			{ "file", required_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ "threads", required_argument, 0, 'i' },
			{ "reference", required_argument, 0, 'r' },
			{ "approximation", no_argument, 0, 'p' },
			{ 0, 0, 0, 0 }
		};
	
		int option_index = 0;
		c = getopt_long(argc, argv, "a:d:c:g:e:s:f:hi:r:p", long_options, &option_index);
		if (c == -1) break;

		switch (c) {
			case 0:
				if (long_options[option_index].flag != 0)
					break;
				std::cout << "Option " << long_options[option_index].name;
				if (optarg)
					std::cout << " with arg " << optarg;
				std::cout << std::endl;
				break;
			case 'a':
				alpha = std::stof(optarg);
				break;
			case 'd':
				delta = std::stof(optarg);
				break;
			case 'c':
				read_cost_matrix = true;
				cost_matrix_file = optarg;
				break;
			case 'g':
				GAP_COST = atoi(optarg);
				break;
			case 'e':
				START_GAP = atoi(optarg);
				break;
			case 's':
				if (strcmp(optarg, "INF") == 0) ignore_special = true;
				else SP = atoi(optarg);
				break;
			case 'f':
				read_file = true;
				file = optarg;
				break;
			case 'h':
				help_flag = true;
				break;
			case 'i':
				threads = atoi(optarg);
				break;
			case 'r':
				reference = optarg;
				break;
			case 'p':
				use_approx = true;
				break;
			case '?': break;
			default: abort();
		}
	}
	
	if (help_flag) {
		return print64_t_usage(argv, 0);
	}
	if (!read_file) {
		return print64_t_usage(argv, 1);
	}

	if (alpha >= 1.0) {
		std::cerr << "Warning: for alpha values >= 1.0, the program will not return any safety windows.\n";
	} else if (alpha < 0.5) {
		std::cerr << "Warning: for alpha values < 0.5, the program will not behave well defined and might crash.\n";
	}


	if (read_cost_matrix) {
		std::ifstream costmat(cost_matrix_file);
		for (int64_t i = 0; i < 20; i++) for (int64_t j = 0; j <= i; j++) {
			costmat >> cost_matrix[i][j];
			cost_matrix[j][i] = cost_matrix[i][j];
		}
	}

	// TODO: Read these from the LTA array instead
	std::vector<char> special_chars = { 'B', 'X', 'Z' };
	auto contains_special = [&](const std::string &a) {
		for (char c: special_chars)
			if (a.find(c) != std::string::npos) return true;
		return false;
	};

	std::ifstream input(file);
	proteins.clear();
	for (std::string line; std::getline(input, line); ) {
		if ((int64_t) line.size() <= 0) continue;
		if (line[0] == '>') {
			if ((int64_t) proteins.size() > 0 && ignore_special && contains_special(proteins.back().sequence))
				proteins.pop_back();
			proteins.push_back(Protein(line));
		} else {
			proteins.back().sequence += line;
		}
	}
	if ((int64_t) proteins.size() > 0 && ignore_special && contains_special(proteins.back().sequence))
		proteins.pop_back();

	int64_t PS = (int64_t) proteins.size();

	if (PS == 0) {
		std::cout << "Protein sequence list is empty.\n";
		return 2;
	}

	char delim = proteins[0].descriptor.find('|') != std::string::npos ? '|' : ':'; // '|' is used by swiss prot, ':' is used by alpha fold
	// find reference protein
	auto split = [&](std::string str, char delimiter) -> std::vector<std::string> {
		size_t pos = 0;
		std::vector<std::string> ret;
		while ((pos = str.find(delimiter)) != std::string::npos) {
			std::string token = str.substr(0, pos);
			ret.push_back(token);
			str.erase(0, pos + 1); // + 1 for delimiter length
		}
		ret.push_back(str);
		return ret;
	};
	if (reference != "-") {
		bool found = false;
		for (int64_t i = 0; i < PS; i++) {
			if (reference == split(proteins[i].descriptor, delim)[1]) {
				if (found)
					std::cerr << "Reference identity found more than once. Using the last found protein as reference.\n";
				ref = i;
				found = true;
			}
		}
		if (!found) std::cerr << "Reference identity not found, using the first protein as reference.\n";
	}

	// reference protein and amount of proteins in the cluster
	std::cout << ref << ' ' << proteins[ref].sequence << '\n';
	std::cout << proteins[ref].descriptor << '\n' << PS << '\n';
	std::vector<std::string> output(PS);
	random_order.clear();
	for (int64_t i = 0; i < PS; i++) if (i != ref) random_order.push_back(i);
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(random_order.begin(), random_order.end(), g);

	#pragma omp parallel for num_threads(threads)
	for (int64_t j = 0; j < PS - 1; j++)
		if (!use_approx)
			run_case<mpz_class, mpq_class>(j, output);
		else
			run_case<double, double>(j, output);

	for (int64_t i = 0; i < PS; i++) if (i != ref) std::cout << output[i];

	std::cout << verbose_flag << std::endl;
}
