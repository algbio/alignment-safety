#include <iostream>
#include <iomanip>
#include <vector>
#include <stack>

#define SRC 0
#define DEST (n-1)

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
std::vector<int> amount_paths(std::vector<std::vector<int>> &dag, int sink) {
	int n = (int) dag.size();
	std::vector<int> sorted = topsort(dag);

	std::vector<int> am(n, 0);
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
	std::vector<int> am = amount_paths(dag, DEST);
	std::vector<int> ram = amount_paths(rdag, SRC);

	std::vector<std::vector<double>> ratios(n);
	for (int i = 0; i < n; i++) for (int v: dag[i]) {
		ratios[i].push_back((double) am[v] * ram[i] / am[SRC]);
	}
	return ratios;
}

std::vector<int> arbitrary_path(std::vector<std::vector<int>> &dag) {
	// dest will be missing from path
	std::vector<int> path;
	int n = (int) dag.size();
	for (int current = SRC; current != DEST; current = dag[current][0])
		path.push_back(current);
	return path;
}

void find_path(int src, int dest, std::vector<int> &path, std::vector<std::vector<int>> &dag,
		std::vector<int> &order) {
	// dest will be missing from path
	if (src == dest) return;

	std::function<bool(int)> dfs = [&](int current) {
		if (current == dest) return true;
		if (order[current] > order[dest]) return false;
		path.push_back(current);
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

	std::vector<std::pair<int, int>> needed;
	for (int i = 0; i < n; i++) for (int j = 0; j < (int) dag[i].size(); j++) {
		int v = dag[i][j];
		double d = ratios[i][j];
		if (d > alpha) needed.emplace_back(i, v);
	}
	if (needed.empty()) return arbitrary_path(dag);

	std::vector<int> sorted = topsort(dag);
	std::vector<int> order(n);
	for (int i = 0; i < n; i++) order[sorted[i]] = i;

	std::vector<int> path;
	int last = SRC;
	for (auto [u,v]: needed) {
		find_path(last, u, path, dag, order);
		path.push_back(u);
		last = v;
	}
	find_path(needed.back().second, DEST, path, dag, order);
	path.push_back(DEST);

	return path;
}

int main() {
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
