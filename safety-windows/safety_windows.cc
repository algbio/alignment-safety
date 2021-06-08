#include <vector>
#include <iostream>

#include <gmpxx.h>

#include "safety_windows.h"

std::vector<mpq_class> find_ratios(std::vector<int> &path, std::vector<std::vector<int>> &dag,
		std::vector<std::vector<mpq_class>> &ratios) {
	int k = (int) path.size();
	std::vector<mpq_class> rat_path(k - 1);

	for (int i = 0; i < k - 1; i++) {
		for (int j = 0; j < (int) dag[path[i]].size(); j++) if (dag[path[i]][j] == path[i + 1]) {
			rat_path[i] = ratios[path[i]][j];
		}
	}

	return rat_path;
}

std::vector<std::pair<int, int>> safety_windows(std::vector<std::vector<int>> &dag,
		std::vector<int> &path, std::vector<mpq_class> &ratios,
		mpq_class alpha) {
	if (path.empty()) return {};
	int k = (int) path.size();

	std::vector<std::pair<int, int>> windows;
	mpq_class a = 1;
	for (int L = 0, R = 0; R < k; a *= ratios[R], R++) {
		while (L < R && a <= alpha) {
			a /= ratios[L];
			L++;
		}
		if (L < R) windows.emplace_back(path[L], path[R]);
	}
	return windows;
}
