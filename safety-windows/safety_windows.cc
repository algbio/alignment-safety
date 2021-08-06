#include <vector>
#include <iostream>
#include <unordered_map>

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

std::vector<std::pair<int, int>> safety_windows(std::vector<mpq_class> &ratios,
		std::vector<int> &path, mpq_class alpha) {
	int k = (int) ratios.size();
	if (k == 0) return {};
	ratios.push_back(1);

	std::unordered_map<int, int> order;
	for (int i = 0; i <= k; i++) order[path[i]] = i;

	std::vector<std::pair<int, int>> windows;
	mpq_class a = 1;
	auto outside = [&](const int &L, const int &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return order[bL] >= order[L] && order[bR] <= order[R];
	};
	auto inside = [&](const int &L, const int &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return order[L] >= order[bL] && order[R] <= order[bR];
	};
	for (int L = 0, R = 0; R <= k; a *= ratios[R], R++) {
		while (L < R && a <= alpha) {
			a /= ratios[L];
			L++;
		}
		while (outside(path[L], path[R])) windows.pop_back();
		if (L < R && !inside(path[L], path[R])) windows.emplace_back(path[L], path[R]);
	}
	return windows;
}
