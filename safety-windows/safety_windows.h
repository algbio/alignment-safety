#include <vector>
#include <gmpxx.h>

// Save the ratios of the path in a vector for the safety_windows function
template<class K>
std::vector<K> find_ratios(std::vector<int64_t> &path, std::vector<std::vector<int64_t>> &dag,
		std::vector<std::vector<K>> &ratios) {
	int64_t k = (int64_t) path.size();
	std::vector<K> rat_path(k - 1);

	for (int64_t i = 0; i < k - 1; i++) {
		for (int64_t j = 0; j < (int64_t) dag[path[i]].size(); j++) if (dag[path[i]][j] == path[i + 1]) {
			rat_path[i] = ratios[path[i]][j];
		}
	}

	return rat_path;
}

// Calculate the safety windows
template<class K>
std::vector<std::pair<int64_t, int64_t>> safety_windows(std::vector<K> &ratios,
		std::vector<int64_t> &path, K alpha) {
	int64_t k = (int64_t) ratios.size();
	if (k == 0) return {};
	ratios.push_back(1);

	std::unordered_map<int64_t, int64_t> order;
	for (int64_t i = 0; i <= k; i++) order[path[i]] = i;

	std::vector<std::pair<int64_t, int64_t>> windows;
	K a = 1;
	auto outside = [&](const int64_t &L, const int64_t &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return order[bL] >= order[L] && order[bR] <= order[R];
	};
	auto inside = [&](const int64_t &L, const int64_t &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return order[L] >= order[bL] && order[R] <= order[bR];
	};
	for (int64_t L = 0, R = 0; R <= k; a *= ratios[R], R++) {
		while (L < R && a <= alpha) {
			a /= ratios[L];
			L++;
		}
		while (outside(path[L], path[R])) windows.pop_back();
		if (L < R && !inside(path[L], path[R])) windows.emplace_back(path[L], path[R]);
	}
	return windows;
}
