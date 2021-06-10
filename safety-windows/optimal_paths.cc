#include <vector>
#include <map>
#include <queue>
#include <string>

#include "optimal_paths.h"

// TODO: Optimize amount of space of the return matrix
std::vector<std::vector<int>> opt_alignment(const std::string &a, const std::string &b) {
	int n = (int) a.size();
	int m = (int) b.size();

	const int SAME_COST = 0;
	const int DIFF_COST = 1;
	const int GAP_COST = 1;

	std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, (1 << 30)));
	for (int i = 0; i <= n; i++) dp[i][0] = i;
	for (int i = 0; i <= m; i++) dp[0][i] = i;
	for (int i = 1; i <= n; i++) for (int j = 1; j <= m; j++) {
		dp[i][j] = std::min(dp[i][j], dp[i - 1][j - 1] +
			(a[i - 1] == b[j - 1] ? SAME_COST : DIFF_COST));
		dp[i][j] = std::min(dp[i][j], dp[i - 1][j] + GAP_COST);
		dp[i][j] = std::min(dp[i][j], dp[i][j - 1] + GAP_COST);
	}
	return dp;
}

Dag gen_dag(const std::vector<std::vector<int>> &dp, const std::string &a,
		const std::string &b) {
	int n = (int) dp.size() - 1;
	assert(n > 0);
	int m = (int) dp[0].size() - 1;
	assert(m > 0);

	const int SAME_COST = 0;
	const int DIFF_COST = 1;
	const int GAP_COST = 1;

	// check if pair of points is contained in the dp matrix
	std::function<bool(int, int)> valid = [&](int x, int y) {
		return x >= 0 && x <= n && y >= 0 && y <= m;
	};

	// build the adjacency list
	std::vector<std::vector<int>> adj(1);
	std::map<std::pair<int, int>, int> trans; // translate to index
	std::map<int, std::pair<int, int>> transr; // translate index to pair
	trans[std::make_pair(n, m)] = 0;
	transr[0] = std::make_pair(n, m);
	int current = 0;
	std::queue<std::pair<int, int>> q;
	q.emplace(n, m);
	
	auto add_node = [&](const int &x, const int &y) {
		if (trans.find(std::make_pair(x, y)) == trans.end()) {
			trans[std::make_pair(x, y)] = ++current;
			transr[current] = std::make_pair(x, y);
			adj.push_back(std::vector<int>());
			q.emplace(x, y);
		}
	};

	while (!q.empty()) {
		auto [x, y] = q.front();
		q.pop();
		if (valid(x - 1, y - 1) && dp[x][y] == dp[x - 1][y - 1] +
				(a[x - 1] == b[y - 1] ? SAME_COST : DIFF_COST)) {
			add_node(x - 1, y - 1);
			adj[trans[std::make_pair(x - 1, y - 1)]].push_back(trans[std::make_pair(x, y)]);
		}
		if (valid(x - 1, y) && dp[x][y] == dp[x - 1][y] + GAP_COST) {
			add_node(x - 1, y);
			adj[trans[std::make_pair(x - 1, y)]].push_back(trans[std::make_pair(x, y)]);
		}
		if (valid(x, y - 1) && dp[x][y] == dp[x][y - 1] + GAP_COST) {
			add_node(x, y - 1);
			adj[trans[std::make_pair(x, y - 1)]].push_back(trans[std::make_pair(x, y)]);
		}
	}

	// source should have index 0 and dest should have index current (it's
	// still the other way round)
	std::swap(trans[std::make_pair(0, 0)], trans[std::make_pair(n, m)]);
	swap(transr[0], transr[current]);
	std::swap(adj[0], adj[current]);
	for (int i = 0; i <= current; i++) {
		for (int j = 0; j < (int) adj[i].size(); j++) {
			if (adj[i][j] == 0) adj[i][j] = current;
			else if (adj[i][j] == current) adj[i][j] = 0;
		}
	}

	return { adj, trans, transr };
}
