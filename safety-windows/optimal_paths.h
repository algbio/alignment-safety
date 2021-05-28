#include <vector>
#include <map>
#include <string>

struct Dag {
	std::vector<std::vector<int>> adj;
	std::map<std::pair<int, int>, int> trans;
};

// Construct optimal alignment matrix
std::vector<std::vector<int>> opt_alignment(const std::string &a, const std::string &b);

// Find the sub-graph of the optimal alignment matrix with optimal paths
Dag gen_dag(const std::vector<std::vector<int>> &dp, const std::string &a, const std::string &b);
