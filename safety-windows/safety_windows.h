#include <vector>

// Safe the ratios of the path in a vector for the safety_windows function
std::vector<double> find_ratios(std::vector<int> &path, std::vector<std::vector<int>> &dag,
		std::vector<std::vector<double>> &ratios);

// Calculate the safety windows
std::vector<std::pair<int, int>> safety_windows(std::vector<std::vector<int>> &dag,
		std::vector<int> &path, std::vector<double> &ratios,
		double alpha);

