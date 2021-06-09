#include <vector>
#include <gmpxx.h>

// Safe the ratios of the path in a vector for the safety_windows function
std::vector<mpq_class> find_ratios(std::vector<int> &path, std::vector<std::vector<int>> &dag,
		std::vector<std::vector<mpq_class>> &ratios);

// Calculate the safety windows
std::vector<std::pair<int, int>> safety_windows(std::vector<std::vector<int>> &dag,
		std::vector<int> &path, std::vector<mpq_class> &ratios,
		mpq_class alpha);

