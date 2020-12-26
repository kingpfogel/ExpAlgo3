#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

template<typename F>
void read_graph_unweighted(std::istream &ins, F fn) {
	std::string line;
	bool seen_header = false;
	while (std::getline(ins, line)) {
		if(line.empty())
			continue;
		if(line.front() == '%')
			continue;

		std::istringstream ls(line);
		unsigned int u, v;
		if (!(ls >> u >> v))
			throw std::runtime_error("Parse error while reading input graph");

		if(!seen_header) {
			seen_header = true;
			continue;
		}

		fn(u, v);
	}
}

struct csr_matrix {
	unsigned int n;
	unsigned int m;
	std::vector<unsigned int> ind;
	std::vector<unsigned int> cols;
	std::vector<float> weights;
};

csr_matrix coordinates_to_csr(unsigned int n,
		std::vector<std::tuple<unsigned int, unsigned int, float>> cv) {
	unsigned int m = cv.size();

	csr_matrix mat;
	mat.n = n;
	mat.m = m;
	mat.ind.resize(n + 1);
	mat.cols.resize(m);
	mat.weights.resize(m);

	// Count the number of neighbors of each node.
	for(auto ct : cv) {
		auto u = std::get<0>(ct);
		++mat.ind[u];
	}

	// Form the prefix sum.
	for(unsigned int i = 1; i <= n; ++i)
		mat.ind[i] += mat.ind[i - 1];
	assert(mat.ind[n] == m);

	// Insert the entries of the matrix in reverse order.
	for(auto it = cv.rbegin(); it != cv.rend(); ++it) {
		auto u = std::get<0>(*it);
		auto v = std::get<1>(*it);
		auto weight = std::get<2>(*it);
		mat.cols[mat.ind[u] - 1] = v;
		mat.weights[mat.ind[u] - 1] = weight;
		--mat.ind[u];
	}

	return mat;
}

int main(int argc, char **argv) {
	std::ifstream ins("foodweb-baydry.konect");
	std::vector<std::tuple<unsigned int, unsigned int, float>> cv;
	std::mt19937 prng{42};
	std::uniform_real_distribution<float> distrib{0.0f, 1.0f};
	read_graph_unweighted(ins, [&] (unsigned int u, unsigned int v) {
		// Generate a random edge weight in [a, b).
		cv.push_back({u, v, distrib(prng)});
	});

	// Determine n as the maximal node ID.
	unsigned int n = 0;
	for(auto ct : cv) {
		auto u = std::get<0>(ct);
		auto v = std::get<1>(ct);
		if(u > n)
			n = u;
		if(v > n)
			n = v;
	}

	auto mat = coordinates_to_csr(n, std::move(cv));

	std::cout << mat.n << " " << mat.m << std::endl;
}
