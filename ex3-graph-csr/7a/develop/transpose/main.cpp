#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "d_ary_addressable_int_heap.hpp"

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

//7a)
bool comp_first(const std::tuple<unsigned int, unsigned int, float> &a, const std::tuple<unsigned int, unsigned int, float> &b)
{
    return (std::get<0>(a) < std::get<0>(b));
}

csr_matrix transpose(csr_matrix matrix){
    auto &ind = matrix.ind;
    auto &cols = matrix.cols;
    auto &weights = matrix.weights;
    std::vector<unsigned int> tmp_new_cols(cols.size());
    #pragma omp parallel for
    for(unsigned i = 0; i < ind.size()-1; ++i)
    {
        for(unsigned j = ind[i]; j < ind[i+1]; ++j)
        {
            tmp_new_cols[j] = i;
        }
    }
    std::vector<std::tuple<unsigned int, unsigned int, float>> target(cols.size());
    #pragma omp parallel for
    for (unsigned i = 0; i < target.size(); ++i){
        std::make_tuple(cols[i], tmp_new_cols[i], weights[i]);
        target[i] = std::make_tuple(cols[i], tmp_new_cols[i], weights[i]);
    }
    std::sort(target.begin(), target.end(), comp_first);
    std::vector<unsigned int> new_cols(target.size());
    std::vector<float> reordered_weights(target.size());
    //// do it in one loop
    #pragma omp parallel for
    for(unsigned int t = 0; t < target.size(); ++t){
        new_cols[t] = std::get<1>(target[t]);
        reordered_weights[t] = std::get<2>(target[t]);
    }

    std::vector<unsigned int> new_ind;
    int real_max = *std::max_element(cols.begin(), cols.end());
    new_ind.resize(real_max+2, 0);
    for(unsigned i = 0; i < cols.size(); ++i){
        std::vector<unsigned int> ones(new_ind.size()-cols[i]+1 ,1);
        std::transform (new_ind.begin()+cols[i]+1, new_ind.end(), ones.begin(), new_ind.begin()+cols[i]+1, std::plus<int>());
    }
    csr_matrix transposed;
    transposed.m = matrix.m;
    transposed.n = real_max;
    transposed.weights = reordered_weights;
    transposed.cols = new_cols;
    transposed.ind = new_ind;
    return transposed;
}


int main(int argc, char **argv) {
    std::ifstream ins("../cit-patent.edges");
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
    //somehow this did not work correctly without the +1...
    auto mat = coordinates_to_csr(n+1, std::move(cv));


    auto start = std::chrono::high_resolution_clock::now();
    csr_matrix transposed = transpose(mat);
    auto t = std::chrono::high_resolution_clock::now() - start;

    std::cout << "threads: " << omp_get_num_threads << std::endl;
    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t).count()
              << " # ms" << std::endl;

    //csr_matrix transposed_again = transpose(transposed);
    //sanity control
    //sort transposed.cols and transposed_again.cols
    //compare transposed.cols == transposed_again.cols
    //compare transposed.ind == transposed_again.ind
    return 0;
}

