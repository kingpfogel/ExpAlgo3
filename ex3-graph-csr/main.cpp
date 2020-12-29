#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <execution>

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



bool comp_first(const std::pair<unsigned int, unsigned int> &a, const std::pair<unsigned int, unsigned int> &b)
{
//    return (std::get<0>(a) < std::get<0>(b));
    return a.first < b.first;
}

csr_matrix transpose(csr_matrix matrix)
{
    auto &ind = matrix.ind;
    auto &cols = matrix.cols;
    auto &weights = matrix.weights;

    std::vector<unsigned int> tmp_new_cols(cols.size());
    for(unsigned i = 0; i < ind.size()-1; ++i)
    {
        for(unsigned c = 0; c < ind[i+1]-ind[i]; ++c)
        {
            tmp_new_cols[ind[i]+c] = i;
        }
    }

    std::vector<std::pair<unsigned int, std::pair<unsigned int, float>>> target(cols.size());
    for (unsigned i = 0; i < target.size(); i++){
        target[i] = std::make_pair(cols[i], std::make_pair(tmp_new_cols[i], weights[i]));
    }

    std::sort(target.begin(), target.end(), comp_first);
    std::vector<unsigned int> new_cols;
    std::vector<float> reordered_weights;

    std::transform(target.begin(), target.end(),
                   std::back_inserter(new_cols),
//                   [](auto const& pair){ return std::get<1>(pair); });
                   [](auto const& pair){ return pair.second.first; });
    std::transform(target.begin(), target.end(),
                   std::back_inserter(reordered_weights),
//                   [](auto const& pair){ return std::get<2>(pair); });
                   [](auto const& pair){ return pair.second.second; });



    std::sort(cols.begin(), cols.end());
    unsigned int uniqueCount = std::unique(cols.begin(), cols.end()) - cols.begin();
    std::vector<unsigned int> new_ind;
    new_ind.resize(uniqueCount, 0);
    for(unsigned i = 0; i < cols.size(); ++i){
        std::vector<unsigned int> ones(ind.size()-cols[i]+1 ,1);
        std::transform (new_ind.begin()+cols[i]+1, new_ind.end(), ones.begin(), new_ind.begin()+cols[i]+1, std::plus<int>());
    }

    matrix.weights = reordered_weights;
    matrix.cols = new_cols;
    matrix.ind = new_ind;
    return matrix;
}

int main() {
    csr_matrix m;
    m.cols = {0,1,1,3,2,3,4,5};
    m.weights = {10.,20.,30.,40.,50.,60.,70.,80.};
    m.ind = {0,2,4,7,8};
    for(auto & i: m.ind){
        std::cout << i << " ";
    }
    std::cout << std::endl;

    for(auto & i: m.cols){
        std::cout << i << " ";
    }
    std::cout << std::endl;

    for(auto & i: m.weights){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    return (int)0;
}

int main2(int argc, char **argv) {
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
	return 0;
}
