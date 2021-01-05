#include <algorithm>
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

//		if(!seen_header) {
//			seen_header = true;
//			continue;
//		}

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



bool comp_first(const std::tuple<unsigned int, unsigned int, float> &a, const std::tuple<unsigned int, unsigned int, float> &b)
{
    return (std::get<0>(a) < std::get<0>(b));
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

    std::vector<std::tuple<unsigned int, unsigned int, float>> target(cols.size());
    for (unsigned i = 0; i < target.size(); ++i){
          std::make_tuple(cols[i], tmp_new_cols[i], weights[i]);
          target[i] = std::make_tuple(cols[i], tmp_new_cols[i], weights[i]);
    }
    std::sort(target.begin(), target.end(), comp_first);
    std::vector<unsigned int> new_cols;
    std::vector<float> reordered_weights;
    //// do it in one loop
    std::transform(target.begin(), target.end(),
                   std::back_inserter(new_cols),
                   [](auto const& pair){ return std::get<1>(pair); });
    std::transform(target.begin(), target.end(),
                   std::back_inserter(reordered_weights),
                   [](auto const& pair){ return std::get<2>(pair); });


    std::vector<unsigned int> cp_cols;
    cp_cols = cols;
    std::sort(cols.begin(), cols.end());

    unsigned int uniqueCount = std::unique(cols.begin(), cols.end()) - cols.begin();
    std::vector<unsigned int> new_ind;
    new_ind.resize(uniqueCount+1+1, 0);
    for(unsigned i = 0; i < cp_cols.size(); ++i){
        std::vector<unsigned int> ones(new_ind.size()-cp_cols[i]+1 ,1);
        std::transform (new_ind.begin()+cp_cols[i]+1, new_ind.end(), ones.begin(), new_ind.begin()+cp_cols[i]+1, std::plus<int>());
    }
//    matrix.weights = reordered_weights;
//    matrix.cols = new_cols;
//    matrix.ind = new_ind;
    return matrix;
}

template<typename T>
struct prio_cmp {
    //pass it by reference - so that the compare class can be updated...
    //not coming from a c++ background - but feels like not very good practice
    std::vector<float> &prios;

    prio_cmp(std::vector<float> &prios) : prios(prios) {
    };
    auto operator()(const T first, const T second) const {
        return prios[first] < prios[second];
    }
};

std::vector<u_int32_t> dijkstra(const csr_matrix& matrix, int source){
    std::vector<float> dist(matrix.n, 500000.0); // Unknown distance from source to v
    std::vector<u_int32_t> prev(matrix.n, UINT32_MAX-1); //Predecessor of v
    auto &ind = matrix.ind;
    auto &cols = matrix.cols;
    auto &weights = matrix.weights;
    //now update source nodes
    dist[source] = 0;
    prev[source] = 0;

    //hand a reference to the distances to the priority comparison struct
    prio_cmp<float> cmp(dist);
    //100% tilt mode... this took so long to find out how to pass my custom comparison...
    DAryAddressableIntHeap<u_int32_t, 2, decltype(cmp)> Q(cmp);
    //now add all nodes to
    //Q - (would normally not name a variable with a capital letter.. but I somehow wanted to stick to pseudocode
    for(int v = 0; v<matrix.n; ++v){
        Q.push(v);
    }
    u_int32_t u;
    while(!Q.empty()){
        u = Q.extract_top();
        auto colStart = ind[u];
        auto colEnd = ind[u+1];
        for(int j = colStart; j < colEnd; j++){
            //j's are basically all outgoing edges from u: [u]-j-[v]
            auto v = cols[j]; //v
            auto w = weights[j]; //dist
            auto alt = dist[u] + w;
            if(alt < dist[v]){
                dist[v] = alt;
                prev[v] = u;
                Q.update(v);
            }

        }

    }

    return prev;
}

int main(int argc, char **argv) {

//    std::ifstream ins("../cit-patent/cit-patent.edges");
    std::ifstream ins("../roadNet-TX/roadNet-TX.mtx");
//    std::ifstream ins("../foodweb-baydry.konect");
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

    auto prev = dijkstra(mat, 1 );
    auto transposedMatrix = transpose(mat);

    //std::cout << mat.n << " " << mat.m << std::endl;
	return 0;
}

