#include <omp.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <sstream>
#include <cstring>
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

template<typename T>
struct prio_cmp {
    std::vector<float> &prios;
    prio_cmp(std::vector<float> &prios) : prios(prios) {};
    auto operator()(const T first, const T second) const {
        return prios[first] < prios[second];
    }
};


//7b implementation
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
    DAryAddressableIntHeap<u_int32_t, 2, decltype(cmp)> Q(cmp);
    //now add all nodes to
    //Q - (would normally not name a variable with a capital letter.. but I somehow wanted to stick to pseudocode
    for(unsigned int v = 0; v<matrix.n; ++v){
        Q.push(v);
    }
    u_int32_t u;
    while(!Q.empty()){
        u = Q.extract_top();
        auto colStart = ind[u];
        auto colEnd = ind[u+1];
        for(unsigned int j = colStart; j < colEnd; j++){
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
    int random_runs = 1;
//    string file = "";
    std::string file;
    const char *arg;
    char **p = argv + 1;

    auto handle_unary_option = [&] (const char *name) -> bool {
        assert(*p);
        if(std::strcmp(*p, name))
            return false;
        ++p;
        if(!(*p))
            std::cerr << "expected argument for unary option" << std::endl;
        arg = *p;
        ++p;
        return true;
    };

    while(*p && !std::strncmp(*p, "--", 2)) {
        if(handle_unary_option("--random-runs")) {
            random_runs = atoi(arg);
        }else if(handle_unary_option("--file")) {
            file = arg;
        }else{
        }
    }
    if(*p)
        std::cerr << "unexpected arguments" << std::endl;

    std::ifstream ins(file);
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

    auto mat = coordinates_to_csr(n+1, std::move(cv));

    std::mt19937 prng2{14};
    std::uniform_int_distribution<int32_t> distrib2{1, 20000000};

    auto start = std::chrono::high_resolution_clock::now();
    for(auto i = 0; i < random_runs; ++i){
        auto prev = dijkstra(mat, distrib2(prng2));
    }
    auto t = std::chrono::high_resolution_clock::now() - start;

    std::cout << "random_runs: " << random_runs << std::endl;
    std::cout << "file: " << file << std::endl;
    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t).count()
              << " # ms" << std::endl;

    return 0;
}
