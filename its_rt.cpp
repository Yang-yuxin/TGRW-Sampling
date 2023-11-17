#include "graph.hpp"
#include "util.hpp"
#include <omp.h>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream> 

// Function to calculate Cumulative Distribution Function (CDF)
std::vector<int> calculateCDF(const std::vector<std::pair<int, int>>& neighbors, int ok_neigh) {
    std::vector<int> cdf;
    // cdf.push_back(0);
    int prefixSum = 0;
    for (int i = ok_neigh ; i < neighbors.size(); i++) {
        prefixSum += neighbors[i].second;
        cdf.push_back(prefixSum);
        // std::cout << neighbors[i].second <<std::endl;
    }
    // std::cout << "Cumulative Distribution Function (CDF): ";
    // for (int value : cdf) {
    //     std::cout << value << " ";
    // }
    // std::cout << std::endl;
    return cdf;
}

int its_sampling(int* t, const std::vector<std::pair<int, int> >& neighbors, 
                    func prob, StdRandNumGenerator& randgen, bool time_decay = true) {
    // find neighbors of u with timestamp > t
    int ok_neigh = binarySearchPair(neighbors, *t, 0, neighbors.size() - 1);
    // std::cout << "Closest available neighbor starts from " << ok_neigh << std::endl;
    if (ok_neigh >= neighbors.size()) {
        return -1;
    }
    std::vector<int> cdf;
    int dart;
    cdf = calculateCDF(neighbors, ok_neigh);
    int upper_bound = cdf.back();
    int rand_num = randgen.gen(upper_bound);
    int candidate = -1;
    for (int i = ok_neigh ; i < cdf.size(); i++) {
        if (rand_num <= cdf[i]){
            dart = i;
        }
    }
    int i_candidate = dart;
    if (DEBUG) std::cout << "Dart " << i_candidate << std::endl;
    if (DEBUG) std::cout << "Neighbors " << neighbors.size() << std::endl;
    assert (i_candidate < neighbors.size());
    int t_candidate = neighbors[i_candidate].second;
    candidate = neighbors[i_candidate].first;
    *t = t_candidate;
    return candidate;
}

float dynamic_prob_comp(int t_u, int t_neigh){
    assert(t_neigh >= t_u);
    // std::cout << t_neigh << " " << t_u << std::endl;
    return std::exp((t_u - t_neigh) / ALPHA);
}


void random_walk(Graph& g, int l, int* walk, std::ostream& output_stream, int& id) {
    StdRandNumGenerator randgen = StdRandNumGenerator();
    int start = randgen.gen(g.v_num);
    walk[0] = start;
    // std::cout << walk[0] << std::endl;
    // std::cout << "ID : " <<id << std::endl;
    int* p_t = new int;
    *p_t = 0;
    for (int i = 1; i < l; i++) {
        // std::cout << "ID : " << id << std::endl;
        // std::cout << "time : " << *p_t << std::endl;
        walk[i] = its_sampling(p_t, g.get_neighbor_list()[walk[i-1]], dynamic_prob_comp, randgen, 1);
        // walk[i] = reject_sampling(p_t, g.neighbor_list[walk[i-1]], dynamic_prob_comp, randgen, 1);
        // std::cout << walk << std::endl << p_t << std::endl << &(g.neighbor_list) << std::endl << &dynamic_prob_comp << std::endl << & randgen << std::endl;
        if (walk[i] == -1) {
            if (LOG) std::cout << "Walk ends in advance" << std::endl;
            if (LOG) output_stream << "Walk ends in advance" << std::endl;
            break;
        }
        if (LOG){
            std::cout << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " | " << *p_t << std::endl;
            output_stream << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " | " << *p_t << std::endl;
        }
    }
    // std::cout<< std::endl;
    delete p_t;

}

int main(int argc, char *argv[]){
    const char* file = "data/out.opsahl-ucsocial";
    Graph g = Graph();
    g.load_graph(1900, 0, file);
    g.check_graph();
    std::ofstream file_stream;                                                  
    file_stream.open ("walks_out.txt"); 
    Timer timer;
    // int num_walk = 10000;
    // int walk_length = 10;
    assert(argc >= 4);
    int num_walk = std::stoi(argv[1]);
    int walk_length = std::stoi(argv[2]);
    int trials = std::stoi(argv[3]);
    int l_walk = walk_length + 1;   //  (l_walk -1) walks - we count the vertices
    int** walks = new int*[num_walk];
    for (int i = 0; i < num_walk; i++){
        walks[i] = new int[l_walk];
    }
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(64); // Use 4 threads for all consecutive parallel regions
    float total_time = 0.0;
    for (int tr = 0; tr < trials; tr++){
        timer.restart();
        #pragma omp parallel for
        for (int k = 0; k < num_walk; k++){
            random_walk(g, l_walk, walks[k], file_stream, k);
        }
        total_time += timer.duration();

    }
    std::cout << "Time lapse: " << total_time / trials << std::endl;
    for (int i = 0; i < num_walk; i++){
        delete[] walks[i];
    }
    delete[] walks;
    
    // delete[] g.precomputed;
    return 0;
}