#include "graph.hpp"
#include "util.hpp"
#include <omp.h>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream> 
#include "constants.hpp"

float dynamic_prob_comp(int t_u, int t_neigh){
    assert(t_neigh >= t_u);
    // std::cout << t_neigh << " " << t_u << std::endl;
    return std::exp((t_u - t_neigh) / ALPHA);
}

// Function to calculate Cumulative Distribution Function (CDF)
std::vector<std::vector<float> >* precompute_CDF(Graph& g, const std::vector<std::pair<int, int>>* sortedNeighbors) {
    std::vector<std::vector<float> >* allCDF = new std::vector<std::vector<float> >;
    for (int i = 0; i < g.v_num; i++) {
        std::vector<float>* cdf = new std::vector<float>;
        float prefixSum = 0;
        for (int j = sortedNeighbors[i].size() - 1; j >= 0; j--) {
            prefixSum += dynamic_prob_comp(0, sortedNeighbors[i][j].second);
            (*cdf).push_back(prefixSum);
        }
        (*allCDF).push_back(*cdf);
    }
    return allCDF;
}



int its_sampling(float* p_t, const std::vector<std::pair<int, int> >& neighbors, std::vector<float>& cdf,
                    StdRandNumGenerator& randgen, bool time_decay = true) {
    float t = *p_t;
    int ok_neigh = binarySearchPair(neighbors, t, 0, neighbors.size()-1);
    // if (DEBUG) std::cout << "Closest available neighbor starts from " << ok_neigh << std::endl;
    if (ok_neigh >= neighbors.size()) {
        return -1;
    }
    int dart;
    float upper_bound = cdf[cdf.size()-ok_neigh-1];
    if (DEBUG) std::cout << "Upper " << upper_bound <<std::endl;
    float rand_num = randgen.gen_float(upper_bound);
    if (DEBUG) std::cout << "Rand " << rand_num <<std::endl;
    dart = binarySearch(cdf, rand_num, ok_neigh, cdf.size() - 1) + 1;
    int i_candidate =  cdf.size() - dart;
    assert (i_candidate < neighbors.size() && i_candidate >= 0);
    int t_candidate = neighbors[i_candidate].second;
    int candidate = neighbors[i_candidate].first;
    assert (candidate >= 0);
    if (DEBUG) std::cout << "Dart " << dart << std::endl;
    if (DEBUG) std::cout << "Candidate " << candidate <<std::endl;
    *p_t = t_candidate;
    return candidate;
}




void random_walk(Graph& g, int l, int* walk, std::ostream& output_stream, int& id) {
    StdRandNumGenerator randgen = StdRandNumGenerator();
    int start = randgen.gen(g.v_num);
    walk[0] = start;
    // std::cout << walk[0] << std::endl;
    // std::cout << "ID : " <<id << std::endl;
    float* p_t  = new float;
    *p_t = 0;
    for (int i = 1; i < l; i++) {
        
        walk[i] = its_sampling(p_t, g.get_neighbor_list()[walk[i-1]], 
        (*((std::vector<std::vector<float> > *)(g.precomputed)))[walk[i-1]], randgen, 1);
        if (walk[i] == -1) {
            if (DEBUG) std::cout << "Walk ends in advance" << std::endl;
            if (LOG) output_stream << "Walk ends in advance" << std::endl;
            break;
        }
        if (LOG) std::cout << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " | " << *p_t << std::endl;
        if (SAVE_WALK) output_stream << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " | " << p_t << std::endl;
        
    }
    delete p_t;
}

int main(int argc, char *argv[]){
    const char* file = "data/out.opsahl-ucsocial";
    Graph g = Graph();
    g.load_graph(1900, 0, file);
    g.precomputed = (void*) precompute_CDF(g, (g.get_neighbor_list()));
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