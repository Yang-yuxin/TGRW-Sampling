#include "graph.hpp"
#include "util.hpp"
#include <omp.h>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream> 

float dynamic_prob_comp(int t_u, int t_neigh){
    assert(t_neigh >= t_u);
    // std::cout << t_neigh << " " << t_u << std::endl;
    return std::exp((t_u - t_neigh) / ALPHA);
}

int reject_sampling(int* t, const std::vector<std::pair<int, int> >& neighbors, 
                    func prob, StdRandNumGenerator& randgen, bool time_decay = true) {
    // find neighbors of u with timestamp > t
    int ok_neigh = binarySearchPair(neighbors, *t, 0, neighbors.size() - 1);
    // std::cout << "Closest available neighbor starts from " << ok_neigh << std::endl;
    if (ok_neigh >= neighbors.size()) {
        return -1;
    }
    int samples = neighbors.size() - ok_neigh;
    // bool succeed = false;
    int t_latest = neighbors[ok_neigh].second;
    
    float prob_upperbound = prob(*t, t_latest);
    // std::cout << "now time and latest neighbor time " << *t << " " << t_latest << std::endl;
    int candidate = -1;
    while (true){
        // throw a dart uniformly
        int dart = randgen.gen(samples);
        int i_candidate = ok_neigh + dart;
        assert (i_candidate < neighbors.size());
        int t_candidate = neighbors[i_candidate].second;
        int candidate = neighbors[i_candidate].first;
        // compute accept probability
        // std::cout << "Candidate neighbor time: " << t_candidate << std::endl;
        assert (t_candidate >= t_latest);
        float prob_ac = prob(*t, t_candidate) / prob_upperbound;
        // std::cout << prob_upperbound << " " << prob(*t, t_candidate, alpha) << std::endl;
        // std::cout << "Accepted prob: " << prob_ac << std::endl;
        if (randgen.gen_float(1.0) < prob_ac) {
            *t = t_candidate;
            return candidate;
        }
    }
    return -1;

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
        walk[i] = reject_sampling(p_t, g.get_neighbor_list()[walk[i-1]], dynamic_prob_comp, randgen, 1);
        // std::cout << walk << std::endl << p_t << std::endl << &(g.neighbor_list) << std::endl << &dynamic_prob_comp << std::endl << & randgen << std::endl;
        if (walk[i] == -1) {
            if (LOG) std::cout << "Walk ends in advance" << std::endl;
            if (LOG) output_stream << "Walk ends in advance" << std::endl;
            break;
        }
        if (LOG) std::cout << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " | " << *p_t << std::endl;
        if (LOG) output_stream << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " | " << *p_t << std::endl;
    }
    // std::cout<< std::endl;
    delete p_t;

}

int main(int argc, char *argv[]){
    // std::string test = "233 this is 244 a test";
    // int end;
    // int num  = get_num(test, 0, &end);
    // std::cout << "test num is " << num << std::endl;
    // std::cout << "end is " << end << std::endl;
    // int start = end;
    // num  = get_num(test, start, &end);
    // std::cout << "test num is " << num << std::endl;
    // std::cout << "end is " << end << std::endl;
const char* file = "data/out.opsahl-ucsocial";
    Graph g = Graph();
    g.load_graph(1900, 0, file);
    g.check_graph();

    std::ofstream file_stream;                                                  
    file_stream.open ("walks_out.txt"); 
    // g.print_info();
    // int tmp = binarySearch(g.neighbor_list[0], 150000);
    // std::cout << "Tmp " << tmp << std::endl;
    Timer timer;
    assert(argc >= 5);
    int num_walk = std::stoi(argv[1]);
    int walk_length = std::stoi(argv[2]);
    int trials = std::stoi(argv[3]);
    int num_threads = std::stoi(argv[4]);
    int l_walk = walk_length + 1;   //  (l_walk -1) walks - we count the vertices
    int** walks = new int*[num_walk];
    for (int i = 0; i < num_walk; i++){
        walks[i] = new int[l_walk];
    }
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(num_threads); // Use 4 threads for all consecutive parallel regions
    float total_time = 0.0;
    for (int tr = 0; tr < trials; tr++){    
        timer.restart();
        #pragma omp parallel for
        for (int k = 0; k < num_walk; k++){
            random_walk(g, l_walk, walks[k], file_stream, k);
        }
    total_time += timer.duration();

    }
    // float time_lapse = timer.duration();
    std::cout << "Time lapse: " << total_time / trials << std::endl;
    for (int i = 0; i < num_walk; i++){
        delete[] walks[i];
    }
    delete[] walks;
    return 0;
}