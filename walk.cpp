#include "graph.hpp"
#include "util.hpp"
#include <omp.h>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream> 

int binarySearch(const std::vector<std::pair<int, int> >& arr, int value) {
    int low = 0;
    int high = arr.size() - 1;
    int result = -1; // Start with -1 to indicate not found

    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (arr[mid].second <= value) {
            result = mid;  // Update result to current mid (potential answer)
            low = mid + 1; // Move to the right half to find larger index
        } else {
            high = mid - 1; // Move to the left half
        }
        // std::cout << "Low: " << low << " High: " << high << std::endl;
    }

    // return result; // Returns -1 if not found
    if (result != -1 && arr[result].second == value){
        return result;
    }
    // std::cout << "Low " << low << " High " << high << std::endl;
    assert(low == arr.size() || arr[low].second > value);
    assert(low == 0 || arr[low-1].second < value);
    // assert(low < arr.size());
    return low;
}

float dynamic_prob_comp(int t_u, int t_neigh, float alpha){
    assert(t_neigh >= t_u);
    // std::cout << t_neigh << " " << t_u << std::endl;
    return std::exp((t_u - t_neigh) / alpha);
}

int reject_sampling(int* t, std::vector<std::pair<int, int> >& neighbors, 
                    func prob, StdRandNumGenerator& randgen, bool time_decay = true) {
    float alpha = 50000.0;
    // find neighbors of u with timestamp > t
    int ok_neigh = binarySearch(neighbors, *t);
    // std::cout << "Closest available neighbor starts from " << ok_neigh << std::endl;
    if (ok_neigh >= neighbors.size()) {
        return -1;
    }
    int samples = neighbors.size() - ok_neigh;
    // bool succeed = false;
    int t_latest = neighbors[ok_neigh].second;
    
    // compute prob upperbound (just assume the latest neighbor is sampled with the largest probability)
    float prob_upperbound = prob(*t, t_latest, alpha);
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
        float prob_ac = prob(*t, t_candidate, alpha) / prob_upperbound;
        // std::cout << prob_upperbound << " " << prob(*t, t_candidate, alpha) << std::endl;
        // std::cout << "Accepted prob: " << prob_ac << std::endl;
        if (randgen.gen_float(1.0) < prob_ac) {
            *t = t_candidate;
            return candidate;
        }
    }
    return -1;

}

void random_walk(Graph& g, int l, int* walk, std::ostream& output_stream, int id) {
    StdRandNumGenerator randgen = StdRandNumGenerator();
    int start = randgen.gen(g.v_num);
    walk[0] = start;
    int* p_t = new int;
    for (int i = 1; i < l; i++) {
        walk[i] = reject_sampling(p_t, g.neighbor_list[walk[i-1]], dynamic_prob_comp, randgen, 1);
        // std::cout << walk << std::endl << p_t << std::endl << &(g.neighbor_list) << std::endl << &dynamic_prob_comp << std::endl << & randgen << std::endl;
        if (walk[i] == -1) {
            std::cout << "Walk ends in advance" << std::endl;
            output_stream << "Walk ends in advance" << std::endl;
            break;
        }
        std::cout << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " " << *p_t << std::endl;
        output_stream << "WALKER ID " << id <<" : Random walk " << i << " is " << walk[i] << " " << *p_t << std::endl;
    }
    // std::cout<< std::endl;
    delete p_t;

}

int main(){
    // std::string test = "233 this is 244 a test";
    // int end;
    // int num  = get_num(test, 0, &end);
    // std::cout << "test num is " << num << std::endl;
    // std::cout << "end is " << end << std::endl;
    // int start = end;
    // num  = get_num(test, start, &end);
    // std::cout << "test num is " << num << std::endl;
    // std::cout << "end is " << end << std::endl;
    const char* file = "data/karate.txt";
    Graph g = Graph();
    g.load_graph(34, 0, file);
    
    std::ofstream file_stream;                                                  
    file_stream.open ("walks_out.txt"); 
    // g.print_info();
    // int tmp = binarySearch(g.neighbor_list[0], 150000);
    // std::cout << "Tmp " << tmp << std::endl;
    Timer timer;
    int num_walk = 8;
    int l_walk = 10;
    int** walks = new int*[num_walk];
    for (int i = 0; i < num_walk; i++){
        walks[i] = new int[l_walk];
    }
    timer.restart();
    #pragma omp parallel for
    for (int k = 0; k < num_walk; k++){
        random_walk(g, l_walk, walks[k], file_stream, k);
    }
    float time_lapse = timer.duration();
    std::cout << "Time lapse: " << time_lapse << std::endl;
    for (int i = 0; i < num_walk; i++){
        delete[] walks[i];
    }
    delete[] walks;
    return 0;
}