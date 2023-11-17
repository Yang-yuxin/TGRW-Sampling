#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <stdint.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include "type.hpp"
#include "constants.hpp"

int get_num(std::string s, int start, int* end){
    int len_s = s.size();
    unsigned int i = start;
    int num = 0;
    while (i < len_s) {
        if (s[i] >= '0' && s[i] <= '9') {
            num *= 10;
            num += s[i] - '0';
            i++;
        }
        else{
            if (i > start) {
                break;
            }
            i++;
            continue;
        }
    }
    *end = i;
    return num;
}

class Graph 
{
private:
    std::vector<std::pair<int, int> >* neighbor_list;
public:
    vertex_id_t v_num;
    edge_id_t e_num;
    
    void* precomputed;
    int* deg_list;
    Graph() {
        this->v_num = 0;
        this->e_num = 0;
        this->neighbor_list = nullptr;
        this->deg_list = nullptr;
    }
    ~Graph() {
        if (this->neighbor_list != nullptr) {
            delete[] this->neighbor_list;
            delete[] this->deg_list;
        }
    }
    
    const std::vector<std::pair<int, int> >* get_neighbor_list(){
        return this->neighbor_list;
    }

    void print_info(){
        std::cout << "Graph node number: " << this->v_num << std::endl;
        std::cout << "Graph edge number: " << this->e_num << std::endl;
        for (int i = 0; i < this->v_num; i++) {
            std::cout << i << std::endl;
            for (int j = 0; j < this->deg_list[i]; j++) {
                std::cout << '\t' << this->neighbor_list[i][j].first << '\t' << this->neighbor_list[i][j].second << std::endl; // neighbor timestamp
            }
        }
    }

    bool check_graph(){
        bool error = 0;
        for (int i = 0; i < this->v_num;i++){
            for (int j = 0; j < neighbor_list[i].size(); j++) {
                if (this->neighbor_list[i][j].first < 0 || this->neighbor_list[i][j].first >= this->v_num){
                    error = 1;
                    if (DEBUG) std::cout << "Node " << i << " has a bad neighbor " << j <<std::endl;
                }
            }
        }
        if (DEBUG && error) std::cout << "Graph error " << error << std::endl;
        return error;
    }

    void load_graph(vertex_id_t v_num_param, bool undirected, const char* graph_path){
        this->v_num = v_num_param;
        this->neighbor_list = new std::vector<std::pair<int, int> >[v_num_param];
        this->deg_list = new int[v_num_param];
        if (DEBUG) std::cout << "Allocated Graph Memory" << std::endl;
        
        std::ifstream f(graph_path);
        std::string line;
        int e_num = 0;
        int start, end, s, t, weight, timestamp;
        
        while (std::getline(f, line)){
            int i_line = line[0] - '0';
            if (!(i_line >= 0 && i_line <= 9)){
                continue;
            }
            e_num++;
            start = 0;
            s = get_num(line, start, &end);
            start = end;
            t = get_num(line, start, &end);
            start = end;
            weight = get_num(line, start, &end);
            start = end;
            timestamp = get_num(line, start, &end);
            if (undirected){
                this->deg_list[t]++;
                this->neighbor_list[t].push_back(std::pair<int, int>(s, timestamp));
            }
            this->deg_list[s]++;
            this->neighbor_list[s].push_back(std::pair<int, int>(t, timestamp));
        }
        if (DEBUG) std::cout << "Loaded Graph " << std::endl;
        this->e_num = e_num;
        struct sort_pred {
            bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
                return left.second < right.second;
            }
        };
        for (int i = 0; i < v_num_param; i++){
            std::sort(this->neighbor_list[i].begin(), this->neighbor_list[i].end(), sort_pred());
        }
        if (DEBUG) std::cout << "Sorted Graph " << std::endl;
        if (DEBUG) this->check_graph();
    }
};
