#pragma once

#include <random>
#include <chrono>
#include "type.hpp"


bool isClose(float a, float b) {
    if ((a > b && a - b < 1e-5) || (a < b && b - a < 1e-5)){
        return true;
    } 
    return false;
}

class RandNumGenerator
{
public:
    virtual vertex_id_t gen(vertex_id_t upper_bound) = 0;
    virtual float gen_float(float upper_bound) = 0;
    virtual ~RandNumGenerator() {}
};

class StdRandNumGenerator : public RandNumGenerator
{
    std::random_device *rd;
    std::mt19937 *mt;
public:
    StdRandNumGenerator()
    {
        rd = new std::random_device();
        mt = new std::mt19937((*rd)());
    }
    ~StdRandNumGenerator()
    {
        delete mt;
        delete rd;
    }
    vertex_id_t gen(vertex_id_t upper_bound)
    {
        std::uniform_int_distribution<vertex_id_t> dis(0, upper_bound - 1);
        return dis(*mt);
    }
    float gen_float(float upper_bound)
    {
        if (isClose(upper_bound, 0.0)) {
            return 0.0;
        }
        std::uniform_real_distribution<float> dis(0.0, upper_bound);
        return dis(*mt);
    }
};

//Timer is used for performance profiling
class Timer
{
    std::chrono::time_point<std::chrono::system_clock> _start = std::chrono::system_clock::now();
public:
    void restart()
    {
        _start = std::chrono::system_clock::now();
    }
    double duration()
    {
        std::chrono::duration<double> diff = std::chrono::system_clock::now() - _start;
        return diff.count();
    }
    static double current_time()
    {
        std::chrono::duration<double> val = std::chrono::system_clock::now().time_since_epoch();
        return val.count();
    }
};

int binarySearchPair(const std::vector<std::pair<int, int>>& arr, int value, int low, int high) {
    // int low = 0;
    // int high = arr.size() - 1;
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

    return result; // Returns -1 if not found
    // if (result != -1 && arr[result].second == value){
    //     return result;
    // }
    // // std::cout << "Low " << low << " High " << high << std::endl;
    // // assert(low == arr.size() || arr[low].second > value);
    // assert(low == 0 || arr[low-1].second < value);
    // // assert(low < arr.size());
    // return low;
}



int binarySearch(const std::vector<float>& arr, float value, int low, int high) {
    // int low = 0;
    // int high = arr.size() - 1;
    int result = -1; // Start with -1 to indicate not found
    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (arr[mid] <= value) {
            result = mid;  // Update result to current mid (potential answer)
            low = mid+1; // Move to the right half to find larger index
        } else {
            high = mid-1; // Move to the left half
        }
        // std::cout << "Low: " << low << " High: " << high << std::endl;
    }

    return result; // Returns -1 if not found
    // if (result != -1){
    //     return result;
    // }
    // return -1;
    // std::cout << "Low " << low << " High " << high << std::endl;
    // std::cout << "arr size " << arr.size() << " value " << value << std::endl;
    // // assert(low == arr.size() || arr[low] > value);
    // assert(low == 0 || arr[low-1] < value);
    // // assert(low < arr.size());
    // return low;
}