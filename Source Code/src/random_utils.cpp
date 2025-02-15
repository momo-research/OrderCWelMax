#include "random_utils.h"

thread_local std::mt19937 rng(std::random_device{}()); 

double randomDouble() {
    return static_cast<double>(rng()) / rng.max();
}
int randomInt(int n) {
    std::uniform_int_distribution<int> dist(0, n-1); 
    return dist(rng); 
}
void resetRNG() {
    rng.seed(std::random_device{}()); 
}