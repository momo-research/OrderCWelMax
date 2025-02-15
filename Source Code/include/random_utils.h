#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <random>

extern thread_local std::mt19937 rng;
double randomDouble();
int randomInt(int n);
void resetRNG();
#endif