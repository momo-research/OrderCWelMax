#pragma once
#ifndef __HEAD_H__
#define __HEAD_H__

#define DEBUG_FLAG 0
#define MC_RUNS 10000
#define INACTIVE 0
#define INFORMED 1
#define SUSPENDED 2
#define ADOPTED 3
#define REJECTED 4
#define POTENTIAL 5
#define LIVE 1
#define BLOCKED -1

#if defined(WIN32)
#include <time.h>
#include <windows.h>
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

#elif defined(__CYGWIN__) // cygwin
#include <sys/time.h>

#else // linux
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif // omp win32

#include <iostream>
#include <set>
#include <list>
#include <sstream>
#include <cmath>
#include <queue>
#include <fstream>
#include <string>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <map>
#include <deque>
#include <ctime>
#include <omp.h>
#include <mutex>
using namespace std;
typedef unsigned int uint;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int, int> ipair;
typedef pair<double, double> dpair;
#define MP make_pair
#define F first
#define S second

#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (2393.910e6)
#endif

typedef char int8;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;

#ifndef WIN32
#ifdef __CYGWIN__
// CYGWIN
inline uint64 rdtsc() {
    uint64 t0;
    asm volatile("rdtsc" : "=A"(t0));
    return t0;
}

#else
// LINUX
inline uint64 rdtsc() {
    unsigned a, d;
    asm volatile("rdtsc" : "=a" (a), "=d" (d));
    return (((uint64)a) | (((uint64)d) << 32));
}
#endif
#endif

inline string double_to_string(double d) {
    std::stringstream s;
    s << d;
    return s.str();
}

#endif // __HEAD_H__