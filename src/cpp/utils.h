#ifndef MPH_UTILS_INCLUDED
#define MPH_UTILS_INCLUDED

#include <algorithm>
#include <thread>
#include <functional>
#include <vector>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include "boost/functional/hash.hpp"
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>

/* // Commented for Linuxed -- MB
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#include <mach/mach.h>
*/

// Hash map implemented with std unorderd map
template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};
template <class Key> class hash_set : public std::unordered_set<Key> {};

// value_t, index_t 64b integers
typedef int64_t value_t;
typedef int64_t index_t;
typedef double_t input_t;

long getMemoryUsage();

void mem_usage(double& vm_usage, double& resident_set);

void get_mem_usage(double& vm_usage, double& resident_set);

#endif // MPH_UTILS_INCLUDED
