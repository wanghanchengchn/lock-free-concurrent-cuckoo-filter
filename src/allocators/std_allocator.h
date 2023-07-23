#pragma once

#include <cstdlib>

struct StdAllocator {
  static void *malloc(size_t size) { return std::malloc(size); }
  static void free(void *ptr) { std::free(ptr); }
};