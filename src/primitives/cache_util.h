#pragma once

// Some utilities to cache align classes and structs.

#include <cstdint>

static const std::size_t S_CACHE_ALIGNMENT = 128;

namespace {
template <class Object>
class DummyCacheAligned : public Object {
 private:
  char padding[S_CACHE_ALIGNMENT - (sizeof(Object) % S_CACHE_ALIGNMENT)];

 public:
  template <typename... _Args>
  DummyCacheAligned(_Args &&...__args) : Object(__args...) {}
};
}  // namespace

template <class Object>
class alignas(alignof(DummyCacheAligned<Object>) > S_CACHE_ALIGNMENT
                  ? alignof(DummyCacheAligned<Object>)
                  : S_CACHE_ALIGNMENT) CacheAligned : public Object {
 private:
  char padding[S_CACHE_ALIGNMENT - (sizeof(Object) % S_CACHE_ALIGNMENT)];

 public:
  template <typename... _Args>
  CacheAligned(_Args &&...__args) : Object(__args...) {}
};