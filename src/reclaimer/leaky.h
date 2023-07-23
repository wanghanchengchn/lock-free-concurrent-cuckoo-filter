#pragma once

// "Leaky" reclaimer. Does nothing but fit interface for reclaimers.

#include <atomic>
#include <cassert>

#include "reclaimer/reclaimer.h"

template <class Allocator>
class LeakyReclaimer : public ReclaimerAllocator<Allocator> {
 public:
  class LeakyBase {};

  class LeakyHandle {
   public:
    LeakyHandle() {}
    ~LeakyHandle() {}
    LeakyHandle(const LeakyHandle &rhs) = delete;
    LeakyHandle &operator=(const LeakyHandle &rhs) = delete;
    LeakyHandle(LeakyHandle &&) {}
    LeakyHandle &operator=(LeakyHandle &&) { return *this; }

    void set([[maybe_unused]] const LeakyBase *ptr) {}

    template <class PtrType, class SourceType, typename Func>
    bool try_protect([[maybe_unused]] PtrType &ptr,
                     [[maybe_unused]] const std::atomic<SourceType> &src,
                     [[maybe_unused]] Func f) {
      return true;
    }
    template <class PtrType>
    bool try_protect(PtrType &ptr, const std::atomic<PtrType> &src) noexcept {
      return this->try_protect(ptr, src, [](PtrType ptr) { return ptr; });
    }

    template <class PtrType, class SourceType, typename Func>
    PtrType get_protected([[maybe_unused]] const std::atomic<SourceType> &src,
                          [[maybe_unused]] Func f) {}
    template <class PtrType>
    PtrType get_protected(const std::atomic<PtrType> &src) noexcept {
      return this->get_protected(src, [](PtrType ptr) { return ptr; });
    }
    friend class LeakyReclaimer;
  };

  typedef LeakyBase RecordBase;
  typedef LeakyHandle RecordHandle;

  LeakyReclaimer(__attribute__((unused)) const std::size_t num_threads) {}
  ~LeakyReclaimer() {}
  bool thread_init([[maybe_unused]] const std::size_t thread_id) {
    return true;
  }
  void enter([[maybe_unused]] const std::size_t thread_id) {}
  void exit([[maybe_unused]] const std::size_t thread_id) {}
  LeakyHandle get_rec([[maybe_unused]] const std::size_t thread_id) {
    return LeakyHandle();
  }
  void retire([[maybe_unused]] const RecordHandle &handle,
              [[maybe_unused]] const std::size_t thread_id) {}
};
