#pragma once

// Enums for reclaimers used.

#include <atomic>
#include <memory>

enum class Reclaimer { Leaky, Epoch };
inline const std::string get_reclaimer_name(const Reclaimer reclaimer) {
  switch (reclaimer) {
    case Reclaimer::Leaky:
      return "Leaky";
    case Reclaimer::Epoch:
      return "Epoch";
    default:
      return "Unknown";
  }
}

template <class Allocator>
class ReclaimerAllocator {
 public:
  void *malloc(size_t size) { return Allocator::malloc(size); }
  void free(void *ptr) { Allocator::free(ptr); }
};

template <class MemReclaimer>
class ReclaimerPin {
  typedef typename MemReclaimer::RecordHandle RecordHandle;

  MemReclaimer *m_reclaimer;
  std::size_t m_thread_id;

 public:
  ReclaimerPin(MemReclaimer *reclaimer, std::size_t thread_id)
      : m_reclaimer(reclaimer), m_thread_id(thread_id) {
    m_reclaimer->enter(m_thread_id);
  }
  ~ReclaimerPin() { m_reclaimer->exit(m_thread_id); }

  RecordHandle get_rec() { return m_reclaimer->get_rec(m_thread_id); }

  void retire(const RecordHandle &handle) {
    m_reclaimer->retire(handle, m_thread_id);
  }
};
