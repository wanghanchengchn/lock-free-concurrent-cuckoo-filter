#pragma once

/*
Implementation of fast K-CAS by
"Reuse, don't Recycle: Transforming Lock-Free Algorithms that Throw Away
Descriptors."
*/

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <type_traits>

#include "primitives/cache_util.h"
#include "reclaimer/reclaimer.h"

template <class Allocator, class MemReclaimer>
class BrownKCAS {
 private:
  // 2 bits
  static const std::size_t S_NO_TAG = 0x0;
  static const std::size_t S_KCAS_TAG = 0x1;
  static const std::size_t S_RDCSS_TAG = 0x2;

  // 8 Bits
  static const std::size_t S_THREAD_ID_SHIFT = 2;
  static const std::size_t S_THREAD_ID_MASK = (1 << 8) - 1;

  // 54 bits
  static const std::size_t S_SEQUENCE_SHIFT = 10;
  static const std::size_t S_SEQUENCE_MASK = (std::size_t(1) << 54) - 1;

  enum class TagType {
    NONE,
    RDCSS,
    KCAS,
  };

  struct KCASDescriptorStatus {
    static const std::uintptr_t UNDECIDED = 0, SUCCEEDED = 1, FAILED = 2;
    std::uintptr_t status : 8;
    std::uintptr_t sequence_number : 54;

    explicit KCASDescriptorStatus() : status(UNDECIDED), sequence_number(0) {}
    explicit KCASDescriptorStatus(const std::uintptr_t status,
                                  const std::uintptr_t sequence_number)
        : status(status), sequence_number(sequence_number) {}
  };

  struct TaggedPointer {
    std::uintptr_t raw_bits;

    static bool is_rdcss(const TaggedPointer tagptr) {
      return (tagptr.raw_bits & S_RDCSS_TAG) == S_RDCSS_TAG;
    }

    static bool is_kcas(const TaggedPointer tagptr) {
      return (tagptr.raw_bits & S_KCAS_TAG) == S_KCAS_TAG;
    }

    static bool is_bits(const TaggedPointer tagptr) {
      return !is_kcas(tagptr) and !is_rdcss(tagptr);
    }

    explicit TaggedPointer() : raw_bits(0) {}

    explicit TaggedPointer(const std::uintptr_t raw_bits) : raw_bits(raw_bits) {
      assert(is_bits(*this));
    }

    explicit TaggedPointer(const std::uintptr_t tag_bits,
                           const std::uintptr_t thread_id,
                           const std::uintptr_t sequence_number)
        : raw_bits(tag_bits | (thread_id << S_THREAD_ID_SHIFT) |
                   (sequence_number << S_SEQUENCE_SHIFT)) {}

    static std::uintptr_t get_thread_id(const TaggedPointer tagptr) {
      assert(!is_bits(tagptr));
      return (tagptr.raw_bits >> S_THREAD_ID_SHIFT) & S_THREAD_ID_MASK;
    }

    static std::uintptr_t get_sequence_number(const TaggedPointer tagptr) {
      assert(!is_bits(tagptr));
      return (tagptr.raw_bits >> S_SEQUENCE_SHIFT) & S_SEQUENCE_MASK;
    }

    static TaggedPointer make_tagged(const std::uintptr_t tag_bits,
                                     const TaggedPointer tagptr) {
      return TaggedPointer{tag_bits, get_thread_id(tagptr),
                           get_sequence_number(tagptr)};
    }

    static TaggedPointer mask_bits(const TaggedPointer tagptr) {
      return make_tagged(S_NO_TAG, tagptr);
    }

    static TaggedPointer make_rdcss(const std::uintptr_t thread_id,
                                    const std::uintptr_t sequence_number) {
      return TaggedPointer{S_RDCSS_TAG, thread_id, sequence_number};
    }

    static TaggedPointer make_rdcss(const TaggedPointer tagptr) {
      return make_tagged(S_RDCSS_TAG, tagptr);
    }

    static TaggedPointer make_kcas(const std::uintptr_t thread_id,
                                   const std::uintptr_t sequence_number) {
      return TaggedPointer{S_KCAS_TAG, thread_id, sequence_number};
    }

    static TaggedPointer make_kcas(const TaggedPointer tagptr) {
      return make_tagged(S_KCAS_TAG, tagptr);
    }

    static TaggedPointer make_bits(const std::uintptr_t raw_bits) {
      return make_bits(TaggedPointer{raw_bits});
    }

    static TaggedPointer make_bits(const TaggedPointer tagptr) {
      return TaggedPointer{tagptr.raw_bits};
    }
  };

  struct RDCSSDescriptor {
    std::atomic_size_t sequence_bits;
    std::atomic<TaggedPointer> *data_location;
    TaggedPointer before, kcas_tagptr;
    std::atomic<KCASDescriptorStatus> *status_location;

    std::size_t increment_sequence() { return sequence_bits.fetch_add(1); }
  };

  bool try_snapshot(RDCSSDescriptor *snapshot, TaggedPointer ptr,
                    [[maybe_unused]] const std::size_t my_thread_id) {
    assert(TaggedPointer::is_rdcss(ptr));
    const std::uintptr_t thread_id = TaggedPointer::get_thread_id(ptr);
    const std::uintptr_t sequence_number =
        TaggedPointer::get_sequence_number(ptr);
    RDCSSDescriptor *snapshot_target = &m_rdcss_descs[thread_id];
    const std::size_t before_sequence = snapshot_target->sequence_bits.load();
    if (before_sequence != sequence_number) {
      assert(my_thread_id != thread_id);
      return false;
    }
    // Snapshot
    snapshot->sequence_bits.store(before_sequence, std::memory_order_relaxed);
    snapshot->data_location = snapshot_target->data_location;
    snapshot->before = snapshot_target->before;
    snapshot->kcas_tagptr = snapshot_target->kcas_tagptr;
    snapshot->status_location = snapshot_target->status_location;
    // Check our sequence number again.
    const std::size_t after_sequence = snapshot_target->sequence_bits.load();
    if (after_sequence != sequence_number) {
      assert(my_thread_id != thread_id);
      return false;
    }
    assert(my_thread_id != thread_id);
    return true;
  }

  TaggedPointer rdcss(TaggedPointer ptr, const std::size_t my_thread_id) {
    assert(TaggedPointer::is_rdcss(ptr));
    const std::size_t thread_id = TaggedPointer::get_thread_id(ptr);
    assert(my_thread_id == thread_id);
    RDCSSDescriptor *rdcss_desc = &m_rdcss_descs[thread_id];
    while (true) {
      TaggedPointer current =
          rdcss_desc->data_location->load(std::memory_order_relaxed);
      if (TaggedPointer::is_rdcss(current)) {
        RDCSSDescriptor snapshot;
        if (try_snapshot(&snapshot, current, my_thread_id)) {
          rdcss_complete(&snapshot, current);
        }
        continue;
      }
      if (current.raw_bits != rdcss_desc->before.raw_bits) {
        assert(!TaggedPointer::is_rdcss(current));
        return current;
      }
      bool success = rdcss_desc->data_location->compare_exchange_strong(
          current, ptr, std::memory_order_relaxed, std::memory_order_relaxed);
      if (success) {
        rdcss_complete(rdcss_desc, ptr);
        return rdcss_desc->before;
      }
    }
  }
  TaggedPointer rdcss_read(
      const std::atomic<TaggedPointer> *location,
      const std::size_t my_thread_id,
      const std::memory_order memory_order = std::memory_order_seq_cst) {
    while (true) {
      TaggedPointer current = location->load(memory_order);
      if (TaggedPointer::is_rdcss(current)) {
        RDCSSDescriptor snapshot;
        if (try_snapshot(&snapshot, current, my_thread_id)) {
          rdcss_complete(&snapshot, current);
        }
      } else {
        return current;
      }
    }
  }

  void rdcss_complete(RDCSSDescriptor *snapshot, TaggedPointer ptr) {
    assert(TaggedPointer::is_rdcss(ptr));
    // Old 1 and address 1.
    const std::uintptr_t sequence_number =
        TaggedPointer::get_sequence_number(snapshot->kcas_tagptr);
    const KCASDescriptorStatus status = snapshot->status_location->load();
    if (status.sequence_number == sequence_number and
        status.status == KCASDescriptorStatus::UNDECIDED) {
      snapshot->data_location->compare_exchange_strong(
          ptr, snapshot->kcas_tagptr, std::memory_order_relaxed,
          std::memory_order_relaxed);
    } else {
      snapshot->data_location->compare_exchange_strong(
          ptr, snapshot->before, std::memory_order_relaxed,
          std::memory_order_relaxed);
    }
  }

  static const std::size_t KCASShift = 2;
  const std::size_t m_num_threads, m_descriptor_size;
  struct DescriptorEntry {
    std::atomic<TaggedPointer> *location;
    TaggedPointer before, desired;
  };

 public:
  class KCASDescriptor;

 private:
  static const std::size_t DescriptorSize = 10;

 public:
  // Class to wrap around the types being KCAS'd
  template <class Type>
  class KCASEntry {
   public:
    KCASEntry() : m_entry(TaggedPointer{}) {}

   private:
    union ItemBitsUnion {
      Type item;
      std::uintptr_t raw_bits;
    };

    std::atomic<TaggedPointer> m_entry;
    static std::size_t from_raw_bits(uintptr_t raw_bits) {
      return (ItemBitsUnion{raw_bits}).item;
    }
    static std::uintptr_t to_raw_bits(const Type &inner) {
      return ItemBitsUnion{inner}.raw_bits;
    }
    friend class BrownKCAS;
    friend class KCASDescriptor;
  };

  class KCASDescriptor {
   private:
    std::size_t m_descriptor_size;
    std::atomic<KCASDescriptorStatus> m_status;
    std::size_t m_num_entries;
    DescriptorEntry m_entries[DescriptorSize];

    KCASDescriptor()
        : m_descriptor_size(DescriptorSize),
          m_status(KCASDescriptorStatus{}),
          m_num_entries(0) {}
    KCASDescriptor(const KCASDescriptor &rhs) = delete;
    KCASDescriptor &operator=(const KCASDescriptor &rhs) = delete;
    KCASDescriptor &operator=(KCASDescriptor &&rhs) = delete;

    void increment_sequence() {
      KCASDescriptorStatus current_status =
          m_status.load(std::memory_order_relaxed);
      m_status.store(KCASDescriptorStatus{KCASDescriptorStatus::UNDECIDED,
                                          current_status.sequence_number + 1});
    }

    bool inactive_sequence() {
      KCASDescriptorStatus current_status =
          m_status.load(std::memory_order_relaxed);
      m_status.store(KCASDescriptorStatus{KCASDescriptorStatus::UNDECIDED,
                                          current_status.sequence_number + 1});
      return true;
    }

    bool active_sequence() {
      KCASDescriptorStatus current_status =
          m_status.load(std::memory_order_relaxed);
      m_status.store(KCASDescriptorStatus{KCASDescriptorStatus::UNDECIDED,
                                          current_status.sequence_number + 1});
      return true;
    }

   public:
    template <class ValType>
    void add_value(KCASEntry<ValType> *location, const ValType &before,
                   const ValType &desired) {
      std::uintptr_t before_raw_bits = KCASEntry<ValType>::to_raw_bits(before);
      TaggedPointer before_desc{before_raw_bits << KCASShift};
      assert(TaggedPointer::is_bits(before_desc));
      std::uintptr_t desired_raw_bits =
          KCASEntry<ValType>::to_raw_bits(desired);
      TaggedPointer desired_desc{desired_raw_bits << KCASShift};
      assert(TaggedPointer::is_bits(desired_desc));
      const std::size_t cur_entry = m_num_entries++;
      assert(cur_entry < m_descriptor_size);
      m_entries[cur_entry].before = before_desc;
      m_entries[cur_entry].desired = desired_desc;
      m_entries[cur_entry].location = &location->m_entry;
    }
    template <class PtrType>
    void add_ptr(const KCASEntry<PtrType> *location, const PtrType &before,
                 const PtrType &desired) {
      std::uintptr_t before_raw_bits = KCASEntry<PtrType>::to_raw_bits(before);
      TaggedPointer before_desc{before_raw_bits};
      assert(TaggedPointer::is_bits(before_desc));
      std::uintptr_t desired_raw_bits =
          KCASEntry<PtrType>::to_raw_bits(desired);
      TaggedPointer desired_desc{desired_raw_bits};
      assert(TaggedPointer::is_bits(desired_desc));
      const std::size_t cur_entry = m_num_entries++;
      assert(cur_entry < m_descriptor_size);
      m_entries[cur_entry].before = before_desc;
      m_entries[cur_entry].desired = desired_desc;
      m_entries[cur_entry].location = &location->m_entry;
    }
    friend class BrownKCAS;
    template <class T>
    friend class CacheAligned;
  };

 private:
  static const std::size_t S_NUM_THREADS = 144;
  CacheAligned<KCASDescriptor> m_kcas_descs[S_NUM_THREADS];
  CacheAligned<RDCSSDescriptor> m_rdcss_descs[S_NUM_THREADS];

  bool try_snapshot(KCASDescriptor *snapshot, TaggedPointer ptr,
                    [[maybe_unused]] const std::size_t my_thread_id) {
    assert(TaggedPointer::is_kcas(ptr));
    const std::uintptr_t thread_id = TaggedPointer::get_thread_id(ptr);
    const std::uintptr_t sequence_number =
        TaggedPointer::get_sequence_number(ptr);
    KCASDescriptor *snapshot_target = &m_kcas_descs[thread_id];
    const KCASDescriptorStatus before_status =
        snapshot_target->m_status.load(std::memory_order_acquire);
    const std::size_t num_entries = snapshot_target->m_num_entries;
    if (before_status.sequence_number != sequence_number) {
      //      assert(my_thread_id != thread_id);
      return false;
    }
    for (std::size_t i = 0; i < num_entries; i++) {
      snapshot->m_entries[i] = snapshot_target->m_entries[i];
    }
    const KCASDescriptorStatus after_status = snapshot_target->m_status.load();
    if (after_status.sequence_number != sequence_number) {
      //      assert(my_thread_id != thread_id);
      return false;
    }
    snapshot->m_num_entries = num_entries;
    //    assert(my_thread_id != thread_id);
    return true;
  }

  bool cas_internal(const std::size_t my_thread_id, const TaggedPointer tagptr,
                    KCASDescriptor *descriptor_snapshot, bool help = false) {
    assert(TaggedPointer::is_kcas(tagptr));
    const std::uintptr_t descriptor_thread_id =
        TaggedPointer::get_thread_id(tagptr);
    const std::size_t num_entries = descriptor_snapshot->m_num_entries;

    // This is the descriptor we're trying to complete.
    // We also have a snapshot of it as an argument.
    KCASDescriptor *original_desc = &m_kcas_descs[descriptor_thread_id];

    // We have two windows into the state:
    // The first is the tagptr which has a state embeded within it.
    // The second is the sequence bits inside the descriptor.
    // They need to match and if they do then the state extracted.
    const std::size_t tagptr_sequence_number =
        TaggedPointer::get_sequence_number(tagptr);
    KCASDescriptorStatus current_status = original_desc->m_status.load();
    if (tagptr_sequence_number != current_status.sequence_number) {
      // This can only fail if we are helping another descriptor
      // since the owner is the only one who can update the sequence bits.
      assert(help);
      return false;
    }
    if (current_status.status == KCASDescriptorStatus::UNDECIDED) {
      std::uintptr_t status = KCASDescriptorStatus::SUCCEEDED;
      for (std::size_t i = help ? 1 : 0;
           i < num_entries and status == KCASDescriptorStatus::SUCCEEDED; i++) {
      retry:
        RDCSSDescriptor *rdcss_desc = &m_rdcss_descs[my_thread_id];
        // Give us a new descriptor.
        rdcss_desc->increment_sequence();
        rdcss_desc->before = descriptor_snapshot->m_entries[i].before;
        rdcss_desc->data_location = descriptor_snapshot->m_entries[i].location;
        rdcss_desc->status_location = &original_desc->m_status;
        rdcss_desc->kcas_tagptr = tagptr;
        // Make it initialised.
        std::size_t new_sequence = rdcss_desc->increment_sequence();

        TaggedPointer rdcss_ptr =
            TaggedPointer::make_rdcss(my_thread_id, new_sequence);
        TaggedPointer value = this->rdcss(rdcss_ptr, my_thread_id);
        if (TaggedPointer::is_kcas(value)) {
          if (value.raw_bits != tagptr.raw_bits) {
            KCASDescriptor descriptor_snapshot;
            if (try_snapshot(&descriptor_snapshot, value, my_thread_id)) {
              cas_internal(my_thread_id, value, &descriptor_snapshot, true);
            }
            goto retry;
          }
        } else if (value.raw_bits != rdcss_desc->before.raw_bits) {
          status = KCASDescriptorStatus::FAILED;
        }
      }
      // Try change descriptor status.
      KCASDescriptorStatus expected_status = current_status;
      KCASDescriptorStatus original_status =
          original_desc->m_status.load(std::memory_order_relaxed);
      if (original_status.sequence_number == expected_status.sequence_number and
          original_status.status == KCASDescriptorStatus::UNDECIDED) {
        assert(status == KCASDescriptorStatus::SUCCEEDED or
               status == KCASDescriptorStatus::FAILED);
        original_desc->m_status.compare_exchange_strong(
            expected_status,
            KCASDescriptorStatus{status, current_status.sequence_number},
            std::memory_order_relaxed, std::memory_order_relaxed);
      }
    }

    // Phase 2.
    KCASDescriptorStatus new_status = original_desc->m_status.load();
    if (new_status.sequence_number != tagptr_sequence_number) {
      assert(help);
      return false;
    }
    bool succeeded = (new_status.status == KCASDescriptorStatus::SUCCEEDED);
    for (std::size_t i = 0; i < num_entries; i++) {
      TaggedPointer new_value = succeeded
                                    ? descriptor_snapshot->m_entries[i].desired
                                    : descriptor_snapshot->m_entries[i].before;
      TaggedPointer expected = tagptr;
      descriptor_snapshot->m_entries[i].location->compare_exchange_strong(
          expected, new_value);
    }
    return succeeded;
  }

 public:
  BrownKCAS(const std::size_t threads, [[maybe_unused]] MemReclaimer *reclaimer)
      : m_num_threads(threads),
        m_descriptor_size(DescriptorSize) /*,
m_kcas_descs(static_cast<CacheAligned<KCASDescriptor> *>(
Allocator::malloc(sizeof(CacheAligned<KCASDescriptor>) * threads))),
m_rdcss_descs(
static_cast<CacheAligned<RDCSSDescriptor> *>(Allocator::malloc(
sizeof(CacheAligned<RDCSSDescriptor>) * threads)))*/
  {
    for (std::size_t i = 0; i < m_num_threads; i++) {
      m_kcas_descs[i].m_status.store(KCASDescriptorStatus{});
      m_kcas_descs[i].m_num_entries = 0;
      m_kcas_descs[i].m_descriptor_size = m_descriptor_size;
    }
  }

  KCASDescriptor *create_descriptor(
      [[maybe_unused]] const std::size_t descriptor_size,
      const std::size_t thread_id) {
    assert(descriptor_size < DescriptorSize);
    // Increment the status sequence number.
    m_kcas_descs[thread_id].increment_sequence();
    m_kcas_descs[thread_id].m_num_entries = 0;
    return &m_kcas_descs[thread_id];
  }

  void free_descriptor([[maybe_unused]] KCASDescriptor *desc) {}

  bool cas(const std::size_t thread_id,
           [[maybe_unused]] ReclaimerPin<MemReclaimer> &pin,
           KCASDescriptor *desc) {
    std::sort(
        desc->m_entries, desc->m_entries + desc->m_num_entries,
        [](const DescriptorEntry &lhs, const DescriptorEntry &rhs) -> bool {
          return lhs.location < rhs.location;
        });
    // Init descriptor...
    desc->increment_sequence();
    TaggedPointer ptr = TaggedPointer::make_kcas(
        thread_id, desc->m_status.load().sequence_number);
    return cas_internal(thread_id, ptr, desc);
  }

  template <class ValType>
  ValType read_value(
      const std::size_t my_thread_id,
      [[maybe_unused]] ReclaimerPin<MemReclaimer> &pin,
      const KCASEntry<ValType> *location,
      const std::memory_order memory_order = std::memory_order_seq_cst) {
    static_assert(!std::is_pointer<ValType>::value, "Type is not a value.");
    while (true) {
      TaggedPointer desc =
          this->rdcss_read(&location->m_entry, my_thread_id, memory_order);
      // Could still be a K-CAS descriptor.
      if (TaggedPointer::is_kcas(desc)) {
        KCASDescriptor descriptor_snapshot;
        if (try_snapshot(&descriptor_snapshot, desc, my_thread_id)) {
          cas_internal(my_thread_id, desc, &descriptor_snapshot, true);
        }
        continue;
      }
      assert(TaggedPointer::is_bits(desc));
      return KCASEntry<ValType>::from_raw_bits(desc.raw_bits) >> KCASShift;
    }
  }

  template <class PtrType>
  PtrType read_ptr(
      const std::size_t thread_id,
      [[maybe_unused]] ReclaimerPin<MemReclaimer> &pin,
      const KCASEntry<PtrType> *location,
      const std::memory_order memory_order = std::memory_order_seq_cst) {
    static_assert(std::is_pointer<PtrType>::value, "Type is not a pointer.");
    while (true) {
      TaggedPointer desc = this->rdcss_read(location, thread_id, memory_order);
      // Could still be a K-CAS descriptor.
      if (TaggedPointer::is_kcas(desc)) {
        KCASDescriptor descriptor_snapshot;
        if (try_snapshot(&descriptor_snapshot, desc, thread_id)) {
          cas_internal(thread_id, desc, &descriptor_snapshot, true);
        }
        continue;
      }
      assert(TaggedPointer::is_bits(desc));
      return KCASEntry<PtrType>::from_raw_bits(desc);
    }
  }

  template <class ValType>
  void write_value(
      KCASEntry<ValType> *location, const ValType &val,
      const std::memory_order memory_order = std::memory_order_seq_cst) {
    static_assert(!std::is_pointer<ValType>::value, "Type is not a value.");
    uintptr_t raw_bits = KCASEntry<ValType>::to_raw_bits(val);
    TaggedPointer desc = TaggedPointer::make_bits(raw_bits << KCASShift);
    assert(TaggedPointer::is_bits(desc));
    location->m_entry.store(desc, memory_order);
  }

  template <class PtrType>
  void write_ptr(
      KCASEntry<PtrType> *location, const PtrType &ptr,
      const std::memory_order memory_order = std::memory_order_seq_cst) {
    static_assert(std::is_pointer<PtrType>::value, "Type is not a pointer.");
    uintptr_t raw_bits = KCASEntry<PtrType>::to_raw_bits(ptr);
    TaggedPointer desc = TaggedPointer::make_bits(raw_bits);
    assert(TaggedPointer::is_bits(desc));
    location->m_entry.store(desc, memory_order);
  }

  template <class ValType>
  bool compare_exchange_strong_value(
      const std::size_t thread_id, KCASEntry<ValType> *location,
      ValType &expected, const ValType &desired,
      std::memory_order success = std::memory_order_seq_cst,
      std::memory_order fail = std::memory_order_seq_cst) {
    std::uintptr_t before_raw_bits = KCASEntry<ValType>::to_raw_bits(expected);
    TaggedPointer before_desc{before_raw_bits << KCASShift};
    std::uintptr_t desired_raw_bits = KCASEntry<ValType>::to_raw_bits(desired);
    TaggedPointer desired_desc{desired_raw_bits << KCASShift};
    bool ret = location->m_entry.compare_exchange_strong(
        before_desc, desired_desc, success, fail);
    if (!ret) {
      while (true) {
        if (TaggedPointer::is_bits(before_desc)) {
          expected = KCASEntry<ValType>::from_raw_bits(before_desc.raw_bits >>
                                                       KCASShift);
          return false;
        }
        if (TaggedPointer::is_rdcss(before_desc)) {
          RDCSSDescriptor snapshot;
          if (try_snapshot(&snapshot, before_desc, thread_id)) {
            rdcss_complete(&snapshot, before_desc);
          }
          before_desc = location->m_entry.load(fail);
          continue;
        }
        // Could still be a K-CAS descriptor.
        if (TaggedPointer::is_kcas(before_desc)) {
          KCASDescriptor descriptor_snapshot;
          if (try_snapshot(&descriptor_snapshot, before_desc, thread_id)) {
            cas_internal(thread_id, before_desc, &descriptor_snapshot, true);
          }
          before_desc = location->m_entry.load(fail);
          continue;
        }
        expected = KCASEntry<ValType>::from_raw_bits(before_desc.raw_bits >>
                                                     KCASShift);
        return false;
      }
    }
    return true;
  }
};
