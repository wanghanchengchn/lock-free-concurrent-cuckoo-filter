#pragma once

#include <cassert>
#include <cstdint>

#include "hashfunction/two_independent.h"
#include "mathutil.h"
#include "reclaimer/reclaimer.h"

template <typename ItemType, typename TagType, typename Reclaimer,
          typename KCAS>
class LockfreeMultislotKCASBfsCuckoo {};

template <typename ItemType, typename Reclaimer, typename KCAS>
class LockfreeMultislotKCASBfsCuckoo<ItemType, uint8_t, Reclaimer, KCAS> {
 public:
  LockfreeMultislotKCASBfsCuckoo(const size_t max_num_keys,
                                 const size_t num_threads)
      : num_buckets_(CalcNumBuckets(max_num_keys)),
        tables_(new
                typename KCAS::template KCASEntry<TableEntry>[num_buckets_]),
        hasher_(),
        reclaimer_(num_threads),
        kcas_(num_threads, &reclaimer_) {}

  bool Lookup(const ItemType &key, const size_t thread_id) {
    int64_t i1, i2;
    uint8_t tag;
    GenerateIndexTagHash(key, &i1, &tag);
    i2 = AltIndex(i1, tag);
    ReclaimerPin<Reclaimer> pin(&reclaimer_, thread_id);

    while (true) {
      TableEntry e1 = kcas_.read_value(thread_id, pin, &tables_[i1]);
      TableEntry e2 = kcas_.read_value(thread_id, pin, &tables_[i2]);
      if (FindTagInBuckets(e1, e2, tag)) {
        return true;
      }

      // second round query
      TableEntry e1x = kcas_.read_value(thread_id, pin, &tables_[i1]);
      TableEntry e2x = kcas_.read_value(thread_id, pin, &tables_[i2]);
      if (FindTagInBuckets(e1x, e2x, tag)) {
        return true;
      }

      if (CheckCounter(GetCounter(e1), GetCounter(e2), GetCounter(e1x),
                       GetCounter(e2x))) {
        continue;
      } else {
        return false;
      }
    }
  }

  bool Delete(const ItemType &key, const size_t thread_id) {
    int64_t i1, i2;
    uint8_t tag;
    GenerateIndexTagHash(key, &i1, &tag);
    i2 = AltIndex(i1, tag);
    ReclaimerPin<Reclaimer> pin(&reclaimer_, thread_id);
    TableEntry dest;

    while (true) {
      TableEntry e1 = kcas_.read_value(thread_id, pin, &tables_[i1]);
      TableEntry e2 = kcas_.read_value(thread_id, pin, &tables_[i2]);
      if (FindTagAndRemove(e1, tag, &dest)) {
        assert(FindEmptySlot(dest) >= 0);
        if (!kcas_.compare_exchange_strong_value(thread_id, &tables_[i1], e1,
                                                 dest)) {
          continue;
        }
        return true;
      }
      if (FindTagAndRemove(e2, tag, &dest)) {
        assert(FindEmptySlot(dest) >= 0);
        if (!kcas_.compare_exchange_strong_value(thread_id, &tables_[i2], e2,
                                                 dest)) {
          continue;
        }
        return true;
      }

      // second round query
      TableEntry e1x = kcas_.read_value(thread_id, pin, &tables_[i1]);
      TableEntry e2x = kcas_.read_value(thread_id, pin, &tables_[i2]);
      if (FindTagAndRemove(e1x, tag, &dest)) {
        assert(FindEmptySlot(dest) >= 0);
        if (!kcas_.compare_exchange_strong_value(thread_id, &tables_[i1], e1x,
                                                 dest)) {
          continue;
        }
        return true;
      }
      if (FindTagAndRemove(e2x, tag, &dest)) {
        assert(FindEmptySlot(dest) >= 0);
        if (!kcas_.compare_exchange_strong_value(thread_id, &tables_[i2], e2x,
                                                 dest)) {
          continue;
        }
        return true;
      }

      if (CheckCounter(GetCounter(e1), GetCounter(e2), GetCounter(e1x),
                       GetCounter(e2x))) {
        continue;
      } else {
        return false;
      }
    }
  }

  bool Insert(const ItemType &key, const size_t thread_id) {
    int64_t i1, i2;
    uint8_t tag;
    GenerateIndexTagHash(key, &i1, &tag);
    i2 = AltIndex(i1, tag);
    ReclaimerPin<Reclaimer> pin(&reclaimer_, thread_id);

    TableEntry e1, e2, dest;
  start_insert:
    e1 = kcas_.read_value(thread_id, pin, &tables_[i1]);
    e2 = kcas_.read_value(thread_id, pin, &tables_[i2]);

    if (FindEmptySlotAndInsertTag(e1, tag, &dest)) {
      assert(FindTagInBuckets(dest, dest, tag));
      if (!kcas_.compare_exchange_strong_value(thread_id, &tables_[i1], e1,
                                               dest)) {
        goto start_insert;
      }
      return true;
    }

    if (FindEmptySlotAndInsertTag(e2, tag, &dest)) {
      assert(FindTagInBuckets(dest, dest, tag));
      if (!kcas_.compare_exchange_strong_value(thread_id, &tables_[i2], e2,
                                               dest)) {
        goto start_insert;
      }
      return true;
    }

    if (Relocate(thread_id, i1, i2, &pin)) {
      goto start_insert;
    }

    return false;
  }

 private:
  // Hashing types and functions
  static inline size_t CalcNumBuckets(const size_t max_num_keys) {
    size_t num_buckets =
        upperpower2(std::max<uint64_t>(1, max_num_keys / kTagsPerBucket));
    double frac = (double)max_num_keys / num_buckets / kTagsPerBucket;
    if (frac > 0.97) {
      num_buckets <<= 1;
    }
    return num_buckets;
  }

  inline int64_t IndexHash(uint32_t hv) const {
    // table_.num_buckets is always a power of two, so modulo can be replaced
    // with bitwise-and:
    return hv & (num_buckets_ - 1);
  }

  inline uint8_t TagHash(uint32_t hv) const {
    uint32_t tag;
    tag = hv & ((1ULL << (sizeof(uint8_t) * 8)) - 1);
    tag += (tag == 0);
    return tag;
  }

  inline int64_t AltIndex(const int64_t index, const uint8_t tag) const {
    // NOTE(binfan): originally we use:
    // index ^ HashUtil::BobHash((const void*) (&tag), 4)) & table_.INDEXMASK;
    // now doing a quick-n-dirty way:
    // 0x5bd1e995 is the hash constant from MurmurHash2
    return IndexHash((uint32_t)(index ^ (tag * 0x5bd1e995)));
  }

  inline void GenerateIndexTagHash(const ItemType &item, int64_t *index,
                                   uint8_t *tag) {
    const uint64_t hash = hasher_(item);
    *index = IndexHash(hash >> 32);
    *tag = TagHash(hash);
  };

  // Hash table entry types and functions

  // table entry has followint components:
  // | empty (2 bit) | 4 tags (32 bit) | relocation counter (30 bit) |
  // leave top 2 bits for KCAS (it will left shift 2 bit for values)
  using TableEntry = uint64_t;
  constexpr static size_t kTagsPerBucket = 4;
  constexpr static uint32_t kTagsShift = 30;
  constexpr static uint32_t kMaxCounterValue = (1 << kTagsShift) - 1;
  constexpr static uint32_t kMaxTagValue = (1 << 8) - 1;

  static uint8_t GetTag(const TableEntry &entry, int i) {
    assert(i >= 0 && i < 4);
    return (entry >> (kTagsShift + i * 8)) & kMaxTagValue;
  }

  static uint32_t GetCounter(const TableEntry &entry) {
    return entry & kMaxCounterValue;
  }

  static TableEntry AddTag(const TableEntry &entry, int i, uint64_t tag) {
    assert(i >= 0 && i < 4);
    assert(GetTag(entry, i) == 0);
    return entry | ((tag & kMaxTagValue) << (kTagsShift + i * 8));
  }

  static TableEntry ClearTag(const TableEntry &entry, int i) {
    assert(i >= 0 && i < 4);
    constexpr uint64_t masks[4]{
        ~(0xfful << kTagsShift), ~(0xfful << (kTagsShift + 8)),
        ~(0xfful << (kTagsShift + 16)), ~(0xfful << (kTagsShift + 24))};
    return entry & masks[i];
  }

  static TableEntry SetCounter(const TableEntry &entry, uint64_t counter) {
    constexpr uint64_t mask = (static_cast<uint64_t>(0xffffffff) << kTagsShift);
    return (entry & mask) | (counter & kMaxCounterValue);
  }

  static bool FindTagInBuckets(const TableEntry &e1, const TableEntry &e2,
                               const uint8_t tag) {
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if ((GetTag(e1, j) == tag) || (GetTag(e2, j) == tag)) {
        return true;
      }
    }
    return false;
  }

  static int8_t FindEmptySlot(const TableEntry &e) {
    static_assert(kTagsPerBucket <= std::numeric_limits<int8_t>::max(),
                  "int8_t can not represent slot index");
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (GetTag(e, j) == 0) {
        return j;
      }
    }
    return -1;
  }

  static bool FindEmptySlotAndInsertTag(const TableEntry &e, const uint8_t tag,
                                        TableEntry *dest) {
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (GetTag(e, j) == 0) {
        *dest = AddTag(e, j, tag);
        return true;
      }
    }
    return false;
  }

  static bool FindTagAndRemove(const TableEntry &e, const uint8_t tag,
                               TableEntry *dest) {
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (GetTag(e, j) == tag) {
        *dest = ClearTag(e, j);
        return true;
      }
    }
    return false;
  }

  // Searching types and functions

  inline bool CheckCounter(uint32_t ts1, uint32_t ts2, uint32_t ts1x,
                           uint32_t ts2x) {
    return (ts1x >= ts1 + 2) && (ts2x >= ts2 + 2) && (ts2x >= ts1 + 3);
  }

  // Insertion types and function

  constexpr static uint8_t kMaxBfsPathLen = 5;
  struct CuckooRecord {
    int64_t bucket_index;
    uint8_t slot_index;
  };
  using CuckooRecords = std::array<CuckooRecord, kMaxBfsPathLen>;

  // A constexpr version of pow that we can use for various compile-time
  // constants and checks.
  static constexpr size_t const_pow(size_t a, size_t b) {
    return (b == 0) ? 1 : a * const_pow(a, b - 1);
  }

  struct BfsSlot {
    // The bucket of the last item in the path.
    int64_t bucket;
    // a compressed representation of the slots for each of the buckets in
    // the path. pathcode is sort of like a base-slot_per_bucket number, and
    // we need to hold at most MAX_BFS_PATH_LEN slots. Thus we need the
    // maximum pathcode to be at least slot_per_bucket()^(MAX_BFS_PATH_LEN).
    uint16_t pathcode;
    static_assert(const_pow(kTagsPerBucket, kMaxBfsPathLen) <
                      std::numeric_limits<decltype(pathcode)>::max(),
                  "pathcode may not be large enough to encode a cuckoo path");
    // The 0-indexed position in the cuckoo path this slot occupies. It must
    // be less than MAX_BFS_PATH_LEN, and also able to hold negative values.
    int8_t depth;
    static_assert(
        kMaxBfsPathLen - 1 <= std::numeric_limits<decltype(depth)>::max(),
        "The depth type must able to hold a value of kMaxBfsPathLen - 1");
    static_assert(-1 >= std::numeric_limits<decltype(depth)>::min(),
                  "The depth type must be able to hold a value of -1");

    BfsSlot() = default;
    BfsSlot(const size_t bucket, const uint16_t pathcode, const int8_t depth)
        : bucket(bucket), pathcode(pathcode), depth(depth) {}
  };

  class BfsQueue {
   public:
    BfsQueue() noexcept : first_(0), last_(0) {}

    void Enqueue(BfsSlot x) { slots_[last_++] = x; }

    BfsSlot Dequeue() {
      assert(!Empty());
      assert(first_ < last_);
      return slots_[first_++];
    }

    BfsSlot Top() {
      assert(!Empty());
      return slots_[first_];
    }

    bool Empty() const { return first_ == last_; }

    bool Full() const { return last_ == kMaxCuckooCount; }

   private:
    constexpr static size_t kMaxCuckooCount =
        2 * ((kTagsPerBucket == 1)
                 ? kMaxBfsPathLen
                 : (const_pow(kTagsPerBucket, kMaxBfsPathLen) - 1) /
                       (kTagsPerBucket - 1));
    // Queue elements, allocate just enough space to complete a full search.
    BfsSlot slots_[kMaxCuckooCount];
    // The index of the head of the queue in the array
    size_t first_;
    // One past the index of the last_ item of the queue in the array.
    size_t last_;
  };

  bool Relocate(const size_t thread_id, int64_t i1, int64_t i2,
                ReclaimerPin<Reclaimer> *pin) {
  path_discovery:
    int8_t depth{-1};
    uint16_t pathcode{};
    BfsQueue q;
    q.Enqueue(BfsSlot(i1, 0, 0));
    q.Enqueue(BfsSlot(i2, 1, 0));
    while (depth == -1 && !q.Empty()) {
      BfsSlot x{q.Dequeue()};
      if (__builtin_expect(!!(!q.Empty()), 1)) {
        __builtin_prefetch(&tables_[q.Top().bucket]);
      }
      // Picks a (sort-of) random slot to start from
      size_t starting_slot{x.pathcode % kTagsPerBucket};
      TableEntry e{kcas_.read_value(thread_id, *pin, &tables_[x.bucket])};
      for (size_t i = 0; i < kTagsPerBucket; ++i) {
        uint16_t slot{
            static_cast<uint16_t>((starting_slot + i) % kTagsPerBucket)};
        uint8_t tag{GetTag(e, slot)};
        if (tag == 0) {
          pathcode = x.pathcode * kTagsPerBucket + slot;
          depth = x.depth;
          break;
        }

        if (x.depth < kMaxBfsPathLen - 1) {
          assert(!q.Full());
          q.Enqueue(BfsSlot(AltIndex(x.bucket, tag),
                            x.pathcode * kTagsPerBucket + slot, x.depth + 1));
        }
      }
    }

    if (depth != -1) {
      CuckooRecords path;
      depth =
          SetPathFromPathCode(thread_id, i1, i2, pathcode, depth, &path, pin);
      if (depth == 0) {
        return true;
      }

      for (; depth > 0; --depth) {
        int64_t bucket_idx{path[depth - 1].bucket_index};
        uint8_t slot_idx{path[depth - 1].slot_index};
        TableEntry e1 = kcas_.read_value(thread_id, *pin, &tables_[bucket_idx]);

        int64_t dest_bucket_idx = AltIndex(bucket_idx, GetTag(e1, slot_idx));
        // avoid a relocation on same bucket (KCAS does not support it).
        if (dest_bucket_idx == bucket_idx) {
          continue;
        }
        TableEntry e2 =
            kcas_.read_value(thread_id, *pin, &tables_[dest_bucket_idx]);
        int8_t dest_slot_idx = FindEmptySlot(e2);
        if (dest_slot_idx < 0) {
          // find another path
          goto path_discovery;
        }
        HelpRelocate(thread_id, bucket_idx, slot_idx, pin);
      }
    }

    return depth != -1;
  }

  void HelpRelocate(const size_t thread_id, int64_t bucket_idx, int8_t slot_idx,
                    ReclaimerPin<Reclaimer> *pin) {
    while (true) {
      TableEntry src = kcas_.read_value(thread_id, *pin, &tables_[bucket_idx]);
      uint8_t tag{GetTag(src, slot_idx)};
      if (tag == 0) {
        return;
      }

      int64_t alt_bucket_index = AltIndex(bucket_idx, tag);
      if (bucket_idx == alt_bucket_index) {
        // KCAS can not perform relocation in same bucket
        return;
      }
      TableEntry dst =
          kcas_.read_value(thread_id, *pin, &tables_[alt_bucket_index]);
      uint32_t ts1{GetCounter(src)}, ts2{GetCounter(dst)};
      int8_t dst_slot_idx = FindEmptySlot(dst);

      if (dst_slot_idx >= 0) {
        uint32_t n_cnt = ts1 > ts2 ? ts1 + 1 : ts2 + 1;
        if (src != kcas_.read_value(thread_id, *pin, &tables_[bucket_idx])) {
          continue;
        }

        auto desc = kcas_.create_descriptor(2, thread_id);
        TableEntry desired_dst{
            SetCounter(AddTag(dst, dst_slot_idx, tag), n_cnt)};
        TableEntry desired_src{SetCounter(ClearTag(src, slot_idx), ts1 + 1)};
        desc->add_value(&tables_[bucket_idx], src, desired_src);
        desc->add_value(&tables_[alt_bucket_index], dst, desired_dst);

        kcas_.cas(thread_id, *pin, desc);
      }
      return;
    }
  }

  int8_t SetPathFromPathCode(const size_t thread_id, size_t i1, size_t i2,
                             uint16_t pathcode, int8_t depth,
                             CuckooRecords *path,
                             ReclaimerPin<Reclaimer> *pin) {
    // Fill in the cuckoo path slots from the end to the beginning.
    for (int i = depth; i >= 0; i--) {
      (*path)[i].slot_index = pathcode % kTagsPerBucket;
      pathcode /= kTagsPerBucket;
    }
    // Fill in the cuckoo_path buckets and keys from the beginning to the
    // end, using the final pathcode to figure out which bucket the path
    // starts on. Since data could have been modified between slot_search
    // and the computation of the cuckoo path, this could be an invalid
    // cuckoo_path.
    CuckooRecord &first = (*path)[0];
    if (pathcode == 0) {
      first.bucket_index = i1;
    } else {
      assert(pathcode == 1);
      first.bucket_index = i2;
    }
    uint8_t tag;
    {
      TableEntry e{
          kcas_.read_value(thread_id, *pin, &tables_[first.bucket_index])};
      tag = GetTag(e, first.slot_index);
      if (tag == 0) {
        return 0;
      }
    }
    for (int i = 1; i <= depth; ++i) {
      CuckooRecord &curr = (*path)[i];
      const CuckooRecord &prev = (*path)[i - 1];
      curr.bucket_index = AltIndex(prev.bucket_index, tag);
      TableEntry e{
          kcas_.read_value(thread_id, *pin, &tables_[curr.bucket_index])};
      tag = GetTag(e, curr.slot_index);
      if (tag == 0) {
        return i;
      }
    }
    return depth;
  }

  // Member variables

  size_t num_buckets_;
  std::unique_ptr<typename KCAS::template KCASEntry<TableEntry>[]> tables_;
  TwoIndependentMultiplyShift hasher_;
  Reclaimer reclaimer_;
  KCAS kcas_;
};