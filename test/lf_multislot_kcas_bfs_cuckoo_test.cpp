#include "filters/lf_multislot_kcas_bfs_cuckoo.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <thread>
#include <unordered_set>

#include "allocators/std_allocator.h"
#include "primitives/brown_kcas.h"
#include "reclaimer/leaky.h"

constexpr uint64_t kNumThreads = 4;
constexpr uint64_t kTotalItemCount = 1000000;
constexpr uint64_t kSingleThreadItemCount = kTotalItemCount / kNumThreads;

TEST(LfMultislotKCASBfsCuckooTest, SingleThreadInsertThenLookup) {
  using Reclaimer = LeakyReclaimer<StdAllocator>;
  using KCAS = BrownKCAS<StdAllocator, LeakyReclaimer<StdAllocator>>;
  LockfreeMultislotKCASBfsCuckoo<uint64_t, uint8_t, Reclaimer, KCAS> filter(
      kTotalItemCount, 1);
  uint64_t num_inserted = 0;
  for (uint64_t i = 0; i < kTotalItemCount; i++, num_inserted++) {
    if (!filter.Insert(i, 0)) {
      break;
    }
  }
  std::cout << "Total items count: " << kTotalItemCount
            << ", inserted items count: " << num_inserted << "\n";

  // Check if previously inserted items are in the filter, expected true for
  // all items
  for (uint64_t i = 0; i < num_inserted; i++) {
    ASSERT_TRUE(filter.Lookup(i, 0)) << "item: " << i;
  }

  // Check non-existing items, a few false positives expected
  uint64_t total_queries = 0;
  uint64_t false_queries = 0;
  for (uint64_t i = kTotalItemCount; i < 2 * kTotalItemCount; i++) {
    if (filter.Lookup(i, 0)) {
      false_queries++;
    }
    total_queries++;
  }

  // Output the measured false positive rate
  std::cout << "false positive rate is "
            << 100.0 * false_queries / total_queries << "%\n";
}

TEST(LfMultislotKCASBfsCuckooTest, MultiThreadInsertThenLookup) {
  static_assert(kTotalItemCount % kNumThreads == 0,
                "total count can not be divided evenly among threads.");

  using Reclaimer = LeakyReclaimer<StdAllocator>;
  using KCAS = BrownKCAS<StdAllocator, LeakyReclaimer<StdAllocator>>;
  using Filter =
      LockfreeMultislotKCASBfsCuckoo<uint64_t, uint8_t, Reclaimer, KCAS>;
  Filter filter(kTotalItemCount, kNumThreads);
  std::array<std::thread, kNumThreads> threads;
  std::array<std::atomic<uint64_t>, kNumThreads> inserted;

  auto thread_insert_func = [](Filter &filter, std::atomic<uint64_t> &inserted,
                               uint64_t id) {
    uint64_t start = kSingleThreadItemCount * id;
    uint64_t end = kSingleThreadItemCount * (id + 1);
    inserted = kSingleThreadItemCount;
    for (uint64_t i = start; i < end; ++i) {
      if (!filter.Insert(i, id)) {
        inserted = i - start;
        return;
      }
    }
  };

  auto thread_lookup_func =
      [](Filter &filter, const std::atomic<uint64_t> &inserted, uint64_t id) {
        uint64_t start = kSingleThreadItemCount * id;
        uint64_t end = start + inserted;
        for (uint64_t i = start; i < end; ++i) {
          ASSERT_TRUE(filter.Lookup(i, id));
        }
      };

  for (uint64_t id = 0; id < kNumThreads; ++id) {
    threads[id] = std::thread(thread_insert_func, std::ref(filter),
                              std::ref(inserted[id]), id);
  }
  for (auto &th : threads) {
    th.join();
  }

  uint64_t sum{};
  for (uint64_t id = 0; id < kNumThreads; ++id) {
    sum += inserted[id];
    threads[id] = std::thread(thread_lookup_func, std::ref(filter),
                              std::cref(inserted[id]), id);
  }
  for (auto &th : threads) {
    th.join();
  }

  std::cout << "Set size: " << sum << '\n';
}

struct SimpleConcurrentSet {
  void Insert(uint64_t item) {
    std::lock_guard<std::mutex> guard(mu_);
    sets_.insert(item);
  }

  bool Lookup(uint64_t item) {
    std::lock_guard<std::mutex> guard(mu_);
    return sets_.find(item) != sets_.end();
  }

  bool Delete(uint64_t item) {
    std::lock_guard<std::mutex> guard(mu_);
    auto it = sets_.find(item);
    if (it == sets_.end()) {
      return false;
    }
    sets_.erase(it);
    return true;
  }

  std::mutex mu_;
  std::unordered_set<uint64_t> sets_;
};

TEST(LfMultislotKCASBfsCuckooTest, MultiThreadInsertAndLookup) {
  SimpleConcurrentSet set;
  using Reclaimer = LeakyReclaimer<StdAllocator>;
  using KCAS = BrownKCAS<StdAllocator, LeakyReclaimer<StdAllocator>>;
  using Filter =
      LockfreeMultislotKCASBfsCuckoo<uint64_t, uint8_t, Reclaimer, KCAS>;
  Filter filter(kTotalItemCount, kNumThreads);
  std::array<std::thread, kNumThreads> threads;

  auto thread_insert_and_lookup_func =
      [](Filter &filter, SimpleConcurrentSet &set, uint64_t id) {
        std::random_device rd;  // obtain a seed for the random number engine
        std::mt19937 gen(
            rd());  // standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<uint64_t> item_gen(0, kTotalItemCount);

        for (uint64_t i = 0; i < kSingleThreadItemCount; ++i) {
          uint64_t item = item_gen(gen);
          if (i % 2 == 0) {
            // try insert
            if (filter.Insert(item, id)) {
              set.Insert(item);
            }
          } else {
            // try lookup
            if (set.Lookup(item)) {
              ASSERT_TRUE(filter.Lookup(item, id));
            }
          }
        }
      };

  for (uint64_t id = 0; id < kNumThreads; ++id) {
    threads[id] = std::thread(thread_insert_and_lookup_func, std::ref(filter),
                              std::ref(set), id);
  }
  for (auto &th : threads) {
    th.join();
  }
}

TEST(LfMultislotKCASBfsCuckooTest, MultiThreadInsertTwiceThenDelete) {
  static_assert(kTotalItemCount % kNumThreads == 0,
                "total count can not be divided evenly among threads.");

  SimpleConcurrentSet set;
  using Reclaimer = LeakyReclaimer<StdAllocator>;
  using KCAS = BrownKCAS<StdAllocator, LeakyReclaimer<StdAllocator>>;
  using Filter =
      LockfreeMultislotKCASBfsCuckoo<size_t, uint8_t, Reclaimer, KCAS>;
  Filter filter(kTotalItemCount * 2, kNumThreads);
  std::array<std::thread, kNumThreads> threads;

  auto thread_insert_delete_func = [](Filter &filter, SimpleConcurrentSet &set,
                                      size_t id) {
    size_t start = kSingleThreadItemCount * id;
    size_t end = kSingleThreadItemCount * (id + 1);
    for (size_t i = start; i < end; ++i) {
      if (filter.Insert(i, id) && filter.Insert(i, id)) {
        set.Insert(i);
      }
    }

    for (size_t i = start; i < end; ++i) {
      if (set.Lookup(i)) {
        ASSERT_TRUE(filter.Delete(i, id));
      }
    }

    for (size_t i = start; i < end; ++i) {
      if (set.Lookup(i)) {
        ASSERT_TRUE(filter.Lookup(i, id));
      }
    }
  };

  for (size_t id = 0; id < kNumThreads; ++id) {
    threads[id] = std::thread(thread_insert_delete_func, std::ref(filter),
                              std::ref(set), id);
  }
  for (auto &th : threads) {
    th.join();
  }
}

TEST(LfMultislotKCASBfsCuckooTest, MultiThreadLookupAndDelete) {
  static_assert(kTotalItemCount % kNumThreads == 0,
                "total count can not be divided evenly among threads.");

  SimpleConcurrentSet set;
  using Reclaimer = LeakyReclaimer<StdAllocator>;
  using KCAS = BrownKCAS<StdAllocator, LeakyReclaimer<StdAllocator>>;
  using Filter =
      LockfreeMultislotKCASBfsCuckoo<size_t, uint8_t, Reclaimer, KCAS>;
  Filter filter(kTotalItemCount * 2, kNumThreads);
  std::array<std::thread, kNumThreads> threads;

  auto thread_lookup_delete_func = [](Filter &filter, SimpleConcurrentSet &set,
                                      size_t id) {
    size_t start = kSingleThreadItemCount * id;
    size_t end = kSingleThreadItemCount * (id + 1);
    for (size_t i = start; i < end; ++i) {
      if (filter.Insert(i, id)) {
        set.Insert(i);
      }
    }

    for (size_t i = start; i < end; ++i) {
      switch (i % 2) {
        case 0: {
          if (set.Lookup(i)) {
            ASSERT_TRUE(filter.Lookup(i, id));
          }
          break;
        }
        case 1: {
          if (set.Delete(i)) {
            ASSERT_TRUE(filter.Delete(i, id));
          }
          break;
        }
      }
    }

    for (size_t i = start; i < end; ++i) {
      if (set.Lookup(i)) {
        ASSERT_TRUE(filter.Lookup(i, id));
      }
    }
  };

  for (size_t id = 0; id < kNumThreads; ++id) {
    threads[id] = std::thread(thread_lookup_delete_func, std::ref(filter),
                              std::ref(set), id);
  }
  for (auto &th : threads) {
    th.join();
  }
}