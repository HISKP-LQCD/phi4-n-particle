#pragma once

#include <omp.h>

#include <random>
#include <iostream>

class GeneratorSingleton {
  public:
    typedef std::mt19937 Generator;

    static GeneratorSingleton &get_instance() {
      static GeneratorSingleton instance;
      return instance;
    }

    static Generator &get() {
      auto const tid = omp_get_thread_num();
      auto &instance = get_instance();
      auto &generator = instance.generators_.at(tid);
      return generator;
    }

  private:
    GeneratorSingleton() : generators_(omp_get_max_threads()) {
      reseed(0);
    }

    void reseed(int const base_seed) {
      for (int i = 0; i < generators_.size(); ++i) {
        generators_[i].seed(base_seed + i);
      }
    }

    std::vector<Generator> generators_;
};
