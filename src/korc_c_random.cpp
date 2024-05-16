#include <chrono>
#include <random>
#include <stdexcept>

#ifdef DOUBLE_PRECISION
typedef double real_type;
#else
typedef float real_type;
#endif

template<class DIST>
class random_dist {
    std::mt19937_64 engine;
    DIST dist;

    public:
    random_dist(const uint64_t offset) :
    engine(static_cast<uint64_t> (std::chrono::system_clock::to_time_t(
                                      std::chrono::system_clock::now())) +
                                  offset),
    dist(0.0, 1.0) {}

    void set_dist(const real_type low, const real_type high) {
        dist = DIST(low, high);
    }
    void set_seed(const uint64_t seed) {
        engine.seed(seed);
    }
    real_type get_number() {
        return dist(engine);
    }
};

typedef random_dist<std::uniform_real_distribution<real_type> > random_U;
typedef random_dist<std::normal_distribution<real_type> > random_N;

extern "C" {
    void *random_construct_U(int seed) {
        return new random_U(static_cast<uint64_t> (seed));
    }

    void *random_construct_N(int seed) {
        return new random_N(static_cast<uint64_t> (seed));
    }
  
    real_type random_get_number_U(void *r) {
        return static_cast<random_U *> (r)->get_number();
    }

    real_type random_get_number_N(void *r) {
        return static_cast<random_N *> (r)->get_number();
    }

    void random_destroy_U(void *r) {
        delete static_cast<random_U *> (r);
    }

    void random_destroy_N(void *r) {
        delete static_cast<random_N *> (r);
    }

    void random_set_dist_U(void *r,
                           const real_type low,
                           const real_type high) {
        static_cast<random_U *> (r)->set_dist(low, high);
    }

    void random_set_dist_N(void *r,
                           const real_type low,
                           const real_type high) {
        static_cast<random_N *> (r)->set_dist(low, high);
    }
    void random_set_seed_U(void *r, int seed) {
        static_cast<random_U *> (r)->set_seed(static_cast<uint64_t> (seed));
    }

    void random_set_seed_N(void *r, int seed) {
        static_cast<random_N *> (r)->set_seed(static_cast<uint64_t> (seed));
    }
}
