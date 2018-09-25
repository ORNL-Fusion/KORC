#include <random>

class random {
    std::mt19937_64 engine;
    std::uniform_real_distribution<double> dist;
    
    public:
    random(uint64_t seed) : engine(seed), dist(0.0, 1.0) {};
    double get_number() {
        return dist(engine);
    }
};

extern "C" {
    void *random_construct(int seed) {
        return new class random(static_cast<uint64_t> (seed));
    }
    
    double random_get_number(void *r) {
        return static_cast<class random *> (r)->get_number();
    }
    
    void random_destroy(void *r) {
        delete  static_cast<class random *> (r);
    }
}