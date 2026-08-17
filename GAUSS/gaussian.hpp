#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace md {

class GaussianGenerator {
public:
    explicit GaussianGenerator(std::uint64_t seed = 80) : engine_(seed) {}

    double operator()() {
        if (has_cached_) {
            has_cached_ = false;
            return cached_;
        }

        double first = 0.0;
        double second = 0.0;
        double radius_squared = 0.0;
        do {
            first = 2.0 * uniform_(engine_) - 1.0;
            second = 2.0 * uniform_(engine_) - 1.0;
            radius_squared = first * first + second * second;
        } while (radius_squared >= 1.0 ||
                 radius_squared <= std::numeric_limits<double>::min());

        const double scale = std::sqrt(-2.0 * std::log(radius_squared) / radius_squared);
        cached_ = first * scale;
        has_cached_ = true;
        return second * scale;
    }

private:
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    bool has_cached_{};
    double cached_{};
};

}  // namespace md
