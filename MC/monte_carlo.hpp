#pragma once

#include "MD/common.hpp"

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace md {

struct MonteCarloParameters {
    double box_length{};
    double cutoff{};
    double sigma{1.0};
    double epsilon{1.0};
    double beta{1.0};
    double max_displacement{0.1};
    double minimum_distance{};
};

class MonteCarloNVT {
public:
    MonteCarloNVT(
        std::vector<Vec3> positions,
        MonteCarloParameters parameters,
        std::uint64_t seed = 80);

    void sweep();
    double tune_displacement(double lower_target = 0.45, double upper_target = 0.55);

    const std::vector<Vec3>& positions() const noexcept;
    double energy() const noexcept;
    double recompute_energy() const;
    double acceptance_ratio() const noexcept;
    std::size_t accepted_moves() const noexcept;
    std::size_t rejected_moves() const noexcept;
    double max_displacement() const noexcept;

private:
    double pair_energy(const Vec3& first, const Vec3& second) const;
    double local_energy(const Vec3& position, std::size_t excluded_index) const;

    std::vector<Vec3> positions_;
    MonteCarloParameters parameters_;
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    double energy_{};
    std::size_t accepted_{};
    std::size_t rejected_{};
};

struct RdfPoint {
    double radius{};
    double value{};
};

class RdfAccumulator {
public:
    RdfAccumulator(std::size_t bin_count, double maximum_radius, double box_length);

    void sample(const std::vector<Vec3>& positions);
    std::vector<RdfPoint> result() const;
    double mean_distance() const noexcept;
    double distance_variance() const noexcept;
    std::size_t sample_count() const noexcept;

private:
    std::vector<double> histogram_;
    double maximum_radius_{};
    double box_length_{};
    std::size_t samples_{};
    std::size_t particle_count_{};
    double distance_sum_{};
    double distance_squared_sum_{};
    std::size_t distance_count_{};
};

}  // namespace md
