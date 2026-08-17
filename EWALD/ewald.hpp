#pragma once

#include "MD/common.hpp"

#include <cstddef>
#include <vector>

namespace md {

struct EwaldParameters {
    double box_length{};
    double cutoff{};
    std::vector<double> charges;
    double lennard_jones_a{};
    double lennard_jones_b{};
    double alpha{};
    int k_max{5};
    int k_squared_limit{20};
};

struct EwaldResult {
    Configuration forces;
    double real_space_energy{};
    double reciprocal_energy{};
    double self_energy{};

    double total_energy() const {
        return real_space_energy + reciprocal_energy + self_energy;
    }
};

class EwaldCalculator {
public:
    explicit EwaldCalculator(EwaldParameters parameters);

    EwaldResult compute(const Configuration& positions) const;
    std::size_t reciprocal_vector_count() const noexcept;
    double alpha() const noexcept;

private:
    struct ReciprocalVector {
        Vec3 wave_vector;
        double weight{};
    };

    EwaldParameters parameters_;
    std::vector<ReciprocalVector> reciprocal_vectors_;

    void validate_configuration(const Configuration& positions) const;
};

}  // namespace md
