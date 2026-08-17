#pragma once

#include "MD/common.hpp"

#include <cstddef>
#include <vector>

namespace md {

struct ShakeResult {
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
    double kinetic_energy{};
    std::size_t iterations{};
};

ShakeResult shake_verlet(
    const std::vector<Vec3>& current_positions,
    const std::vector<Vec3>& previous_positions,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance = 1.0e-6,
    std::size_t max_iterations = 100);

}  // namespace md
