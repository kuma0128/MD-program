#pragma once

#include "MD/common.hpp"

#include <cstddef>
#include <vector>

namespace md {

struct RattleState {
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
};

std::size_t rattle_move_a(
    RattleState& state,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance = 1.0e-8,
    std::size_t max_iterations = 100);

std::size_t rattle_move_b(
    RattleState& state,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance = 1.0e-8,
    std::size_t max_iterations = 100);

double kinetic_energy(const RattleState& state, const std::vector<double>& masses);

}  // namespace md
