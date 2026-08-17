#include "RATTLE/rattle.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace md {
namespace {

void validate_state(
    const RattleState& state,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance,
    std::size_t max_iterations) {
    if (state.positions.size() != state.velocities.size()) {
        throw std::invalid_argument("RATTLE positions and velocities must have equal sizes");
    }
    validate_particle_data(state.positions.size(), forces, masses);
    validate_constraints(state.positions.size(), constraints);
    if (!(time_step > 0.0) || !(box_length > 0.0) || !(relative_tolerance > 0.0) ||
        max_iterations == 0) {
        throw std::invalid_argument("invalid RATTLE integration parameters");
    }
}

}  // namespace

std::size_t rattle_move_a(
    RattleState& state,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance,
    std::size_t max_iterations) {
    validate_state(state, forces, masses, constraints, time_step, box_length,
                   relative_tolerance, max_iterations);
    const auto reference_positions = state.positions;
    for (std::size_t particle = 0; particle < state.positions.size(); ++particle) {
        const Vec3 half_acceleration =
            forces[particle] * (0.5 * time_step / masses[particle]);
        state.velocities[particle] += half_acceleration;
        state.positions[particle] += state.velocities[particle] * time_step;
    }

    if (constraints.empty()) {
        return 0;
    }
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        bool converged = true;
        for (const auto& constraint : constraints) {
            const Vec3 predicted_displacement = minimum_image(
                state.positions[constraint.first] - state.positions[constraint.second],
                box_length);
            const double target_squared = constraint.distance * constraint.distance;
            const double error = norm_squared(predicted_displacement) - target_squared;
            if (std::abs(error) <= relative_tolerance * target_squared) {
                continue;
            }

            const Vec3 reference_displacement = minimum_image(
                reference_positions[constraint.first] -
                    reference_positions[constraint.second],
                box_length);
            const double inverse_first_mass = 1.0 / masses[constraint.first];
            const double inverse_second_mass = 1.0 / masses[constraint.second];
            const double denominator = 2.0 * (inverse_first_mass + inverse_second_mass) *
                                       dot(reference_displacement, predicted_displacement);
            if (std::abs(denominator) <= std::numeric_limits<double>::epsilon() *
                                             target_squared) {
                throw std::runtime_error("RATTLE position correction is singular");
            }

            const double multiplier = error / denominator;
            const Vec3 first_correction =
                -reference_displacement * (multiplier * inverse_first_mass);
            const Vec3 second_correction =
                reference_displacement * (multiplier * inverse_second_mass);
            state.positions[constraint.first] += first_correction;
            state.positions[constraint.second] += second_correction;
            state.velocities[constraint.first] += first_correction / time_step;
            state.velocities[constraint.second] += second_correction / time_step;
            converged = false;
        }
        if (converged) {
            return iteration;
        }
    }
    throw std::runtime_error("RATTLE move A did not converge within max_iterations");
}

std::size_t rattle_move_b(
    RattleState& state,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance,
    std::size_t max_iterations) {
    validate_state(state, forces, masses, constraints, time_step, box_length,
                   relative_tolerance, max_iterations);
    for (std::size_t particle = 0; particle < state.positions.size(); ++particle) {
        state.velocities[particle] +=
            forces[particle] * (0.5 * time_step / masses[particle]);
    }

    if (constraints.empty()) {
        return 0;
    }
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        bool converged = true;
        for (const auto& constraint : constraints) {
            const Vec3 displacement = minimum_image(
                state.positions[constraint.first] - state.positions[constraint.second],
                box_length);
            const Vec3 relative_velocity =
                state.velocities[constraint.first] - state.velocities[constraint.second];
            const double projection = dot(displacement, relative_velocity);
            const double target_squared = constraint.distance * constraint.distance;
            if (std::abs(projection) <=
                relative_tolerance * target_squared / time_step) {
                continue;
            }

            const double inverse_first_mass = 1.0 / masses[constraint.first];
            const double inverse_second_mass = 1.0 / masses[constraint.second];
            const double denominator = (inverse_first_mass + inverse_second_mass) *
                                       norm_squared(displacement);
            if (denominator <= std::numeric_limits<double>::min()) {
                throw std::runtime_error("RATTLE velocity correction is singular");
            }
            const double multiplier = projection / denominator;
            state.velocities[constraint.first] -=
                displacement * (multiplier * inverse_first_mass);
            state.velocities[constraint.second] +=
                displacement * (multiplier * inverse_second_mass);
            converged = false;
        }
        if (converged) {
            return iteration;
        }
    }
    throw std::runtime_error("RATTLE move B did not converge within max_iterations");
}

double kinetic_energy(const RattleState& state, const std::vector<double>& masses) {
    if (state.velocities.size() != masses.size()) {
        throw std::invalid_argument("RATTLE velocities and masses must have equal sizes");
    }
    double result = 0.0;
    for (std::size_t particle = 0; particle < masses.size(); ++particle) {
        if (!(masses[particle] > 0.0)) {
            throw std::invalid_argument("all masses must be positive");
        }
        result += 0.5 * masses[particle] * norm_squared(state.velocities[particle]);
    }
    return result;
}

}  // namespace md
