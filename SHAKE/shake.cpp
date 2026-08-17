#include "SHAKE/shake.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace md {

ShakeResult shake_verlet(
    const std::vector<Vec3>& current_positions,
    const std::vector<Vec3>& previous_positions,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses,
    const std::vector<Constraint>& constraints,
    double time_step,
    double box_length,
    double relative_tolerance,
    std::size_t max_iterations) {
    if (current_positions.size() != previous_positions.size()) {
        throw std::invalid_argument("current and previous positions must have equal sizes");
    }
    validate_particle_data(current_positions.size(), forces, masses);
    validate_constraints(current_positions.size(), constraints);
    if (!(time_step > 0.0) || !(box_length > 0.0) || !(relative_tolerance > 0.0) ||
        max_iterations == 0) {
        throw std::invalid_argument("invalid SHAKE integration parameters");
    }

    ShakeResult result;
    result.positions.resize(current_positions.size());
    for (std::size_t particle = 0; particle < current_positions.size(); ++particle) {
        result.positions[particle] =
            2.0 * current_positions[particle] - previous_positions[particle] +
            forces[particle] * (time_step * time_step / masses[particle]);
    }

    bool converged = constraints.empty();
    for (std::size_t iteration = 1; iteration <= max_iterations && !converged; ++iteration) {
        converged = true;
        for (const auto& constraint : constraints) {
            const Vec3 predicted_displacement = minimum_image(
                result.positions[constraint.first] - result.positions[constraint.second],
                box_length);
            const double target_squared = constraint.distance * constraint.distance;
            const double error = norm_squared(predicted_displacement) - target_squared;
            if (std::abs(error) <= relative_tolerance * target_squared) {
                continue;
            }

            const Vec3 reference_displacement = minimum_image(
                current_positions[constraint.first] - current_positions[constraint.second],
                box_length);
            const double inverse_first_mass = 1.0 / masses[constraint.first];
            const double inverse_second_mass = 1.0 / masses[constraint.second];
            const double denominator = 2.0 * (inverse_first_mass + inverse_second_mass) *
                                       dot(reference_displacement, predicted_displacement);
            if (std::abs(denominator) <= std::numeric_limits<double>::epsilon() *
                                             target_squared) {
                throw std::runtime_error("SHAKE correction is singular");
            }

            const double multiplier = error / denominator;
            result.positions[constraint.first] -=
                reference_displacement * (multiplier * inverse_first_mass);
            result.positions[constraint.second] +=
                reference_displacement * (multiplier * inverse_second_mass);
            converged = false;
        }
        result.iterations = iteration;
    }
    if (!converged) {
        throw std::runtime_error("SHAKE did not converge within max_iterations");
    }

    result.velocities.resize(current_positions.size());
    for (std::size_t particle = 0; particle < current_positions.size(); ++particle) {
        const Vec3 displacement = minimum_image(
            result.positions[particle] - previous_positions[particle], box_length);
        result.velocities[particle] = displacement / (2.0 * time_step);
        result.kinetic_energy +=
            0.5 * masses[particle] * norm_squared(result.velocities[particle]);
    }
    return result;
}

}  // namespace md
