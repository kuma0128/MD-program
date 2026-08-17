#include "EWALD/ewald.hpp"
#include "GAUSS/gaussian.hpp"
#include "MC/monte_carlo.hpp"
#include "RATTLE/rattle.hpp"
#include "ROTATE/rotate.hpp"
#include "SHAKE/shake.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_near(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            message + ": expected " + std::to_string(expected) + ", got " +
            std::to_string(actual));
    }
}

void test_gaussian() {
    md::GaussianGenerator generator(1234);
    constexpr std::size_t samples = 200'000;
    double sum = 0.0;
    double squared_sum = 0.0;
    for (std::size_t index = 0; index < samples; ++index) {
        const double value = generator();
        sum += value;
        squared_sum += value * value;
    }
    const double mean = sum / static_cast<double>(samples);
    const double variance = squared_sum / static_cast<double>(samples) - mean * mean;
    require(std::abs(mean) < 0.01, "Gaussian mean is outside deterministic tolerance");
    require(std::abs(variance - 1.0) < 0.02,
            "Gaussian variance is outside deterministic tolerance");
}

void test_rotations() {
    const auto identity = md::euler_zyz(0.0, 0.0, 0.0);
    const md::Vec3 point{1.2, -3.4, 5.6};
    const auto unchanged = md::rotate(identity, point);
    require_near(md::norm(unchanged - point), 0.0, 1.0e-12, "zero Euler rotation");

    std::mt19937_64 engine(99);
    const auto matrix = md::quaternion_matrix(md::random_quaternion(engine));
    const md::Vec3 first{1.0, 2.0, 3.0};
    const md::Vec3 second{-2.0, 0.5, 4.0};
    const auto rotated_first = md::rotate(matrix, first);
    const auto rotated_second = md::rotate(matrix, second);
    require_near(md::norm(rotated_first), md::norm(first), 1.0e-12,
                 "quaternion rotation length");
    require_near(md::dot(rotated_first, rotated_second), md::dot(first, second), 1.0e-12,
                 "quaternion rotation inner product");
}

void test_ewald_force_gradient() {
    md::EwaldParameters parameters;
    parameters.box_length = 10.0;
    parameters.cutoff = 4.9;
    parameters.charges = {1.0};
    parameters.lennard_jones_a = 0.01;
    parameters.lennard_jones_b = 0.02;
    md::EwaldCalculator calculator(parameters);
    require(calculator.reciprocal_vector_count() > 0, "Ewald reciprocal vectors are empty");

    md::Configuration positions{{{2.0, 2.0, 2.0}}, {{5.0, 2.3, 2.0}}};
    const auto result = calculator.compute(positions);
    const md::Vec3 net_force = result.forces[0][0] + result.forces[1][0];
    require(md::norm(net_force) < 1.0e-12, "Ewald net force is not zero");
    require(std::isfinite(result.total_energy()), "Ewald energy is not finite");

    constexpr double step = 1.0e-5;
    auto plus = positions;
    auto minus = positions;
    plus[0][0].x += step;
    minus[0][0].x -= step;
    const double numerical_force =
        -(calculator.compute(plus).total_energy() - calculator.compute(minus).total_energy()) /
        (2.0 * step);
    require_near(result.forces[0][0].x, numerical_force, 2.0e-7,
                 "Ewald force/energy gradient");
}

void test_monte_carlo_and_rdf() {
    std::vector<md::Vec3> positions{
        {1.0, 1.0, 1.0}, {3.0, 1.0, 1.0}, {1.0, 3.0, 1.0}, {3.0, 3.0, 1.0}};
    md::MonteCarloParameters parameters;
    parameters.box_length = 10.0;
    parameters.cutoff = 4.0;
    parameters.max_displacement = 0.2;
    parameters.minimum_distance = 0.5;
    md::MonteCarloNVT simulation(positions, parameters, 42);
    for (int sweep = 0; sweep < 100; ++sweep) {
        simulation.sweep();
    }
    require(simulation.accepted_moves() + simulation.rejected_moves() == 400,
            "Monte Carlo move count");
    require_near(simulation.energy(), simulation.recompute_energy(), 1.0e-11,
                 "Monte Carlo incremental energy");
    require(simulation.acceptance_ratio() >= 0.0 && simulation.acceptance_ratio() <= 1.0,
            "Monte Carlo acceptance ratio");

    md::RdfAccumulator rdf(10, 5.0, 10.0);
    rdf.sample({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
    require_near(rdf.mean_distance(), 1.0, 1.0e-12, "RDF mean distance");
    require_near(rdf.distance_variance(), 0.0, 1.0e-12, "RDF variance");
    const auto points = rdf.result();
    require(points.size() == 10, "RDF bin count");
    require(points[2].value > 0.0, "RDF expected bin is empty");
}

void test_shake() {
    const std::vector<md::Vec3> current{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    const std::vector<md::Vec3> previous{{0.01, 0.0, 0.0}, {0.99, 0.0, 0.0}};
    const std::vector<md::Vec3> forces(2);
    const std::vector<double> masses{1.0, 2.0};
    const std::vector<md::Constraint> constraints{{0, 1, 1.0}};
    const auto result = md::shake_verlet(
        current, previous, forces, masses, constraints, 0.01, 20.0, 1.0e-12, 100);
    const double distance = md::norm(md::minimum_image(
        result.positions[0] - result.positions[1], 20.0));
    require_near(distance, 1.0, 1.0e-10, "SHAKE bond distance");
    require(result.iterations > 0, "SHAKE did not iterate");
    require(std::isfinite(result.kinetic_energy), "SHAKE kinetic energy is not finite");
}

void test_rattle() {
    const double height = std::sqrt(3.0) / 2.0;
    md::RattleState state{
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, height, 0.0}},
        {{0.2, 0.1, 0.0}, {-0.1, 0.05, 0.0}, {0.0, -0.2, 0.0}},
    };
    const std::vector<md::Vec3> forces{{0.1, 0.0, 0.0}, {0.0, 0.2, 0.0},
                                        {-0.1, -0.2, 0.0}};
    const std::vector<double> masses{16.0, 1.0, 1.0};
    const std::vector<md::Constraint> constraints{{0, 1, 1.0}, {1, 2, 1.0},
                                                   {2, 0, 1.0}};
    md::rattle_move_a(state, forces, masses, constraints, 0.002, 20.0, 1.0e-10, 200);
    md::rattle_move_b(state, forces, masses, constraints, 0.002, 20.0, 1.0e-10, 200);
    for (const auto& constraint : constraints) {
        const auto displacement = md::minimum_image(
            state.positions[constraint.first] - state.positions[constraint.second], 20.0);
        const auto relative_velocity =
            state.velocities[constraint.first] - state.velocities[constraint.second];
        require_near(md::norm(displacement), constraint.distance, 1.0e-9,
                     "RATTLE bond distance");
        require(std::abs(md::dot(displacement, relative_velocity)) < 1.0e-7,
                "RATTLE velocity constraint");
    }
    require(std::isfinite(md::kinetic_energy(state, masses)),
            "RATTLE kinetic energy is not finite");
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"Gaussian", test_gaussian},
        {"rotations", test_rotations},
        {"Ewald", test_ewald_force_gradient},
        {"Monte Carlo and RDF", test_monte_carlo_and_rdf},
        {"SHAKE", test_shake},
        {"RATTLE", test_rattle},
    };

    std::size_t failures = 0;
    for (const auto& [name, test] : tests) {
        try {
            test();
            std::cout << "[PASS] " << name << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[FAIL] " << name << ": " << error.what() << '\n';
        }
    }
    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All tests passed\n";
    return 0;
}
