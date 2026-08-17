#include "EWALD/ewald.hpp"

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>

namespace md {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

}  // namespace

EwaldCalculator::EwaldCalculator(EwaldParameters parameters)
    : parameters_(std::move(parameters)) {
    if (!(parameters_.box_length > 0.0) || !(parameters_.cutoff > 0.0) ||
        parameters_.cutoff > 0.5 * parameters_.box_length) {
        throw std::invalid_argument("Ewald box_length and cutoff are invalid");
    }
    if (parameters_.charges.empty()) {
        throw std::invalid_argument("at least one site charge is required");
    }
    if (parameters_.k_max < 1 || parameters_.k_squared_limit < 1) {
        throw std::invalid_argument("reciprocal-space limits must be positive");
    }
    if (parameters_.alpha == 0.0) {
        parameters_.alpha = 5.0 / parameters_.box_length;
    }
    if (!(parameters_.alpha > 0.0)) {
        throw std::invalid_argument("alpha must be positive");
    }

    const double wave_scale = 2.0 * pi / parameters_.box_length;
    const double attenuation = 1.0 / (4.0 * parameters_.alpha * parameters_.alpha);
    for (int kx = -parameters_.k_max; kx <= parameters_.k_max; ++kx) {
        for (int ky = -parameters_.k_max; ky <= parameters_.k_max; ++ky) {
            for (int kz = -parameters_.k_max; kz <= parameters_.k_max; ++kz) {
                const int integer_norm = kx * kx + ky * ky + kz * kz;
                if (integer_norm == 0 || integer_norm >= parameters_.k_squared_limit) {
                    continue;
                }
                const Vec3 wave{
                    wave_scale * static_cast<double>(kx),
                    wave_scale * static_cast<double>(ky),
                    wave_scale * static_cast<double>(kz),
                };
                const double wave_norm = norm_squared(wave);
                reciprocal_vectors_.push_back(
                    {wave, std::exp(-wave_norm * attenuation) / wave_norm});
            }
        }
    }
}

void EwaldCalculator::validate_configuration(const Configuration& positions) const {
    if (positions.empty()) {
        throw std::invalid_argument("the Ewald configuration must contain a molecule");
    }
    for (const auto& molecule : positions) {
        if (molecule.size() != parameters_.charges.size()) {
            throw std::invalid_argument("every molecule must contain one position per charge");
        }
    }
}

EwaldResult EwaldCalculator::compute(const Configuration& positions) const {
    validate_configuration(positions);

    EwaldResult result;
    result.forces.reserve(positions.size());
    for (const auto& molecule : positions) {
        result.forces.emplace_back(molecule.size());
    }

    const double cutoff_squared = parameters_.cutoff * parameters_.cutoff;
    for (std::size_t first_molecule = 0; first_molecule + 1 < positions.size();
         ++first_molecule) {
        for (std::size_t second_molecule = first_molecule + 1;
             second_molecule < positions.size(); ++second_molecule) {
            for (std::size_t first_site = 0; first_site < parameters_.charges.size();
                 ++first_site) {
                for (std::size_t second_site = 0; second_site < parameters_.charges.size();
                     ++second_site) {
                    const Vec3 displacement = minimum_image(
                        positions[first_molecule][first_site] -
                            positions[second_molecule][second_site],
                        parameters_.box_length);
                    const double distance_squared = norm_squared(displacement);
                    if (distance_squared >= cutoff_squared) {
                        continue;
                    }
                    if (distance_squared <= std::numeric_limits<double>::min()) {
                        throw std::domain_error("overlapping Ewald sites have singular energy");
                    }

                    const double inverse_distance_squared = 1.0 / distance_squared;
                    const double distance = std::sqrt(distance_squared);
                    const double alpha_distance = parameters_.alpha * distance;
                    const double charge_product = parameters_.charges[first_site] *
                                                  parameters_.charges[second_site];
                    const double screened = std::erfc(alpha_distance);
                    result.real_space_energy += charge_product * screened / distance;

                    double force_coefficient = charge_product *
                        (screened / (distance_squared * distance) +
                         2.0 * parameters_.alpha / std::sqrt(pi) *
                             std::exp(-alpha_distance * alpha_distance) /
                             distance_squared);

                    if (first_site == 0 && second_site == 0) {
                        const double inverse_sixth =
                            inverse_distance_squared * inverse_distance_squared *
                            inverse_distance_squared;
                        const double repulsive = parameters_.lennard_jones_a *
                                                 inverse_sixth * inverse_sixth;
                        const double attractive = parameters_.lennard_jones_b * inverse_sixth;
                        result.real_space_energy += repulsive - attractive;
                        force_coefficient +=
                            (12.0 * repulsive - 6.0 * attractive) * inverse_distance_squared;
                    }

                    const Vec3 pair_force = displacement * force_coefficient;
                    result.forces[first_molecule][first_site] += pair_force;
                    result.forces[second_molecule][second_site] -= pair_force;
                }
            }
        }
    }

    const double volume = parameters_.box_length * parameters_.box_length *
                          parameters_.box_length;
    for (const auto& reciprocal : reciprocal_vectors_) {
        std::complex<double> structure_factor{};
        for (const auto& molecule : positions) {
            for (std::size_t site = 0; site < molecule.size(); ++site) {
                const double phase = dot(reciprocal.wave_vector, molecule[site]);
                structure_factor += parameters_.charges[site] *
                                    std::complex<double>(std::cos(phase), std::sin(phase));
            }
        }

        result.reciprocal_energy +=
            (2.0 * pi / volume) * reciprocal.weight * std::norm(structure_factor);

        for (std::size_t molecule = 0; molecule < positions.size(); ++molecule) {
            for (std::size_t site = 0; site < positions[molecule].size(); ++site) {
                const double phase = dot(reciprocal.wave_vector, positions[molecule][site]);
                const std::complex<double> phase_factor(std::cos(phase), std::sin(phase));
                const double imaginary_product =
                    std::imag(phase_factor * std::conj(structure_factor));
                const double coefficient = 4.0 * pi / volume * reciprocal.weight *
                                           parameters_.charges[site] * imaginary_product;
                result.forces[molecule][site] += reciprocal.wave_vector * coefficient;
            }
        }
    }

    double squared_charge_sum = 0.0;
    for (double charge : parameters_.charges) {
        squared_charge_sum += charge * charge;
    }
    result.self_energy = -parameters_.alpha / std::sqrt(pi) * squared_charge_sum *
                         static_cast<double>(positions.size());
    return result;
}

std::size_t EwaldCalculator::reciprocal_vector_count() const noexcept {
    return reciprocal_vectors_.size();
}

double EwaldCalculator::alpha() const noexcept {
    return parameters_.alpha;
}

}  // namespace md
