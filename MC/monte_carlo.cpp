#include "MC/monte_carlo.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace md {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

}  // namespace

MonteCarloNVT::MonteCarloNVT(
    std::vector<Vec3> positions,
    MonteCarloParameters parameters,
    std::uint64_t seed)
    : positions_(std::move(positions)), parameters_(parameters), engine_(seed) {
    if (positions_.empty()) {
        throw std::invalid_argument("Monte Carlo requires at least one particle");
    }
    if (!(parameters_.box_length > 0.0) || !(parameters_.cutoff > 0.0) ||
        parameters_.cutoff > 0.5 * parameters_.box_length || !(parameters_.sigma > 0.0) ||
        !(parameters_.epsilon > 0.0) || !(parameters_.beta > 0.0) ||
        !(parameters_.max_displacement > 0.0) || parameters_.minimum_distance < 0.0) {
        throw std::invalid_argument("invalid Monte Carlo parameters");
    }
    for (auto& position : positions_) {
        position = wrap_position(position, parameters_.box_length);
    }
    energy_ = recompute_energy();
    if (!std::isfinite(energy_)) {
        throw std::invalid_argument("initial configuration contains overlapping particles");
    }
}

double MonteCarloNVT::pair_energy(const Vec3& first, const Vec3& second) const {
    const Vec3 displacement = minimum_image(first - second, parameters_.box_length);
    const double distance_squared = norm_squared(displacement);
    const double minimum_squared = parameters_.minimum_distance * parameters_.minimum_distance;
    if (distance_squared <= std::max(minimum_squared, std::numeric_limits<double>::min())) {
        return std::numeric_limits<double>::infinity();
    }
    if (distance_squared >= parameters_.cutoff * parameters_.cutoff) {
        return 0.0;
    }

    const double sigma_squared_over_distance_squared =
        parameters_.sigma * parameters_.sigma / distance_squared;
    const double sigma_sixth = sigma_squared_over_distance_squared *
                               sigma_squared_over_distance_squared *
                               sigma_squared_over_distance_squared;
    return 4.0 * parameters_.epsilon * (sigma_sixth * sigma_sixth - sigma_sixth);
}

double MonteCarloNVT::local_energy(const Vec3& position, std::size_t excluded_index) const {
    double result = 0.0;
    for (std::size_t index = 0; index < positions_.size(); ++index) {
        if (index != excluded_index) {
            result += pair_energy(position, positions_[index]);
        }
    }
    return result;
}

void MonteCarloNVT::sweep() {
    for (std::size_t index = 0; index < positions_.size(); ++index) {
        const Vec3 current = positions_[index];
        const double old_energy = local_energy(current, index);
        Vec3 trial = current;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            trial[axis] += (2.0 * uniform_(engine_) - 1.0) * parameters_.max_displacement;
        }
        trial = wrap_position(trial, parameters_.box_length);

        const double new_energy = local_energy(trial, index);
        const double energy_change = new_energy - old_energy;
        const bool accepted = std::isfinite(new_energy) &&
                              (energy_change <= 0.0 ||
                               std::exp(-parameters_.beta * energy_change) > uniform_(engine_));
        if (accepted) {
            positions_[index] = trial;
            energy_ += energy_change;
            ++accepted_;
        } else {
            ++rejected_;
        }
    }
}

double MonteCarloNVT::tune_displacement(double lower_target, double upper_target) {
    if (!(lower_target >= 0.0) || !(upper_target <= 1.0) ||
        !(lower_target < upper_target)) {
        throw std::invalid_argument("acceptance targets must satisfy 0 <= lower < upper <= 1");
    }
    const double ratio = acceptance_ratio();
    if (accepted_ + rejected_ != 0) {
        if (ratio > upper_target) {
            parameters_.max_displacement *= 1.05;
        } else if (ratio < lower_target) {
            parameters_.max_displacement *= 0.95;
        }
    }
    accepted_ = 0;
    rejected_ = 0;
    return ratio;
}

const std::vector<Vec3>& MonteCarloNVT::positions() const noexcept {
    return positions_;
}

double MonteCarloNVT::energy() const noexcept {
    return energy_;
}

double MonteCarloNVT::recompute_energy() const {
    double result = 0.0;
    for (std::size_t first = 0; first + 1 < positions_.size(); ++first) {
        for (std::size_t second = first + 1; second < positions_.size(); ++second) {
            result += pair_energy(positions_[first], positions_[second]);
        }
    }
    return result;
}

double MonteCarloNVT::acceptance_ratio() const noexcept {
    const std::size_t attempted = accepted_ + rejected_;
    return attempted == 0 ? 0.0
                          : static_cast<double>(accepted_) / static_cast<double>(attempted);
}

std::size_t MonteCarloNVT::accepted_moves() const noexcept {
    return accepted_;
}

std::size_t MonteCarloNVT::rejected_moves() const noexcept {
    return rejected_;
}

double MonteCarloNVT::max_displacement() const noexcept {
    return parameters_.max_displacement;
}

RdfAccumulator::RdfAccumulator(
    std::size_t bin_count,
    double maximum_radius,
    double box_length)
    : histogram_(bin_count), maximum_radius_(maximum_radius), box_length_(box_length) {
    if (bin_count == 0 || !(maximum_radius > 0.0) || !(box_length > 0.0) ||
        maximum_radius > 0.5 * box_length) {
        throw std::invalid_argument("invalid RDF dimensions");
    }
}

void RdfAccumulator::sample(const std::vector<Vec3>& positions) {
    if (positions.empty()) {
        throw std::invalid_argument("cannot sample an empty RDF configuration");
    }
    if (samples_ == 0) {
        particle_count_ = positions.size();
    } else if (positions.size() != particle_count_) {
        throw std::invalid_argument("RDF particle count cannot change between samples");
    }

    const double bin_width = maximum_radius_ / static_cast<double>(histogram_.size());
    for (std::size_t first = 0; first + 1 < positions.size(); ++first) {
        for (std::size_t second = first + 1; second < positions.size(); ++second) {
            const double distance = norm(minimum_image(
                positions[first] - positions[second], box_length_));
            distance_sum_ += distance;
            distance_squared_sum_ += distance * distance;
            ++distance_count_;
            if (distance < maximum_radius_) {
                const auto bin = static_cast<std::size_t>(distance / bin_width);
                histogram_[std::min(bin, histogram_.size() - 1)] += 1.0;
            }
        }
    }
    ++samples_;
}

std::vector<RdfPoint> RdfAccumulator::result() const {
    if (samples_ == 0) {
        throw std::logic_error("RDF has no samples");
    }
    std::vector<RdfPoint> points;
    points.reserve(histogram_.size());
    const double volume = box_length_ * box_length_ * box_length_;
    const double bin_width = maximum_radius_ / static_cast<double>(histogram_.size());
    const double pair_count = 0.5 * static_cast<double>(particle_count_) *
                              static_cast<double>(particle_count_ - 1);

    for (std::size_t bin = 0; bin < histogram_.size(); ++bin) {
        const double inner = static_cast<double>(bin) * bin_width;
        const double outer = inner + bin_width;
        const double shell_volume = 4.0 * pi / 3.0 *
                                    (outer * outer * outer - inner * inner * inner);
        const double ideal_count = pair_count * shell_volume / volume;
        const double value = ideal_count == 0.0
                                 ? 0.0
                                 : histogram_[bin] /
                                       (static_cast<double>(samples_) * ideal_count);
        points.push_back({inner + 0.5 * bin_width, value});
    }
    return points;
}

double RdfAccumulator::mean_distance() const noexcept {
    return distance_count_ == 0 ? 0.0
                                : distance_sum_ / static_cast<double>(distance_count_);
}

double RdfAccumulator::distance_variance() const noexcept {
    if (distance_count_ == 0) {
        return 0.0;
    }
    const double mean = mean_distance();
    return distance_squared_sum_ / static_cast<double>(distance_count_) - mean * mean;
}

std::size_t RdfAccumulator::sample_count() const noexcept {
    return samples_;
}

}  // namespace md
