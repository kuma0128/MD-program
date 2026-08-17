#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace md {

struct Vec3 {
    double x{};
    double y{};
    double z{};

    double& operator[](std::size_t index) {
        if (index > 2) {
            throw std::out_of_range("Vec3 index must be 0, 1, or 2");
        }
        return index == 0 ? x : (index == 1 ? y : z);
    }

    const double& operator[](std::size_t index) const {
        if (index > 2) {
            throw std::out_of_range("Vec3 index must be 0, 1, or 2");
        }
        return index == 0 ? x : (index == 1 ? y : z);
    }
};

inline Vec3& operator+=(Vec3& lhs, const Vec3& rhs) {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
}

inline Vec3& operator-=(Vec3& lhs, const Vec3& rhs) {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    return lhs;
}

inline Vec3 operator+(Vec3 lhs, const Vec3& rhs) {
    return lhs += rhs;
}

inline Vec3 operator-(Vec3 lhs, const Vec3& rhs) {
    return lhs -= rhs;
}

inline Vec3 operator-(const Vec3& value) {
    return {-value.x, -value.y, -value.z};
}

inline Vec3 operator*(const Vec3& value, double scalar) {
    return {value.x * scalar, value.y * scalar, value.z * scalar};
}

inline Vec3 operator*(double scalar, const Vec3& value) {
    return value * scalar;
}

inline Vec3 operator/(const Vec3& value, double scalar) {
    return {value.x / scalar, value.y / scalar, value.z / scalar};
}

inline double dot(const Vec3& lhs, const Vec3& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline double norm_squared(const Vec3& value) {
    return dot(value, value);
}

inline double norm(const Vec3& value) {
    return std::sqrt(norm_squared(value));
}

inline Vec3 minimum_image(Vec3 displacement, double box_length) {
    if (!(box_length > 0.0)) {
        throw std::invalid_argument("box_length must be positive");
    }
    displacement.x -= std::nearbyint(displacement.x / box_length) * box_length;
    displacement.y -= std::nearbyint(displacement.y / box_length) * box_length;
    displacement.z -= std::nearbyint(displacement.z / box_length) * box_length;
    return displacement;
}

inline double wrap_coordinate(double coordinate, double box_length) {
    coordinate = std::fmod(coordinate, box_length);
    return coordinate < 0.0 ? coordinate + box_length : coordinate;
}

inline Vec3 wrap_position(Vec3 position, double box_length) {
    position.x = wrap_coordinate(position.x, box_length);
    position.y = wrap_coordinate(position.y, box_length);
    position.z = wrap_coordinate(position.z, box_length);
    return position;
}

using Molecule = std::vector<Vec3>;
using Configuration = std::vector<Molecule>;

struct Constraint {
    std::size_t first{};
    std::size_t second{};
    double distance{};
};

inline void validate_particle_data(
    std::size_t particle_count,
    const std::vector<Vec3>& forces,
    const std::vector<double>& masses) {
    if (forces.size() != particle_count || masses.size() != particle_count) {
        throw std::invalid_argument("positions, forces, and masses must have equal sizes");
    }
    for (double mass : masses) {
        if (!(mass > 0.0)) {
            throw std::invalid_argument("all masses must be positive");
        }
    }
}

inline void validate_constraints(
    std::size_t particle_count,
    const std::vector<Constraint>& constraints) {
    for (const auto& constraint : constraints) {
        if (constraint.first >= particle_count || constraint.second >= particle_count ||
            constraint.first == constraint.second || !(constraint.distance > 0.0)) {
            throw std::invalid_argument("invalid distance constraint");
        }
    }
}

}  // namespace md
