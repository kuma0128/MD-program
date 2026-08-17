#include "ROTATE/rotate.hpp"

#include <cmath>
#include <stdexcept>

namespace md {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

}  // namespace

Matrix3 euler_zyz(double phi, double theta, double psi) {
    const double cosine_phi = std::cos(phi);
    const double sine_phi = std::sin(phi);
    const double cosine_theta = std::cos(theta);
    const double sine_theta = std::sin(theta);
    const double cosine_psi = std::cos(psi);
    const double sine_psi = std::sin(psi);

    return {{
        {{cosine_phi * cosine_theta * cosine_psi - sine_phi * sine_psi,
          -cosine_phi * cosine_theta * sine_psi - sine_phi * cosine_psi,
          cosine_phi * sine_theta}},
        {{sine_phi * cosine_theta * cosine_psi + cosine_phi * sine_psi,
          -sine_phi * cosine_theta * sine_psi + cosine_phi * cosine_psi,
          sine_phi * sine_theta}},
        {{-sine_theta * cosine_psi, sine_theta * sine_psi, cosine_theta}},
    }};
}

Matrix3 quaternion_matrix(const Quaternion& quaternion) {
    const double magnitude_squared = quaternion.w * quaternion.w +
                                     quaternion.x * quaternion.x +
                                     quaternion.y * quaternion.y +
                                     quaternion.z * quaternion.z;
    if (!(magnitude_squared > 0.0)) {
        throw std::invalid_argument("a zero quaternion cannot define a rotation");
    }

    const double inverse_magnitude = 1.0 / std::sqrt(magnitude_squared);
    const double w = quaternion.w * inverse_magnitude;
    const double x = quaternion.x * inverse_magnitude;
    const double y = quaternion.y * inverse_magnitude;
    const double z = quaternion.z * inverse_magnitude;

    return {{
        {{1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),
          2.0 * (x * z + y * w)}},
        {{2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z),
          2.0 * (y * z - x * w)}},
        {{2.0 * (x * z - y * w), 2.0 * (y * z + x * w),
          1.0 - 2.0 * (x * x + y * y)}},
    }};
}

Quaternion random_quaternion(std::mt19937_64& engine) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double first = uniform(engine);
    const double second = uniform(engine);
    const double third = uniform(engine);
    const double root_one_minus_first = std::sqrt(1.0 - first);
    const double root_first = std::sqrt(first);

    return {
        root_first * std::cos(2.0 * pi * third),
        root_one_minus_first * std::sin(2.0 * pi * second),
        root_one_minus_first * std::cos(2.0 * pi * second),
        root_first * std::sin(2.0 * pi * third),
    };
}

Vec3 rotate(const Matrix3& matrix, const Vec3& point) {
    return {
        matrix[0][0] * point.x + matrix[0][1] * point.y + matrix[0][2] * point.z,
        matrix[1][0] * point.x + matrix[1][1] * point.y + matrix[1][2] * point.z,
        matrix[2][0] * point.x + matrix[2][1] * point.y + matrix[2][2] * point.z,
    };
}

std::vector<Vec3> rotate(const Matrix3& matrix, const std::vector<Vec3>& points) {
    std::vector<Vec3> result;
    result.reserve(points.size());
    for (const auto& point : points) {
        result.push_back(rotate(matrix, point));
    }
    return result;
}

}  // namespace md
