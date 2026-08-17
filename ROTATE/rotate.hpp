#pragma once

#include "MD/common.hpp"

#include <array>
#include <random>
#include <vector>

namespace md {

struct Quaternion {
    double w{1.0};
    double x{};
    double y{};
    double z{};
};

using Matrix3 = std::array<std::array<double, 3>, 3>;

Matrix3 euler_zyz(double phi, double theta, double psi);
Matrix3 quaternion_matrix(const Quaternion& quaternion);
Quaternion random_quaternion(std::mt19937_64& engine);
Vec3 rotate(const Matrix3& matrix, const Vec3& point);
std::vector<Vec3> rotate(const Matrix3& matrix, const std::vector<Vec3>& points);

}  // namespace md
